"""
Memory-bounded introspection of a large pandas DataFrame pickle.

A plain pd.read_pickle() on this file gets OOM-killed on a ~3.9GB-RAM
sandbox (the file is 4.8GB on disk). This script intercepts numpy array
reconstruction *during* unpickling and, for any array whose full byte
size exceeds a threshold, keeps only its dtype/shape (and a small sample
of the first N elements it can decode cheaply) instead of retaining the
full buffer -- so peak memory stays bounded by the single largest
in-flight array's raw bytes, not the sum of the whole DataFrame.
"""
from __future__ import annotations

import io
import pickle
import resource
import struct
import sys
import numpy as np

MAX_KEEP_BYTES = 20_000_000  # keep full data for arrays under ~20MB
SAMPLE_N = 20
# Large BINBYTES/BINUNICODE payloads (>= pickle's 64KB frame-size target)
# are written UNFRAMED and unbuffered by CPython's pickler
# (_Framer.write_large_bytes), which means we can skip over them on the
# raw file with .seek() instead of ever reading them into memory. This is
# what makes a >4GB pickle introspectable in <4GB of RAM: only small
# opcodes and arrays under this threshold ever get materialized.
SKIP_THRESHOLD = 5_000_000
# Hard memory ceiling. This sandbox has ~3.9GB total RAM; stop well
# before that rather than let the OOM killer do it -- a controlled abort
# preserves everything printed/discovered so far, a SIGKILL preserves
# nothing. ru_maxrss is KB on Linux.
RSS_ABORT_KB = 2_200_000


class MemoryBudgetExceeded(Exception):
    """Raised when RSS approaches the sandbox ceiling. Not a bug -- an
    intentional early stop so we can report partial-but-real results
    instead of getting SIGKILLed with nothing to show."""


def _check_memory_budget():
    rss_kb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    if rss_kb >= RSS_ABORT_KB:
        raise MemoryBudgetExceeded(f"RSS reached {rss_kb / 1e6:.2f} GB")


class ArrayShell:
    """Stand-in for a numpy ndarray during unpickling. Records dtype/shape
    always; keeps full data only if small, else a decoded head sample."""

    def __init__(self, subtype, shape, dtype):
        self.subtype = subtype
        self.shape = shape
        self.dtype = dtype
        self.full = None
        self.sample = None
        self.nbytes_est = None

    def __setstate__(self, state):
        # ndarray.__reduce__ state tuple: (version, shape, dtype, fortran_order, data)
        version, shape, dtype, forder, data = state
        self.shape = shape
        self.dtype = dtype
        _check_memory_budget()
        try:
            nbytes = int(np.prod(shape)) * np.dtype(dtype).itemsize if shape else np.dtype(dtype).itemsize
        except Exception:
            nbytes = None
        self.nbytes_est = nbytes

        if isinstance(data, SkipMarker):
            # The BINBYTES payload itself was never read (raw file seek);
            # nothing to sample from here without a second, seek-back pass.
            self.sample = f"<elided {data.nbytes:,} bytes, not read>"
        elif isinstance(data, (bytes, bytearray)):
            if nbytes is not None and nbytes <= MAX_KEEP_BYTES:
                try:
                    arr = np.frombuffer(data, dtype=dtype)
                    arr = arr.reshape(shape) if shape else arr
                    self.full = arr.copy()
                except Exception as e:
                    self.sample = f"<decode error: {e}>"
            else:
                try:
                    itemsize = np.dtype(dtype).itemsize
                    head = np.frombuffer(data[: itemsize * SAMPLE_N], dtype=dtype)
                    self.sample = head.copy()
                except Exception as e:
                    self.sample = f"<decode error: {e}>"
            del data  # drop the (potentially huge) raw bytes ASAP
        elif isinstance(data, list):
            # object-dtype arrays pickle their elements as a plain list
            # (possibly a CappedList if it was large -- see CappedList).
            total = getattr(data, "total_appended", len(data))
            if total > len(data):
                self.sample = f"<{total:,} elements, kept first {len(data)}: {data[:SAMPLE_N]!r}>"
            else:
                self.full = data if len(data) <= SAMPLE_N * 50 else None
                self.sample = data[:SAMPLE_N]
        else:
            self.sample = repr(data)[:500]

        print(f"[array] shape={shape} dtype={dtype} sample={self._sample_repr()}",
              file=sys.stderr, flush=True)

    def _sample_repr(self):
        if self.full is not None:
            vals = self.full[:SAMPLE_N] if hasattr(self.full, "__getitem__") else self.full
            return f"FULL kept: {list(vals)!r}"
        return repr(self.sample)[:400]

    def describe(self):
        if self.full is not None:
            return f"dtype={self.dtype} shape={self.shape} FULL kept ({self.nbytes_est} bytes)"
        return f"dtype={self.dtype} shape={self.shape} ~{self.nbytes_est} bytes -- sample={self.sample!r}"


def fake_reconstruct(subtype, shape, dtype):
    return ArrayShell(subtype, shape, dtype)


class Recorder(dict):
    """Generic stand-in for pandas internal objects we don't need fully
    reconstructed (BlockManager, blocks, index wrappers, categoricals,
    extension arrays). Captures constructor args positionally."""

    def __init__(self, cls_path, args):
        super().__init__()
        self.cls_path = cls_path
        self.args = args
        _check_memory_budget()
        if cls_path.endswith("_new_Index"):
            print(f"[index] _new_Index args={args!r}", file=sys.stderr, flush=True)

    def __setstate__(self, state):
        self["__state__"] = state

    def __repr__(self):
        return f"Recorder({self.cls_path}, nargs={len(self.args)})"


INTERCEPT = {
    ("numpy.core.multiarray", "_reconstruct"): fake_reconstruct,
    ("numpy._core.multiarray", "_reconstruct"): fake_reconstruct,
}

PASSTHROUGH_RECORD = {
    ("pandas.core.internals.managers", "_unpickle_block"),
    ("pandas._libs.internals", "_unpickle_block"),
    ("pandas._libs.arrays", "__pyx_unpickle_NDArrayBacked"),
    ("pandas.core.arrays.categorical", "Categorical"),
    ("pandas.core.arrays.string_", "StringArray"),
    ("pandas.core.indexes.base", "_new_Index"),
}


class SkipMarker:
    """Stands in for a large bytes/str payload we deliberately never read
    into memory -- we seek past it on the raw file instead."""
    __slots__ = ("nbytes",)

    def __init__(self, nbytes):
        self.nbytes = nbytes

    def __repr__(self):
        return f"<elided {self.nbytes:,} bytes>"


class CappedList(list):
    """A list that silently stops retaining elements past MAX_KEEP but
    still reports how many were actually appended. Object-dtype numpy
    arrays (e.g. a column of many small Python strings) pickle as an
    EMPTY_LIST + batched APPENDS opcodes (batches of <=1000 items each,
    per pickle's internal _BATCH_SIZE) growing the SAME list object across
    the whole array -- for an array with millions of elements, that one
    list is the actual OOM culprit, not any single APPENDS call. Capping
    the list itself (while every individual per-batch items list pickle
    hands us stays small) is what bounds this."""

    MAX_KEEP = 200

    def __init__(self):
        super().__init__()
        self.total_appended = 0

    def append(self, x):
        self.total_appended += 1
        if len(self) < self.MAX_KEEP:
            super().append(x)
        if self.total_appended % 50_000 == 0:
            print(f"[capped-list] {self.total_appended:,} elements streamed so far "
                  f"(keeping first {self.MAX_KEEP})", file=sys.stderr, flush=True)
            _check_memory_budget()

    def extend(self, xs):
        for x in xs:
            self.append(x)


class BoundedUnpickler(pickle._Unpickler):
    """Pure-Python pickle.Unpickler subclass (not the C-accelerated one --
    we need overridable per-opcode methods). Two independent memory-bound
    mechanisms:

    1. find_class() intercepts numpy's array-reconstruction callable so
       small arrays get fully decoded but large ones only keep
       dtype/shape/a head sample (see ArrayShell above).
    2. load_binbytes/binbytes8/binunicode/binunicode8 intercept the raw
       opcode handlers themselves: for any payload at/above
       SKIP_THRESHOLD, seek the underlying file past it instead of
       reading -- these payloads are guaranteed unframed at this size
       (CPython's pickler writes anything >= 64KB via
       _Framer.write_large_bytes, unbuffered), so a raw seek is safe.
       This is what actually prevents the OOM: mechanism (1) alone still
       requires the full bytes object to exist momentarily before our
       __setstate__ can discard it, which is exactly the multi-GB spike
       that killed the plain pd.read_pickle() attempt.
    """

    def __init__(self, file, **kw):
        super().__init__(file, **kw)
        self._raw_file = file

    def _skip_or_read(self, n, as_text):
        if n >= SKIP_THRESHOLD:
            self._raw_file.seek(n, io.SEEK_CUR)
            self.append(SkipMarker(n))
        else:
            data = self.read(n)
            self.append(str(data, "utf-8", "surrogatepass") if as_text else data)

    def load_binbytes(self):
        (n,) = struct.unpack("<I", self.read(4))
        self._skip_or_read(n, as_text=False)

    def load_binbytes8(self):
        (n,) = struct.unpack("<Q", self.read(8))
        self._skip_or_read(n, as_text=False)

    def load_binunicode(self):
        (n,) = struct.unpack("<I", self.read(4))
        self._skip_or_read(n, as_text=True)

    def load_binunicode8(self):
        (n,) = struct.unpack("<Q", self.read(8))
        self._skip_or_read(n, as_text=True)

    def load_empty_list(self):
        self.append(CappedList())

    # NOTE: an earlier version of this also tried to cap self.memo's
    # growth (pickle's MEMOIZE opcode fires for every individual object,
    # including each element of a multi-million-row object column, well
    # before CappedList ever gets a chance to bound the *list*). That
    # approach was WRONG: memo indices are positional and backreferenced
    # by later opcodes (STACK_GLOBAL, further BINGETs) for perfectly
    # ordinary structural reasons (shared dtype objects, repeated class
    # references) that have nothing to do with which array they happen
    # to trail -- dropping "old" entries by a count-based heuristic
    # corrupted unrelated, still-needed backreferences and crashed with
    # `STACK_GLOBAL requires str`. There is no safe way to know at
    # MEMOIZE-time whether a given object will be referenced again later
    # without full lookahead, so memo is left untouched (default,
    # correct behavior) and the *actual* memory bound comes from the
    # watchdog below instead: stop before the ceiling, rather than try
    # to selectively forget.
    dispatch = pickle._Unpickler.dispatch.copy()
    dispatch[pickle.BINBYTES[0]] = load_binbytes
    dispatch[pickle.BINBYTES8[0]] = load_binbytes8
    dispatch[pickle.BINUNICODE[0]] = load_binunicode
    dispatch[pickle.BINUNICODE8[0]] = load_binunicode8
    dispatch[pickle.EMPTY_LIST[0]] = load_empty_list

    def find_class(self, module, name):
        key = (module, name)
        if key in INTERCEPT:
            return INTERCEPT[key]
        if key in PASSTHROUGH_RECORD:
            def factory(*args, __cls_path=f"{module}.{name}"):
                return Recorder(__cls_path, args)
            return factory
        try:
            return super().find_class(module, name)
        except Exception:
            def generic_factory(*args, __cls_path=f"{module}.{name}"):
                return Recorder(f"{module}.{name}", args)
            return generic_factory


def main(path):
    with open(path, "rb") as f:
        try:
            obj = BoundedUnpickler(f).load()
            print("TOP-LEVEL TYPE:", type(obj), file=sys.stderr)
            return obj
        except MemoryBudgetExceeded as e:
            print(f"\n[STOPPED EARLY: {e}] -- this is a controlled abort, not a "
                  f"crash. Everything printed above ([array]/[capped-list] lines) "
                  f"is real, fully-decoded schema info seen before the cutoff; "
                  f"columns/data after this point in the pickle stream were not "
                  f"reached.", file=sys.stderr, flush=True)
            return None


def resource_report():
    ru = resource.getrusage(resource.RUSAGE_SELF)
    print(f"[peak RSS: {ru.ru_maxrss / 1e6:.2f} GB]", file=sys.stderr)


if __name__ == "__main__":
    result = main(sys.argv[1])
    if result is not None:
        import pprint
        print("=== repr ===")
        pprint.pprint(result, depth=4)
    resource_report()
