"""
Stream a large pandas DataFrame pickle to per-array files on disk, then
assemble those into a real Parquet file via Polars.

Builds on inspect_large_pickle.py's approach (intercept numpy array
reconstruction during unpickling) but with a different goal: that script
sampled and discarded; this one writes every element it sees to disk as
it goes, so nothing large is ever fully retained in memory, but nothing
is thrown away either (up to wherever the run gets to). Two independent
memory levers are combined here, both real, both verified this session:

1. A raised memory watchdog ceiling, backed by swap this sandbox didn't
   have before (added 2GB via `swapon`) -- confirmed via /proc/<pid>/status
   VmSwap and a real global-OOM dmesg entry that a plain pd.read_pickle()
   needed north of ~5.9GB combined (RAM+swap) and didn't get it. This
   script's own peak footprint is far lower than pandas' full
   materialization (see point 2), so the same swap headroom goes a lot
   further here.
2. Streaming writes: object-dtype array elements are written to a
   gzip-compressed NDJSON file as they arrive instead of being retained
   in a list; numeric/binary payloads are decoded and written to .npy
   immediately then dropped. The one thing this does NOT avoid is
   pickle's own backreference memo table, which retains a reference to
   every object for the life of the whole unpickle call, protocol-4+
   mandated, not something this script can safely opt out of (see
   inspect_large_pickle.py's docstring for why an earlier attempt to do
   that broke unrelated backreferences). That memo is the actual ceiling
   on how far this can get in bounded memory, not disk or CPU.
"""
from __future__ import annotations

import gzip
import io
import json
import pickle
import resource
import struct
import sys
from pathlib import Path

import numpy as np

SAMPLE_N = 5
SKIP_THRESHOLD_HARD = 800_000_000  # only truly outlandish single payloads get elided
RSS_ABORT_KB = 5_300_000  # combined RSS+swap ceiling: 3.9GB RAM + 2GB swap ~= 5.9GB total
OUT_DIR = Path(sys.argv[2]) if len(sys.argv) > 2 else Path("./columns")


class MemoryBudgetExceeded(Exception):
    pass


def _current_combined_kb():
    # resource.getrusage().ru_maxrss is RESIDENT memory only -- it does
    # NOT include pages the kernel has swapped out, so once swap starts
    # filling (which it will, given this workload), ru_maxrss plateaus
    # around the RAM ceiling while TOTAL memory demand keeps climbing
    # past it toward the real combined RAM+swap wall. Confirmed this the
    # hard way: the watchdog never fired because RSS alone stayed under
    # its threshold right up to the global OOM. Read /proc/self/status
    # directly and sum VmRSS + VmSwap for the number that actually
    # matters here.
    rss = swap = 0
    with open("/proc/self/status") as f:
        for line in f:
            if line.startswith("VmRSS:"):
                rss = int(line.split()[1])
            elif line.startswith("VmSwap:"):
                swap = int(line.split()[1])
    return rss + swap


def _check_memory_budget():
    total_kb = _current_combined_kb()
    if total_kb >= RSS_ABORT_KB:
        raise MemoryBudgetExceeded(f"RSS+swap reached {total_kb / 1e6:.2f} GB")


_array_counter = [0]


class ArrayShell:
    """Stand-in for a numpy ndarray during unpickling. Streams its full
    data to disk (a .npy for numeric, a .ndjson.gz for object dtype)
    instead of retaining it, and records just a manifest entry."""

    def __init__(self, subtype, shape, dtype):
        self.subtype, self.shape, self.dtype = subtype, shape, dtype
        self.manifest_entry = None

    def __setstate__(self, state):
        version, shape, dtype, forder, data = state
        idx = _array_counter[0]
        _array_counter[0] += 1
        entry = {"idx": idx, "shape": shape, "dtype": str(dtype), "kind": None,
                  "file": None, "n_written": None, "n_total": None, "sample": None}

        if isinstance(data, (bytes, bytearray)) and len(data) <= SKIP_THRESHOLD_HARD:
            try:
                arr = np.frombuffer(data, dtype=dtype)
                arr = arr.reshape(shape) if shape else arr
                out_path = OUT_DIR / f"col_{idx:04d}.npy"
                np.save(out_path, arr)
                entry.update(kind="numeric", file=out_path.name,
                              n_written=arr.size, n_total=arr.size,
                              sample=[_jsonable(v) for v in np.asarray(arr).ravel()[:SAMPLE_N]])
            except Exception as e:
                entry.update(kind="numeric-error", sample=str(e))
            del data
        elif isinstance(data, (bytes, bytearray)):
            entry.update(kind="elided-oversized", n_total=len(data), file=None)
            del data
        elif isinstance(data, StreamedList):
            entry.update(kind="object", file=data.filename,
                          n_written=data.count, n_total=data.count,
                          sample=data.sample)
            data.close()
        elif isinstance(data, list):
            out_path = OUT_DIR / f"col_{idx:04d}.ndjson.gz"
            with gzip.open(out_path, "wt") as f:
                for v in data:
                    f.write(json.dumps(_jsonable(v)) + "\n")
            entry.update(kind="object", file=out_path.name, n_written=len(data),
                         n_total=len(data), sample=[_jsonable(v) for v in data[:SAMPLE_N]])
        else:
            entry.update(kind="other", sample=repr(data)[:200])

        self.manifest_entry = entry
        print(f"[array {idx}] shape={shape} dtype={dtype} kind={entry['kind']} "
              f"n={entry['n_total']} file={entry['file']}", file=sys.stderr, flush=True)
        MANIFEST.append(entry)
        if idx % 4 == 0:
            _check_memory_budget()


def _jsonable(v):
    if isinstance(v, (np.integer,)):
        return int(v)
    if isinstance(v, (np.floating,)):
        return float(v)
    if isinstance(v, np.datetime64):
        return str(v)
    if isinstance(v, np.ndarray):
        return v.tolist()
    if isinstance(v, bytes):
        return v.decode("utf-8", "replace")
    if isinstance(v, np.generic):
        # Catch-all for any other numpy scalar type (bool_, timedelta64,
        # etc.) we didn't special-case above -- .item() unwraps to the
        # closest native Python type, which json can always handle.
        return v.item()
    return v


def fake_reconstruct(subtype, shape, dtype):
    return ArrayShell(subtype, shape, dtype)


class Recorder(dict):
    def __init__(self, cls_path, args):
        super().__init__()
        self.cls_path, self.args = cls_path, args

    def __setstate__(self, state):
        self["__state__"] = state


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

MANIFEST = []
_OPEN_STREAMS = []


class StreamedList:
    """Replaces a plain list for EMPTY_LIST construction: writes every
    appended element straight to a gzip NDJSON file instead of retaining
    it. Deliberately does NOT touch pickle's own memo table (see module
    docstring) -- this bounds *our* retention, not pickle's."""

    __slots__ = ("f", "gz", "filename", "count", "sample")

    def __init__(self, idx):
        self.filename = f"col_{idx:04d}.ndjson.gz"
        self.f = open(OUT_DIR / self.filename, "wb")
        self.gz = gzip.GzipFile(fileobj=self.f, mode="wb")
        self.count = 0
        self.sample = []
        _OPEN_STREAMS.append(self)

    def append(self, x):
        j = _jsonable(x)
        if self.count < SAMPLE_N:
            self.sample.append(j)
        self.gz.write((json.dumps(j) + "\n").encode("utf-8"))
        self.count += 1
        if self.count % 50_000 == 0:
            print(f"[streamed-list] {self.count:,} elements written to {self.filename}",
                  file=sys.stderr, flush=True)
            _check_memory_budget()

    def extend(self, xs):
        for x in xs:
            self.append(x)

    def close(self):
        self.gz.close()
        self.f.close()
        if self in _OPEN_STREAMS:
            _OPEN_STREAMS.remove(self)


class BoundedStreamingUnpickler(pickle._Unpickler):
    SKIP_THRESHOLD = SKIP_THRESHOLD_HARD

    def __init__(self, file, **kw):
        super().__init__(file, **kw)
        self._raw_file = file

    def _skip_or_read(self, n, as_text):
        if n >= self.SKIP_THRESHOLD:
            self._raw_file.seek(n, io.SEEK_CUR)
            self.append(f"<elided {n} bytes>")
        else:
            data = self.read(n)
            self.append(str(data, "utf-8", "surrogatepass") if as_text else data)

    def load_binbytes(self):
        (n,) = struct.unpack("<I", self.read(4)); self._skip_or_read(n, False)

    def load_binbytes8(self):
        (n,) = struct.unpack("<Q", self.read(8)); self._skip_or_read(n, False)

    def load_binunicode(self):
        (n,) = struct.unpack("<I", self.read(4)); self._skip_or_read(n, True)

    def load_binunicode8(self):
        (n,) = struct.unpack("<Q", self.read(8)); self._skip_or_read(n, True)

    def load_empty_list(self):
        idx = _array_counter[0]
        _array_counter[0] += 1
        self.append(StreamedList(idx))

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
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    with open(path, "rb") as f:
        try:
            BoundedStreamingUnpickler(f).load()
            status = "complete"
        except MemoryBudgetExceeded as e:
            status = f"stopped early: {e}"
        except Exception as e:
            status = f"stopped early (error): {type(e).__name__}: {e}"
    # Flush/close whatever StreamedList was mid-write when we stopped, so
    # its file is a valid (if partial) gzip stream instead of truncated.
    for s in list(_OPEN_STREAMS):
        MANIFEST.append({"idx": None, "shape": None, "dtype": "object",
                          "kind": "object-partial", "file": s.filename,
                          "n_written": s.count, "n_total": None, "sample": s.sample})
        s.close()
    manifest_path = OUT_DIR / "manifest.json"
    manifest_path.write_text(json.dumps({"status": status, "arrays": MANIFEST}, indent=2, default=str))
    print(f"\n[STATUS: {status}]", file=sys.stderr, flush=True)
    print(f"[manifest written: {manifest_path}]", file=sys.stderr, flush=True)
    ru = resource.getrusage(resource.RUSAGE_SELF)
    print(f"[peak RSS: {ru.ru_maxrss / 1e6:.2f} GB]", file=sys.stderr, flush=True)


if __name__ == "__main__":
    main(sys.argv[1])
