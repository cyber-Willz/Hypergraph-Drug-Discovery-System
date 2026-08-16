//! `krylov_ds`: matrix-free Krylov subspace data structures shared across
//! this workspace's spectral-graph tools.
//!
//! - [`operator::LinearOperator`] -- the matrix-free operator trait.
//! - [`Arnoldi`] -- Arnoldi iteration ([`arnoldi::Arnoldi`], re-exported
//!   here since it's the crate's main entry point).
//! - [`eig::arnoldi_real_ritz_pairs`] -- read real Ritz (eigenvalue,
//!   eigenvector) pairs off an [`arnoldi::ArnoldiResult`].

pub mod arnoldi;
pub mod eig;
pub mod operator;

pub use arnoldi::Arnoldi;
