# Heston Model Solver

A high-performance C solver for the Heston stochastic-volatility PDE. It prices European and American calls and puts on a 2-D `(S, v)` grid via the Crank–Nicolson scheme, computes six Greeks, and writes a large, compressed HDF5 dataset across thousands of parameter combinations — suitable as training data for ML/AI models.

## Overview

For each parameter set `(κ, θ, ξ, ρ, τ, r)` and strike `K`, the solver:

1. Builds non-uniform `S` and `v` grids that concentrate nodes where they matter most (near `S = K` and `v = 0`).
2. Assembles the spatial Heston operator
   `L = ½ v S² ∂ₛₛ + ρ ξ v S ∂ₛᵥ + ½ ξ² v ∂ᵥᵥ + r S ∂ₛ + κ(θ−v) ∂ᵥ − r I`
   as a CSR sparse matrix with a 9-point non-uniform stencil (boundary rows handled as Dirichlet identities; `v = 0` row reduced per the Feller-degeneracy of the PDE).
3. Forms `A = I − (Δt/2) L` and `B = I + (Δt/2) L` for Crank–Nicolson and factorises `A` **once** via UMFPACK.
4. Steps backward in time from the payoff at maturity for all four option types (European/American × Call/Put), enforcing the early-exercise constraint for American options.
5. Computes Δ, Γ, Θ, Vega, Vanna and Volga with 2nd-order non-uniform finite-difference stencils.
6. Streams each result block into `dataset_Options.h5`.

## Output

A single HDF5 file `dataset_Options.h5` containing:

```
/european/call    (extendible 2-D dataset, float32, gzip-4)
/european/put
/american/call
/american/put
/grid/S           non-uniform asset-price grid
attributes:       N_S, N_v, S_min, S_max, v_min
```

Each row of every dataset has 18 columns:

```
kappa, theta, xi, rho, r, tau, K, feller, feller_ok,
S, v,
prix, delta, gamma, theta_gr, vega, vanna, volga
```

`feller = 2κθ/ξ²` and the `feller_ok` flag (1 if `feller ≥ 1`) are stored alongside each row; sets with `feller < 0.2` are skipped entirely.

## Build

### Dependencies

- A C compiler with OpenMP (`-fopenmp`)
- HDF5 (serial), discovered via `pkg-config hdf5-serial` or `pkg-config hdf5`
- BLAS and LAPACK (`-lblas -llapack`)
- SuiteSparse — at least UMFPACK and its transitive deps (AMD, CHOLMOD, COLAMD, CCOLAMD, CAMD, `SuiteSparse_config`)

SuiteSparse is **not** bundled — install it through your system package manager or from source:

```bash
# Debian / Ubuntu
sudo apt install libsuitesparse-dev libhdf5-dev libblas-dev liblapack-dev

# Fedora / RHEL
sudo dnf install suitesparse-devel hdf5-devel blas-devel lapack-devel

# macOS (Homebrew)
brew install suite-sparse hdf5 openblas
```

Upstream source: <https://people.engr.tamu.edu/davis/suitesparse.html>.

The current [Makefile](Makefile) links against a bundled `SuiteSparse/` tree via static archives and local include paths. If you install SuiteSparse system-wide, replace the `LIB` and `IFLAGS` lines with something like:

```make
LIB    = -lumfpack -lamd -lcholmod -lcolamd -lccolamd -lcamd \
         -lsuitesparseconfig -lm -lblas -llapack -fopenmp \
         $(shell pkg-config --libs hdf5-serial 2>/dev/null || pkg-config --libs hdf5)
IFLAGS = $(shell pkg-config --cflags suitesparse 2>/dev/null)
```

### Targets

```bash
make           # release build (-O3, -fopenmp)
make release   # same, after a clean
make debug     # -O0 -g for debugging
make clean
```

`-O3` is essential — it typically yields a 2–3× speed-up over `-O0` on this numerical workload.

## Run

```bash
./main
```

Grid sizes and parameter ranges are configured in [main.c](main.c) and [dataset.c](dataset.c). Defaults:

- `N_S = 100`, `N_v = 81`, `N_t = 100`
- `S ∈ [0, 3K]`, `K = 1`
- Tavella–Randall density `α_S = 0.2`, relative `v`-density `α_v_rel = 0.1`

Constraints on grid sizes (enforced by convention):

- `(N_S − 1)` divisible by 3
- `(N_v − 1)` divisible by 8
- `N_t ≥ N_S / 10`

The built-in parameter sweep covers `5 × 4 × 5 × 4 × 4 × 5 = 8000` parameter sets, each producing four option surfaces (European/American × Call/Put).

## Parallelism

The outer sweep over parameter combinations is parallelised with OpenMP (`schedule(dynamic, 1)`). Each thread performs an independent matrix build, factorisation and back-sweep; HDF5 writes are serialised in a single critical section. Set thread count with `OMP_NUM_THREADS`.

## File layout

| File | Purpose |
|---|---|
| [main.c](main.c) | Entry point, grid configuration, timing |
| [dataset.c](dataset.c) | Parameter sweep, OpenMP loop, HDF5 orchestration |
| [Crank_Nicolson.c](Crank_Nicolson.c) | Crank–Nicolson time loop, all four option types |
| [prob.c](prob.c) | Assembly of the CSR Heston operator `L` |
| [grid.c](grid.c) | Tavella–Randall (`S`) and sinh (`v`) non-uniform grids |
| [grid_index.c](grid_index.c) | 2-D → 1-D index helper |
| [greeks.c](greeks.c) | Non-uniform finite-difference Greeks |
| [umfpack.c](umfpack.c) | UMFPACK factor/solve/free wrappers |
| [hdf5_writer.c](hdf5_writer.c) | Chunked, compressed HDF5 output |
| [norme.c](norme.c) | Vector-norm utilities |
| [time.c](time.c) | CPU / wall-clock timers |
| [Option_Characteristic.h](Option_Characteristic.h) | `CALL`/`PUT`, `EUROPEAN`/`AMERICAN` macros |
| [Makefile](Makefile) | Build configuration |
