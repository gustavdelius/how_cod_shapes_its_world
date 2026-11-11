# Background model: scripts and results

This folder holds the computations behind the model in which cod and sprat
are embedded in a power-law background: a fixed background mortality, and
background prey that is depleted by what cod and sprat eat of it. The write-up is
[`background_model.qmd`](../background_model.qmd) in the repository root,
rendered to [`docs/background_model.html`](../docs/background_model.html).

| File | Contents |
|:--|:--|
| [`run_scan.R`](run_scan.R) | The scan over the strengths of the two predation links, with and without cod larvae feeding on plankton: invasion fitness of cod, its decomposition, and the stability of the coexistence state |
| `scan.rds` | Its results |
| [`run_cases.R`](run_cases.R) | Four points of the regime map: size spectra, cod started below and above the coexistence state, cod introductions, sensitivity, and density-dependent cod reproduction. Two slow points of the scan over 3000 years |
| `cases.rds` | Its results |
| [`run_fishing.R`](run_fishing.R) | The equilibrium branches under cod fishing, and the effect of a sprat fishery |
| `fishing.rds` | Its results |

The model itself is in [`R/background_model.R`](../R/background_model.R). It
needs the [mizerEcopath](https://github.com/gustavdelius/mizerEcopath)
package, which brings in mizerFromEcopath.

## Reproducing

Run from the repository root, with BLAS limited to one thread per worker:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 Rscript background/run_scan.R
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 Rscript background/run_cases.R
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 Rscript background/run_fishing.R
quarto render background_model.qmd
```

The scan takes about an hour and a half on 4 cores, and the other two scripts
about 35 and 20 minutes. Set
the environment variable `MC_CORES` to use fewer workers. The document reads
the three results files and renders in under a minute.
