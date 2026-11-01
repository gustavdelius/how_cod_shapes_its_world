# Combined model: scripts and results

This folder holds the computations behind the prototype that combines the
emergent Allee effect of Van Leeuwen et al. (2008) with
cultivation–depensation in one mizer model. The write-up is
[`combined_model.qmd`](../combined_model.qmd) in the repository root,
rendered to [`docs/combined_model.html`](../docs/combined_model.html) and
published at
<https://gustavdelius.github.io/how_cod_shapes_its_world/combined_model.html>.

| File | Contents |
|:--|:--|
| [`run_prototype.R`](run_prototype.R) | The scan over the diet-shift size `w_shift` and the larval predation θ, and the cod introductions |
| [`results.rds`](results.rds) | Its results |
| [`fishing_branches.R`](fishing_branches.R) | The equilibrium branches under cod fishing, with and without a sprat fishery |
| [`fishing_branches.rds`](fishing_branches.rds) | Its results |

The model itself is in [`R/combined_model.R`](../R/combined_model.R).

## Reproducing

Run from the repository root, with BLAS limited to one thread per worker:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 Rscript combined/run_prototype.R
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 Rscript combined/fishing_branches.R
quarto render combined_model.qmd
```

The scan takes about an hour on 8 cores, and the fishing branches about 15
minutes on 4 cores. The document reads the two results files and renders in
seconds.
