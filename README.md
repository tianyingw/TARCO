# TARCO

R implementation of tree-aggregated regression with measurement error correction.

Requires R and CVXR with the ECOS and SCS solvers:

```r
install.packages(c("CVXR", "ECOSolveR", "scs"))
```

From this directory, run:

```sh
Rscript --vanilla example.R
```

The example generates 100 observations and an eight-leaf tree, fits TARCO with
known measurement error covariance and five-fold cross-validation, and prints
the coefficients and groups. Results are also saved to `output/coefficients.csv`
and `output/groups.csv`. Each training fold is centered before fitting.
Groups are maximal subtrees with coefficient range at most 0.01.
