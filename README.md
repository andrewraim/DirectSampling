This package contains code to support the paper "Direct Sampling with a Step
Function" (<https://doi.org/10.1007/s11222-022-10188-x>). Several examples
in the paper demonstrate use of the sampler.

To install the package from Github, first ensure that you have a recent
installation of R and the packages `datastructures`  and `R6`. At the time of
this writing, `datastructures` is not actively maintained on CRAN but is
available from the CRAN archive.

The `DirectSampling` package can be installed and built using `devtools`.
The tag `v0.2.0` represents the code corresponding to the paper. 

```
# Install from the master branch
R> devtools::install_github("andrewraim/DirectSampling")

# Install the code tagged v0.2.0
R> devtools::install_github("andrewraim/DirectSampling", ref = "v0.2.0")
```
