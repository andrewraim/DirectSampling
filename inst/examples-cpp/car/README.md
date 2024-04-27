Code to reproduce conditional autoregression (CAR) mixed model example in
<https://link.springer.com/article/10.1007/s11222-022-10188-x>. The paper uses
a chain of 100,000 draws, burn = 20,000, and thin = 10.

- pricedata.R: Prepare data and run Gibbs sampler.
- gibbs.R: Gibbs sampler to fit mixed model with direct sampling step for rho.
- sampler.cpp: C++ version of direct sampling step.

This example uses older versions of CARBayes (v1.6) and CARBayesdata (v2.2),
which can be obtained using the following.

```{r}
devtools::install_version("CARBayes", version = "1.6")
devtools::install_version("CARBayesdata", version = "2.2")
```

