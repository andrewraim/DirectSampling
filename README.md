This package contains code to support the manuscript "Direct Sampling with a
Step Function". Replication materials for the manuscript contain several
examples which demonstrate use of the package.

To install the package from Github, first ensure that you have a recent
installation of R and the packages datastructures and R6. Now the DirectSampling
package can be installed and built using devtools.

```
R> devtools::install_github("andrewraim/DirectSampling")
```
That will install whatever is currently in `HEAD`. The tag `v0.2.0` represents
the code which corresponds to the latest version of the manuscript. To install
this instead, use the following.

```
R> devtools::install_github("andrewraim/DirectSampling", ref = "v0.2.0")
```

