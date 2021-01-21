It appears difficult to export the C++ libraries in a way that can easily be
linked outside of the package. Rcpp provides some tools to do this, like the
`Rcpp::interfaces` attribute, but that seems to require function arguments to
be objects that Rcpp knows how to work with. Here we have some custom classes
that are part of the interface.

For now, we'll just construct samplers for specific cases inside the package,
and export those. This is okay for now, but is not particularly great for a
public package.

This is only an issue for linking C++ code. For R code, it should be much
easier to make the interface stuff available outside the package.

