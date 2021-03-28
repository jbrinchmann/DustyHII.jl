# Dusty HII models

This package implements the dusty HII models with radiation pressure
presented in
[Draine (2011, D11 hereafter)](https://ui.adsabs.harvard.edu/abs/2011ApJ...732..100D/abstract).
However I have extended this to also include self-gravity.

*NOTE*: This is currently under development. It works but is not fully tested.

# Basic usage

The usage is fairly straightforward. You need to define the model
parameter which is a `Drain2011Parameters` type, and the model results
are encapsulated in the  `Draine2011` type. 


The model is specified in a dimensionless form in terms of three
variables:

- y is a dimensionless radius
- u(y) is a dimensionless inverse density
- ϕ(y) is a dimensionless variable tracing the ionizing photons as a
  function of radius.


# Governing equations

These are defined through a set of differential equations covering the
physics of the system. The derivation of these equations is given
briefly in D11. A more extended version and the expansion to
self-gravity can be found in the `docs` folder.

