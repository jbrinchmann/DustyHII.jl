# Dusty HII models

This package implements the dusty HII models with radiation pressure
presented in
[Draine (2011)](https://ui.adsabs.harvard.edu/abs/2011ApJ...732..100D/abstract).
However I have extended this to also include self-gravity.

# Basic usage

The usage is fairly straightforward. You need to define the model
parameter which is a `Drain2011Parameters` type, and the model results
are encapsulated in the  `Draine2011` type. 



# Governing equations

These are defined through a set of differential equations covering the physics of the system. I want to extend these with the Poisson equation and including self-gravity, but this notebook is for the vanilla Drain (2011) model.

In this model the defining equations are as follows:

$$
2k_BT \frac{dn_H}{dr} = n_H\sigma_d \frac{L_n e^{-\tau} + L_i\phi}{4\pi r^2 c} + \alpha_B n_H^2 \frac{\langle h\nu\rangle_i}{c}, $$

$$\frac{d\phi}{dr} = -\frac{\alpha_B n_H^2 4 \pi r^2}{Q_0} - n_H \sigma_d \phi,$$
and
$$ \frac{d\tau}{dr} = n_H \sigma_d.$$

