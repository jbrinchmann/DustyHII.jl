### Reproduce Figure 2c from Draine (2011)

using DustyHII
using Plots

# Define the fixed model parameters for that panel
γ = 10.0
β = 3.0
T4=0.94
hnui = 18.0 # eV

# The model grid used for the lines.
grid = (
    a=(τ=0.21, Qn=1),
    b=(τ=0.45, Qn=1e1),
    c=(τ=0.96, Qn=1e2),
    d=(τ=2.07, Qn=1e3),
    e=(τ=4.46, Qn=1e4),
    f=(τ=9.60, Qn=1e5),
    e=(τ=20.7, Qn=1e6),
    f=(τ=44.6, Qn=1e7))


keys = grid.keys()

# Go through and calculate the models
# 
# Since the solver expects u0 as an initial value and the
# model grid gives τ and Q_0 n_rms, we here need to convert the Qn 
# value to u0. For this we need the n_rms value which however is a
# value calculated by the model. The way to tackle this is to create
# a table of n_rms vs u0 values.
#



for k in keys

    
m = DustyHII.draine2011model
