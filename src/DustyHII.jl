#
# The main file focuses on defining types and the dimensionless form
# of the equations and parameters.
#
# The functions to convert the model results to physical units are all
# given in the physical.jl file, while utility functions to calculate
# stellar parameters for a spectrum (in particular black body ones) are
# in stars.jl
#


module DustyHII

using DifferentialEquations
using QuadGK
using NumericalIntegration
using Unitful, UnitfulRecipes, UnitfulAstro
using Dierckx

abstract type AbstractDustyHII end
abstract type AbstractModelParameters end
abstract type AbstractPhysicalModel end

"A struct that contains the parameters (γ, β) of the Draine 2011 model."
struct Draine2011Parameters <: AbstractModelParameters
γ::Real
β::Real
end

"""
Draine2011(solution, parameterss, u0, ϕ0, τ0, ymax, ymin)

A type containing the results of solving the differential equations for the Draine 2011 model

The solution is stored in `solution`, which is a `ODESolution`.
The parameters in `parameters`, which is a `Draine2011Parameters`
The initial conditions in `u0`, `ϕ0`, `τ0`
The range of validity in `ymax` and `ymin`.
"""
struct Draine2011 <: AbstractDustyHII
solution::ODESolution
parameters::Draine2011Parameters
u0::Real
ϕ0::Real
τ0::Real
ymin::Real
ymax::Real
end


"""
Draine2011Physical(model, T, Q0, hnui, σ)

A type containing the necessary information to give a physical realisation
of a model. It encapsulates the abstract unitless model solution as well
as the physical characteristics needed to create a model.

TODO: Support self-gravity
TODO (maybe): allow for type-checking for the extra parameters. 
"""
mutable struct Draine2011Physical <: AbstractPhysicalModel
model::Draine2011
T
Q0
hnui
σ
end



# I would like a pretty printing method as well.
#
# This is far from being ready, it is here mostly as a placeholder. 
function Base.show(io::IO, m::Draine2011)
    println(io, "Dusty HII region model")
    ymin=m.ymin ; ymax=m.ymax
    println(io, "    Defined for $ymin < y < $ymax")
end



# The actual equations are defined here
include("equations.jl")
# Physical conversion routines are kept in this file.
include("physical.jl")
# Stellar quantities are calculated here
include("stars.jl")





end # module
