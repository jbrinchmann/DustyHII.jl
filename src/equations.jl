using Unitful, UnitfulRecipes, UnitfulAstro
using Dierckx

#
# The file containing the various differential equation systems to solve.
#

"The differential equation system derived by Draine 2011"
function d2011system!(dv, v, p, y)
    u, ϕ, τ = v
    γ, β = p
    
    dv[1] = du = -1 - γ*(β*exp(-τ) + ϕ)*u/y^2
    dv[2] = dϕ = -(y/u)^2 - γ * ϕ/u
    dv[3] = dτ = γ/u
end


"""
    draine2011model(γ, β; u0=454, ystart=10, yend=100)

Solve the Draine (2011) model with model parameters γ and β. 
This formulates the ODEProblem and solves it

"""
function draine2011model(γ, β; u0=454.0, ystart=10.0, yend=100.0)
    
    starting_point = [u0, 1.0, 0.0]
    yspan = [ystart, yend]

    p = [γ, β]

    # The system is only defined for ϕ>=0 so we interrupt on ϕ=u[2]=0
    condition(u, y, integrator) = u[2]
    affect!(integrator) = terminate!(integrator)
    cb = ContinuousCallback(condition, affect!);
    prob = ODEProblem(d2011system!, starting_point, yspan, p)
    
    sol = solve(prob, Tsit5(), callback=cb);

    result = Draine2011(sol, Draine2011Parameters(γ, β), u0,
                        starting_point[2], starting_point[3],
                        ystart, ymax(sol))
    
    return result
    
end


draine2011model(d::NamedTuple; kwargs...) = draine2011model(d.γ, d.β; u0=d.u0, kwargs...)


function d2011wgravity!(dv, v, h, p, y)
    u, ϕ, τ = v
    γ, β, η, ymin = p # Ymin is indirectly a parameter of the model
    
    rint(x) = x^2/h(p, x)[1]
    Ry, err = quadgk(rint, ymin, y)
            
    dv[1] = du = u^2*η*Ry/y^2 -1 - γ*(β*exp(-τ) + ϕ)*u/y^2
    dv[2] = dϕ = -(y/u)^2 - γ * ϕ/u
    dv[3] = dτ = γ/u
end



function d2011wgravitymodel(γ, β, η; u0=454, ystart=10, yend=100)
    
    starting_point = [u0, 1, 0]
    yspan = [ystart, yend]

    p = [γ, β, η, ystart]

    # The system is only defined for ϕ>=0 so we interrupt on ϕ=u[2]=0
    condition(u, y, integrator) = u[2]
    affect!(integrator) = terminate!(integrator)
    cb = ContinuousCallback(condition, affect!);
    
    # This specifies the original starting point. I think it is probably
    # necessary to rethink this  but haven't done this yet
    history(p, t) = [u0, 1, 0]
    prob = DDEProblem(d2011wgravity!, starting_point, history, yspan, p)
    
    # The choice for stepper is following https://diffeq.sciml.ai/stable/solvers/dde_solve/
    # see also https://github.com/SciML/DifferentialEquations.jl/issues/383
    sol = solve(prob, MethodOfSteps(RK4()), callback=cb);

    return sol
    
end



#
# Handling of solutions themselves
#

ymax(solution) = maximum(solution.t)
ymin(solution) = minimum(solution.t)


# Return the functions of the differential equation
tau(y, model::AbstractDustyHII) = model.solution(y; idxs=3)
phi(y, model::AbstractDustyHII) = model.solution(y; idxs=2)
ufunc(y, model::AbstractDustyHII) = model.solution(y; idxs=1)



#
# To-from initial conditions.
#
# In this case we have two options:
# 
#    1. The user might give physical parameters of a model and want
#       the parameters to use for initial conditions and model parameters γ and β.
#
#    2. The user might have a solution for a given γ and β and wants physical
#       parameters that are consistent with this. This is the more complex case
#       because of the various extra parameters that might have to be passed.
#

"""
    physical_to_model(; Teff=45000u"K", Tgas=1e4u"K", ninner=1e3u"1/cm^3")

Calculates β, Q0, hν_i for a star with the given Teff and calculates γ and
η for the other properties. 

"""
function physical_to_model(;Teff=45000u"K", Tgas=1e4u"K", ninner=1e3u"1/cm^3",
                           σ=1e-21*u"1/cm^2")

    # Calculate first a Black-body spectrum - I fix wavelength ranges to
    # be reasonable for a hot star.
    lrange = collect(LinRange(1., 6e4, 10000)) 
    flux = [PlanckAstro(l*u"angstrom", Teff) for l in lrange];
    itp = interpolate_spectrum(lrange, flux)
    resStar = specQuantities(lrange, itp)

    # Then calculate γ
    γ = gamma(;T=Tgas, Q0=resStar.Q0, hnui=resStar.hnui, σ=σ)

    #  η
    η = upreferred(Eg(Tgas; Q0=resStar.Q0, hnui=resStar.hnui, σ=σ)/Et(Tgas))

    # u0
    u0 = ustrip(uval(ninner; T=Tgas, Q0=resStar.Q0, hnui=resStar.hnui))

    return (Ln=resStar.Ln, Li=resStar.Li, hnui=resStar.hnui,
            Q0=resStar.Q0, β=resStar.β, γ=γ, η=η, u0=u0, Teff=Teff,
            Tgas=Tgas, ninner=ninner, σ=σ)

end


"""
    model_to_physical(m)

Takes a model which has been calculated in unitless form and converts quantities
to physical units

NOT YET FUNCTIONAL
"""
function model_to_physical(m, T, Q0, hnui, σ)
    return Draine2011Physical(m, T, Q0, hnui, σ)

end
