#
# This file contains functions that deal with physical ingredients
# in the model.
#
using Unitful, UnitfulRecipes, UnitfulAstro
using PhysicalConstants.CODATA2018
G = NewtonianConstantOfGravitation


T4(T) = ustrip(T/(1e4*u"K"))
Q0p49(Q0) = ustrip(Q0/(1e49*u"1/s"))
hnui18(hnui) = ustrip(hnui/(18.0*u"eV"))
σ21(σ) = ustrip(σ/(1e-21*u"1/cm^2"))

"The Case B recombination coefficient from Draine 2011"
alpha_B(T) = 2.56e-13*T4(T)^(-0.83)*u"cm^3/s"



#
# Scaling parameters of the model
#

"The characteristic density used in defining u = n0/n"
n0(;T=1e4*u"K", Q0=1e49*u"1/s", hnui=18*u"eV", σ=nothing) = 4.54e5*(T4(T)^4.66/Q0p49(Q0))*hnui18(hnui)^(-3)*u"1/cm^3"
"The characteristic length scale used to define y = r/λ0"
λ0(;T=1e4*u"K", Q0=1e49*u"1/s", hnui=18*u"eV", σ=nothing) = 2.47e16*(Q0p49(Q0)/T4(T)^2.83)*hnui18(hnui)^2*u"cm"
"The γ parameter used in the equation definition"
gamma(;T=1e4*u"K", Q0=nothing, hnui=18*u"eV", σ=1e-21*u"1/cm^2") = 11.2*T4(T)^1.83*σ21(σ)/hnui18(hnui)

"The characteristic mass"
M0(T; kwargs...) = 4π*λ0(;T=T4(T), kwargs...)^3*n0(;T=T4(T), kwargs...)*ProtonMass/3
"The characteristic gravitational energy"
Eg(T; kwargs...) = G*ProtonMass*M0(T4(T); kwargs...)/λ0(;T=T4(T), kwargs...)
"The characteristic thermal energy"
Et(T) = 2*BoltzmannConstant*T4(T)*1e4u"K"/3


#
# To-from (n, r) <-> (u, y)
#
"Convert density to u"
uval(n; kwargs...) = n0(; kwargs...)/n
"Calculate radius from y"
radius(y; kwargs...) = λ0(; kwargs...)*y
"Calculate y from radius"
yval(r; kwargs...) = r/λ0(; kwargs...)
"Calculate density from the solution"
density(y, model::AbstractDustyHII) = n0(;kwargs...)/ufunc(y, model)


#
# Physical properties of the solution
#

"RMS density given a solution"
function n_rms(m; kwargs...)
    integrand_rms(y) = y^2/ufunc(y, m)^2 
    int, dint = quadgk(integrand_rms, m.ymin, m.ymax)
    return n0(;kwargs...)*sqrt(3.0*int/m.ymax^3)
end

"Mean density given a solution"
function n_mean(m; kwargs...)
    integrand_mean(y) = y^2/ufunc(y, m)
    int, dint = quadgk(integrand_mean, m.ymin, m.ymax)
    return 3.0*n0(;kwargs...)*int/m.ymax^3
end
    
"The Stromgren radius for n(RMS)"
function Rs(sol; Q0=1e49u"1/s", T=1e4u"K", kwargs...)
    nrms_3 = n_rms(sol; Q0=Q0, T=T, kwargs...)/(1e3u"1/cm^3")
    return 2.10e18*Q0p49(Q0)^(1/3.)*T4(T)^0.25/nrms_3^(2/3.)*u"cm"
end

# RMS optical depth
taud0(sol; σ=1e-21*u"1/cm^2", kwargs...) = n_rms(sol; kwargs...)*Rs(sol; kwargs...)*σ

"Pressure on the edge of the model"
p_edge(model::AbstractDustyHII; T=1e4u"K", kwargs...) = 2.0*n0(;kwargs...)*Unitful.k*T/ufunc(ymax(model), model)

"The radius of the ionised zone"
Rion(sol;kwargs...) = ymax(sol)*λ0(;kwargs...)

"Fraction of absorbed E>13.6 eV photons absorbed by H"
fion(sol;kwargs...) = (Rion(sol; kwargs...)/Rs(sol; kwargs...))^2

 
