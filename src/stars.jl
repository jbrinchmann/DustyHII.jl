# Functions to calculate basic quantities for stellar spectra.
# For these models this is only the Q_0, <hν>_i, L_i, L_n, β quantities. 
#
using Unitful, UnitfulRecipes, UnitfulAstro
using Dierckx


"Interpolate a spectrum. The routine strips units so the user needs to take care with that"
function interpolate_spectrum(λ, flux)
    itp = Spline1D(ustrip(λ), ustrip(flux));
    return itp
end
"""
  res = specQuantities(λ, itp)

Calculate various quantities for an interpolated spectrum.

    λ: Wavelength array, assumed to be in Angstrom
  itp: Interpolated spectrum returned from the interpolate_spectrum function above. 

Returns:
  res: Named tuple with the following keys: 
     Ln: Luminosity in photons with E<13.6 eV  [erg/s]
     Li: Luminosity in photons with E>13.6 eV  [erg/s]
     Q0: Number of photons with E>13.6 eV      [1/s]
   hnui: Mean energy in photons with E>13.6eV  [eV]
      β: Ln/Li
"""
function specQuantities(λ, itp)
        
    max_l = maximum(λ)
    min_l = minimum(λ)
    lambda_i = 911.267

    norm, err = quadgk(itp, min_l, max_l)
    Ln, err = quadgk(itp, lambda_i, max_l)
    Li, err = quadgk(itp, min_l, lambda_i)
    
    β = Ln/Li
    # The units are expected to be erg/s/angstrom for the flux so
    # this gives the units for the luminosities
    Ln = Ln*u"erg/s"
    Li = Li*u"erg/s"
    
    m1_int(x) = itp(x)/x
    mm1_int(x) = x*itp(x)
    norm, err = quadgk(itp, min_l, lambda_i)
    pre = Unitful.h*Unitful.c
    hnui, err = quadgk(m1_int, min_l, lambda_i)
    hnui = pre*hnui/norm
    Q0, err = quadgk(mm1_int, min_l, lambda_i)
    Q0 = Q0*u"erg/s"/pre
    
    # Finally adjust units because the integrals are done unitless so
    # we need to multiply my Angstrom (except for Ln/Li which are done above)
    hnui = uconvert(u"eV", hnui*1u"1/angstrom")
    Q0 = upreferred(Q0*1u"angstrom")
    
    return (Ln=Ln, Li=Li, hnui=hnui, Q0=Q0, β=β)
end

#
# Black-body functions below
#

function Bnu(ν, T)
    z = Unitful.h.*ν./(Unitful.k.*T)
    2 .*Unitful.h.*ν^3/(Unitful.c^2 .*(exp(z).-1))
end

function Blambda(λ, T)
    z = Unitful.h.*Unitful.c ./(Unitful.k.*T*λ)
    2 .*Unitful.h.*Unitful.c^2/(λ^5 .*(exp(z).-1))
end

"""
   PlanckAstro(λ, T; convert_to_L=true)

The Planck function in astronomical units and following planck.pro in IDL

Arguments:
    λ: Wavelength (must have units)
    T: Temperature (must have units)

Keywords:
    convert_to_L: Whether the result should be converted to a luminosity.
                If set to true, the radius-T_eff relationship from 
                Sternberg, Hoffmann & Pauldrach (2003) is used

"""
function PlanckAstro(λ, T; convert_to_L=true)
    Bl = Blambda(λ, T)
    
    # Now convert to erg/s/cm^3/A
    Bl = π*uconvert(u"erg/s/cm^2/angstrom", Bl)
    
    # And then to luminosities if requested
    if (convert_to_L)
        Bl=uconvert(u"erg/s/angstrom", Bl*4*π*radius_for_T(T)^2)
    end
    return ustrip(Bl)
end

#
# Finally we need a T_eff-radius relation to calculate luminosities.
#
"The radius relation for dwarfs from Sternberg, Hoffmann & Pauldrach (2003) "
radius_for_T(T) = (0.2*ustrip(T)/1e3+2)*u"Rsun"
