using DustyHII
using Plots

# Define the fixed model parameters
γ = 10.0
β = 3.0
T4=0.94
hnui = 18.0 # eV


# Next, run various u0 values and tabulate n_rms
# We use different Q0 values but it also depends on the inner
# density. Since I want 

N = 100
Q0 = 10.0.^LinRange(48, 51, N)
ninner = 1e3u"1/cm^3"
u0s = [ustrip(DustyHII.uval(ninner, Q0=q*u"1/s", hnui=hnui)) for q in Q0];

Qnrms = zeros(N)
nrms = zeros(N)
for i=1:N
    m = DustyHII.draine2011model(γ, β; u0=u0s[i], ystart=10.0, yend=1e3)
    nrms[i] = ustrip(DustyHII.n_rms(m; Q0=Q0[i], T=1e4*T4*u"K"))
    Qnrms[i] = ustrip(Q0[i]*nrms[i])
end



    
