#=
Reference:
Chiba H, Schneider NS, Matsuoka S, Noma A. A Simulation Study on the Activation
of Cardiac CaMKII δ-Isoform and Its Regulation by Phosphatases.
Biophysical Journal. 2008;95(5):2139-2149. doi:10.1529/biophysj.107.118505.
DOI: 10.1529/biophysj.107.118505, PMID: 18502812, PMCID: PMC2517018
=#

using DifferentialEquations, Parameters, LabelledArrays, Plots

# Michaelis-Menton and Hill function (from my CMC model)
_mm(x, k=one(x)) = x / (x + k)
_mm_reci(x, k=one(x)) = (x + k) / x
_hill(x, k, n) = _mm(x^n, k^n)
_hill_reci(x, k, n) = _mm_reci(x^n, k^n)
_comp(x) = one(x) - x

@with_kw struct CaMParams
    ΣCAM = 6.0e-3
    k1_p = 2.5
    k1_m = 0.05
    k2_p = 88.25
    k2_m = 0.05
    k3_p = 12.5
    k3_m = 1.25
    k4_p = 250
    k4_m = 1.25
end

# Calmodulin system
function cam_sys(camca1, camca2, camca3, camca4, ca, pCaM::CaMParams)
    cam = pCaM.ΣCAM - (camca1 + camca2 + camca3 + camca4)
    @unpack k1_p, k1_m, k2_p, k2_m, k3_p, k3_m, k4_p, k4_m = pCaM
    v1 = k1_p * cam * ca - k1_m * camca1
    v2 = k2_p * camca1 * ca - k2_m * camca2
    v3 = k3_p * camca2 * ca - k3_m * camca3
    v4 = k4_p * camca3 * ca - k4_m * camca4
    dcamca1 = v1 - v2
    dcamca2 = v2 - v3
    dcamca3 = v3 - v4
    dcamca4 = v4
    return (dcamca1, dcamca2, dcamca3, dcamca4)
end

function cam_sys!(du, u, pCaM::CaMParams)
    du.camca1, du.camca2, du.camca3, du.camca4 = cam_sys(u.camca1, u.camca2, u.camca3, u.camca4, u.ca, pCaM)
end

# δ Isoform by default (All concentrations in mM and time in ms)
@with_kw struct CAMKIIParams
    ΣCAMKII = 0.1e-3  # CaMKII concentrations
    PP1 = 0.1e-3     # Phosphatase 1 concentrations
    # rate constantS
    KA = 2.1
    KD = 0.7E-4
    KD_CA = 0.95e-3
    KD_2 = 1E-3 * KD
    KD_CA_2 = 1E-3 * KD_CA
    KM_CAM = 3E-5
    KCAT = 5.4E-3  # 310K
    KM_ATP = 19.1E-3
    KCAT_PP1 = 1.72E-3
    KM_PP1 = 11E-3
    VMAX_PP1 = PP1 * KCAT_PP1
end

pCAMKIIδ = CAMKIIParams()
pCAMKIIα = CAMKIIParams(pCAMKIIδ;
                        KD = 2 * pCAMKIIδ.KD,
                        KD_CA = 2 * pCAMKIIδ.KD_CA,
                        KD_2 = 2 * pCAMKIIδ.KD_2,
                        KD_CA_2 = 2 * pCAMKIIδ.KD_CA_2,
                        KCAT = pCAMKIIδ.KCAT / 6)

_camkii(camkii_camca4, camkii_p_camca4, camkii_p, pCaMKII::CAMKIIParams) = pCaMKII.ΣCAMKII - camkii_camca4 - camkii_p_camca4 - camkii_p
_camkii(u, pCaMKII::CAMKIIParams) = _camkii(u.camkii_camca4, u.camkii_p_camca4, u.camkii_p, pCaMKII)
inactive_camkii(camkii, p::CAMKIIParams) = camkii / p.ΣCAMKII
active_camkii(camkii, p::CAMKIIParams) = _comp(inactive_camkii(camkii, p))
inactive_camkii(camkii_camca4, camkii_p_camca4, camkii_p, p::CAMKIIParams) = inactive_camkii(_camkii(camkii_camca4, camkii_p_camca4, camkii_p, p), p)
active_camkii(camkii_camca4, camkii_p_camca4, camkii_p, p::CAMKIIParams) = (camkii_camca4 + camkii_p_camca4 + camkii_p) / p.ΣCAMKII

function camkii_sys(camkii_camca4, camkii_p_camca4, camkii_p, camca4, ca, atp, pCaMKII::CAMKIIParams)
    camkii = _camkii(camkii_camca4, camkii_p_camca4, camkii_p, pCaMKII)
    @unpack KA, KD, KD_CA, KM_CAM, KD_2, KD_CA_2 = pCaMKII
    a1 = KA * camkii * camca4
    c2 = KA * camkii_p * camca4
    ϕcam = _hill(ca, KM_CAM, 3)
    a2 = (KD * ϕcam + KD_CA * _comp(ϕcam)) * camkii_camca4
    c1 = (KD_2 * ϕcam + KD_CA_2 * _comp(ϕcam)) * camkii_p_camca4

    @unpack KCAT, KM_ATP, VMAX_PP1, KM_PP1 = p
    P = _comp(inactive_camkii(camkii, p)^2)
    b1 = KCAT * P * _mm(atp, KM_ATP) * camkii_camca4
    b2 = VMAX_PP1 * _mm(camkii_p_camca4, KM_PP1)
    d1 = VMAX_PP1 * _mm(camkii_p, KM_PP1)
    v1 = a1 - a2
    v2 = b1 - b2
    v3 = c1 - c2
    v4 = d1
    dcamkii_p = v3 - v4
    dcamkii_camca4 = v1 - v2
    dcamkii_p_camca4 = v2 - v3
    return  dcamkii_camca4, dcamkii_p_camca4, dcamkii_p
end

# Parameters for Calcium transient
@with_kw struct CaTransientParams
    INTERVAL = 1000  # BCL = 1s by defaul
    KP = 12e-6
    KM_PUMP = 0.1e-3
    QPUMP_REST = 40e-6
    QM = 60E-6
    TP = 25.0
end

function ca_transient(ca, t, p::CaTransientParams)
    @unpack KP, KM_PUMP, QM, TP, QPUMP_REST, INTERVAL = p
    localTime = rem(t, INTERVAL)
    tt = localTime/TP
    q_pump = KP * _hill(ca, KM_PUMP, 2)
    q_rel = QM * tt^4 * exp(4*(1 - tt)) + QPUMP_REST
    dca = q_rel - q_pump
    return dca
end

caTProb = ODEProblem((u,p,t)->ca_transient(u, t, p), 1e-3, (0.0, 10000.0),  CaTransientParams())
sol = solve(caTProb, tstops=0.0:1000.0:10000.0)

plot(sol)

plot(t->ca_transient(1e-3, t, CaTransientParams()), 0.0, 1e4)

# Phosphorylatio of targets (mtCK, ion channels) by active camkii
@with_kw struct TargetParams
    KmCaMK = 0.15    # From OR dmodel
end

phosphorylation_portion(camkii, pCaMKII::CAMKIIParams, pTarget::TargetParams) = _mm(active_camkii(camkii, pCaMKII), pTarget.KmCaMK)

# The model parameters
@with_kw struct ModelParams
    pCaM = CaMParams()
    pCAMK = CAMKIIParams()
    pCaT = CaTransientParams()
    pTarget = TargetParams()
    ATP = 0.0
end

p_ref = ModelParams()

# The RHS of the model
function rhs!(du, u, p::ModelParams, t)
    @unpack ATP, pCAMK, pCaM, pTarget, pCaT = p
    du.camkii_camca4, du.camkii_p_camca4, du.camkii_p = camkii_sys(u.camkii_camca4, u.camkii_p_camca4, u.camkii_p, u.camca4, u.ca, ATP, pCAMK)
    du.camca1, du.camca2, du.camca3, du.camca4 = cam_sys(u.camca1, u.camca2, u.camca3, u.camca4, u.ca, pCaM)
    du.ca = ca_transient(u.ca, t, pCaT)
end

# Reference Initial conditions
u0_ref = LVector(camca1 = 0.0, camca2 = 0.0, camca3 = 0.0, camca4 = 0.0, camkii_p = 0.0, camkii_camca4 = 0.0, camkii_p_camca4 = 0.0, ca = 0.1e-3)

u0 = copy(u0_ref)
u0.ca = 0.5

p = ModelParams(pCaM=CaMParams(ΣCAM = 1e-3))
prob = SteadyStateProblem(rhs!, u0, p)
sol = solve(prob, DynamicSS(Tsit5()))
uend = sol.u[end]
active_camkii(uend.camkii_camca4, uend.camkii_p_camca4, uend.camkii_p, p.pCAMK)

# Generate figure 2 A in the paper
function fig2A()
    function get_camkii_act(cam_tot)
        u0 = copy(u0_ref)
        u0.ca = 0.5
        p = ModelParams(pCaM=CaMParams(ΣCAM = cam_tot), pCAMK = CAMKIIParams(pCAMKIIδ))
        prob = SteadyStateProblem(rhs!, u0, p)
        sol = solve(prob, SSRootfind())
        # uend = sol.u[end]
        uend = sol.u
        active_camkii(uend.camkii_camca4, uend.camkii_p_camca4, uend.camkii_p, p.pCAMK)
    end
    cam_tots = 10 .^ (-6:0.01:-2)
    camkii_acts = get_camkii_act.(cam_tots)
    plot(cam_tots, camkii_acts, xscale = :log10, xlabel = "Total CaM", ylabel="Active CaMKII")
end

fig2A()
