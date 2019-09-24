#=
Reference:
Chiba H, Schneider NS, Matsuoka S, Noma A. A Simulation Study on the Activation
of Cardiac CaMKII δ-Isoform and Its Regulation by Phosphatases.
Biophysical Journal. 2008;95(5):2139-2149. doi:10.1529/biophysj.107.118505.
DOI: 10.1529/biophysj.107.118505, PMID: 18502812, PMCID: PMC2517018
=#

using DifferentialEquations, Parameters, LabelledArrays, Plots, DiffEqBiological
using Latexify

const QPUMP_REST = 12e-6 * (40e-6^2 / (0.1e-3^2 + 40e-6^2))

@reaction_func function qrel(t, period, tp, Qm)
    tt = (t % period) / tp
    Qrel = Qm * (tt * exp(1 - tt))^4 + QPUMP_REST
    return Qrel
end

camkIIRn = @reaction_network begin
    qrel(t, period, tp, Qm), 0 → ca
    sign(ca) * hill(ca, Kp, KmPump, 2), ca ⇒ 0
    (k1p * ca, k1m), CaM ↔ CaMCa1
    (k2p * ca, k2m), CaMCa1 ↔ CaMCa2
    (k3p * ca, k3m), CaMCa2 ↔ CaMCa3
    (k4p * ca, k4m), CaMCa3 ↔ CaMCa4
    kasso * CaMCa4, CaMKII → CaMKII_CaMCa4
    hill(ca, kdisso, kmCaM, 3), CaMKII_CaMCa4 → CaMKII
    hillr(ca, kdissoCa, kmCaM, 3), CaMKII_CaMCa4 → CaMKII
    hill(ca, kdisso2, kmCaM, 3), CaMKIIp_CaMCa4 → CaMKIIp
    hillr(ca, kdissoCa2, kmCaM, 3), CaMKIIp_CaMCa4 → CaMKIIp
    mm(ATP, kcat, KmATP) * (1 - (CaMKII / (CaMKII + CaMKII_CaMCa4 + CaMKIIp_CaMCa4 + CaMKIIp))^2 ), CaMKII_CaMCa4 → CaMKIIp_CaMCa4
    mm(CaMKIIp_CaMCa4, kcatPP1, KmPP1) * PP1, CaMKIIp_CaMCa4 ⇒ CaMKII_CaMCa4
    mm(CaMKIIp, kcatPP1, KmPP1) * PP1, CaMKIIp ⇒ CaMKII
end period tp Qm Kp KmPump k1p k1m k2p k2m k3p k3m k4p k4m kasso kdisso kmCaM kdissoCa kdisso2 kdissoCa2 ATP kcat KmATP kcatPP1 KmPP1 PP1

speciesmap(camkIIRn)

paramsmap(camkIIRn)

odeexprs(camkIIRn)

latexify(odeexprs(camkIIRn))

# Parameters for δ type CaMKII
pδ = LVector(
      ca0 = 40e-6,
      period = 1000,
      tp = 25,
      Qm = 60e-6,
      Kp = 12e-6,
      KmPump = 0.1e-3,
      k1p = 2.5,
      k1m = 0.05,
      k2p = 88.25,
      k2m = 0.05,
      k3p = 12.5,
      k3m = 1.25,
      k4p =  250,
      k4m = 1.25,
      kasso = 2.1,
      kdisso = 0.7E-4,
      kmCaM = 3e-5,
      kdissoCa = 0.95e-3,
      kdisso2 = 1E-3 * 0.7E-4,
      kdissoCa2 = 1E-3 * 0.95e-3,
      ATP = 1.0,
      kcat = 5.4E-3,
      KmATP = 19.1E-3,
      kcatPP1 = 1.72E-3,
      KmPP1 = 11E-3,
      PP1 = 0.1e-3,
      ΣCaM = 6.0e-3,
      ΣCaMKII = 0.1e-3)

function make_α(pδ)
    pα = LVector(pδ)
    pα.kdisso *= 2
    pα.kdissoCa *= 2
    pα.kdisso2 *= 2
    pα.kdissoCa2 *= 2
    pα.kcat /= 6
    return pα
end

active_camkii(u) = (u.camkii_camca4 + u.camkii_p_camca4 + u.camkii_p) / (u.camkii + u.camkii_camca4 + u.camkii_p_camca4 + u.camkii_p)
inactive_camkii(u) = 1 - active_camkii(u)


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
