#=
Reference:
Chiba H, Schneider NS, Matsuoka S, Noma A. A Simulation Study on the Activation
of Cardiac CaMKII δ-Isoform and Its Regulation by Phosphatases.
Biophysical Journal. 2008;95(5):2139-2149. doi:10.1529/biophysj.107.118505.
DOI: 10.1529/biophysj.107.118505, PMID: 18502812, PMCID: PMC2517018
=#


using DifferentialEquations, Parameters, LabelledArrays

#=
Default params
TOT_CAM = 6e-3
K1 = 2.5
Km1 = 0.05
K2 = 88.25
Km2 = 0.05
K3 = 12.5
Km3 = 1.25
K4 = 250
Km4 = 1.25
Ca
K_ASSO = 1.4e-4  (α isoform), 0.7e-4 (δ isoform)
K_DISSO = 1.9e-3  (α isoform), 0.95e-3 (δ isoform)
K_DISSO_2 = 1E-3 * K_DISSO
K_DISSO_CA = 1.9e-3  (α isoform), 0.95e-3 (δ isoform)
K_DISSO_CA_2 = 1E-3 * K_DISSO_CA
KM_CAM = 3E-5
KCAT = 9E-4  (α isoform), 54E-4 (δ isoform)
KM_ATP = 19.1E-3
KCAT_PP1 = 1.72E-3
KM_PP1 = 11E-3
=#

_hills(x, k, n) = x^n / (x^n + k^n)
_mm(x, k) = x / (x + k)

# δ Isoform by default (All concentrations in mM and time in ms)
@with_kw_noshow struct CAMKIIParams
    # concentrations
    TOT_CAM = 6.0e-3
    TOT_CAMKII = 0.1e-3
    ATP = 1.0
    PP1 = 1e-6
    # CaM association and idssociation constant
    K_P1 = 2.5
    K_M1 = 0.05
    K_P2 = 88.25
    K_M2 = 0.05
    K_P3 = 12.5
    K_M3 = 1.25
    K_P4 = 250
    K_M4 = 1.25
    # CaMKII association and idssociation constant
    K_ASSO = 2.1
    K_DISSO = 0.7E-4
    K_DISSO_2 = 1E-3 * K_DISSO
    K_DISSO_CA = 0.95e-3
    K_DISSO_CA_2 = 1E-3 * K_DISSO_CA
    KM_CAM = 3E-5
    KCAT = 5.4E-3
    KM_ATP = 19.1E-3

    ATP_FAC = _hills(ATP, KM_ATP, 1)
    KCAT_PP1 = 1.72E-3
    KM_PP1 = 11E-3
    PP1 = 1e-6
    # Calcium Pulse parameters
    INTERVAL = 0
    KP = 12e-6
    KM_PUMP = 0.1e-3
    QPUMP_REST = 40e-6
    QM = 60E-6
    TP = 25.0
    # CaMKII target paramaters
    KCAT_TARGET = 1E-3
    MAX_SS_INHIB = 0.1
    KD_TARGET = KCAT_TARGET * MAX_SS_INHIB / (1 - MAX_SS_INHIB)
end

u0 = @LVector zeros(8) (:ca, :camca1, :camca2, :camca3, :camca4, :camk2_camca4, :camk2_p, :camk2_p_camca4)

function pulse_ca!(du, u, p, t)
    @unpack KP, KM_PUMP, QM, TP, QPUMP_REST, INTERVAL = p
    if INTERVAL > 0
        localTime = rem(t, INTERVAL)
        ca = u.ca
        q_pump = KP * _hills(ca, KM_PUMP, 2)
        q_rel = QM * (localTime/TP * exp(1 - localTime/TP))^4 + QPUMP_REST
        du.ca = q_rel - q_pump
    else
        du.ca = 0.0
    end
    return du
end

function cam_system!(du, u, p, t)
    @unpack TOT_CAM, KP1, KM1, KP2, KM2, KP3, KM3, KP4, KM4 = p
    cam = TOT_CAM - (u.camca1 + u.camca2 + u.camca3 + u.camca4)
    v1 = u.ca * cam * KP1 - u.camca1 * KM1
    v2 = u.ca * u.camca1 * KP2 - u.camca2 * KM2
    v1 = u.ca * u.camca1 * KP3 - u.camca3 * KM3
    v1 = u.ca * u.camca1 * KP4 - u.camca4 * KM4
end

function camkk2_system!(du, u, p, t)

end

function cam!(du, u, p, t)
    # Rate of CaM-Ca association
    _v_cam(ca, cam1, cam2, k⁺, k⁻) = k⁺ * ca * cam1 - k⁻ * cam2
    @unpack TOT_CAM, K1, Km1, K2, Km2, K3, Km3, K4, Km4 = p
    @views cam_ca1 = u[:, 1]
    @views cam_ca2 = u[:, 2]
    @views cam_ca3 = u[:, 3]
    @views cam_ca4 = u[:, 4]
    @views Ca = u[:, 8]
    cam = TOT_CAM - (cam_ca1 + cam_ca2 + cam_ca3 + cam_ca4)
    v_ca1 = _v_cam(Ca, cam, cam_ca1, K1, Km1)
    v_ca2 = _v_cam(Ca, cam_ca1, cam_ca2, K2, Km2)
    v_ca3 = _v_cam(Ca, cam_ca2, cam_ca3, K3, Km3)
    v_ca4 = _v_cam(Ca, cam_ca3, cam_ca4, K4, Km4)
    @. @views du[:, 1] = _v_cam(Ca, TOT_CAM - (cam_ca1 + cam_ca2 + cam_ca3 + cam_ca4), cam_ca1, K1, Km1) - _v_cam(Ca, cam_ca1, cam_ca2, K2, Km2)
    @. @views du[:, 2] = _v_cam(Ca, cam_ca1, cam_ca2, K2, Km2) - _v_cam(Ca, cam_ca2, cam_ca3, K3, Km3)

    du[1] = d_cam_ca1 = v_ca1 - v_ca2
    du[2] = d_cam_ca2 = v_ca2 - v_ca3
    du[3] = d_cam_ca3 = v_ca3 - v_ca4
    du[4] = d_cam_ca4 = v_ca4
    return du
end


function ca_pulse!(du, u, p, t)
    @unpack KP, KM_PUMP, QM, TP, QPUMP_REST, INTERVAL = p
    Ca = u[8]
    q_pump = KP /(1 + (KM_PUMP / Ca)^2)
    q_rel = QM * (t/TP * exp(1 - t/TP))^4 + QPUMP_REST
    du[8] = dCa = q_rel - q_pump
    return du
end


function target!(du, u, p, t)
    @unpack KCAT_TARGET, KD_TARGET, TOT_CAMKII = p
    camkii = TOT_CAMKII - u[5] - u[6] - u[7]
    target = u[9]
    du[9] = -KCAT_TARGET * target * (TOT_CAMKII - camkii) + KD_TARGET * (1 - target)
    return du
end


function camkii!(du, u, p, t)
    @unpack (ATP_FAC, PP1, TOT_CAMKII,
            K_ASSO, K_DISSO, K_DISSO_2,K_DISSO_CA, K_DISSO_CA_2, KM_CAM, KCAT,
            KM_ATP, KCAT_PP1, KM_PP1) = p

    cam_ca4 = u[4]
    camkii_camca4 = u[5]
    camkii_p = u[6]
    camkii_p_camca4 = u[7]
    Ca = u[8]
    camkii = TOT_CAMKII - camkii_camca4 - camkii_p - camkii_p_camca4
    # Association of CaMKII and CaMCa4
    a1 = K_ASSO * camkii * cam_ca4

    # Dissociation of CaMKII and CaMCa4
    h = _hills(Ca, KM_CAM, 3)
    a2 = (K_DISSO * (1 - h) + K_DISSO_CA * h) * camkii_camca4

    # Autophosphorylation of CaMKII subunits
    p = 1 - (camkii / TOT_CAMKII)^2
    b1 = KCAT * p * ATP_FAC * camkii_camca4

    # Dephosphorylation by phosphorylase (PP1)
    b2 = KCAT_PP1 * PP1 * _mm(camkii_p_camca4, KM_PP1)
    d1 = KCAT_PP1 * PP1 *_mm(camkii_p, KM_PP1)

    # Dissociation of CaMKII-p and CaMCa4 (1000x slower)
    c1 = (K_DISSO_2 * (1 - h) + K_DISSO_CA_2 * h) * camkii_p_camca4

    # Association of CaMKII-p and CaMCa4
    c2 = K_ASSO * camkii_p * cam_ca4

    du[5] = d_camkii_camca4 = a1 - a2 - b1 + b2
    du[6] = d_camkii_p = c1 - c2 - d1
    du[7] = d_camkii_p_camca4 = b1 - b2 - c1 + c2
    return du
end

param
u₀ = zeros(length(labels))
u₀[labels[:Ca]] = 0.5
u₀[labels[:target]] = 1.0
tspan = (0.0, 1000.0)
prob = ODEProblem(camkii, u₀, tspan, CAMKIIParams())

sol = solve(prob, Tsit5())


using Plots
gr()

sol[3,:]

savefig("plot.png")

scatter(rand(100000))
