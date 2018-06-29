import math

# Calmodulin model

def get_v_ca(kf, kb, ca, a ,b):
    """ a + Ca <-> b """
    return -kf * ca * a + kb * b


def get_rates_cam_ca(cam_cas, ca, cam_total, k1, km1, k2, km2, k3, km3, k4, km4):
    """ returns dCaMCai /dt """
    cam = cam_total - sum(cam_cas)
    ca1 = get_v_ca(k1, km1, ca, cam, cam_cas[0])
    ca2 = get_v_ca(k2, km2, ca, cam_cas[0], cam_cas[1])
    ca3 = get_v_ca(k3, km3, ca, cam_cas[1], cam_cas[2])
    ca4 = get_v_ca(k4, km4, ca, cam_cas[2], cam_cas[3])
    
    return [-ca1 + ca2, -ca2 + ca3, -ca3 + ca4, -ca4]

def get_rates_camkii(camkiis, camca4, ca, camkii_total, atp, pp1, k_asso, k_disso, k_disso_ca, k_disso2, k_disso_ca2, km_cam, kcat, km_atp, kcat_pp1, km_pp1):
    """CaMKII model"""
    
    camkii_camca4, camkii_p_camca4, camkii_p = camkiis
    camkii = camkii_total - camkii_camca4 - camkii_p_camca4 - camkii_p
    
    # Association of CaMKII and CaMCa4
    a1 = k_asso * camkii * camca4
    
    # Dissociation of CaMKII and CaMCa4
    km_cam_portion = 1 / (1 + (ca / km_cam)**3)
    a2 = (k_disso * (1 - km_cam_portion) + k_disso_ca * km_cam_portion) * camkii_camca4
    
    # Autophosphorylation of CaMKII subunits
    p = 1 - (camkii / camkii_total)**2
    b1 = kcat * p * (1 / (1 + (km_atp / atp))) * camkii_p_camca4
    
    # Dephosphorylation by phosphorylase (PP1)
    b2 = kcat_pp1 * pp1 * (1 / (1 + (km_pp1 / camkii_p_camca4)))
    d1 = kcat_pp1 * pp1 * (1 / (1 + (km_pp1 / camkii_p)))
    
    # Dissociation of CaMKII-p and CaMCa4 (1000x slower)
    c1 = (k_disso2 * (1 - km_cam_portion) + k_disso_ca2 * km_cam_portion) * camkii_p_camca4
    
    # Association of CaMKII-p and CaMCa4
    c2 = k_asso * camkii_p * camca4
    
    return [a1 - a2 - b1 + b2, b1 - b2 - c1 + c2, c1 - c2 - d1]


def get_rates_phospho(camkiis, camkii_total, target, kcat, kd):
    """Get phosphorylation rates by CaMKII  """
    active_camkii = 1 - camkiis[0] / camkii_total
    return [-kcat * target * active_camkii + kd * (1 - target)]


def full_model(t, y, ca, cam_total, k1, km1, k2, km2, k3, km3, k4, km4, 
               camkii_total, atp, pp1, k_asso, k_disso, k_disso_ca, k_disso2, k_disso_ca2, 
               km_cam, kcat, km_atp, kcat_pp1, km_pp1, kcat_complex1, kd_complex1, kcat_ck_mito, kd_ck_mito):
    cam_cas = y[:4]
    camca4 = cam_cas[3]
    camkiis = y[4:7]
    complex1 = y[7]
    ck_mito = y[8]
    v_cam_ca = get_rates_cam_ca(cam_cas, ca, cam_total, k1, km1, k2, km2, k3, km3, k4, km4)
    v_camkii = get_rates_camkii(camkiis, camca4, ca, camkii_total, atp, pp1, k_asso, k_disso, k_disso_ca, k_disso2, k_disso_ca2, km_cam, kcat, km_atp, kcat_pp1, km_pp1)
    v_complex1 = get_rates_phospho(camkiis, camkii_total, complex1, kcat_complex1, kd_complex1)
    v_ck_mito = get_rates_phospho(camkiis, camkii_total, ck_mito, kcat_complex1, kd_complex1)
    return v_cam_ca + v_camkii + v_complex1 + v_ck_mito # Merge two lists


# Reference constants
K1 = 2.5
KM1 = KM2 = 0.05
K2 = 88.25
K3 = 12.5
KM3 = KM4 = 1.25
K4 = 250

K_ASSO = 2.1
K_DISSO_ALPHA = 1.4E-4
K_DISSO_DELTA = 0.5 * K_DISSO_ALPHA
K_DISSO_CA_ALPHA = 1.9E-3
K_DISSO_CA_DELTA = 0.5 * K_DISSO_CA_ALPHA
K_DISSO2_ALPHA = K_DISSO_ALPHA * 1e-3
K_DISSO2_DELTA = K_DISSO_DELTA * 1e-3
K_DISSO2_CA_ALPHA = K_DISSO_CA_ALPHA * 1e-3
K_DISSO2_CA_DELTA = K_DISSO_CA_DELTA * 1e-3

KM_CAM = 3E-5
KCAT_0_ALPHA = 1E-5
KCAT_0_DELTA = 6 * KCAT_0_ALPHA
KCAT_30_ALPHA = 10 * KCAT_0_ALPHA
KCAT_30_DELTA = 6 * KCAT_30_ALPHA
KCAT_37_ALPHA = 3 * KCAT_30_ALPHA
KCAT_37_DELTA = 6 * KCAT_37_ALPHA
KM_ATP = 19.1E-3
KCAT_PP1 = 1.72E-3
KM_PP1 = 11.0E-3

CAM_TOTAL = 6e-3
PP1 = 0.01E-3
CAMKII_TOTAL = 0.12E-3


