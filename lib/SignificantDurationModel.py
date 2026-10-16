# Based on the script by Baker Research Group: https://www.jackwbaker.com/GMMs_archive.html 

import numpy as np

def abrahamson_silva_ds_1999(magnitude=7.0,distance=10.0,soil=True,duration_type='DS575H'):
    
    """
    Significant duration model by Abrahamson and Silva (1996) Empirical ground motion 
    models, report prepared for Brookhaven National Laboratory.
    Input
    magnitude: earthquake magnitude
    distance: site-rupture distance
    soil: True for foil prediction, False for rock prediction
    duration_type: 'DS575H': Ds5-75 (horizontal), 'DS575V': Ds5-75 (vertical)
                   'DS595H': Ds5-95 (horizontal), 'DS595V': Ds5-95 (vertical)
    Output:
    log(ds_median): log(median) significant duration prediction
    ds_sigma: logarithmic standard deviation of the prediction
    """


    # map the duration_type to integer key
    dur_map = {'DS575H': 0,
               'DS575V': 1,
               'DS595H': 2,
               'DS595V': 3}
    dur_tag = dur_map.get(duration_type.upper(),None)
    if dur_tag is None:
        print("SignificantDurationModel.abrahamson_silva_ds_1999: duration_type='DS575H','DS575V','DS595H','DS595V'?")
        return None, None
    # modeling coefficients
    beta = [3.2, 3.2, 3.2, 3.2]
    b1 = [5.204, 4.610, 5.204, 4.610]
    b2 = [0.851, 1.536, 0.851, 1.536]
    m_star = [6, 6, 6, 6]
    c1 = [0.805, 1.076, 0.805, 1.076]
    c2 = [0.063, 0.107, 0.063, 0.107]
    rc = [10, 10, 10, 10]
    Drat = [0.000, 0.000, 0.845, 0.646]
    sigma = [0.55, 0.46, 0.49, 0.45]
    # median
    if distance > rc[dur_tag]:
        ds_median = np.exp(np.log((np.exp(b1[dur_tag]+b2[dur_tag]* \
                    (magnitude-m_star[dur_tag]))/(10**(1.5*magnitude+ \
                    16.05)))**(-1/3)/(4.9e6*beta[dur_tag])+soil* \
                    c1[dur_tag]+c2[dur_tag]*(distance-rc[dur_tag]))+Drat[dur_tag])
    else:
        ds_median = np.exp(np.log((np.exp(b1[dur_tag]+b2[dur_tag]* \
                    (magnitude-m_star[dur_tag]))/(10**(1.5*magnitude+ \
                    16.05)))**(-1/3)/(4.9e6*beta[dur_tag])+soil* \
                    c1[dur_tag])+Drat[dur_tag])
    # sigma
    ds_sigma = sigma[dur_tag]

    # return
    return np.log(ds_median), ds_sigma


def bommer_stafford_alarcon_ds_2009(magnitude=7.0, distance=10.0, vs30=760.0, ztor=0.0, duration_type='DS575H'):
    """
    Singificant duration model by Bommer, Stafford, Alarcon (2009) Empirical 
    Equations for the Prediction of the Significant, Bracketed, and Uniform 
    Duration of Earthquake Ground Motion
    Input
    magnitude: earthquake magnitude
    distance: site-rupture distance
    vs30: average soil shear-wave velocity over the top 30 meters
    ztor: depth to top of rupture (km)
    duration_type: 'DS575H': Ds5-75 (horizontal), 'DS595H': Ds5-95 (horizontal)
    Output:
    log(ds_median): log(median) significant duration prediction
    ds_sigma: logarithmic standard deviation of the prediction
    ds_tau: within-event logarithmic standard deviation
    ds_phi: between-event logarithmic standard deviation
    """

    # duration type map
    dur_map = {'DS575H':0, 'DS595H': 1}
    dur_tag = dur_map.get(duration_type.upper(), None)
    if dur_tag is None:
        print("SignificantDurationModel.bommer_stafford_alarcon_ds_2009: duration_type='DS575H','DS595H'?")
        return None, None, None, None
    
    # modeling coefficients
    c0 = [-5.6298, -2.2393]
    m1 = [1.2619, 0.9368]
    r1 = [2.0063, 1.5686]
    r2 = [-0.2520, -0.1953]
    h1 = [-2.3316, 2.5000]
    v1 = [-0.2900, -0.3478]
    z1 = [-0.0522, -0.0365]
    tauCoeff = [0.3527, 0.3252]
    phiCoeff = [0.4304, 0.3460]
    sigma_c = [0.1729, 0.1114]
    sigma_Tgm = [0.5289, 0.4616]

    # median
    ds_median = np.exp(c0[dur_tag]+m1[dur_tag]*magnitude+(r1[dur_tag]+ \
                r2[dur_tag]*magnitude)*np.log(np.sqrt(distance**2+h1[dur_tag]**2))+ \
                v1[dur_tag]*np.log(vs30)+z1[dur_tag]*ztor)
    # standard deviations
    ds_sigma = sigma_Tgm[dur_tag]
    ds_tau = tauCoeff[dur_tag]
    ds_phi = phiCoeff[dur_tag]

    # return
    return np.log(ds_median), ds_sigma, ds_tau, ds_phi


def afshari_stewart_ds_2016(magnitude=7.0, distance=10.0, vs30=760.0, mechanism='unknown', 
                            z1=None, region='california', duration_type='DS575H'):
    """
    Significant duration model by Afshari and Stewart (2016) hysically Parameterized
    Prediction Equations for Significant Duration in Active Crustal Regions
    Input
    magnitude: earthquake magnitude
    distance: site-rupture distance
    vs30: average soil shear-wave velocity over the top 30 meters
    mechanism: 'unknown', 'normal', 'reverse', 'strike-slip'
    z1: depth to shear velocity of 1 km/s isosurface (m)
    region: 'california', 'japan', 'other'
    duration_type: 'DS575H': Ds5-75 (horizontal), 'DS595H': Ds5-95 (horizontal), 'DS2080H': Ds20-80 (horizontal)
    Output:
    log(ds_median): log(median) significant duration prediction
    ds_sigma: logarithmic standard deviation of the prediction
    ds_tau: within-event logarithmic standard deviation
    ds_phi: between-event logarithmic standard deviation
    """

    # mechanism map
    mech_map = {'unknown':0, 'normal': 1, 'reverse': 2, 'strike-slip': 3}
    mech_tag = mech_map.get(mechanism.lower(), None)
    if mech_tag is None:
        print("SignificantDurationModel.afshari_stewart_ds_2016: mechanism='unknown','normal','reverse','strike-slip'?")
        return None, None, None, None
    # region map
    reg_map = {'california':0, 'japan': 1, 'other': 2}
    reg_tag = reg_map.get(region.lower(), None)
    if reg_tag is None:
        print("SignificantDurationModel.afshari_stewart_ds_2016: region='california', 'japan', 'other'?")
        return None, None, None, None
    # duration type map
    dur_map = {'DS575H':0, 'DS595H': 1, 'DS2080H': 2}
    dur_tag = dur_map.get(duration_type.upper(), None)
    if dur_tag is None:
        print("SignificantDurationModel.afshari_stewart_ds_2016: duration_type='DS575H','DS595H','DS2080H'?")
        return None, None, None, None
    
    # source coefficients
    M1 = [5.35, 5.20, 5.20]
    M2 = [7.15, 7.40, 7.40]
    b0 = [[1.2800, 2.1820, 0.8822],
          [1.5550, 2.5410, 1.4090],
          [0.7806, 1.6120, 0.7729],
          [1.2790, 2.3020, 0.8804]]
    b1 = [[5.576, 3.628, 6.182],
          [4.992, 3.170, 4.778],
          [7.061, 4.536, 6.579],
          [5.578, 3.467, 6.188]]
    b2 = [0.9011, 0.9443, 0.7414]
    b3 = [-1.684, -3.911, -3.164]
    Mstar = [6, 6, 6]
    # path coefficients
    c1 = [0.1159, 0.3165, 0.0646]
    RR1 = [10, 10, 10]
    RR2 = [50, 50, 50]
    c2 = [0.1065, 0.2539, 0.0865]
    c3 = [0.0682, 0.0932, 0.0373]
    # site coefficients
    c4 = [-0.2246, -0.3183, -0.4237]
    Vref = [368.2, 369.9, 369.6]
    V1 = [600, 600, 600]
    c5 = [0.0006, 0.0006, 0.0005]
    dz1ref = [200, 200, 200]
    # standard deviation coefficients
    phi1 = [0.54, 0.43, 0.56]
    phi2 = [0.41, 0.35, 0.45]
    tau1 = [0.28, 0.25, 0.30]
    tau2 = [0.25, 0.19, 0.19]

    # basin depth
    if reg_tag == 0:
        mu_z1 = np.exp(-7.15/4*np.log((vs30**4+570.94**4)/(1360**4+570.94**4)))
    else:
        mu_z1 = np.exp(-5.23/4*np.log((vs30**4+412.39**4)/(1360**4+412.39**4)))
    # differential basin depth
    if z1 is None or z1 < 0 or reg_tag == 2:
        dz1 = 0
    else:
        dz1 = z1-mu_z1

    # source term
    if magnitude < M1[dur_tag]:
        F_E = b0[mech_tag][dur_tag]
    else:
        if magnitude < M2[dur_tag]:
            deltaSigma = np.exp(b1[mech_tag][dur_tag]+b2[dur_tag]*(magnitude-Mstar[dur_tag]))
        else:
            deltaSigma = np.exp(b1[mech_tag][dur_tag]+b2[dur_tag]*(M2[dur_tag]-Mstar[dur_tag])+ \
                         b3[dur_tag]*(magnitude-M2[dur_tag]))
        
        M_0 = 10**(1.5*magnitude+16.05)
        f_0 = 4.9e6*3.2*(deltaSigma/M_0)**(1/3)
        F_E = 1/f_0
    # path term
    if distance < RR1[dur_tag]:
        F_P = c1[dur_tag]*distance
    elif distance < RR2[dur_tag]:
        F_P = c1[dur_tag]*RR1[dur_tag]+c2[dur_tag]*(distance-RR1[dur_tag])
    else:
        F_P = c1[dur_tag]*RR1[dur_tag]+c2[dur_tag]*(RR2[dur_tag]-RR1[dur_tag])+c3[dur_tag]*(distance-RR2[dur_tag])    
    # F_deltaz term
    if dz1 <= dz1ref[dur_tag]:
        F_deltaz = c5[dur_tag]*dz1
    else:
        F_deltaz = c5[dur_tag]*dz1ref[dur_tag]
    # site term
    if vs30 < V1[dur_tag]:
        F_S = c4[dur_tag]*np.log(vs30/Vref[dur_tag])+F_deltaz
    else:
        F_S = c4[dur_tag]*np.log(V1[dur_tag]/Vref[dur_tag])+F_deltaz

    # median
    ds_median = np.exp(np.log(F_E+F_P)+F_S)
    # standard deviations
    # between event
    if magnitude < 5.5:
        ds_phi = phi1[dur_tag]
    elif magnitude < 5.75:
        ds_phi = phi1[dur_tag]+(phi2[dur_tag]-phi1[dur_tag])*(magnitude-5.5)/(5.75-5.5)
    else:
        ds_phi = phi2[dur_tag]
    # within event
    if magnitude < 6.5:
        ds_tau = tau1[dur_tag]
    elif magnitude < 7:
        ds_tau = tau1[dur_tag]+(tau2[dur_tag]-tau1[dur_tag])*(magnitude-6.5)/(7-6.5)
    else:
        ds_tau = tau2[dur_tag]
    # total
    ds_sigma = np.sqrt(ds_phi**2+ds_tau**2)

    # return
    return np.log(ds_median), ds_sigma, ds_tau, ds_phi


def pinilla_ramos_abrahamson_ds_2025(M, Rrup, Vs30, event_type="interface", region="JP",
                                     epistemic_other_region=False):
    """
    Predict D5-75 significant duration using Pinilla-Ramos & Abrahamson (2025).

    Parameters
    ----------
    M           : float
        Moment magnitude
    Rrup        : float
        Closest rupture distance (km)
    Vs30        : float
        Time-averaged shear-wave velocity in top 30 m (m/s)
    event_type  : str
        "interface" or "slab" (intraslab)
    region      : str
        "JP" (Japan), "NZ" (New Zealand), "SA" (South America), "TW" (Taiwan, slab only)
        For Cascadia/US/Alaska or other unlisted regions → use "JP"
    epistemic_other_region : bool
        If True, adds 0.10 epistemic uncertainty (ln units) to sigma_ln,
        as recommended by the authors for regions outside JP/NZ/SA/TW.

    Returns
    -------
    dict with keys:
        D575_median_s   : median D5-75 in seconds
        sigma_total     : total std dev in s^n2 units (power-normal, NOT ln units)
        sigma_ln        : equivalent std dev of ln(D5-75), (ln D84 - ln D16)/2,
                          including the epistemic term if epistemic_other_region=True
        D575_16th_s     : 16th percentile D5-75 (seconds)
        D575_84th_s     : 84th percentile D5-75 (seconds)
        n2              : power-normal exponent used for sigma
        components      : dict of Dsource, Dpath, Dsite (in seconds)
        rho_PGA_dur     : correlation coefficient between lnPGA and D5-75 epsilon
    """
    # ── Table 1: Region-independent model coefficients ────────────────────────
    coeffs_global = {
        "slab": {
            "c1":  0.340,   # source amplitude at Mth
            "c2":  0.500,   # exponential magnitude scaling exponent
            "d1":  2.770,   # linear scaling slope above Mth
            "M0":  4.250,   # reference magnitude for source term
            "Mth": 6.000,   # magnitude hinge point (exponential -> linear)
            "c31": 0.013,   # bilinear path slope (R <= R1)
            "c32": 0.020,   # bilinear path slope (R > R1)
            "R1":  250.0,   # path breakpoint distance (km)
            "V1":  3100.0,  # Vs30 above which site effect = 0 (m/s)
            "n1":  0.7,     # power-normal exponent (within-event residuals)
            "n2":  0.25,    # power-normal exponent (total residuals)
            "phi_S": 0.744,
            "phi_P": 0.687,
            "tau":   0.688,
            "sigma": 1.466,
        },
        "interface": {
            "c1":  5.268,   # source amplitude at Mth
            "c2":  0.275,   # magnitude scaling (fixed to anchor Maule & Tohoku)
            "d1":  None,    # not used (single exponential for interface)
            "M0":  None,    # not used
            "Mth": 7.000,   # reference magnitude for source term
            "c31": 0.031,   # bilinear path slope (R <= R1)
            "c32": 0.029,   # bilinear path slope (R > R1)
            "R1":  250.0,   # path breakpoint distance (km)
            "V1":  3100.0,  # Vs30 above which site effect = 0 (m/s)
            "n1":  0.7,
            "n2":  0.15,
            "phi_S": 0.756,
            "phi_P": 0.671,
            "tau":   0.827,
            "sigma": 1.665,
        },
    }

    # ── Table 2: Regionalized coefficients (c3: distance slope, c4: site scaling)
    # Note: interface Taiwan excluded. For other regions use Japan + 10% epistemic.
    coeffs_regional = {
        "slab": {
            "JP": {"c3": 0.046, "c4": -1.480},
            "NZ": {"c3": 0.032, "c4": -1.252},
            "SA": {"c3": 0.042, "c4": -1.510},
            "TW": {"c3": 0.038, "c4": -1.151},
        },
        "interface": {
            "JP": {"c3": 0.030, "c4": -2.485},
            "NZ": {"c3": 0.024, "c4": -2.420},
            "SA": {"c3": 0.024, "c4": -2.248},
            # TW excluded for interface
        },
    }

    # ── Table 3: Standard deviation model coefficients (Eq. 9)
    # sigma(s^n2) = a0 + a1*(R/100) + a2*(R/100)^2 + b1*M + b2*M^2 + c1*ln(Vs30)
    sigma_coeffs = {
        "slab": {
            "a0":  0.3035,
            "a1": -0.0188,
            "a2":  0.0040,
            "b1": -0.0288,
            "b2":  0.0019,
            "c1":  0.0038,
        },
        "interface": {
            "a0":  0.3107,
            "a1": -0.0393,
            "a2":  0.0062,
            "b1": -0.0519,
            "b2":  0.0041,
            "c1":  0.0028,
        },
    }

    # Correlation between lnPGA and D5-75 epsilons
    rho_pga_dur = {
        "slab":      -0.376,
        "interface": -0.539,
    }

    # V3: reference Vs30 for site term (zero-crossing point, ~800 m/s rock)
    V3 = 800.0

    # ── Input validation ──────────────────────────────────────────────────────
    event_type = event_type.lower()
    region = region.upper()

    if event_type not in ("interface", "slab"):
        raise ValueError("event_type must be 'interface' or 'slab'")
    if event_type == "interface" and region == "TW":
        raise ValueError("Taiwan (TW) is excluded from the interface model. Use JP, NZ, or SA.")
    if region not in coeffs_regional[event_type]:
        raise ValueError(
            f"Region '{region}' not available for {event_type}. "
            f"Valid: {list(coeffs_regional[event_type].keys())}. "
            f"For Cascadia/US use region='JP' with epistemic_other_region=True."
        )

    c   = coeffs_global[event_type]
    c3  = coeffs_regional[event_type][region]["c3"]
    c4  = coeffs_regional[event_type][region]["c4"]
    n2  = c["n2"]
    V1  = c["V1"]
    R1  = c["R1"]
    sc  = sigma_coeffs[event_type]

    # ── Source term Dsource(M) ────────────────────────────────────────────────
    # Interface: single exponential  Dsource = c1 * 10^(c2 * (M - Mth))
    # Slab:      two-part model
    #   M <= Mth: Dsource = c1 * 10^(c2 * (M - M0))
    #   M >  Mth: Dsource = c1 * 10^(c2*(Mth-M0)) + d1*(M-Mth)/(8.5-Mth)
    if event_type == "interface":
        D_source = c["c1"] * 10 ** (c["c2"] * (M - c["Mth"]))
    else:
        if M <= c["Mth"]:
            D_source = c["c1"] * 10 ** (c["c2"] * (M - c["M0"]))
        else:
            base     = c["c1"] * 10 ** (c["c2"] * (c["Mth"] - c["M0"]))
            D_source = base + c["d1"] * (M - c["Mth"]) / (8.5 - c["Mth"])

    # ── Path term Dpath(Rrup) ─────────────────────────────────────────────────
    # Dpath = c3(reg)*Rrup + f_path(Rrup)   [Eq. 5, median ΔP=1]
    # f_path = c31*Rrup           for Rrup <= R1
    #        = c31*R1 + (Rrup-R1)*c32   for Rrup > R1
    if Rrup <= R1:
        f_path = c["c31"] * Rrup
    else:
        f_path = c["c31"] * R1 + (Rrup - R1) * c["c32"]
    D_path = c3 * Rrup + f_path

    # ── Site term Dsite(Vs30) ─────────────────────────────────────────────────
    # Dsite = c4(reg) * ln(Vs30 / V3)   for Vs30 <= V1
    #       = 0                           for Vs30 > V1
    D_site = 0.0 if Vs30 > V1 else c4 * np.log(Vs30 / V3)

    # ── Median duration (Eq. 11) ──────────────────────────────────────────────
    mu_median_s = D_source + D_path + D_site

    # ── Standard deviation sigma(s^n2) (Eq. 9, Table 3) ──────────────────────
    R100  = Rrup / 100.0
    sigma = (sc["a0"]
             + sc["a1"] * R100
             + sc["a2"] * R100**2
             + sc["b1"] * M
             + sc["b2"] * M**2
             + sc["c1"] * np.log(Vs30))

    # ── Percentiles using power-normal distribution (Eq. 10) ─────────────────
    mu_n2     = mu_median_s ** n2
    D575_16th = max(0.0, mu_n2 - sigma) ** (1.0 / n2)
    D575_84th = (mu_n2 + sigma) ** (1.0 / n2)

    # ── Equivalent lognormal standard deviation ──────────────────────────────
    # sigma is defined for D^n2, so it is converted to ln units via the 16th/84th
    # percentiles; the epistemic term (ln units) is then added in ln space
    sigma_ln = (np.log(D575_84th) - np.log(D575_16th)) / 2.0
    if epistemic_other_region:
        sigma_ln = np.sqrt(sigma_ln**2 + 0.10**2)

    return {
        "D575_median_s":  mu_median_s,
        "sigma_total":    sigma,
        "sigma_ln":       sigma_ln,
        "D575_16th_s":    D575_16th,
        "D575_84th_s":    D575_84th,
        "n2":             n2,
        "components": {
            "Dsource_s": D_source,
            "Dpath_s":   D_path,
            "Dsite_s":   D_site,
        },
        "rho_PGA_dur": rho_pga_dur[event_type],
    }
