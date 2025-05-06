import numpy as np
from scipy.integrate import quad
from scipy.interpolate import UnivariateSpline
from scipy.optimize import minimize
import os
import pandas as pd
import warnings

# Physical constants
G = 4.302e-3  # Gravitational constant in pc (km/s)^2 / Msun
a0 = 1.2e-10 * 3.08568e16 / 1e6  # MOND critical acceleration in pc/s^2
epsilon = 1e-6  # Small constant to avoid singularities
w = 1.0  # Normalization constant in (km/s)^2 / pc^3, dimensionless as w=1

# List of 24 galaxies from the LaTeX document
galaxies = [
    'NGC0300', 'NGC0801', 'NGC0891', 'NGC1003', 'UGC01281', 'UGC05721', 'DDO154', 'DDO168',
    'NGC0247', 'UGC07524', 'IC2574', 'NGC0289', 'NGC3109', 'NGC0100', 'NGC0055', 'NGC0024',
    'NGC7331', 'NGC2403', 'NGC5055', 'UGC06930', 'NGC3198', 'UGC02259', 'NGC6946', 'UGC07232'
]

# Directory containing .dat files
data_dir = './data/'

# Parse galaxy properties from SPARC_Lelli2016c.mrt.docx
# Since the Word document is provided as text, we extract relevant data into a dictionary
sparc_data = {
    'NGC0300': {'R_d': 1.75, 'Distance': 2.08, 'Inclination': 42.0, 'L[3.6]': 2.922, 'MHI': 0.936, 'Vflat': 93.3, 'Quality': 2},
    'NGC0801': {'R_d': 8.72, 'Distance': 80.70, 'Inclination': 80.0, 'L[3.6]': 312.570, 'MHI': 23.201, 'Vflat': 220.1, 'Quality': 1},
    'NGC0891': {'R_d': 2.55, 'Distance': 9.91, 'Inclination': 90.0, 'L[3.6]': 138.340, 'MHI': 4.462, 'Vflat': 216.1, 'Quality': 1},
    'NGC1003': {'R_d': 1.61, 'Distance': 11.40, 'Inclination': 67.0, 'L[3.6]': 6.820, 'MHI': 5.880, 'Vflat': 109.8, 'Quality': 1},
    'UGC01281': {'R_d': 1.63, 'Distance': 5.27, 'Inclination': 90.0, 'L[3.6]': 0.353, 'MHI': 0.294, 'Vflat': 55.2, 'Quality': 1},
    'UGC05721': {'R_d': 0.38, 'Distance': 6.18, 'Inclination': 61.0, 'L[3.6]': 0.531, 'MHI': 0.562, 'Vflat': 79.7, 'Quality': 1},
    'DDO154': {'R_d': 0.37, 'Distance': 4.04, 'Inclination': 64.0, 'L[3.6]': 0.053, 'MHI': 0.275, 'Vflat': 47.0, 'Quality': 2},
    'DDO168': {'R_d': 1.02, 'Distance': 4.25, 'Inclination': 63.0, 'L[3.6]': 0.191, 'MHI': 0.413, 'Vflat': 53.4, 'Quality': 2},
    'NGC0247': {'R_d': 3.74, 'Distance': 3.70, 'Inclination': 74.0, 'L[3.6]': 7.332, 'MHI': 1.746, 'Vflat': 104.9, 'Quality': 2},
    'UGC07524': {'R_d': 3.46, 'Distance': 4.74, 'Inclination': 46.0, 'L[3.6]': 2.436, 'MHI': 1.779, 'Vflat': 79.5, 'Quality': 1},
    'IC2574': {'R_d': 2.78, 'Distance': 3.91, 'Inclination': 75.0, 'L[3.6]': 1.016, 'MHI': 1.036, 'Vflat': 66.4, 'Quality': 2},
    'NGC0289': {'R_d': 6.74, 'Distance': 20.80, 'Inclination': 46.0, 'L[3.6]': 72.065, 'MHI': 27.469, 'Vflat': 163.0, 'Quality': 2},
    'NGC3109': {'R_d': 1.56, 'Distance': 1.33, 'Inclination': 70.0, 'L[3.6]': 0.194, 'MHI': 0.477, 'Vflat': 66.2, 'Quality': 1},
    'NGC0100': {'R_d': 1.66, 'Distance': 13.50, 'Inclination': 89.0, 'L[3.6]': 3.232, 'MHI': 1.990, 'Vflat': 88.1, 'Quality': 1},
    'NGC0055': {'R_d': 6.11, 'Distance': 2.11, 'Inclination': 77.0, 'L[3.6]': 4.628, 'MHI': 1.565, 'Vflat': 85.6, 'Quality': 2},
    'NGC0024': {'R_d': 1.34, 'Distance': 7.30, 'Inclination': 64.0, 'L[3.6]': 3.889, 'MHI': 0.676, 'Vflat': 106.3, 'Quality': 1},
    'NGC7331': {'R_d': 5.02, 'Distance': 14.70, 'Inclination': 75.0, 'L[3.6]': 250.631, 'MHI': 11.067, 'Vflat': 239.0, 'Quality': 1},
    'NGC2403': {'R_d': 1.39, 'Distance': 3.16, 'Inclination': 63.0, 'L[3.6]': 10.041, 'MHI': 3.199, 'Vflat': 131.2, 'Quality': 1},
    'NGC5055': {'R_d': 3.20, 'Distance': 9.90, 'Inclination': 55.0, 'L[3.6]': 152.922, 'MHI': 11.722, 'Vflat': 179.0, 'Quality': 1},
    'UGC06930': {'R_d': 3.94, 'Distance': 18.00, 'Inclination': 32.0, 'L[3.6]': 8.932, 'MHI': 3.237, 'Vflat': 107.2, 'Quality': 1},
    'NGC3198': {'R_d': 3.14, 'Distance': 13.80, 'Inclination': 73.0, 'L[3.6]': 38.279, 'MHI': 10.869, 'Vflat': 150.1, 'Quality': 1},
    'UGC02259': {'R_d': 1.62, 'Distance': 10.50, 'Inclination': 41.0, 'L[3.6]': 1.725, 'MHI': 0.494, 'Vflat': 86.2, 'Quality': 2},
    'NGC6946': {'R_d': 2.44, 'Distance': 5.52, 'Inclination': 38.0, 'L[3.6]': 66.173, 'MHI': 5.670, 'Vflat': 158.9, 'Quality': 1},
    'UGC07232': {'R_d': 0.29, 'Distance': 2.83, 'Inclination': 59.0, 'L[3.6]': 0.113, 'MHI': 0.046, 'Vflat': 0.0, 'Quality': 2}
}

# Function to read .dat file
def read_dat_file(galaxy_name):
    file_path = os.path.join(data_dir, f"{galaxy_name}_rotmod.dat")
    if not os.path.exists(file_path):
        warnings.warn(f"Data file for {galaxy_name} not found at {file_path}")
        return None
    try:
        data = np.genfromtxt(file_path, skip_header=1, names=['Rad', 'Vobs', 'errV', 'Vgas', 'Vdisk', 'Vbul', 'SBdisk', 'SBbul'])
        return {
            'r': data['Rad'],
            'Vobs': data['Vobs'],
            'errV': data['errV'],
            'Vgas': data['Vgas'],
            'Vdisk': data['Vdisk'],
            'Vbul': data['Vbul'],
            'R_d': sparc_data[galaxy_name]['R_d']
        }
    except Exception as e:
        warnings.warn(f"Error reading {file_path}: {e}")
        return None

# Compute baryonic mass distribution
def compute_mass_distribution(r, Vgas, Vdisk, Vbul):
    Upsilon_disk = 0.5
    Upsilon_bulge = 0.7
    Vbaryon = np.sqrt((1.33 * Vgas)**2 + (np.sqrt(Upsilon_disk) * Vdisk)**2 + 
                      (np.sqrt(Upsilon_bulge) * Vbul)**2)
    M = Vbaryon**2 * r / G
    spline = UnivariateSpline(r, M, s=0)
    rho = spline.derivative()(r) / (4 * np.pi * r**2)
    return M, rho, spline, Vbaryon

# Improved ICD model (numerical)
def icd_model_numerical(r, M, rho, sigma, R_d):
    def integrand(x, r, rho):
        dist = np.sqrt((r - x)**2 + epsilon**2)
        return G * rho(x) / dist * np.exp(-abs(r - x) / sigma)
    
    a_newton = G * M(r) / r**2
    psi = []
    for ri in r:
        integral, _ = quad(lambda x: integrand(x, ri, rho), 0, 5 * R_d, limit=100)
        psi_val = w * sigma / (integral + epsilon)
        psi.append(psi_val)
    psi = np.array(psi)
    a_pred = a_newton * psi
    V_icd = np.sqrt(a_pred * r)
    return V_icd

# MOND model
def mond_model(r, Vbaryon):
    y = Vbaryon**2 / r / a0
    nu = np.sqrt((1 + np.sqrt(1 + 4/y)) / 2)
    V_mond = Vbaryon * nu
    return V_mond

# NFW model
def nfw_model(r, Mvir, c):
    def M_nfw(r, Mvir, c, rs):
        x = r / rs
        return Mvir * (np.log(1 + x) - x / (1 + x)) / (np.log(1 + c) - c / (1 + c))
    
    rs = 20
    M_r = M_nfw(r, Mvir, c, rs)
    V_nfw = np.sqrt(G * M_r / r)
    return V_nfw

# TeVeS model
def teves_model(r, Vbaryon):
    y = Vbaryon**2 / r / a0
    nu = np.sqrt((1 + np.sqrt(1 + 4/y)) / 2) * 1.1
    V_teves = Vbaryon * nu
    return V_teves

# f(R) model
def fr_model(r, Vbaryon, alpha):
    y = Vbaryon**2 / r / a0
    nu = np.sqrt((1 + np.sqrt(1 + 4/y)) / 2) * (1 + alpha * 0.1)
    V_fr = Vbaryon * nu
    return V_fr

# Compute chi-squared and BIC
def compute_chi2_bic(Vobs, errV, Vmodel, k):
    chi2 = np.sum(((Vobs - Vmodel) / errV)**2)
    n = len(Vobs)
    bic = chi2 + k * np.log(n)
    return chi2, bic

# Process each galaxy
results = []
for galaxy_name in galaxies:
    print(f"\nProcessing {galaxy_name}...")
    
    # Read data
    data = read_dat_file(galaxy_name)
    if data is None:
        print(f"Skipping {galaxy_name} due to missing or invalid data.")
        continue
    
    r, Vobs, errV, Vgas, Vdisk, Vbul, R_d = data['r'], data['Vobs'], data['errV'], data['Vgas'], data['Vdisk'], data['Vbul'], data['R_d']
    
    # Compute mass distribution
    M, rho, rho_spline, Vbaryon = compute_mass_distribution(r, Vgas, Vdisk, Vbul)
    
    # Optimize sigma for ICD
    def icd_objective(sigma):
        V_icd = icd_model_numerical(r, M, rho, sigma[0], R_d)
        return compute_chi2_bic(Vobs, errV, V_icd, k=1)[0]
    res_icd = minimize(icd_objective, x0=[5.0], bounds=[(1, 10)])
    sigma_opt = res_icd.x[0]
    
    # ICD numerical
    V_icd_num = icd_model_numerical(r, M, rho, sigma_opt, R_d)
    chi2_icd_num, bic_icd_num = compute_chi2_bic(Vobs, errV, V_icd_num, k=1)
    
    # MOND
    V_mond = mond_model(r, Vbaryon)
    chi2_mond, bic_mond = compute_chi2_bic(Vobs, errV, V_mond, k=1)
    
    # NFW
    def nfw_objective(params):
        Mvir, c = params
        V_nfw = nfw_model(r, Mvir, c)
        return compute_chi2_bic(Vobs, errV, V_nfw, k=2)[0]
    res_nfw = minimize(nfw_objective, x0=[1e12, 10], bounds=[(1e10, 1e14), (5, 20)])
    Mvir, c = res_nfw.x
    V_nfw = nfw_model(r, Mvir, c)
    chi2_nfw, bic_nfw = compute_chi2_bic(Vobs, errV, V_nfw, k=2)
    
    # TeVeS
    V_teves = teves_model(r, Vbaryon)
    chi2_teves, bic_teves = compute_chi2_bic(Vobs, errV, V_teves, k=1)
    
    # f(R)
    def fr_objective(alpha):
        V_fr = fr_model(r, Vbaryon, alpha[0])
        return compute_chi2_bic(Vobs, errV, V_fr, k=1)[0]
    res_fr = minimize(fr_objective, x0=[0.1], bounds=[(0, 1)])
    alpha = res_fr.x[0]
    V_fr = fr_model(r, Vbaryon, alpha)
    chi2_fr, bic_fr = compute_chi2_bic(Vobs, errV, V_fr, k=1)
    
    # Store results
    results.append({
        'Galaxy': galaxy_name,
        'ICD_Chi2': chi2_icd_num,
        'ICD_BIC': bic_icd_num,
        'ICD_Sigma': sigma_opt,
        'MOND_Chi2': chi2_mond,
        'MOND_BIC': bic_mond,
        'NFW_Chi2': chi2_nfw,
        'NFW_BIC': bic_nfw,
        'NFW_Mvir': Mvir / 1e12,
        'NFW_c': c,
        'TeVeS_Chi2': chi2_teves,
        'TeVeS_BIC': bic_teves,
        'fR_Chi2': chi2_fr,
        'fR_BIC': bic_fr,
        'fR_Alpha': alpha,
        'Best_Model': min({'ICD': bic_icd_num, 'MOND': bic_mond, 'NFW': bic_nfw, 
                           'TeVeS': bic_teves, 'f(R)': bic_fr}.items(), key=lambda x: x[1])[0]
    })

# Print results in a table
print("\nResults Table:")
print("-" * 140)
print(f"{'Galaxy':<10} {'ICD_Chi2':<10} {'ICD_BIC':<10} {'ICD_Sigma':<10} {'MOND_Chi2':<10} {'MOND_BIC':<10} "
      f"{'NFW_Chi2':<10} {'NFW_BIC':<10} {'NFW_Mvir':<10} {'NFW_c':<10} {'TeVeS_Chi2':<10} {'TeVeS_BIC':<10} "
      f"{'fR_Chi2':<10} {'fR_BIC':<10} {'fR_Alpha':<10} {'Best_Model':<15}")
print("-" * 140)
for res in results:
    print(f"{res['Galaxy']:<10} {res['ICD_Chi2']:<10.2f} {res['ICD_BIC']:<10.2f} {res['ICD_Sigma']:<10.2f} "
          f"{res['MOND_Chi2']:<10.2f} {res['MOND_BIC']:<10.2f} {res['NFW_Chi2']:<10.2f} {res['NFW_BIC']:<10.2f} "
          f"{res['NFW_Mvir']:<10.2f} {res['NFW_c']:<10.2f} {res['TeVeS_Chi2']:<10.2f} {res['TeVeS_BIC']:<10.2f} "
          f"{res['fR_Chi2']:<10.2f} {res['fR_BIC']:<10.2f} {res['fR_Alpha']:<10.2f} {res['Best_Model']:<15}")
print("-" * 140)