Integrated Collective Dynamics (ICD) for Galactic Rotation Curves
This repository contains the Python code, data, and supplementary materials for the research paper "A Novel Framework for Modeling Galactic Rotation Curves: Integrated Collective Dynamics" by Esmail Dabbaghsaz and Alireza Mohebalhojeh. The study introduces the Integrated Collective Dynamics (ICD) framework to model galactic rotation curves as an alternative to dark matter and modified gravity models, applied to 24 galaxies from the SPARC database.
Repository Overview
The repository includes:

Python Code: A script (Icd_Code.py) to compute rotation curves using the ICD framework, alongside comparisons with MOND, NFW, TeVeS, and $f(R)$ models.
Data: Rotation curve data for 24 galaxies in .dat format, sourced from the SPARC database, and an additional metadata file (SPARC_Lelli2016c.mrt.docx) containing galaxy properties.
Documentation: This README provides instructions for using the code and understanding the data.

Data Description
Galaxy Rotation Curve Data
The data/ directory contains 24 .dat files, one for each galaxy analyzed in the study. Each file (e.g., NGC0247_rotmod.dat) provides rotation curve data with the following columns:

Rad: Radial distance from the galactic center (kpc).
Vobs: Observed rotation velocity (km/s).
errV: Uncertainty in the observed velocity (km/s).
Vgas: Velocity contribution from gas (km/s).
Vdisk: Velocity contribution from the disk (km/s).
Vbul: Velocity contribution from the bulge (km/s).
SBdisk: Surface brightness of the disk.
SBbul: Surface brightness of the bulge.

The galaxies included are:

NGC0300, NGC0801, NGC0891, NGC1003, UGC01281, UGC05721, DDO154, DDO168, NGC0247, UGC07524, IC2574, NGC0289, NGC3109, NGC0100, NGC0055, NGC0024, NGC7331, NGC2403, NGC5055, UGC06930, NGC3198, UGC02259, NGC6946, UGC07232.

Galaxy Properties
The file SPARC_Lelli2016c.mrt.docx contains metadata for the 24 galaxies, extracted from the SPARC database. It includes:

Distance (Mpc)
Inclination (degrees)
Luminosity at 3.6 μm ($L_{[3.6]}$, in $10^9 L_\odot$)
Disk scale length ($R_d$, in kpc)
Neutral hydrogen mass ($M_{\text{HI}}$, in $10^9 M_\odot$)
Flat rotation velocity ($V_{\text{flat}}$, in km/s)
Data quality (Q)

This metadata is embedded in the Python script for convenience.
Code Description
The Icd_Code.py script implements the ICD framework to model galactic rotation curves and compares its performance with alternative models (MOND, NFW, TeVeS, and $f(R)$). Key features include:

Data Parsing: Reads .dat files for each galaxy.
Mass Distribution: Computes baryonic mass and density using observed velocity components.
ICD Model: Numerically evaluates the ICD correction function and predicts rotation velocities.
Comparison Models: Implements MOND, NFW, TeVeS, and $f(R)$ models for benchmarking.
Optimization: Dynamically optimizes parameters (e.g., $\sigma$ for ICD, $M_{\text{vir}}$ and $c$ for NFW) using scipy.optimize.minimize.
Metrics: Calculates $\chi^2$ and Bayesian Information Criterion (BIC) to evaluate model fits.

Dependencies

Python 3.x
NumPy
SciPy
Pandas
Warnings (for error handling)

Install dependencies using:
pip install numpy scipy pandas

Usage Instructions

Clone the Repository:
git clone https://github.com/ICDGalacticDynamics.git
cd ICDGalacticDynamics


Prepare the Data:

Ensure the data/ directory contains the 24 .dat files.
The SPARC_Lelli2016c.mrt.docx file is not directly used by the script but is provided for reference. Its data is embedded in Icd_Code.py.


Run the Script:Execute the Python script to process all galaxies and generate results:
python Icd_Code.py

The script will:

Read data for each galaxy.
Compute rotation curves using ICD, MOND, NFW, TeVeS, and $f(R)$ models.
Optimize model parameters.
Output a table with $\chi^2$, BIC, and best-fit parameters for each model.


Output:The script prints a table summarizing the results for each galaxy, including:

$\chi^2$ and BIC for each model.
Optimized parameters (e.g., $\sigma$ for ICD, $M_{\text{vir}}$ and $c$ for NFW).
The best model based on the lowest BIC.



Reproducibility
The code is designed to ensure reproducibility:

All data sources are publicly available or provided in the repository.
The script uses deterministic numerical methods (e.g., Clenshaw-Curtis quadrature for integration).
Parameter optimization is performed with bounded constraints to ensure physically meaningful results.

To reproduce the results from the paper:

Use the provided .dat files or download the latest SPARC data from http://astroweb.cwru.edu/SPARC/.
Run the script as described above.
Compare the output table with Table 5 in the paper.

Citation
If you use this code or data in your research, please cite the following paper:
Dabbaghsaz, E., & Mohebalhojeh, A. (2025). A Novel Framework for Modeling Galactic Rotation Curves: Integrated Collective Dynamics. [Include journal or arXiv link if available]
Acknowledgments

The SPARC team for providing high-quality rotation curve data.
Stacy McGaugh and Federico Lelli for their insights into galactic dynamics.
Iran Meteorological Organization and Institute of Geophysics, University of Tehran, for their support.

License
This project is licensed under the MIT License. See the LICENSE file for details.

For any questions or issues, please open an issue on this GitHub repository.
