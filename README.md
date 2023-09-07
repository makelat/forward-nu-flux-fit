# forward-nu-flux-fit

A numerical package for estimating the discovery potential of existing and proposed forward neutrino experiments running concurrently with the Large Hadron Collider. 
A broad selection of neutrino samples derived from existing predictions for the parent hadron spectra are included, and can be combined in a model parametrizing the expected correlations in the neutrino spectra. 
This allows the determination of the highest achievable precision for their observation based on Fisher information.

Developed and maintained by Toni Mäkelä (toni.makela@cern.ch).
The provided physics results examples reproduce the work by TM, Felix Kling (felix.kling@desy.de), and Sebastian Trojanowski (sebastian.trojanowski@ncbj.gov.pl).

# installation requirements and usage
The program is implemented into the notebook main.ipynb. 
The requirements are working installations of
- `jupyter notebook`
- `python3`

Once these are there, the usage is as simple as opening the notebook main.ipynb e.g. in your favorite webbrowser and pressing the "fast-forward" button ("restart and run all").

# citations

Upon using this program and/or the predictions provided with it, please cite the following papers:

`SIBYLL 2.3d`

* E.-J. Ahn, R. Engel, T. K. Gaisser, P. Lipari, and T. Stanev, 
  “Cosmic ray interaction event generator SIBYLL 2.1,” 
  Phys. Rev. D 80 (2009) 094003, 
  arXiv:0906.4113 [hep-ph].

* E.-J. Ahn, R. Engel, T. K. Gaisser, P. Lipari, and T. Stanev, 
  “Sibyll with charm,” 
  in 16th International Symposium on Very High Energy Cosmic Ray Interactions. 2, 2011. 
  arXiv:1102.5705 [astro-ph.HE].

* F. Riehn, R. Engel, A. Fedynitch, T. K. Gaisser, and T. Stanev, 
  “A new version of the event generator Sibyll,” 
  PoS ICRC2015 (2016) 558, 
  arXiv:1510.00568 [hep-ph].

* A. Fedynitch, F. Riehn, R. Engel, T. K. Gaisser, and T. Stanev, 
  “Hadronic interaction model sibyll 2.3c and inclusive lepton fluxes,” 
  Phys. Rev. D 100 (2019) no. 10, 103018, 
  arXiv:1806.04140 [hep-ph].


`EPOS-LHC`

* T. Pierog, I. Karpenko, J. M. Katzy, E. Yatsenko, and K. Werner, 
  “EPOS LHC: Test of collective hadronization with data measured at the CERN Large Hadron Collider,” 
  Phys. Rev. C 92 (2015) no. 3, 034906, 
  arXiv:1306.0121 [hep-ph].


`DPMJET 3.2019.1`

* S. Roesler, R. Engel, and J. Ranft, 
  “The Monte Carlo event generator DPMJET-III,” 
  in International Conference on Advanced Monte Carlo for Radiation Physics, Particle Transport Simulation and Applications (MC 2000), pp. 1033–1038. 12, 2000. 
  arXiv:hep-ph/0012252.

* A. Fedynitch, 
  “Cascade equations and hadronic interactions at very high energies”. 
  PhD thesis, KIT, Karlsruhe,
  Dept. Phys., 11, 2015.


`QGSJET II-04`

* S. Ostapchenko, 
  “Monte Carlo treatment of hadronic interactions in enhanced Pomeron scheme: I. QGSJET-II model,” 
  Phys. Rev. D 83 (2011) 014018, 
  arXiv:1010.1869 [hep-ph].


`Pythia 8.2 (forward tune)`

* M. Fieg, F. Kling, H. Shulz, and T. Sjostrand, 
  “Tuning pythia for forward physics experiments,”. 
  In preparation


`BKRS`

* L. Buonocore, F. Kling, L. Rottoli, and J. Sominka. 
  In preparation.


`BDGJKR`

* W. Bai, M. Diwan, M. V. Garzelli, Y. S. Jeong, and M. H. Reno, 
  “Far-forward neutrinos at the Large Hadron Collider,” 
  JHEP 06 (2020) 032, 
  arXiv:2002.03012 [hep-ph].

* W. Bai, M. Diwan, M. V. Garzelli, Y. S. Jeong, F. K. Kumar, and M. H. Reno, 
  “Parton distribution function uncertainties in theoretical predictions for far-forward tau neutrinos at the Large Hadron Collider,” 
  JHEP 06 (2022) 148, 
  arXiv:2112.11605 [hep-ph].

* W. Bai, M. Diwan, M. V. Garzelli, Y. S. Jeong, K. Kumar, and M. H. Reno, 
  “Forward production of prompt neutrinos from charm in the atmosphere and at high energy colliders,” 
  arXiv:2212.07865 [hep-ph].


`BKSS kT`

* A. Bhattacharya, F. Kling, I. Sarcevic, and A. M. Stasto, 
  “Forward Neutrinos from Charm at Large Hadron Collider,” 
  arXiv:2306.01578 [hep-ph].


`MS kT`

* R. Maciula and A. Szczurek, 
  “Far-forward production of charm mesons and neutrinos at forward physics facilities at the LHC and the intrinsic charm in the proton,” 
  Phys. Rev. D 107 (2023) no. 3, 034002,
  arXiv:2210.08890 [hep-ph].


