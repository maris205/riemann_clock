# ESPRESSO_HE0515-4414

<a href="https://doi.org/10.5281/zenodo.5512490"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.5512490.svg" alt="DOI"></a>

VLT/ESPRESSO spectra of HE0515-4414 and absorption profile fits associated with Murphy et al. (2021, A&amp;A, accepted, <a href="https://arxiv.org/abs/2112.05819">arXiv:2112.05819</a>).

Read this README, and view/download/use the files within this repository, in conjunction with a careful read of the paper itself.

If you use any of the materials in this repository, please cite the paper. If you want to cite only the data and/or fits (for some reason), please use the DOI:
<a href="https://doi.org/10.5281/zenodo.5512490"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.5512490.svg" alt="DOI"></a>


<!---
The paper is available at the following websites, in published or pre-print form:<br>
&ndash; arXiv.org:  <a href="https://arxiv.org/abs/1708.00014">arxiv.org/abs/1708.00014</a><br>
&ndash; MNRAS (via DOI): <a href="https://doi.org/10.1093/mnras/stx1949">10.1093/mnras/stx1949</a><br>
&ndash; NASA/ADS: <a href="http://adsabs.harvard.edu/cgi-bin/bib_query?arXiv:1708.00014">2017arXiv170800014M</a> (to be updated once final version is published)<br>
--->

This repository contains the following folders:<br><br>

&ndash; <b>DRS_extracted_spectra:</b><br>
FITS files produced by the ESPRESSO Data Reduction Software (DRS) containing the extracted, wavelength calibrated spectrum of both traces, of each echelle order, for each of the 17 quasar exposures.<br><br>

&ndash; <b>Final_spectrum:</b><br>
FITS file of the final, combined spectrum (hes0515m44114.fits) produced by the command in the UVES_popler.command file using the hes0515m4414.upl UVES_popler Log (UPL) file and the atmospheric line mask in atmomask.dat.<br><br>

&ndash; <b>Fitting_results:</b>
<br>Absorption profile parameter files (fit\_[lmr]\_iso.f13) containing the input parameter values used to produce the final fitting sequence recorded in the corresponding fit\_[lmr]\_iso.f18 files. The absorption profile fit was run in VPFIT using the commented-out command in the first line of each fit*.f13 file. There are 3 fits, one for each of the "l"eft, "m"iddle and "r"ight regions. MM_VPFIT_2013-11-10.dat is a copy of the atomic data file used by VPFIT from the <a href="https://github.com/MTMurphy77/MMatomdat">repository</a> associated with <a href="http://adsabs.harvard.edu/abs/2014MNRAS.438..388M">Murphy & Berengut (2014, MNRAS, 438, 388, arXiv:1311.2949)</a>. The fit_[lmr].f1[38] files are the corresponding fits when no isotopic structures were used in fitting the various transitions (i.e. using the MM_VPFIT_2013-11-10_noiso.dat atomic data file). vp_setup.dat is the setup file for VPFIT used for all fits.
