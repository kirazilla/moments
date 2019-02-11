Compute spectral line moments
Language: Python3

Python packages:
matplotlib
docopt
numpy
astropy

Compute the equivalent width and first three moments (radial velocity, variance and 
skewness) for a line profile or CCF profile with xrange in km/s.
The moments are computed following Aerts et al. (1992, A&A 226, 294, 
http://adsabs.harvard.edu/abs/1992A%26A...266..294A).
The errors are computed using statistical uncertainties as in the FAMIAS software 
(http://www.ster.kuleuven.be/~zima/famias/famias_manual/Data_Manager.html).
