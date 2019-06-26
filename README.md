# Pulse-Identification

Inputs of the function:

* st - Obspy stream

* scale - Magnifier for the trace (default = 1)

* zero_adding - Adding 100 data point to beginning and the end of the data (default = False)

The code uses the modified version of the article Torrence, Christopher, and Gilbert P. Compo. "A practical guide to wavelet analysis." Bulletin of the American Meteorological society 79.1 (1998): 61-78.

Github: [Link](https://github.com/chris-torrence/wavelets)

waveletAnalysis.py and waveletFunctions.py are related to above-mentioned article.

Following python packages must be installed on the system:

* scipy
* numpy
* obspy

# Citation

Ertuncay, D. & Costa, G. J Seismol (2019). https://doi.org/10.1007/s10950-019-09845-y
