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

Ertuncay, Deniz, and Giovanni Costa. "An alternative pulse classification algorithm based on multiple wavelet analysis." Journal of Seismology 23.4 (2019): 929-942. https://doi.org/10.1007/s10950-019-09845-y
