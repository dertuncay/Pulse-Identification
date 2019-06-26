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

Bibtex:
@Article{Ertuncay2019,
author="Ertuncay, Deniz
and Costa, Giovanni",
title="An alternative pulse classification algorithm based on multiple wavelet analysis",
journal="Journal of Seismology",
year="2019",
month="Jun",
day="26",
abstract="Near fault ground motions may contain impulse behavior on velocity records. Such signals have a particular indicator which makes it possible to distinguish them from non-impulsive signals. These signals have significant effects on structures; therefore, they have been investigated for more than 20 years. In this study, we used Ricker and Morlet wavelets in order to analyze the wavelet power spectrum of the strong motion signals to investigate the impulsiveness. Both the area around the PGV and the area that exceeds the minimum threshold for the energy function are used in order to determine the position of the pulse. On both of these cases, particular criteria are used in order to characterize the signal. Then, we calculate the pulse period of the pulse region. Ricker and Morlet wavelets are also used to mimic the pulse signal. This method provides advanced information about the position of the maximum energy of the pulse part of the signal. We found that the impulsive part of the signal is frequently at the position where PGV occurs and the Ricker wavelet is better than the Morlet wavelet on mimicking the pulse part of the waveform. Spectral responses of strong motion waveform and the wavelets have strong correlation at around pulse period. Results show consistency with previous studies; hence, it can be used as a solid alternative on pulse shape signal investigations.",
issn="1573-157X",
doi="10.1007/s10950-019-09845-y",
url="https://doi.org/10.1007/s10950-019-09845-y"
}

