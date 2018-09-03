import numpy as np
from waveletFunctions import wavelet, wave_signif
import matplotlib.pylab as plt
from matplotlib.gridspec import GridSpec
import matplotlib.ticker as ticker
from matplotlib.pyplot import savefig
from mpl_toolkits.axes_grid1 import make_axes_locatable

# WAVETEST Example Python script for WAVELET, using NINO3 SST dataset
#
# See "http://paos.colorado.edu/research/wavelets/"
# The Matlab code written January 1998 by C. Torrence is modified to Python by Evgeniya Predybaylo, December 2014
#
# Modified Oct 1999, changed Global Wavelet Spectrum (GWS) to be sideways,
#   changed all "log" to "log2", changed logarithmic axis on GWS to
#   a normal axis.
# ------------------------------------------------------------------------------------------------------------------

## READ THE DATA
#sst = np.loadtxt('sst_nino3.dat')  # input SST time series
#sst = sst - np.mean(sst)
#variance = np.std(sst, ddof=1) ** 2
##print("variance = ", variance)

#----------C-O-M-P-U-T-A-T-I-O-N------S-T-A-R-T-S------H-E-R-E------------------------------------------------------

# normalize by standard deviation (not necessary, but makes it easier
# to compare with plot on Interactive Wavelet page, at
# "http://paos.colorado.edu/research/wavelets/plot/"
def wave_analy(sst,time,dt,mother):
  variance = np.std(sst, ddof=1) ** 2
  if 0:
      variance = 1.0
      sst = sst / np.std(sst, ddof=1)
  n = len(sst)
  #dt = 0.25
  #time = np.arange(len(sst)) * dt + 1871.0  # construct time array
  #xlim = ([1870, 2000])  # plotting range
  pad = 1  # pad the time series with zeroes (recommended)
  dj = 0.25  # this will do 4 sub-octaves per octave
  s0 = 2 * dt  # this says start at a scale of 6 months
  j1 = 4 / dj  # this says do 7 powers-of-two with dj sub-octaves each
  alpha = np.corrcoef(sst[0:-1], sst[1:])[0,1]; 
  lag1 = alpha
  #lag1 = 0.72  # lag-1 autocorrelation for red noise background
  #print("lag1 = ", lag1)
  #mother = 'DOG'

  # Wavelet transform: #, J1=j1
  wave, period, scale, coi = wavelet(Y=sst, dt=dt, pad=pad, dj=dj, s0=s0, mother=mother)
  power = (np.abs(wave)) ** 2  # compute wavelet power spectrum
  global_ws = (np.sum(power, axis=1) / n)  # time-average over all times

  # Significance levels:
  signif = wave_signif(([variance]), dt=dt, sigtest=0, scale=scale,
      lag1=lag1, mother=mother)
  sig95 = signif[:, np.newaxis].dot(np.ones(n)[np.newaxis, :])  # expand signif --> (J+1)x(N) array
  sig95 = power / sig95  # where ratio > 1, power is significant

  # Global wavelet spectrum & significance levels:
  dof = n - scale  # the -scale corrects for padding at edges
  global_signif = wave_signif(variance, dt=dt, scale=scale, sigtest=1,
      lag1=lag1, dof=dof, mother=mother)

  ## Scale-average between El Nino periods of 2--8 years
  #avg = np.logical_and(scale >= 2, scale < 8)
  #Cdelta = 0.776  # this is for the MORLET wavelet
  #scale_avg = scale[:, np.newaxis].dot(np.ones(n)[np.newaxis, :])  # expand scale --> (J+1)x(N) array
  #scale_avg = power / scale_avg  # [Eqn(24)]
  #scale_avg = dj * dt / Cdelta * sum(scale_avg[avg, :])  # [Eqn(24)]
  #scaleavg_signif = wave_signif(variance, dt=dt, scale=scale, sigtest=2,
      #lag1=lag1, dof=([2, 7.9]), mother=mother)
      
  #print('MOTHER ', mother)
  #wave_plot(time,sst,power,period, sig95, coi, global_ws, global_signif,' ',plotting = True)
  return power, period, sig95, coi, global_ws, global_signif

#------------------------------------------------------ Plotting

def wave_plot(time,sst,power,period, sig95, coi, global_ws, global_signif,title,plotting = False):
  if plotting:
    #--- Plot time series
    fig = plt.figure(figsize=(9, 10))
    gs = GridSpec(2, 4, hspace=0.4, wspace=0.75)
    plt.subplots_adjust(left=0.1, bottom=0.05, right=0.9, top=0.95, wspace=0, hspace=0)
    plt.subplot(gs[0, 0:4])
    plt.plot(time, sst, 'k')
    #plt.xlim(xlim[:])
    plt.xlabel('Time (s)')
    plt.ylabel('Amplitude (cm/s)')
    plt.title(title)

    #plt.text(time[-1] + 35, 0.5,'Wavelet Analysis\nC. Torrence & G.P. Compo\n' +
	#'http://paos.colorado.edu/\nresearch/wavelets/',
	#horizontalalignment='center', verticalalignment='center')

    #--- Contour plot wavelet power spectrum
    # plt3 = plt.subplot(3, 1, 2)
    plt3 = plt.subplot(gs[1, 0:4])
    levels = [0, 0.5, 1, 2, 4, 999]
    CS = plt.contourf(time, period, power, len(levels))  #*** or use 'contour'
    im = plt.contourf(CS, levels=levels, colors=['white','bisque','orange','orangered','darkred'])
    plt.xlabel('Time (s)')
    plt.ylabel('Wavelet Power Spectrum')
    #plt.title(title)
    
    #plt.xlim(xlim[:])
    # 95# significance contour, levels at -99 (fake) and 1 (95# signif)
    #plt.contour(time, period, sig95, [-99, 1], colors='k')
    # cone-of-influence, anything "below" is dubious
    #plt.plot(time, coi, 'k')
    # format y-scale
    plt3.set_yscale('log', basey=2, subsy=None)
    plt.ylim([np.min(period), np.max(period)])
    ax = plt.gca().yaxis
    ax.set_major_formatter(ticker.ScalarFormatter())
    plt3.ticklabel_format(axis='y', style='plain')
    plt3.invert_yaxis()
    # set up the size and location of the colorbar
    # position=fig.add_axes([0.5,0.36,0.2,0.01]) 
    # plt.colorbar(im, cax=position, orientation='horizontal') #, fraction=0.05, pad=0.5)

    # plt.subplots_adjust(right=0.7, top=0.9)

    #--- Plot global wavelet spectrum
    #plt4 = plt.subplot(gs[1, -1])
    #plt.plot(global_ws, period)
    #plt.plot(global_signif, period, '--')
    #plt.xlabel('Power')
    #plt.title('Global Wavelet Spectrum')
    ##plt.xlim([0, 1.25 * np.max(global_ws)])
    ## format y-scale
    #plt4.set_yscale('log', basey=2, subsy=None)
    #plt.ylim([np.min(period), np.max(period)])
    #ax = plt.gca().yaxis
    #ax.set_major_formatter(ticker.ScalarFormatter())
    #plt4.ticklabel_format(axis='y', style='plain')
    #plt4.invert_yaxis()

    # --- Plot 2--8 yr scale-average time series
    #plt.subplot(gs[2, 0:3])
    #plt.plot(time, scale_avg, 'k')
    ##plt.xlim(xlim[:])
    #plt.xlabel('Time (s)')
    #plt.ylabel('Avg variance')
    #plt.title('d) 2-8 yr Scale-average Time Series')
    #plt.plot(xlim, scaleavg_signif + [0, 0], '--')
    
    #fig.savefig('/home/dertuncay/Wavelet/' + title + '.eps', dpi=300)
    #plt.show()
    #fig.clf()
    #plt.close('all')
# end of code

