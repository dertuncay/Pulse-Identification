import numpy as np
from obspy import read
from scipy.signal import hilbert
import matplotlib.pyplot as plt
from waveletAnalysis import * 
from pulse_identifier import wavelet, shift, ricker, pgv_pos_finder

def corrcoef(rec1,rec2):
  return np.corrcoef(rec1,rec2)


def ricker_with_phase(f, dt,data, phase):
  # PGV pos
  _, PGV_pos = pgv_pos_finder(data)

  # Create Ricker Wavelet
  t, y = ricker(f, dt,data)
  
  # Hilbert Transformation
  analytic_signal = hilbert(y)

  # Imaginary Part
  im = np.imag(analytic_signal)
  # Real Part
  re = np.real(analytic_signal)
  
  ricker_phased = np.cos(phase)*re - np.sin(phase)*im;
  
  # Find Time Shift
  af = scipy.fft(data)
  bf = scipy.fft(ricker_phased)

  c = scipy.ifft(af[:len(af)//2] * scipy.conj(bf[:len(bf)//2]))
  time_shift = np.argmax(abs(c))
  
  # Shift the wavelet
  shifted_wavelet = shift(ricker_phased, time_shift, fill_value=np.nan)
  shifted_wavelet = shifted_wavelet[~np.isnan(shifted_wavelet)]
  
  _, shifted_ricker_PGV_pos = pgv_pos_finder(shifted_wavelet)
  shifted_wavelet = shift(shifted_wavelet, PGV_pos-shifted_ricker_PGV_pos, fill_value=np.nan)
  # Add zeros to the end
  shifted_wavelet = np.pad(shifted_wavelet, (0, time_shift), 'constant')
  shifted_wavelet = np.nan_to_num(shifted_wavelet)
  
  return shifted_wavelet

st = read('NGA.RSN1482..HGR.SAC',format='SAC')
st = read(sac,format='SAC')
st[0].filter('bandpass',freqmin = 0.05, freqmax = 10)
st[0].integrate()
scale = scale_finder(st)
data = st[0].data*scale
dt = st[0].stats.delta
times = st[0].times()

  # Energy drop
  en_drop = []
  
  # Phase angles
  phase_angles = np.arange(0,356,5)
  for angle in phase_angles:
    # Phase Shift to Ricker Wavelet
    phase = np.deg2rad(angle)
    
    # Apply phase angle to ricker wavelet
    shifted_wavelet = ricker_with_phase(1/Tp, dt, data, phase)
    
    # Correlation Coefficient
    corr_coef = corrcoef(data,shifted_wavelet)
    #abs(corr_coef[0][1])
    en_drop.append(corr_coef[0][1])
    
    max_corr = np.where(en_drop == np.array(en_drop).max())[0][0]
    phased_ricker = ricker_with_phase(1/Tp, dt,data, np.deg2rad(phase_angles[max_corr]))
    
    
  plt.plot(times,data,linewidth=1,color='k')
  plt.plot(times,phased_ricker,linewidth=2,color='r',label=str(phase_angles[max_corr]))
  plt.xlabel('Time (s)')
  plt.ylabel('Amplitude (cm/s)')
  plt.xlim(0,max(times))
  plt.show()
