import scipy.optimize
import numpy as np
from obspy import read, Stream, Trace
from waveletAnalysis import *  

def wavelet(data,times,dt,PGV_pos,method):
  ''' 
  Use The Matlab code written January 1998 by C. Torrence is modified to Python by Evgeniya Predybaylo, December 2014
  '''
  power,period, sig95, coi, global_ws, global_signif = wave_analy(data,times,dt,method)
    
  # PGV Energy
  tp_pos = np.where(power[:,PGV_pos] == np.max(power[:,PGV_pos]))[0][0]
  tp = period[tp_pos]
  
  # Maximum energy region of the waveform
  powers = power.sum(axis=0)
  max_power_pos = np.where(powers == max(powers))[0][0]
  max_power_tp_pos = np.where(power[:,max_power_pos] == np.max(power[:,max_power_pos]))[0][0]
  max_power_tp = period[max_power_tp_pos]
  
  return tp,max_power_tp,max_power_pos,power,period, sig95, coi, global_ws, global_signif

def morlet(f,dt,data,order = 4):
  '''
  This function will produce morlet wavelet with a given order
  '''
  time = np.linspace(-len(data)/2, (len(data)-dt)/2, int(len(data)/dt))
  # Sinus Wave
  sine_wave = np.exp(2*np.pi*complex(1j)*f*time) 
  s = order/(2*np.pi*f)
  # Gaussian Wave
  gaussian = np.exp(-1*(time)**2/(2*(s**2)))
  # Morlet Wavelet
  wavelet = sine_wave*gaussian
    
  # Shifting Morlet wavelet to match PGV positions
  PGV, PGV_pos = pgv_pos_finder(data)
  _, morlet_PGV_pos = pgv_pos_finder(wavelet)
  wavelet = shift(wavelet, PGV_pos-morlet_PGV_pos, fill_value=np.nan)
  wavelet = wavelet[~np.isnan(wavelet)]
  wavelet = wavelet*PGV 
  if morlet_PGV_pos > PGV_pos:
    wavelet = np.append(wavelet,np.zeros(morlet_PGV_pos - PGV_pos))
    if len(data) != len(wavelet):
      wavelet = wavelet[0:len(data)]
  elif morlet_PGV_pos < PGV_pos:
    wavelet = np.append(np.zeros(abs(morlet_PGV_pos - PGV_pos)),wavelet)
    if len(data) != len(wavelet):
      wavelet = wavelet[0:len(data)]
  return wavelet

def morlet_maxpos(f, dt,data,max_pos,EMAX_GV,order = 4):
  '''
  This function will produce morlet wavelet with a given order
  '''
  time = np.linspace(-len(data)/2, (len(data)-dt)/2, int(len(data)/dt))
  # Sinus Wave
  sine_wave = np.exp(2*np.pi*complex(1j)*f*time) 
  s = order/(2*np.pi*f)
  # Gaussian Wave
  gaussian = np.exp(-1*(time)**2/(2*(s**2)))
  # Morlet Wavelet
  wavelet = sine_wave*gaussian
  if np.sign(data[max_pos]) != np.sign(EMAX_GV):
    wavelet = wavelet*EMAX_GV*-1
  else:
    wavelet = wavelet*EMAX_GV
    
  # Shifting Ricker wavelet to match PGV positions
  _, morlet_PGV_pos = pgv_pos_finder(wavelet)
  wavelet = shift(wavelet, max_pos-morlet_PGV_pos, fill_value=np.nan)
  wavelet = wavelet[~np.isnan(wavelet)]
  wavelet = wavelet*data[max_pos] 
  if morlet_PGV_pos > max_pos:
    wavelet = np.append(wavelet,np.zeros(morlet_PGV_pos - max_pos))
    if len(data) != len(wavelet):
      wavelet = wavelet[0:len(data)]
  elif morlet_PGV_pos < max_pos:
    wavelet = np.append(np.zeros(abs(morlet_PGV_pos - max_pos)),wavelet)
    if len(data) != len(wavelet):
      wavelet = wavelet[0:len(data)]
  return wavelet

def shift(arr, num, fill_value=np.nan):
  '''
  Shifts waveform. Adds NaN in both ends if necessery.
  '''
  result = np.empty_like(arr)
  if num > 0:
      result[:num] = fill_value
      result[num:] = arr[:-num]
  elif num < 0:
      result[num:] = fill_value
      result[:num] = arr[-num:]
  else:
      result = arr
  return result

def ricker(f, dt,data):
  ''' ricker wavelet '''
  t = np.linspace(-len(data)/2, (len(data)-dt)/2, int(len(data)/dt))
  y = (1.-2.*(np.pi**2)*(f**2)*(t**2))*np.exp(-(np.pi**2)*(f**2)*(t**2))
  # Shifting Ricker wavelet to match PGV positions
  PGV, PGV_pos = pgv_pos_finder(data)
  _, ricker_PGV_pos = pgv_pos_finder(y)
  y = shift(y, PGV_pos-ricker_PGV_pos, fill_value=np.nan)
  y = y[~np.isnan(y)]
  y = y*PGV 
  if ricker_PGV_pos > PGV_pos:
    y = np.append(y,np.zeros(ricker_PGV_pos - PGV_pos))
    if len(data) != len(y):
      y = y[0:len(data)]
  elif ricker_PGV_pos < PGV_pos:
    y = np.append(np.zeros(abs(ricker_PGV_pos - PGV_pos)),y)
    if len(data) != len(y):
      y = y[0:len(data)]
  return t, y

def ricker_maxpos(f, dt,data,max_pos,EMAX_GV):
  ''' ricker wavelet '''
  t = np.linspace(-len(data)/2, (len(data)-dt)/2, int(len(data)/dt))
  y = (1.-2.*(np.pi**2)*(f**2)*(t**2))*np.exp(-(np.pi**2)*(f**2)*(t**2))
  # Shifting Ricker wavelet to match PGV positions
  _, ricker_PGV_pos = pgv_pos_finder(y)
  y = shift(y, max_pos-ricker_PGV_pos, fill_value=np.nan)
  y = y[~np.isnan(y)]
  #y = y*data[max_pos]
  if np.sign(data[max_pos]) != np.sign(EMAX_GV):
    y = y*EMAX_GV*-1
  else:
    y = y*EMAX_GV
    
  if ricker_PGV_pos > max_pos:
    y = np.append(y,np.zeros(ricker_PGV_pos - max_pos))
    if len(data) != len(y):
      y = y[0:len(data)]
  elif ricker_PGV_pos < max_pos:
    y = np.append(np.zeros(abs(ricker_PGV_pos - max_pos)),y)
    if len(data) != len(y):
      y = y[0:len(data)]  
  return t, y

def pgv_pos_finder(data):
  '''
  PGV position finder'''
  if abs(data.min()) > data.max():
    PGV = max(data.min(), data.max(), key=abs)
    PGV_pos = np.where(data == PGV)
    PGV_pos= PGV_pos[0]
  else:
      PGV = data.max()
      PGV_pos = np.where(data == PGV)
      PGV_pos = PGV_pos[0]     
  return PGV, PGV_pos[0]

def analysis(data,times,dt,PGV,PGV_pos,method,order = -1):
  '''
  Analyze data with wavelet analysis for a given method and order if necessery
  '''
  
  if order != -1:
    order = order
  
  # Total Waveform Energy
  total_e = np.sum(data ** 2)
  
  if method == 'Ricker':
    # RICKER ANALYSIS
    # Wavelet Analysis
    pgv_tp, max_tp,max_tp_pos, power,period, sig95, coi, global_ws, global_signif = wavelet(data,times,dt,PGV_pos,'DOG')
    
    EMAX_GV = max(data[max_tp_pos-int(max_tp/2*1/dt):max_tp_pos+int(max_tp/2*1/dt)], key=abs)
      
    # PGV Pulse Waveform Energy
    pgv_e = np.sum(data[PGV_pos - int((pgv_tp/2)*(1/dt)):PGV_pos + int((pgv_tp/2)*(1/dt))] ** 2)
    # Maximum Period Waveform Energy
    max_e = np.sum(data[max_tp_pos - int((max_tp/2)*(1/dt)):max_tp_pos + int((max_tp/2)*(1/dt))] ** 2)
    
    # PGV Pulse Spec Energy
    pgv_spec_e = np.sum(power[:,PGV_pos - int((pgv_tp/2)*(1/dt)):PGV_pos + int((pgv_tp/2)*(1/dt))])
    # Maximum Period Spec Energy
    max_spec_e = np.sum(power[:,max_tp_pos - int((max_tp/2)*(1/dt)):max_tp_pos + int((max_tp/2)*(1/dt))])
    
    # Total Spec Energy
    spec_e = np.sum(power)


    # Determination of Pulse type
    if ((max_e/total_e)+(max_spec_e/spec_e))/2 >= 0.30 and max_e >= pgv_e*1.1 and max_spec_e >= pgv_spec_e*1.1 and abs(EMAX_GV) >= 25 and abs(PGV_pos - max_tp_pos)*dt > pgv_tp/4:
      late = 1
      is_pulse_dz = 1
      Tp = max_tp
      p_e = (max_e/total_e)*100
      spec_p_e = (max_spec_e/spec_e)*100
      _, data_ricker = ricker_maxpos(1/Tp, dt , data,max_tp_pos,EMAX_GV)
      # Energy Change
      e_change = ((total_e-np.sum(data_ricker ** 2))/total_e)
    elif ((pgv_e/total_e) + (pgv_spec_e/spec_e))/2 >= 0.30 and abs(PGV) >= 30:
      is_pulse_dz = 1
      Tp = pgv_tp
      p_e = (pgv_e/total_e)*100
      spec_p_e = (pgv_spec_e/spec_e)*100
      _, data_ricker = ricker(1/Tp, dt , data)
      # Energy Change
      e_change = ((total_e-np.sum(data_ricker ** 2))/total_e)
      late = 0
    else:
      is_pulse_dz = 0
      Tp = max(pgv_tp,max_tp)
      p_e = (pgv_e/total_e)*100
      spec_p_e = (pgv_spec_e/spec_e)*100
      late = float('NaN')
      e_change = 0
      data_ricker = []; power = [];period = []; sig95 = []; coi = []; global_ws = []; global_signif = [];
  elif method == 'Morlet':
    # MORLET ANALYSIS
    # Wavelet Analysis
    pgv_tp, max_tp,max_tp_pos, power,period, sig95, coi, global_ws, global_signif = wavelet(data,times,dt,PGV_pos,'MORLET')
    
    EMAX_GV = max(data[max_tp_pos-int(max_tp/2*1/dt):max_tp_pos+int(max_tp/2*1/dt)], key=abs)
    
    # PGV Pulse Waveform Energy
    pgv_e = np.sum(data[PGV_pos - int((pgv_tp/2)*(1/dt)):PGV_pos + int((pgv_tp/2)*(1/dt))] ** 2)
    # Maximum Period Waveform Energy
    max_e = np.sum(data[max_tp_pos - int((max_tp/2)*(1/dt)):max_tp_pos + int((max_tp/2)*(1/dt))] ** 2)
    
    # PGV Pulse Spec Energy
    pgv_spec_e = np.sum(power[:,PGV_pos - int((pgv_tp/2)*(1/dt)):PGV_pos + int((pgv_tp/2)*(1/dt))])
    # Maximum Period Spec Energy
    max_spec_e = np.sum(power[:,max_tp_pos - int((max_tp/2)*(1/dt)):max_tp_pos + int((max_tp/2)*(1/dt))])
    # Total Spec Energy
    spec_e = np.sum(power)
        
    # Determination of Pulse type
    if ((max_e/total_e)+(max_spec_e/spec_e))/2 >= 0.30 and max_e >= pgv_e*1.1 and max_spec_e >= pgv_spec_e*1.1 and abs(EMAX_GV) >= 25 and abs(PGV_pos - max_tp_pos)*dt > pgv_tp/4:
      late = 1
      is_pulse_dz = 1
      Tp = max_tp
      p_e = (max_e/total_e)*100
      spec_p_e = (max_spec_e/spec_e)*100
      data_ricker = morlet_maxpos(1/Tp, dt ,data,max_tp_pos,EMAX_GV,order)
      # Energy Change
      e_change = ((total_e-np.sum(np.real(data_ricker) ** 2))/total_e)
    elif ((pgv_e/total_e) + (pgv_spec_e/spec_e))/2 >= 0.30 and abs(PGV) >= 30:
      is_pulse_dz = 1
      Tp = pgv_tp
      p_e = (pgv_e/total_e)*100
      spec_p_e = (pgv_spec_e/spec_e)*100
      data_ricker = morlet(1/Tp, dt , data,order)
      # Energy Change
      e_change = ((total_e-np.sum(np.real(data_ricker) ** 2))/total_e)
      late = 0
    else:
      is_pulse_dz = 0
      Tp = max(pgv_tp,max_tp)
      p_e = (pgv_e/total_e)*100
      spec_p_e = (pgv_spec_e/spec_e)*100
      late = float('NaN')
      e_change = float('NaN')
      data_ricker = []; power = [];period = []; sig95 = []; coi = []; global_ws = []; global_signif = [];      
  return is_pulse_dz, Tp, p_e, spec_p_e, late, e_change, data_ricker, power,period, sig95, coi, global_ws, global_signif

def vel_wf_detector(st, scale = 1, zero_adding = False):
  ''' Main Code that does all the work '''
  # Adding points at the beginning and at the end to overcome problems in the calculation of pulse waveform areas.
  if zero_adding == False:
    data = st[0].data*scale
    times = st[0].times()
  elif zero_adding == True:
    random_add = np.random.uniform(low=-0.0000000000000001, high=0.0000000000000001, size=(100,))
    data = np.concatenate([random_add, st[0].data*scale, random_add])
    times = np.arange(0,len(data)*st[0].stats.delta,st[0].stats.delta)
    if len(data) > len(times):
      data = data[:-1*(len(data)-len(times))]
    elif len(times) > len(data):
      times = times[:-1*(len(times)-len(data))]
  dt = 1/st[0].stats.sampling_rate
  
  # PGV Position
  PGV, PGV_pos = pgv_pos_finder(data)

  if abs(PGV) >= 30:
    r1 = []; r2 = []; r3 = []; r4 = []; r5 = []; r6 = []; r7 = []; r8 = []; r9 = []; r10 = []; r11 = []; r12 = []; r13 = [];
    is_pulse_dz, Tp, p_e, spec_p_e, late, e_change, data_ricker, power,period, sig95, coi, global_ws, global_signif = analysis(data,times,dt,PGV,PGV_pos,method = 'Ricker')
    r1.append(is_pulse_dz);r2.append(Tp);r3.append(p_e);r4.append(spec_p_e);r5.append(late);r6.append(e_change);
    r7.append(data_ricker);r8.append(power);r9.append(period);r10.append(sig95);r11.append(coi);r12.append(global_ws);r13.append(global_signif)
    
    is_pulse_dz, Tp, p_e, spec_p_e, late, e_change, data_ricker, power,period, sig95, coi, global_ws, global_signif = analysis(data,times,dt,PGV,PGV_pos,method = 'Morlet',order = 3)
    r1.append(is_pulse_dz);r2.append(Tp);r3.append(p_e);r4.append(spec_p_e);r5.append(late);r6.append(e_change);
    r7.append(data_ricker);r8.append(power);r9.append(period);r10.append(sig95);r11.append(coi);r12.append(global_ws);r13.append(global_signif)
    
    is_pulse_dz, Tp, p_e, spec_p_e, late, e_change, data_ricker, power,period, sig95, coi, global_ws, global_signif = analysis(data,times,dt,PGV,PGV_pos,method = 'Morlet',order = 4)
    r1.append(is_pulse_dz);r2.append(Tp);r3.append(p_e);r4.append(spec_p_e);r5.append(late);r6.append(e_change);
    r7.append(data_ricker);r8.append(power);r9.append(period);r10.append(sig95);r11.append(coi);r12.append(global_ws);r13.append(global_signif)
    # Finding the most proper wavelet type and getting its info
    amax = np.argmax(np.asarray(r6))
    p_e = r3[amax]
    spec_p_e = r4[amax]
    Tp = r2[amax]
    is_pulse_dz = r1[amax]
    late = r5[amax]
    fit_wavelet = r7[amax]
    if is_pulse_dz == 1:
      wavelet_type = amax
    elif is_pulse_dz == 0:
     wavelet_type = 0
  else:
#     _, _, p_e, spec_p_e, _, _, _, _,_, _, _, _, _ = analysis(data,times,dt,PGV,PGV_pos,method = 'Ricker')
    is_pulse_dz = 0
    p_e = float('NaN')
    spec_p_e = float('NaN')
    Tp = 0
    late = float('NaN')
    wavelet_type = float('NaN')
    fit_wavelet = False
  return p_e, spec_p_e, Tp, is_pulse_dz, late, wavelet_type, fit_wavelet
