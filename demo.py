from obspy import read
from pulse_identifier import *

st = read('NGA.RSN1482..HGT.SAC',format = 'SAC')
p_e, spec_p_e, Tp, is_pulse_dz, late, wavelet_type = vel_wf_detector(st, scale = 1, zero_adding = False)
print(st[0].stats)
print(p_e, spec_p_e, Tp, is_pulse_dz, late, wavelet_type)

st = read('NGA.RSN1482..HGR.SAC',format = 'SAC')
p_e, spec_p_e, Tp, is_pulse_dz, late, wavelet_type = vel_wf_detector(st, scale = 1, zero_adding = False)
print(st[0].stats)
print(p_e, spec_p_e, Tp, is_pulse_dz, late, wavelet_type)
