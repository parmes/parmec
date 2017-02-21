# PARMEC test --> acc_sweep module
import sys
sys.path.append('python')
from acc_sweep import *

step = 1E-4  # time step
stop = 5.0   # duration of the simulation
lofq = 1     # low frequency for the sweep
hifq = 10    # high frequency for the sweep
amag = 10.0  # acceleration magnitude

(vt, vd, vv, va) = acc_sweep (step, stop, lofq, hifq, amag,
                              acc_plot = 'tests/asweep_acc.png',
			      vel_plot = 'tests/asweep_vel.png',
			      dsp_plot = 'tests/asweep_dsp.png',
			      dsp_envelope = 'tests/asweep_dsp_envelope.png')
