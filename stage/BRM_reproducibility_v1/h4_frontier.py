# Local execution copy of the canonical reference-model generative constants/functions.
import numpy as np
CYCLE_LEN=28.; DA_OPT=.50; U_WIDTH=.35; K_GAIN=.15; SIGMA_STATE=.085; N_SUBJ=39
_grid=np.arange(0,CYCLE_LEN,.25)
def _e2_raw(d): return .15+np.exp(-((d-13.)**2)/(2*2.**2))+.55*np.exp(-((d-21.)**2)/(2*3.5**2))
_E2_MAX=_e2_raw(_grid).max(); _E2_MEAN=(_e2_raw(_grid)/_E2_MAX).mean()
def e2(days): return _e2_raw(np.asarray(days,float)%CYCLE_LEN)/_E2_MAX
def inverted_u(da): return np.exp(-((da-DA_OPT)**2)/(2*U_WIDTH**2))
