import numpy as np
CYCLE_LEN=28.0; DA_OPT=.50; U_WIDTH=.35; K_GAIN=.15; SIGMA_STATE=.085; N_SUBJ=39
_grid=np.arange(0,CYCLE_LEN,.25)
def _e2_raw(d, peak_shift=0.0):
    d=np.asarray(d,float)
    return .15+np.exp(-((d-(13.+peak_shift))**2)/(2*2.**2))+.55*np.exp(-((d-21.)**2)/(2*3.5**2))
_E2_MAX=_e2_raw(_grid).max(); _E2_MEAN=(_e2_raw(_grid)/_E2_MAX).mean()
def e2(days, peak_shift=0.0): return _e2_raw(np.asarray(days,float)%CYCLE_LEN,peak_shift)/_E2_MAX
def inverted_u(da): return np.exp(-((da-DA_OPT)**2)/(2*U_WIDTH**2))
def within_r(Y,x):
    Y=np.asarray(Y,float); x=np.asarray(x,float)
    xc=x-np.nanmean(x)
    if Y.ndim==1: Y=Y[None,:]
    out=[]
    for y in Y:
        m=np.isfinite(y)&np.isfinite(x)
        if m.sum()<4 or np.nanstd(y[m])<1e-12 or np.nanstd(x[m])<1e-12: out.append(np.nan); continue
        out.append(float(np.corrcoef(y[m],x[m])[0,1]))
    return np.array(out)
def true_coupling(b, K=None, form='inverted_u'):
    xd=e2(np.linspace(0,CYCLE_LEN,400,endpoint=False))
    b=np.asarray(b,float); K=np.full_like(b,K_GAIN) if K is None else np.asarray(K,float)
    da=b[:,None]+K[:,None]*xd[None,:]
    if form=='inverted_u': sig=inverted_u(da)
    elif form=='sigmoid': sig=1/(1+np.exp(-10*(da-DA_OPT)))
    elif form=='asym_linear':
        z=da-DA_OPT; sig=np.where(z<0,1.1*z,-0.65*z)
    elif form=='piecewise':
        z=da-DA_OPT; sig=np.where(z<-0.08,0.7*z,np.where(z<0.10,1.0*z,-0.5*z+0.15))
    else: raise ValueError(form)
    return within_r(sig,xd)
