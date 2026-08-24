"""Fast ML likelihood-ratio tests for between-person slope variance."""
import numpy as np
from scipy import stats
from scipy.linalg import cho_factor, cho_solve
from scipy.optimize import minimize
_EPS=1e-8

def _neg_ll(params,Y,X,random_slope):
    N,n=Y.shape; sig2=np.exp(2*params[0])
    if random_slope:
        a,b,c=params[1],params[2],params[3]; L=np.array([[a,0.],[b,c]]); G=L@L.T
    else:
        a=params[1]; G=np.array([[a*a,0.],[0.,0.]])
    V=X@G@X.T+sig2*np.eye(n)
    try:
        cf=cho_factor(V,lower=True,check_finite=False); logdetV=2*np.sum(np.log(np.diag(cf[0])))
        ViX=cho_solve(cf,X,check_finite=False); A=X.T@ViX; ybar=Y.mean(axis=0)
        beta=np.linalg.solve(A,ViX.T@ybar); R=Y-(X@beta)[None,:]; ViR=cho_solve(cf,R.T,check_finite=False); quad=float(np.sum(R.T*ViR))
    except Exception: return 1e12
    if not np.isfinite(logdetV) or not np.isfinite(quad): return 1e12
    return .5*(N*n*np.log(2*np.pi)+N*logdetV+quad)

def _mom_start(Y,X):
    N,n=Y.shape; XtXi=np.linalg.pinv(X.T@X); B=Y@(XtXi@X.T).T; R=Y-B@X.T; dof=max(n-2,1)
    s2=float((R**2).sum()/(N*dof)); se2=s2*np.diag(XtXi)
    g00=max(B[:,0].var(ddof=1)-se2[0],_EPS); g11=max(B[:,1].var(ddof=1)-se2[1],_EPS)
    return np.sqrt(max(s2,_EPS)),np.sqrt(g00),np.sqrt(g11)

def lrt_random_slope_fast(Y,x):
    Y=np.asarray(Y,float); x=np.asarray(x,float); n=len(x); X=np.column_stack([np.ones(n),x]); sd_e,sd0,sd1=_mom_start(Y,X); ls=np.log(max(sd_e,1e-4))
    try:
        r0=minimize(_neg_ll,np.array([ls,sd0]),args=(Y,X,False),method='L-BFGS-B',options=dict(maxiter=500)); best1=None
        for st in ([ls,sd0,0.,sd1],[ls,sd0,0.,sd1*3+1e-3],[ls,sd0*.5,0.,1e-3]):
            r=minimize(_neg_ll,np.array(st),args=(Y,X,True),method='L-BFGS-B',options=dict(maxiter=500)); best1=r if best1 is None or r.fun<best1.fun else best1
    except Exception: return np.nan
    if r0.fun>=1e11 or best1.fun>=1e11: return np.nan
    stat=max(2*(r0.fun-best1.fun),0.); return float(.5*stats.chi2.sf(stat,1)+.5*stats.chi2.sf(stat,2)) if np.isfinite(stat) else np.nan

def _neg_ll_varying(params,Y,Xpred,random_slope):
    Y=np.asarray(Y,float); Xpred=np.asarray(Xpred,float); N,n=Y.shape; sig2=np.exp(2*params[0])
    if random_slope:
        a,b,c=params[1],params[2],params[3]; L=np.array([[a,0.],[b,c]]); G=L@L.T
    else:
        a=params[1]; G=np.array([[a*a,0.],[0.,0.]])
    A=np.zeros((2,2)); rhs=np.zeros(2); cached=[]; logdet=0.
    try:
        for i in range(N):
            Xi=np.column_stack([np.ones(n),Xpred[i]]); Vi=Xi@G@Xi.T+sig2*np.eye(n); cf=cho_factor(Vi,lower=True,check_finite=False)
            logdet+=2*np.sum(np.log(np.diag(cf[0]))); ViXi=cho_solve(cf,Xi,check_finite=False); Viy=cho_solve(cf,Y[i],check_finite=False)
            A+=Xi.T@ViXi; rhs+=Xi.T@Viy; cached.append((cf,Xi))
        beta=np.linalg.solve(A,rhs); quad=0.
        for i,(cf,Xi) in enumerate(cached):
            ri=Y[i]-Xi@beta; quad+=float(ri@cho_solve(cf,ri,check_finite=False))
    except Exception: return 1e12
    if not np.isfinite(logdet) or not np.isfinite(quad): return 1e12
    return .5*(N*n*np.log(2*np.pi)+logdet+quad)

def _mom_start_varying(Y,Xpred):
    coefs=[]; rss=0.; dof=0
    for i in range(Y.shape[0]):
        Xi=np.column_stack([np.ones(Y.shape[1]),Xpred[i]]); bi=np.linalg.pinv(Xi)@Y[i]; ri=Y[i]-Xi@bi
        coefs.append(bi); rss+=float(ri@ri); dof+=max(len(ri)-2,1)
    B=np.asarray(coefs); s2=max(rss/max(dof,1),_EPS)
    return np.sqrt(s2),np.sqrt(max(B[:,0].var(ddof=1),_EPS)),np.sqrt(max(B[:,1].var(ddof=1),_EPS))

def lrt_random_slope_varying_x(Y,Xpred):
    """Random-slope LRT when each participant has her own predictor trajectory."""
    Y=np.asarray(Y,float); Xpred=np.asarray(Xpred,float)
    if Y.shape!=Xpred.shape or Y.ndim!=2: raise ValueError('Y and Xpred must have the same 2D shape')
    sd_e,sd0,sd1=_mom_start_varying(Y,Xpred); ls=np.log(max(sd_e,1e-4))
    try:
        r0=minimize(_neg_ll_varying,np.array([ls,sd0]),args=(Y,Xpred,False),method='L-BFGS-B',options=dict(maxiter=500)); best1=None
        for st in ([ls,sd0,0.,sd1],[ls,sd0,0.,sd1*3+1e-3],[ls,sd0*.5,0.,1e-3]):
            r=minimize(_neg_ll_varying,np.array(st),args=(Y,Xpred,True),method='L-BFGS-B',options=dict(maxiter=500)); best1=r if best1 is None or r.fun<best1.fun else best1
    except Exception: return np.nan
    if r0.fun>=1e11 or best1.fun>=1e11: return np.nan
    stat=max(2*(r0.fun-best1.fun),0.); return float(.5*stats.chi2.sf(stat,1)+.5*stats.chi2.sf(stat,2)) if np.isfinite(stat) else np.nan
