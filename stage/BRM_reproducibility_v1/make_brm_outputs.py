from pathlib import Path
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
import h4_frontier as H

import argparse
ap=argparse.ArgumentParser()
ap.add_argument('--root', default='.')
ap.add_argument('--results-dir', default='results')
ap.add_argument('--figures-dir', default='figures')
ap.add_argument('--tables-dir', default='tables')
a=ap.parse_args()
ROOT=Path(a.root).resolve()
RES=ROOT/a.results_dir; FIG=ROOT/a.figures_dir; TAB=ROOT/a.tables_dir
FIG.mkdir(exist_ok=True); TAB.mkdir(exist_ok=True)
rob=pd.read_csv(RES/'brm_robustness_metrics.csv')
samp=pd.read_csv(RES/'brm_sampling_design.csv')
ssf=pd.read_csv(RES/'brm_ssf_benchmark.csv')
sur=pd.read_csv(RES/'brm_surrogate_calibration_power.csv')

# ---------- Tables ----------
ref=rob[(rob.factor=='reference')&(rob.level=='reference')].copy()
ref.to_csv(TAB/'table1_reference_recovery.csv',index=False)

rows=[]
for (fac,lev),g in rob[rob.n_cycles==2].groupby(['factor','level']):
    g=g.sort_values('n_obs'); hit=g[g.fidelity>=.70]
    get=lambda n,c: float(g.loc[g.n_obs==n,c].iloc[0])
    rows.append(dict(factor=fac,level=lev,min_n_fidelity_070=(int(hit.n_obs.iloc[0]) if len(hit) else np.nan),
                     fidelity_n20=get(20,'fidelity'),fidelity_n56=get(56,'fidelity'),
                     rmse_n20=get(20,'rmse'),sign_accuracy_n20=get(20,'sign_accuracy'),bias_n20=get(20,'bias')))
t2=pd.DataFrame(rows).sort_values(['factor','level']);t2.to_csv(TAB/'table2_robustness_summary.csv',index=False)

ssfagg=ssf.groupby('method').agg(mean_abs_bias=('bias',lambda x:float(np.mean(np.abs(x)))),mean_mse=('mse','mean'),
                                 mean_stability_var=('stability_var','mean'),mean_failure_rate=('failure_rate','mean')).reset_index().sort_values('mean_mse')
ssfagg.to_csv(TAB/'table3_ssf_benchmark_summary.csv',index=False)

sw=samp.groupby(['schedule','obs_per_cycle']).agg(worst_fidelity=('fidelity','min'),median_fidelity=('fidelity','median'),
                                                  worst_sign_accuracy=('sign_accuracy','min'),worst_rmse=('rmse','max'),min_mean_observed=('mean_observed','min')).reset_index()
sw.to_csv(TAB/'table4_sampling_worstcase.csv',index=False)
sur.to_csv(TAB/'table5_surrogate_calibration_power.csv',index=False)

# ---------- Figure 1: mechanism + reference recovery ----------
fig,ax=plt.subplots(1,3,figsize=(14,4.1))
x=np.linspace(.1,.9,400); ax[0].plot(x,H.inverted_u(x),lw=2)
mu=H.DA_OPT-H.K_GAIN*H._E2_MEAN
for b in (mu-.12,mu+.12):
    da=b+H.K_GAIN*np.linspace(.13,1,70); ax[0].plot(da,H.inverted_u(da),lw=4)
ax[0].axvline(H.DA_OPT,ls=':'); ax[0].set_xlabel('Dopaminergic tone'); ax[0].set_ylabel('Performance'); ax[0].set_title('(a) Bidirectional mechanism')
for cyc,g in ref.groupby('n_cycles'):
    g=g.sort_values('n_obs'); ax[1].plot(g.n_obs,g.fidelity,marker='o',label=f'{cyc} cycle(s)')
ax[1].axhline(.70,ls='--'); ax[1].set_ylim(.45,.86); ax[1].set_xlabel('Total observations/person'); ax[1].set_ylabel('Recovery fidelity'); ax[1].set_title('(b) Fidelity under reference model'); ax[1].legend(fontsize=8)
for cyc,g in ref.groupby('n_cycles'):
    g=g.sort_values('n_obs'); ax[2].plot(g.n_obs,g.sign_accuracy,marker='o',label=f'{cyc} cycle(s)')
ax[2].set_ylim(.65,.90); ax[2].set_xlabel('Total observations/person'); ax[2].set_ylabel('Correct sign proportion'); ax[2].set_title('(c) Directional accuracy')
fig.tight_layout(); fig.savefig(FIG/'figure1_reference_recovery.png',dpi=220,bbox_inches='tight'); fig.savefig(FIG/'figure1_reference_recovery.svg',bbox_inches='tight'); plt.close(fig)

# ---------- Figure 2: robustness heatmaps at 2 cycles ----------
g=rob[rob.n_cycles==2].copy(); order=t2.sort_values('fidelity_n20').level.tolist(); ns=[10,20,30,40,56]
def mat(metric):
    p=g.pivot(index='level',columns='n_obs',values=metric).reindex(index=order,columns=ns); return p
fig,ax=plt.subplots(1,2,figsize=(15,8.5))
for a,metric,title,vmin,vmax in [(ax[0],'fidelity','(a) Recovery fidelity',0,1),(ax[1],'sign_accuracy','(b) Directional accuracy',.5,1)]:
    p=mat(metric); im=a.imshow(p.values,aspect='auto',vmin=vmin,vmax=vmax)
    a.set_xticks(range(len(ns)));a.set_xticklabels(ns);a.set_yticks(range(len(order)));a.set_yticklabels([s.replace('_',' ') for s in order],fontsize=8)
    a.set_xlabel('Total observations/person (2 cycles)');a.set_title(title);fig.colorbar(im,ax=a,fraction=.035,pad=.02)
fig.tight_layout();fig.savefig(FIG/'figure2_robustness_heatmaps.png',dpi=220,bbox_inches='tight');fig.savefig(FIG/'figure2_robustness_heatmaps.svg',bbox_inches='tight');plt.close(fig)

# ---------- Figure 3: SSF benchmark ----------
methods=['ssf_020','ssf_025','ssf_adaptive','acf_linear','spline_gcv','state_space']
labels={'ssf_020':'SSF cutoff .20','ssf_025':'SSF cutoff .25','ssf_adaptive':'SSF adaptive','acf_linear':'ACF-linear','spline_gcv':'Spline CV','state_space':'State-space','haar_wavelet':'Haar wavelet'}
fig,ax=plt.subplots(1,3,figsize=(15,4.4))
a=ssf[ssf.method.isin(methods)].groupby(['n','method']).mse.mean().reset_index()
for m in methods:
    z=a[a.method==m];ax[0].plot(z.n,z.mse,marker='o',label=labels[m])
ax[0].set_xlabel('Series length');ax[0].set_ylabel('Mean squared error');ax[0].set_title('(a) MSE across shapes and truth levels');ax[0].legend(fontsize=7)
b=ssf[ssf.method.isin(methods)].groupby(['shape','method']).mse.mean().unstack('method').reindex(columns=methods)
b.plot(kind='bar',ax=ax[1],legend=False);ax[1].set_ylabel('Mean squared error');ax[1].set_xlabel('Signal shape');ax[1].set_title('(b) Shape sensitivity');ax[1].tick_params(axis='x',rotation=0)
c=ssf[ssf.method.isin(['ssf_020','ssf_025','ssf_030','ssf_adaptive'])].groupby(['n','method']).failure_rate.mean().reset_index()
for m in ['ssf_020','ssf_025','ssf_030','ssf_adaptive']:
    z=c[c.method==m];ax[2].plot(z.n,z.failure_rate,marker='o',label=labels.get(m,m))
ax[2].set_xlabel('Series length');ax[2].set_ylabel('Failure rate');ax[2].set_title('(c) Spectral estimator feasibility');ax[2].legend(fontsize=7)
fig.tight_layout();fig.savefig(FIG/'figure3_ssf_benchmark.png',dpi=220,bbox_inches='tight');fig.savefig(FIG/'figure3_ssf_benchmark.svg',bbox_inches='tight');plt.close(fig)

# ---------- Figure 4: sampling design worst across peak error ----------
fig,ax=plt.subplots(2,2,figsize=(12,8.2),sharex=True,sharey=True)
for a,miss in zip(ax.ravel(),['none','mcar10','mcar25','weekend']):
    z=samp[samp.missingness==miss].groupby(['schedule','obs_per_cycle']).fidelity.min().reset_index()
    for sch in ['uniform','phase_targeted','adaptive']:
        q=z[z.schedule==sch];a.plot(q.obs_per_cycle,q.fidelity,marker='o',label=sch.replace('_',' '))
    a.axhline(.70,ls='--');a.set_title(miss.replace('mcar','MCAR ').replace('weekend','Weekend non-response'))
    a.set_xlabel('Planned observations/cycle');a.set_ylabel('Worst-case fidelity across peak error')
ax[0,0].legend(fontsize=8);fig.tight_layout();fig.savefig(FIG/'figure4_sampling_design.png',dpi=220,bbox_inches='tight');fig.savefig(FIG/'figure4_sampling_design.svg',bbox_inches='tight');plt.close(fig)

# ---------- Figure 5: surrogate calibration and power ----------
fig=plt.figure(figsize=(14,7));gs=fig.add_gridspec(2,3,height_ratios=[1,2])
a=fig.add_subplot(gs[0,:]);cal=sur[sur.analysis.str.startswith('calibration')].sort_values(['beta_mean','rho']);xpos=np.arange(len(cal));
a.bar(xpos,cal.rejection_rate,yerr=1.96*cal.mcse,capsize=4);a.axhline(.05,ls='--');a.set_xticks(xpos);a.set_xticklabels([f'beta={b:.1f}, rho={r:.1f}' for b,r in zip(cal.beta_mean,cal.rho)]);a.set_ylabel('Type I error');a.set_ylim(0,.09);a.set_title('(a) Surrogate-test calibration (1,000 Monte Carlo replications)')
for j,sd in enumerate((.1,.2,.3)):
    a=fig.add_subplot(gs[1,j]);p=sur[(sur.analysis=='power')&(sur.beta_sd==sd)].pivot(index='N',columns='nobs',values='rejection_rate').reindex(index=[20,40,60],columns=[30,60,90]);im=a.imshow(p.values,origin='lower',vmin=0,vmax=1,aspect='auto')
    a.set_xticks(range(3));a.set_xticklabels([30,60,90]);a.set_yticks(range(3));a.set_yticklabels([20,40,60]);a.set_xlabel('Observations/person');a.set_ylabel('Participants');a.set_title(f'(b{j+1}) Slope SD = {sd:.1f}')
    for ii in range(3):
        for jj in range(3): a.text(jj,ii,f'{p.values[ii,jj]:.2f}',ha='center',va='center',fontsize=9)
fig.colorbar(im,ax=fig.axes[1:],fraction=.015,pad=.02,label='Power');fig.subplots_adjust(hspace=.42,wspace=.34);fig.savefig(FIG/'figure5_surrogate_power.png',dpi=220,bbox_inches='tight');fig.savefig(FIG/'figure5_surrogate_power.svg',bbox_inches='tight');plt.close(fig)

# ---------- Figure 6: empirical illustration reconstructed from audited outputs ----------
emp=pd.DataFrame([
    ['fatigue',.1387,.1395,.533],['mood swing',.1301,.1468,.880],['resting heart rate',.2161,.1950,.156],['skin temperature',.1155,.1324,.874]
],columns=['outcome','observed_sd','null_median','p'])
ssfemp=pd.DataFrame([
    ['E3G',.4271],['fatigue',.3618],['mood swing',.3783],['cramps',.5724],['bloating',.3491],['resting heart rate',.8914],['skin temperature',.0820]
],columns=['measure','ssf'])
fig,ax=plt.subplots(1,3,figsize=(15,4.5));xx=np.arange(len(emp));w=.36
ax[0].bar(xx-w/2,emp.observed_sd,w,label='Observed SD(r_i)');ax[0].bar(xx+w/2,emp.null_median,w,label='Null median');ax[0].set_xticks(xx);ax[0].set_xticklabels(emp.outcome,rotation=25,ha='right');ax[0].set_ylabel('Between-person dispersion');ax[0].set_title('(a) Coupling heterogeneity');ax[0].legend(fontsize=8)
for i,p in enumerate(emp.p):ax[0].text(i,max(emp.observed_sd.iloc[i],emp.null_median.iloc[i])+.006,f'p={p:.3f}',ha='center',fontsize=7)
xx=np.arange(len(ssfemp));ax[1].bar(xx,ssfemp.ssf);ax[1].set_xticks(xx);ax[1].set_xticklabels(ssfemp.measure,rotation=30,ha='right',fontsize=8);ax[1].set_ylim(0,1);ax[1].set_ylabel('SSF');ax[1].set_title('(b) Empirical smooth-signal fraction')
ax[2].bar([0,1],[.0782,.1744]);ax[2].set_xticks([0,1]);ax[2].set_xticklabels(['|r| with E3G level','Phase eta-squared']);ax[2].set_ylim(0,.22);ax[2].set_title('(c) Cramps: predictor alignment')
fig.tight_layout();fig.savefig(FIG/'figure6_empirical_illustration.png',dpi=220,bbox_inches='tight');fig.savefig(FIG/'figure6_empirical_illustration.svg',bbox_inches='tight');plt.close(fig)

print('tables',len(list(TAB.glob('*.csv'))),'figures',len(list(FIG.glob('*.png'))))
