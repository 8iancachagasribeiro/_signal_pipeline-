# mcPHASES empirical re-audit - 2026-08-24

This directory records the empirical re-run performed from the original credentialed mcPHASES 1.0.0 archive after the Article 1 pipeline audit.

## Structural corrections

- The dataset contains 42 unique participants and 62 participant-intervals. A study interval is now the temporal segment for FFT and surrogate operations. Intervals separated by the long 2022-to-2024 gap are never bridged.
- Repeated intervals are aggregated back to one participant-level coupling estimate with Fisher-z weighting, so intervals are not treated as independent participants.
- Pooled within-person associations are centered within participant-interval.
- Daily Fitbit files are collapsed to one participant-interval-day before joining with E3G. Resting-heart-rate values <= 0 are treated as missing.
- The smooth-signal fraction (SSF) is not treated as classical reliability. Attenuation based on SSF is retained only as an identification sensitivity analysis.

## Confirmatory surrogate results

| Outcome | N participants | SD observed r_i | Null median | p |
|---|---:|---:|---:|---:|
| Fatigue | 41 | 0.1387 | 0.1395 | 0.5329 |
| Mood swing | 41 | 0.1301 | 0.1468 | 0.8802 |

The confirmatory conclusion remains null: neither fatigue nor mood swing shows excess between-person heterogeneity beyond the calendar-preserving phase-randomized null.

## Exploratory objective outcomes

| Outcome | N participants | intervals | SD observed r_i | Null median | p |
|---|---:|---:|---:|---:|---:|
| Resting heart rate | 41 | 61 | 0.2161 | 0.1950 | 0.1557 |
| Skin temperature | 40 | 59 | 0.1155 | 0.1324 | 0.8743 |

The previously reported skin-temperature p approximately .027 does not survive the audit.

## Why the old skin-temperature result appeared significant

The legacy objective-analysis merge joined raw Fitbit rows directly to hormone rows by participant/day. The Fitbit files contain repeated rows for the same nominal daily measure. This multiplied observations before the correlation and surrogate calculations. Re-running that legacy path reproduces the apparent skin-temperature result: SD(r_i)=0.1782, null median=0.1486, p=0.0240, with a median of 137.5 merged rows per participant. After collapsing to one daily value and respecting study intervals, p=0.8743.

This is therefore a data-join artifact, not evidence of a temperature-specific heterogeneous E3G coupling.

## Corrected participant-balanced spectral SSF

- E3G: 0.4271
- Fatigue: 0.3618
- Mood swing: 0.3783
- Cramps: 0.5724
- Bloating: 0.3491
- Resting heart rate: 0.8914
- Skin temperature: 0.0820, but only 13 participants had an eligible >=25-sample contiguous segment under the strict gap policy.

The large changes from older manuscript values are expected because the audited estimator no longer compacts missing-row gaps and no longer treats duplicate daily Fitbit rows as consecutive samples.

## Phase-locked result

Cramping remains the strongest phase-locked result. Its absolute pooled within-episode correlation with E3G level is 0.0782, while cycle phase explains eta-squared=0.1744 of within-episode standardized variance and the menstrual-phase mean is +0.836 within-episode SD. This supports predictor misspecification for cramping: cycle phase carries substantially more information than linear E3G level.

## Manuscript consequences

Keep:
- the simulated masking mechanism and falsifiability boundary;
- the distinction between detecting heterogeneity and recovering person-specific effects;
- the phase-misspecification demonstration for cramping;
- the statement that SSF is not reliability.

Remove or rewrite:
- the exploratory skin-temperature p approximately .027;
- claims that objective temperature detects heterogeneity where self-report does not;
- the old SSF values 0.469/0.323 as if they were fixed reliability parameters;
- empirical-power claims derived directly from SSF;
- statements that instrument repair is proven to recover exact power values.

The raw mcPHASES archive is not committed because its access and redistribution conditions belong to the source dataset. Only derived non-identifying audit outputs are included here.
