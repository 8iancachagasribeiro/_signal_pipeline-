#!/usr/bin/env bash
set -euo pipefail
SIGMAS="${SIGMAS:-0.05 0.10 0.15 0.20}"
NCHUNKS="${NCHUNKS:-8}"
NSIMS="${NSIMS:-500}"
OUT="${OUT:-./results/chunks}"
mkdir -p "$OUT"
echo "sigma_b: $SIGMAS | chunks: $NCHUNKS | reps/cell: $NSIMS | out: $OUT"
for SB in $SIGMAS; do
  for CH in $(seq 0 $((NCHUNKS-1))); do
    TARGET="$OUT/h4_500_sb${SB}_c${CH}.csv"
    if [ -f "$TARGET" ]; then echo "exists: $TARGET"; continue; fi
    python run_h4_500.py "$SB" "$CH" "$NCHUNKS" --n-sims "$NSIMS" --out-dir "$OUT" > "$OUT/progress_${SB}_${CH}.log" 2>&1 &
  done
  wait
done
python - "$OUT" <<'PY'
import glob,os,sys,pandas as pd
out=sys.argv[1]; files=sorted(glob.glob(os.path.join(out,'h4_500_sb*_c*.csv')))
if not files: raise SystemExit('no chunk CSVs found')
df=pd.concat([pd.read_csv(f) for f in files],ignore_index=True).drop_duplicates(['sigma_b','obs_per_cycle','reliability'],keep='last').sort_values(['sigma_b','obs_per_cycle','reliability'])
path=os.path.join(os.path.dirname(out.rstrip('/')) or '.','h4_frontier_500reps.csv'); df.to_csv(path,index=False); print(f'[saved] {path} ({len(df)} cells)')
PY
