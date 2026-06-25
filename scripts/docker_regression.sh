#!/usr/bin/env bash
# Docker-based DEA regression for prolfquapp.
#
# Runs a fixture through several released Docker image versions and models,
# then correlates the protein-level fold-change / p.value / FDR estimates.
#
# Reference chain (see README "Docker regression"):
#   - 0.1.8 : old reference, only `lm` and `lm_missing` (model via config
#             `model_missing`; no `-m` flag, older config schema)
#   - 2.2.6 : model-complete reference (model via `-m <facade>`)
#
# Comparisons:
#   0.1.8/lm         vs 2.2.6/lm          (core linear model stable across span)
#   0.1.8/lm_missing vs 2.2.6/lm_missing  (legacy missing-group model stable)
#   2.2.6/lm_impute  vs 2.2.6/lm_missing  (new rescue ~ legacy rescue)
#
# Usage: scripts/docker_regression.sh [fixture] [software_2x] [dataset] [min_cor]
set -euo pipefail

FIXTURE="${1:-maxquant_ionstar}"
SW_2X="${2:-prolfquapp.MAXQUANT}"
DATASET="${3:-dataset.csv}"
MIN_COR="${4:-0.95}"
SW_018="${SW_2X#prolfquapp.}"   # 0.1.x uses the unprefixed software key

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
DOCKER_SH="$ROOT/prolfquapp_docker.sh"
OUT="$ROOT/test-outputs/regression/$FIXTURE"
mkdir -p "$OUT"

# run_dea <tag> <run_label> <extra dea args...>
# Sets up an isolated workdir, runs prolfqua_dea.sh in the given image, and
# echoes the produced DE_*.xlsx path.
run_dea() {
  local tag="$1" label="$2"; shift 2
  local wd="$OUT/$label"
  rm -rf "$wd"; mkdir -p "$wd"
  cp -r "$ROOT/fixtures/$FIXTURE/." "$wd/"
  ( cd "$wd" && "$DOCKER_SH" --image-version "$tag" prolfqua_dea.sh \
      -i . -d "$DATASET" -y config.yaml -o . -w "$label" "$@" ) >"$wd/dea.log" 2>&1
  find "$wd" -name 'DE_*.xlsx' | head -1
}

echo ">>> 0.1.8 (native config; model via model_missing)"
# 0.1.8 rejects the 2.x config schema, so generate a native one, then toggle
# model_missing to pick lm (no) vs lm_missing (yes).
gen_018_config() {  # <workdir> <model_missing yes|no>
  ( cd "$1" && "$DOCKER_SH" --image-version 0.1.8 prolfqua_yaml.sh \
      -n robscale -y config.yaml -w r -s "$SW_018" -o . ) >/dev/null 2>&1
  sed -i.bak "s/^  model_missing: .*/  model_missing: $2/" "$1/config.yaml"
}
for combo in "lm:no" "lm_missing:yes"; do
  model="${combo%%:*}"; mm="${combo##*:}"
  wd="$OUT/018_$model"; rm -rf "$wd"; mkdir -p "$wd"
  cp -r "$ROOT/fixtures/$FIXTURE/." "$wd/"
  gen_018_config "$wd" "$mm"
  ( cd "$wd" && "$DOCKER_SH" --image-version 0.1.8 prolfqua_dea.sh \
      -i . -d "$DATASET" -y config.yaml -s "$SW_018" -o . -w "018_$model" ) \
      >"$wd/dea.log" 2>&1
  declare "XLSX_018_$model=$(find "$wd" -name 'DE_*.xlsx' | head -1)"
done

echo ">>> 2.2.6 (config schema 2.x; model via -m)"
XLSX_226_lm=$(run_dea 2.2.6 226_lm -s "$SW_2X" -m lm)
XLSX_226_lm_missing=$(run_dea 2.2.6 226_lm_missing -s "$SW_2X" -m lm_missing)
XLSX_226_lm_impute=$(run_dea 2.2.6 226_lm_impute -s "$SW_2X" -m lm_impute)

echo ">>> comparing"
Rscript "$ROOT/scripts/compare_docker_versions.R" --min-cor="$MIN_COR" \
  --pairs="0.1.8/lm:2.2.6/lm,0.1.8/lm_missing:2.2.6/lm_missing,2.2.6/lm_impute:2.2.6/lm_missing" \
  "0.1.8/lm=$XLSX_018_lm" \
  "0.1.8/lm_missing=$XLSX_018_lm_missing" \
  "2.2.6/lm=$XLSX_226_lm" \
  "2.2.6/lm_missing=$XLSX_226_lm_missing" \
  "2.2.6/lm_impute=$XLSX_226_lm_impute"
