#!/usr/bin/env bash
set -u

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
FIXTURE_DIR="${ROOT_DIR}/fixtures/diann_wu345302"
OUT_ROOT="${ROOT_DIR}/test-outputs/wu345302_facades"
LOG_DIR="${OUT_ROOT}/logs"
CONFIG_DIR="${OUT_ROOT}/configs"
STATUS_FILE="${OUT_ROOT}/status.tsv"

# Same-needs facades (protein-aggregated input → protein FC).
SAME_MODELS=(
  deqms
  deqms_voom
  firth
  limma
  limma_impute
  limma_voom
  limma_voom_impute
  limpa
  lm
  lm_impute
  lm_missing
  rlm
  saint
)

# Nested-needs facades (peptide-level input → protein FC).
NESTED_MODELS=(
  firth_nested
  limpa_nested
  lmer_nested
  ropeca_nested
)

MODELS=("${SAME_MODELS[@]}" "${NESTED_MODELS[@]}")

mkdir -p "${LOG_DIR}" "${CONFIG_DIR}"
printf 'model\trun_status\texit_code\tsoftware\tconfig\toutput_dir\n' > "${STATUS_FILE}"

PROLFQUA_DEA_SH="$(Rscript --vanilla -e "cat(system.file('application/bin/prolfqua_dea.sh', package = 'prolfquapp'))")"
SOURCE_PROLFQUA_DEA_SH="${ROOT_DIR}/../prolfquapp/inst/application/bin/prolfqua_dea.sh"

if [[ ! -x "${PROLFQUA_DEA_SH}" && -x "${SOURCE_PROLFQUA_DEA_SH}" ]]; then
  PROLFQUA_DEA_SH="${SOURCE_PROLFQUA_DEA_SH}"
  export PROLFQUAPP_APP_DIR="${ROOT_DIR}/../prolfquapp/inst/application"
  export PROLFQUAPP_SOURCE_DIR="${ROOT_DIR}/../prolfquapp"
  export PROLFQUA_SOURCE_DIR="${ROOT_DIR}/../prolfqua"
fi

if [[ ! -x "${PROLFQUA_DEA_SH}" ]]; then
  echo "Cannot find executable prolfqua_dea.sh at: ${PROLFQUA_DEA_SH}" >&2
  exit 1
fi

software_for_model() {
  for nested in "${NESTED_MODELS[@]}"; do
    if [[ "${nested}" == "$1" ]]; then
      printf 'prolfquapp.DIANN_PEPTIDE'
      return
    fi
  done
  printf 'prolfquapp.DIANN'
}

make_config() {
  local model="$1"
  local config="$2"

  MODEL="${model}" \
    TEMPLATE="${FIXTURE_DIR}/configs/config_template.yaml" \
    CONFIG="${config}" \
    Rscript --vanilla -e '
      x <- readLines(Sys.getenv("TEMPLATE"))
      x[grepl("^[[:space:]]*model:", x)] <- paste0("  model: ", Sys.getenv("MODEL"))
      writeLines(x, Sys.getenv("CONFIG"))
    '
}

for model in "${MODELS[@]}"; do
  software="$(software_for_model "${model}")"
  config="${CONFIG_DIR}/config_${model}.yaml"
  outdir="${OUT_ROOT}/${model}"
  mkdir -p "${outdir}"
  make_config "${model}" "${config}"

  echo "=== running ${model} (${software}) ==="
  bash "${PROLFQUA_DEA_SH}" \
    --indir "${FIXTURE_DIR}/out-DIANN" \
    --dataset "${FIXTURE_DIR}/dataset_saint.csv" \
    --yaml "${config}" \
    --workunit 345302 \
    --software "${software}" \
    --outdir "${outdir}" \
    > "${LOG_DIR}/${model}.stdout.log" 2>&1
  exit_code=$?

  if [[ "${exit_code}" -eq 0 ]]; then
    run_status="ok"
  else
    run_status="fail"
  fi

  printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
    "${model}" "${run_status}" "${exit_code}" "${software}" "${config}" "${outdir}" \
    >> "${STATUS_FILE}"
  echo "=== ${model} ${run_status} (${exit_code}) ==="
done

Rscript --vanilla "${SCRIPT_DIR}/summarize_wu345302_facades.R" "${OUT_ROOT}"
