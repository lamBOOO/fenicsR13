#!/usr/bin/env sh
set -eu

# Build a ParaView-readable XDMF file series from Knudsen sweep folders.
#
# Usage:
#   ./build_paraview_series.sh [field|all] [source_glob] [output_dir_or_prefix]
#
# Examples:
#   ./build_paraview_series.sh
#   ./build_paraview_series.sh all
#   ./build_paraview_series.sh u
#   ./build_paraview_series.sh theta 'results_lid_driven_cavity_kn_*'
#   ./build_paraview_series.sh all 'results_lid_driven_cavity_kn_*' results_lid_driven_cavity_series
#
# Behavior:
# - field=all (default): builds 5 series folders:
#   <prefix>_theta, <prefix>_s, <prefix>_p, <prefix>_u, <prefix>_sigma
# - single field: builds one folder. If output_dir_or_prefix is set, it is used
#   as the folder name; otherwise defaults to results_lid_driven_cavity_series_<field>

field="${1:-all}"
source_glob="${2:-results_lid_driven_cavity_kn_*}"
output_arg="${3:-results_lid_driven_cavity_series}"

# Enter script directory so relative paths match this example folder.
script_dir="$(cd "$(dirname "$0")" && pwd)"
cd "${script_dir}"

tmp_dirs_file="$(mktemp)"
tmp_sort_file="$(mktemp)"
cleanup() {
  rm -f "${tmp_dirs_file}"
  rm -f "${tmp_sort_file}"
}
trap cleanup EXIT HUP INT TERM

replace_inplace() {
  file_path="$1"
  from_pattern="$2"
  to_pattern="$3"
  if sed --version >/dev/null 2>&1; then
    sed -i -e "s|${from_pattern}|${to_pattern}|g" "${file_path}"
  else
    sed -i '' -e "s|${from_pattern}|${to_pattern}|g" "${file_path}"
  fi
}

for d in ${source_glob}; do
  if [ -d "${d}" ]; then
    kn_value="${d##*_kn_}"
    # Keep only directories whose suffix can be parsed as a number.
    if printf '%s\n' "${kn_value}" | grep -Eq '^[0-9]+([.][0-9]+)?([eE][-+]?[0-9]+)?$'; then
      printf '%s\t%s\n' "${kn_value}" "${d}"
    fi
  fi
done | sort -g -k1,1 > "${tmp_sort_file}"

cut -f2- "${tmp_sort_file}" > "${tmp_dirs_file}"

if [ ! -s "${tmp_dirs_file}" ]; then
    echo "No folders found for pattern: ${source_glob}" >&2
    echo "Run the parameter study first or provide a matching glob." >&2
    exit 1
fi

build_one_field() {
  one_field="$1"
  one_output_dir="$2"
  map_file="${one_output_dir}/index_kn_map.csv"
  mkdir -p "${one_output_dir}"
  printf 'index,kn,source_dir\n' > "${map_file}"

  count=0
  while IFS= read -r d; do
    xdmf_src="${d}/${one_field}_0.xdmf"
    h5_src="${d}/${one_field}_0.h5"

    if [ ! -f "${xdmf_src}" ] || [ ! -f "${h5_src}" ]; then
      echo "Skip ${d}: missing ${one_field}_0.xdmf or ${one_field}_0.h5"
      continue
    fi

    xdmf_dst="${one_output_dir}/${one_field}_${count}.xdmf"
    cp "${xdmf_src}" "${xdmf_dst}"
    cp "${h5_src}" "${one_output_dir}/${one_field}_${count}.h5"

    # ParaView needs internal references/time to match renamed files.
    replace_inplace "${xdmf_dst}" "${one_field}_0" "${one_field}_${count}"

    kn_value="${d##*_kn_}"
    replace_inplace "${xdmf_dst}" '<Time Value="0" />' "<Time Value=\"${kn_value}\" />"
    printf '%s,%s,%s\n' "${count}" "${kn_value}" "${d}" >> "${map_file}"
    count=$((count + 1))
  done < "${tmp_dirs_file}"

  if [ "${count}" -eq 0 ]; then
    echo "No matching field files were copied for ${one_field}." >&2
    return 1
  fi

  echo "Created ${count} steps in ${one_output_dir}"
  echo "Wrote index map: ${map_file}"
  echo "Open ${one_output_dir}/${one_field}_0.xdmf in ParaView and enable File Series."
  return 0
}

if [ "${field}" = "all" ]; then
  failed=0
  for one_field in theta s p u sigma; do
    if ! build_one_field "${one_field}" "${output_arg}_${one_field}"; then
      failed=1
    fi
  done
  echo "Then use Animation View -> Save Animation to export the video."
  if [ "${failed}" -ne 0 ]; then
    exit 1
  fi
else
  case "${field}" in
    theta|s|p|u|sigma)
      ;;
    *)
      echo "Unknown field: ${field}" >&2
      echo "Use one of: theta, s, p, u, sigma, all" >&2
      exit 1
      ;;
  esac

  if [ "${3:-}" != "" ]; then
    single_output_dir="${output_arg}"
  else
    single_output_dir="results_lid_driven_cavity_series_${field}"
  fi

  build_one_field "${field}" "${single_output_dir}"
  echo "Then use Animation View -> Save Animation to export the video."
fi
