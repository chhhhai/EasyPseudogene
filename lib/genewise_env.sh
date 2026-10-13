# Shared GeneWise executable/configuration discovery. Source from Bash.
# Explicit WISEPATH/WISECONFIGDIR values are validated, never silently replaced.

epg_wise_config_ok() {
    local directory="$1" file
    [[ -d "${directory}" ]] || return 1
    for file in codon.table BLOSUM62.bla gene.stat; do
        [[ -r "${directory}/${file}" && -s "${directory}/${file}" ]] || return 1
    done
}

epg_setup_genewise() {
    local binary="" resolved="" prefix candidate file
    local -a prefixes=() candidates=()

    if [[ -n "${WISEPATH:-}" ]]; then
        binary="${WISEPATH}/genewise"
        if [[ ! -f "${binary}" || ! -x "${binary}" ]]; then
            echo "ERROR: --wisepath/WISEPATH does not contain an executable genewise: ${binary}" >&2
            echo "  --wisepath expects the directory containing genewise, not wisecfg." >&2
            return 1
        fi
    else
        binary="$(type -P genewise || true)"
        if [[ -z "${binary}" && -n "${CONDA_PREFIX:-}" && -x "${CONDA_PREFIX}/bin/genewise" ]]; then
            binary="${CONDA_PREFIX}/bin/genewise"
        fi
    fi
    if [[ -z "${binary}" || ! -f "${binary}" || ! -x "${binary}" ]]; then
        echo "ERROR: GeneWise executable (genewise) was not found in PATH or the active Conda environment." >&2
        echo "  Activate the environment containing wise2: conda activate <environment>" >&2
        echo "  Or specify: --wisepath /path/to/genewise/bin" >&2
        echo "  The package is named wise2; the executable is named genewise." >&2
        return 1
    fi

    WISEPATH="$(cd "$(dirname "${binary}")" && pwd -P)" || return 1
    GENEWISE_BIN="${WISEPATH}/genewise"
    export WISEPATH GENEWISE_BIN
    export PATH="${WISEPATH}:${PATH}"

    if [[ -n "${WISECONFIGDIR:-}" ]]; then
        if ! epg_wise_config_ok "${WISECONFIGDIR}"; then
            echo "ERROR: --wiseconfig/WISECONFIGDIR is not a usable GeneWise configuration directory: ${WISECONFIGDIR}" >&2
            for file in codon.table BLOSUM62.bla gene.stat; do
                if [[ ! -r "${WISECONFIGDIR}/${file}" || ! -s "${WISECONFIGDIR}/${file}" ]]; then
                    echo "  Missing, empty or unreadable: ${WISECONFIGDIR}/${file}" >&2
                fi
            done
            return 1
        fi
    else
        # Inspect both the invoked path and its symlink target. The second
        # parent also covers source installations with src/bin/genewise.
        resolved="$(python3 -c 'import pathlib,sys; print(pathlib.Path(sys.argv[1]).resolve())' "${GENEWISE_BIN}")" || return 1
        for binary in "${GENEWISE_BIN}" "${resolved}"; do
            prefix="$(dirname "$(dirname "${binary}")")"
            prefixes+=("${prefix}" "$(dirname "${prefix}")")
        done
        [[ -z "${CONDA_PREFIX:-}" ]] || prefixes+=("${CONDA_PREFIX}")
        prefixes+=(/usr/local /usr)
        for prefix in "${prefixes[@]}"; do
            candidates+=(
                "${prefix}/share/wise2/wisecfg"
                "${prefix}/share/wise/wisecfg"
                "${prefix}/wisecfg"
                "${prefix}"/share/wise*/wisecfg
                "${prefix}"/share/wise*
                "${prefix}"/lib/wise*/wisecfg
                "${prefix}"/etc/wise*
            )
        done
        WISECONFIGDIR=""
        for candidate in "${candidates[@]}"; do
            if epg_wise_config_ok "${candidate}"; then
                WISECONFIGDIR="${candidate}"
                break
            fi
        done
        if [[ -z "${WISECONFIGDIR}" ]]; then
            echo "ERROR: Found ${GENEWISE_BIN}, but no usable GeneWise configuration directory." >&2
            echo "  Configuration needs readable, nonempty codon.table, BLOSUM62.bla and gene.stat." >&2
            echo "  Specify: --wiseconfig /path/to/wisecfg (or export WISECONFIGDIR)." >&2
            echo "  Checked directories (including patterns with no matching directory):" >&2
            printf '    %s\n' "${candidates[@]}" >&2
            return 1
        fi
    fi
    WISECONFIGDIR="$(cd "${WISECONFIGDIR}" && pwd -P)" || return 1
    export WISECONFIGDIR
    echo "GeneWise executable: ${GENEWISE_BIN}"
    echo "GeneWise configuration: ${WISECONFIGDIR}"
}
