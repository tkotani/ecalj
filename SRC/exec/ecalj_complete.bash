# ecalj_complete.bash -- bash completion for ecalj scripts.
# Source this in ~/.bashrc:    source /home/takao/bin/ecalj_complete.bash
#
# Adds tab-completion for the <sname> positional of every ecalj script.
# What is offered depends on the input file each script consumes:
#
#   Legacy2toml.py <TAB>   -> ctrl.*       (legacy ctrl files)
#   ctrlgenToml.py <TAB>   -> ctrls.*      (ctrls.<sname> seed files)
#   lmf / lmfa / lmchk /
#   gwsc / gw_lmfh /
#   eps_lmfh / epsPP_lmfh /
#   genMLWF / genMLWFx     -> ctrlg.*.toml (TOML inputs)
#
# All globs are scoped to the current working directory; completion only
# shows candidates that actually exist for the script you are about to run.

_ecalj_complete_from_glob() {
    local cur="${COMP_WORDS[COMP_CWORD]}"
    local pattern="$1"  prefix="$2"  suffix="$3"
    local f stripped snames=()
    shopt -q nullglob; local saved_nullglob=$?
    shopt -s nullglob
    for f in $pattern; do
        [ -f "$f" ] || continue
        stripped="${f#$prefix}"
        stripped="${stripped%$suffix}"
        snames+=("$stripped")
    done
    [ $saved_nullglob -ne 0 ] && shopt -u nullglob
    COMPREPLY=( $(compgen -W "${snames[*]}" -- "$cur") )
}

_ecalj_legacy_sname()  { _ecalj_complete_from_glob 'ctrl.*'        'ctrl.'  ''      ; }
_ecalj_ctrls_sname()   { _ecalj_complete_from_glob 'ctrls.*'       'ctrls.' ''      ; }
# Offer the full filename (`ctrlg.nio.toml`) rather than the bare sname
# (`nio`) so the user can see what's actually in cwd. m_ext_init strips
# both the `ctrlg.` prefix and the `.toml` suffix at startup, so the
# binary accepts either form.
_ecalj_toml_sname()    { _ecalj_complete_from_glob 'ctrlg.*.toml'  ''       ''      ; }

# Tab-complete --<flag> and --ctrlg:<dotted.path>= for the Fortran
# binaries (lmf, lmfa, lmchk).
#
#   <TAB>          -> sname from ctrlg.*.toml
#   --<TAB>        -> all registered cmdopt0/cmdopt2 flags + --ctrlg:
#                     (dumped at install time into BINDIR/ecalj_cmdopts.list
#                      by InstallAll.py via `lmf --listcmdopt`)
#   --ctrlg:<TAB>  -> every dotted-path key actually present in the
#                     cwd's ctrlg.*.toml ([[spec]] / [[site]] arrays are
#                     expanded into spec.1.r, spec.2.r, ...).
_ecalj_fortran_complete() {
    local cur="${COMP_WORDS[COMP_CWORD]}"
    local prev="${COMP_WORDS[COMP_CWORD-1]:-}"
    # Reconstruct the full word when bash split it on a `:`
    # (COMP_WORDBREAKS default), e.g. "--ctrlg:bz." -> ("--ctrlg",":","bz.").
    # Two cursor positions are possible:
    #   `lmf --ctrlg:<TAB>`     -> COMP_CWORD points at ":"  (cur=":")
    #   `lmf --ctrlg:bz<TAB>`   -> COMP_CWORD points at "bz" (cur="bz", prev=":")
    # Both need the leading "--ctrlg" stitched back on. For a plain sname
    # argument ("lmf nio<TAB>") leave cur as the raw current token so the
    # positional branch below ($cur not matching -*) falls into
    # _ecalj_toml_sname.
    if [ "$prev" = ":" ] || [ "$cur" = ":" ]; then
        local i=$COMP_CWORD
        cur=""
        while [ $i -ge 0 ]; do
            local w="${COMP_WORDS[$i]}"
            cur="$w$cur"
            [[ "$w" == -* ]] && break
            i=$((i - 1))
        done
    fi

    local script_dir="$(dirname "${BASH_SOURCE[0]}")"
    local cmdopt_list="$script_dir/ecalj_cmdopts.list"
    case "$cur" in
        --ctrlg:*)
            local toml
            toml=$(ls ctrlg.*.toml 2>/dev/null | head -1)
            [ -z "$toml" ] && return 0
            local prefix="${cur#--ctrlg:}"
            # Bypass compgen here so that values containing TOML/glob
            # metacharacters ([ ] etc.) are not subjected to filename
            # expansion. Python emits `key=<current-value>` lines
            # filtered by $prefix; mapfile reads them verbatim.
            mapfile -t COMPREPLY < <(python3 - "$toml" "$prefix" 2>/dev/null <<'PYEOF'
import sys, tomllib, shlex
# to_toml: render the value as a valid TOML literal so the user can
# accept the completion and run the command without further editing.
#   bool -> true / false
#   num  -> as repr (canonical Python form is TOML-compatible)
#   str  -> "..."   (basic string, inner " and \ escaped)
#   list -> [a,b,c] / nested
# Multi-line strings have no single-line TOML literal that fits in
# one CLI arg; return None so the key prints as `key=` with no value.
def to_toml(v):
    if isinstance(v, bool):           return 'true' if v else 'false'
    if isinstance(v, (int, float)):   return repr(v)
    if isinstance(v, str):
        if '\n' in v: return None
        return '"' + v.replace('\\', '\\\\').replace('"', '\\"') + '"'
    if isinstance(v, list):
        parts = [to_toml(x) for x in v]
        if any(p is None for p in parts): return None
        return '[' + ','.join(parts) + ']'
    return None
def walk(d, p=''):
    for k, v in d.items():
        if isinstance(v, dict):
            yield from walk(v, p + k + '.')
        elif isinstance(v, list) and v and isinstance(v[0], dict):
            for i, it in enumerate(v, 1):
                yield from walk(it, f'{p}{k}.{i}.')
        else:
            t = to_toml(v)
            # shlex.quote wraps the TOML value in single quotes when
            # it contains bash-special chars ([ ] " etc.), so bash
            # quote-removal at execute time gives the binary the
            # unescaped TOML form (`bz.nkabc=[8,8,8]`, `symgrp="find"`).
            yield f'{p}{k}=' if t is None else f'{p}{k}={shlex.quote(t)}'
prefix = sys.argv[2]
for line in walk(tomllib.load(open(sys.argv[1], 'rb'))):
    if line.startswith(prefix):
        print(line)
PYEOF
)
            # User may want to edit the value, so keep the cursor at
            # the end of the inserted text rather than padding a space.
            compopt -o nospace 2>/dev/null
            ;;
        -*)
            local flags=""
            [ -r "$cmdopt_list" ] && flags=$(cat "$cmdopt_list")
            COMPREPLY=( $(compgen -W "$flags --ctrlg:" -- "$cur") )
            # If the unique completion ends with ':' (--ctrlg:) or '='
            # (a cmdopt2 like --jobgw=), the user still needs to type
            # the value -- suppress the trailing space.
            if [ ${#COMPREPLY[@]} -eq 1 ]; then
                case "${COMPREPLY[0]}" in
                    *:|*=) compopt -o nospace 2>/dev/null ;;
                esac
            fi
            ;;
        *)
            _ecalj_toml_sname
            ;;
    esac
}

complete -F _ecalj_legacy_sname    Legacy2toml.py
complete -F _ecalj_ctrls_sname     ctrlgenToml.py
complete -F _ecalj_toml_sname \
    gwsc gw_lmfh eps_lmfh epsPP_lmfh \
    genMLWF genMLWFx
complete -F _ecalj_fortran_complete lmf lmfa lmchk

# Delegate to _ecalj_fortran_complete when the user types
# `mpirun [-np N] lmf|lmfa|lmchk ...`. Without this, tab-completion
# is governed by mpirun's own completion (or the bash default) and
# misses our --<flag> / --ctrlg:<path> handling.
_ecalj_mpirun_passthrough() {
    local i=1
    while [ $i -lt $COMP_CWORD ]; do
        case "${COMP_WORDS[$i]##*/}" in
            lmf|lmfa|lmchk)
                _ecalj_fortran_complete
                return $?
                ;;
        esac
        i=$((i + 1))
    done
    # No ecalj binary seen yet -- let bash's default (filename) kick in.
    COMPREPLY=()
}
complete -F _ecalj_mpirun_passthrough -o default mpirun mpiexec
