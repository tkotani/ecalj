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
_ecalj_toml_sname()    { _ecalj_complete_from_glob 'ctrlg.*.toml'  'ctrlg.' '.toml' ; }

complete -F _ecalj_legacy_sname  Legacy2toml.py
complete -F _ecalj_ctrls_sname   ctrlgenToml.py
complete -F _ecalj_toml_sname \
    lmf lmfa lmchk \
    gwsc gw_lmfh eps_lmfh epsPP_lmfh \
    genMLWF genMLWFx
