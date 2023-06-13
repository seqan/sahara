#pragma once

#include <cassert>
#include <fmt/format.h>

inline void printCompletion(std::string gen) {
    if (gen == std::string{"bash"}) {
        fmt::print("{}\ncomplete -F clice_GetOpts {}",
R"(function clice_GetOpts ()
{
    WORD_TO_COMPLETE=${COMP_WORDS[COMP_CWORD]}
    LINE="${COMP_LINE[@]:0:COMP_POINT}"
    if [ "${COMP_LINE[@]:COMP_POINT-1:1}" = ' ' ]; then
        LINE="${LINE} ''"
        WORD_TO_COMPLETE=""
    fi

    readarray -t HINTS <<< $(CLICE_COMPLETION= eval ${LINE})

    if [ "${HINTS[0]}" == " -d " ]; then
        compopt -o filenames
        local IFS=$'\n'
        COMPREPLY=( $(compgen -d -- ${2}) )
    elif [ "${HINTS[0]}" == " -f " ]; then
        compopt -o filenames
        local IFS=$'\n'
        COMPREPLY=( $(compgen -f -- ${2}) )
    elif [[ "${HINTS[0]}" =~ " -f " ]]; then
        compopt -o filenames
        parts=(${HINTS[0]})
        EXT="${parts[1]}"

        HINTS=$(printf "%q:" "${HINTS[@]}")
        unset HINTS[0]
        unset HINTS[1]
        local IFS=":"
        if [ ${#HINTS[@]} -gt 0 ]; then
            readarray -t COMPREPLY <<< $(compgen -W '${HINTS[@]}' -- "${WORD_TO_COMPLETE}")
        fi

        local IFS=$'\n'
        COMPREPLY+=( $(compgen -d -- ${2}) )
        COMPREPLY+=( $(compgen -f -X '!*'${EXT} -- ${2}) )

    elif [ ${#HINTS} -eq 0 ]; then
        COMPREPLY=()
    else
        HINTS=$(printf "%q:" "${HINTS[@]}")
        local IFS=":"
        mapfile -t COMPREPLY <<< $(compgen -W '${HINTS[@]}' -- "${WORD_TO_COMPLETE}")
    fi
})", argv0);
    } else if (gen == std::string{"zsh"}) {
        fmt::print("{}\ncompdef clice_GetOpts -P $(basename {}) -N",
R"abc(clice_GetOpts () {
    completions=$(CLICE_COMPLETION= ${(@q)words[@]:0:$CURRENT})
    completions=(${(@f)completions})

    if [ "${completions}" = " -d " ]; then
        _files -/
    elif [ "${completions}" = " -f " ]; then
        _files
    elif [[ "${completions}" =~ " -f " ]]; then
        EXT="$(echo "${completions}" | cut -d ' ' -f 3)"
        _files -g "*${EXT}"
    else
        compadd ${completions[@]}
    fi
})abc", argv0);
    } else {
        fmt::print("unknown generator for shell '{}'\n", gen);
        exit(1);
    }
}
