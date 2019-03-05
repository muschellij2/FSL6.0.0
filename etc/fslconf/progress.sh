#!/usr/bin/env bash

if [ -z "$1" ]; then
    length=150
else
    length=$1
fi

quiet=0
if [ ! -z "$2" ]; then
    quiet=$2
fi
lineno=0
while read -r line; do
    p=$(( ( ${lineno} * 100 ) / ${length} ))
    if [ $p -gt 100 ]; then
        p=100
    fi
    if [ $quiet -eq 0 ]; then
        printf '\r%s%%' ${p} >&2
    fi
    echo ${line} >&1
    lineno=$(( ${lineno} + 1 ))
done
if [ $quiet -eq 0 ]; then
    printf '\r100%%\n' >&2
fi