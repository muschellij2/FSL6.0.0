#!/bin/sh

#Only unlink if scripts are in subdirectory of FSLDIR
if [ ! -z ${FSLDIR} ]; then
    case "$PREFIX" in "${FSLDIR}"*)
	for script in $@; do
	    targetScript="${FSLDIR}/bin/$script"
	    if [ -L $targetScript ]; then
		rm "$targetScript"
	    fi
	done
    esac
fi