#!/bin/sh

#Only link if scripts are in subdirectory of FSLDIR
if [ ! -z ${FSLDIR} ]; then 
    case "$PREFIX" in "${FSLDIR}"*)
	for script in $@; do
	    sourceScript="${PREFIX}/bin/$script"
	    if [ -e $sourceScript ]; then
		targetScript="${FSLDIR}/bin/$script"
		if [ -e $targetScript ]; then
		    rm "$targetScript"
		fi
		ln -s "$sourceScript" "$targetScript"
	    fi
	done
    esac
fi