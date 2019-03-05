#!/usr/bin/env bash
# Script to link GUI applications into /Applications

# Call with -f <FSLDIR path>, e.g. /usr/local/fsl (with use FSLDIR if given
# no arguments)

arch=`uname -s`
if [ "${arch}" != "Darwin" ]; then
    exit 0
fi

# Where is this script?
script_dir=$( cd $(dirname $0) ; pwd)

# Set some defaults
OPTIND=1
fsl_dir=""
quiet=0
app_list="bin/fslview.app bin/assistant.app bin/fsleyes.app"

function syntax {
    echo "make_application_links.sh [-f <FSLDIR>] [-q]"
    echo "  -f <FSLDIR> Location of installed FSL, e.g. /usr/local/fsl"
    echo "                  if not provided looks for FSLDIR in environment"
    echo "  -q          Don't report progress'"
}

while getopts "h?qf:" opt; do
    case "${opt}" in
    h|\?)
        syntax
        exit 0
        ;;
    q)  quiet=1
        ;;
    f)  fsl_dir=${OPTARG}
        ;;
    esac
done

shift $((OPTIND-1))

[ "$1" = "--" ] && shift

if [ -z "${fsl_dir}" ]; then
    if [ -z "${FSLDIR}" ]; then
        echo "Error - FSLDIR not given as an argument and \$FSLDIR not set!" >&2
        exit 1
    else
        fsl_dir=${FSLDIR}
    fi
fi

if [ ! -e "${fsl_dir}/bin" ]; then
    echo "Error - ${fsl_dir}/bin does not exist!" >&2
    exit 1
fi

if [ $quiet -eq 0 ]; then
    echo "Creating links to applications."
fi

if [ ! -z "${SUDO_UID}" ]; then
    su_do=''
else
    echo "We will require your password to authorise this."
    su_do="sudo"
fi

echo $app_list
for app in ${app_list}; do
    target="${fsl_dir}/${app}"
    if [ ! -e "${target}" ]; then
        continue
    fi
    link="/Applications/`basename ${app}`"
    if [ "${link}" = "/Applications/" ]; then
        continue
    fi
    if [ -L ${link} ]; then
        if [ ! "`readlink ${link}`" = "${target}" ]; then
            ${su_do} rm -f "${link}"
        else
            continue
        fi
    elif [ -e "${link}" ]; then
        echo "${link} exists and isn't a symbolic link - not changing"
        continue
    fi
    ${su_do} ln -sf "${target}" "${link}"
done
