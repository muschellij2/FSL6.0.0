#!/usr/bin/env bash
#
# This script installs miniconda, and configures a
# python environment with all of the dependencies
# required by FSL. This involves:
#
#  1. Downloading the miniconda install script from https://repo.continuum.io
#
#  2. Installing miniconda to $FSLDIR/fslpython/
#
#  3. Creating a miniconda environment called 'fslpython', with all of the
#     packages listed in fslpython_environment.yml (this file is assumed
#     to be present in the same location as this script).
#
#  4. Creating a symlink from the fslpython environment binary to
#     $FSLDIR/bin/fslpython.
#
# Call with -f <FSLDIR path>, e.g. /usr/local/fsl (with use FSLDIR if given
# no arguments)

# Where is this script?
set -o pipefail
script_dir=$( cd $(dirname $0) ; pwd)

# Set some defaults
OPTIND=1
fsl_dir=""
quiet=0
debug=0
dropprivileges=0

# Have we been called by sudo?
if [ ! -z "${SUDO_UID}" ]; then
    dropprivileges=1
fi

function syntax {
    echo "fslpython_install.sh [-f <FSLDIR>] [-q]"
    echo "  -f <FSLDIR> Location of installed FSL, e.g. /usr/local/fsl"
    echo "                  if not provided looks for FSLDIR in environment"
    echo "  -q          Accept Miniconda license automatically"
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
    d)  debug=1
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

if [ ! -w "${fsl_dir}" ]; then
    echo "Error - cannot write to ${fsl_dir}!" >&2
    exit 1
fi

if [ ! -w "${fsl_dir}/bin" ]; then
    echo "Error - cannot write to ${fsl_dir}/bin!" >&2
    exit 1
fi

if [  -e "${fsl_dir}/fslpython" ]; then
    echo "Error - ${fsl_dir}/fslpython already exists!" >&2
    exit 1
fi


function drop_sudo {
    if [ ${dropprivileges} -eq 1 ]; then
        sudo -u \#${SUDO_UID} "$@"
        if [ $? -eq 1 ]; then
            sudo -u \#${SUDO_UID} -g \#${SUDO_GID} "$@"
        fi
    else
        "$@"
    fi
}

#####################################
# Download miniconda installer script
#####################################


platform=`uname -s`
miniconda_url="https://repo.continuum.io/miniconda"
miniconda_tmp=`drop_sudo mktemp -d -t fslpythonXXXX`
if [ $? -ne 0 ]; then
    echo "Failed to create temporary directory"
    exit 2
fi
if [ $debug -eq 1 ]; then
    echo "Debug log: $miniconda_tmp" >&2
fi
miniconda_installer="${miniconda_tmp}/fslpython_miniconda_installer.sh"
miniconda_install_log="${miniconda_tmp}/fslpython_miniconda_installer.log"
miniconda_root_dir="${fsl_dir}/fslpython"
miniconda_bin_dir="${miniconda_root_dir}/bin"
fslpython_env_dir="${miniconda_root_dir}/envs/fslpython/"
drop_sudo echo "Installing FSL conda distribution into ${miniconda_root_dir}" >> "${miniconda_install_log}"

if [ "$platform" = "Linux" ]; then
    if [ `getconf LONG_BIT` -ne 64 ]; then
        echo "We only support 64 bit Linux" >&2
        exit 2
    fi

    miniconda_script="Miniconda3-latest-Linux-x86_64.sh"
elif [ "$platform" = "Darwin" ]; then
    miniconda_script="Miniconda3-latest-MacOSX-x86_64.sh"
else
    echo "Unknown platform" >&2
    exit 2
fi

# Find a downloader...
dl_cmd=$(which wget 2>> "${miniconda_install_log}")
if [ $? -eq 1 ]; then
    dl_cmd=$(which curl 2>> "${miniconda_install_log}")
    if [ $? -eq 1 ]; then
        echo "Unable to locate wget or curl" >&2
        cat "${miniconda_install_log}" >&2
        exit 2
    fi
    dl_cmd_opts="--fail"
    dl_out="-o"
    dl_quiet="-s"
else
    dl_cmd_opts=""
    dl_out="-O"
    dl_quiet="--quiet"
fi

if [ ${quiet} -eq 1 ]; then
    dl_cmd_opts="${dl_cmd_opts} ${dl_quiet}"
fi

drop_sudo ${dl_cmd} ${dl_out} "${miniconda_installer}" ${dl_cmd_opts} \
    ${miniconda_url}/${miniconda_script} 2>> "${miniconda_install_log}"
status=$?
if [ ${status} -ne 0 ]; then
    echo "Failed to download Miniconda - see ${miniconda_install_log} for details" >&2
    exit ${status}
fi

###################
# Install miniconda
###################
if [ ${quiet} -ne 1 ]; then
    echo "Stage 1"
    echo "By installing this python distribution you agree to the license terms in"
    echo "${miniconda_root_dir}/LICENSE.txt"
fi
/usr/bin/env bash ${miniconda_installer} -b -p "${miniconda_root_dir}" \
    2>> "${miniconda_install_log}" | \
    ${script_dir}/progress.sh 26 ${quiet} >> "${miniconda_install_log}"

if [ $? -ne 0 ]; then
    echo "Failed to install Miniconda - see ${miniconda_install_log} for details" >&2
    exit 3
fi
rm "${miniconda_installer}"


##############################
# Create fslpython environment
##############################


# Figure out the location of this script
# so we can figure out the location of the
# conda environment specification file.
if [ ${quiet} -ne 1 ]; then
    echo "Stage 2"
fi
FSLDIR=$fsl_dir "${miniconda_bin_dir}/conda" env create \
    -f "${script_dir}/fslpython_environment.yml" \
    2>> "${miniconda_install_log}" | \
    ${script_dir}/progress.sh 133 ${quiet} 1>> "${miniconda_install_log}"
if [ $? -ne 0 ]; then
    echo "Failed to create FSL Python environment - see ${miniconda_install_log} for details" >&2
    exit 4
fi

# Symlink the environment python
# binary into $FSLDIR/bin/
ln -sf "${fslpython_env_dir}/bin/python" "${fsl_dir}/bin/fslpython"
if [ -e "${fslpython_env_dir}/bin/ipython" ]; then
    ln -sf "${fslpython_env_dir}/bin/ipython" "${fsl_dir}/bin/fslipython"
fi

if [ -e "${fslpython_env_dir}/bin/pythonw" ]; then
  ln -sf "${fslpython_env_dir}/bin/pythonw" "${fsl_dir}/bin/fslpythonw"
else
  ln -sf "${fslpython_env_dir}/bin/python" "${fsl_dir}/bin/fslpythonw"
fi

if [ ${debug} -ne 1 ]; then
    echo "${miniconda_install_log}"
fi

rm "${miniconda_install_log}"
drop_sudo rmdir "${miniconda_tmp}"
