#!/usr/bin/env bash
#
# This script installs the current version of FSLeyes
#
# Call with -f <FSLDIR path>, e.g. /usr/local/fsl (will use FSLDIR if given
# no arguments)

# Where is this script?
script_dir=$( cd $(dirname $0) ; pwd)

# Set some defaults
OPTIND=1
fsl_dir=""
quiet=0
dropprivileges=0

# Have we been called by sudo?
if [ ! -z "${SUDO_UID}" ]; then
    dropprivileges=1
fi    

function syntax {
    echo "install_fsleyes.sh [-f <FSLDIR>] [-q]"
    echo "  -f <FSLDIR> Location of installed FSL, e.g. /usr/local/fsl"
    echo "                  if not provided looks for FSLDIR in environment"
    echo "  -q          Install quietly"
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

if [ ! -w "${fsl_dir}" ]; then
    echo "Error - cannot write to ${fsl_dir}!" >&2
    exit 1
fi

if [ ! -w "${fsl_dir}/bin" ]; then
    echo "Error - cannot write to ${fsl_dir}/bin!" >&2
    exit 1
fi

if [ ! -e "${fsl_dir}/etc/fslversion" ]; then
    echo "${fsl_dir} doesn't look like an FSL installation folder!" >&2
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
# Download FSLeyes
#####################################

function debian_ubuntu_lts_map {
    if [ $1 -eq 6 ]; then
        echo "10.04"
    elif [ $1 -eq 7 ]; then
        echo "12.04"
    elif [ $1 -eq 8 ]; then
        echo "14.04"
    elif [ $1 -eq 9 ]; then
        echo "16.04"
    fi
}

platform=`uname -s`
if [ "$platform" = "Linux" ]; then
    dl_cmd_opts=""
    dl_out="-O"
    dl_quiet="--quiet"
    dl_cmd="/usr/bin/wget"
    
    release_file="/etc/os-release"
    if [ -e "${release_file}" ]; then
        version=`grep 'VERSION_ID=' ${release_file} | cut -d'=' -f 2 | tr -d '"'`
        os_name=`grep '^ID=' ${release_file} | cut -d'=' -f2 | tr -d '""'`
        if [ "${os_name}" = "rhel" ] || [ "${os_name}" = "fedora" ]; then
            os_name="centos"
            exit 0
        elif [ "${os_name}" = "ubuntu" ]; then
            version=`echo ${version} | tr -d '.'`
        elif [ "${os_name}" = "debian" ]; then
            os_name="ubuntu"
            version=`debian_ubuntu_lts_map ${version} | tr -d '.'`
        fi
    else
        # OS doesn't have a os-release!
        echo "Failed to find OS version" >&2
        exit 2
    fi
elif [ "$platform" = "Darwin" ]; then
    dl_cmd="/usr/bin/curl"
    dl_out="-o"
    dl_cmd_opts="--fail" # Returns 22 on error
    dl_quiet="-s"

    os_name="macos"
    version=""
    exit 0
fi
extension="tar.gz"
fsleyes_dist="FSLeyes-latest-${os_name}${version}.${extension}"
fsleyes_url="https://fsl.fmrib.ox.ac.uk/fsldownloads/fsleyes/${fsleyes_dist}"
fsleyes_tmp=`drop_sudo mktemp -d -t fsleyesXXXX`
fsleyes_archive="${fsleyes_tmp}/${fsleyes_dist}"

if [ $? -ne 0 ]; then
    echo "Failed to create temporary directory" >&2
    exit 2
fi
fsleyes_install_log="${fsleyes_tmp}/fsleyes_installer.log"
fsleyes_root_dir="${fsl_dir}/bin"
fsleyes_dir="${fsleyes_root_dir}/FSLeyes"

if [ ${quiet} -eq 1 ]; then
    dl_cmd_opts="${dl_cmd_opts} ${dl_quiet}"
fi

drop_sudo echo "Downloading FSLeyes archive" >> "${fsleyes_install_log}"
drop_sudo ${dl_cmd} ${dl_out} "${fsleyes_archive}" ${dl_cmd_opts} ${fsleyes_url} 2>> "${fsleyes_install_log}"
status=$?
if [ ${status} -ne 0 ]; then
    echo "Failed to download FSLeyes - see ${fsleyes_install_log} for details" >&2
    exit ${status}
fi

###################
# Install FSLeyes
###################
drop_sudo echo "Installing FSLeyes into ${fsleyes_root_dir}" >> "${fsleyes_install_log}"
if [ -d "${fsleyes_dir}" ]; then
    drop_sudo echo "Removing old FSLeyes" >> "${fsleyes_install_log}"
    /bin/rm -rf "${fsleyes_dir}"
fi
nfiles=`drop_sudo /bin/tar -ztvf "${fsleyes_archive}" | /usr/bin/wc -l`
if [ -z "$nfiles" ]; then
    nfiles=3500
fi

/bin/tar -zxvf "${fsleyes_archive}" -C "${fsleyes_root_dir}" \
    2>> "${fsleyes_install_log}" | \
    ${script_dir}/progress.sh ${nfiles} ${quiet} >> "${fsleyes_install_log}"
 
if [ $? -ne 0 ]; then
    echo "Failed to install FSLeyes - see ${fsleyes_install_log} for details" >&2
    exit 3
fi
drop_sudo rm -f "${fsleyes_archive}"

drop_sudo rm -f "${fsleyes_install_log}"
drop_sudo rmdir "${fsleyes_tmp}"