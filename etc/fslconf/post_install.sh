#!/usr/bin/env bash
#
# Post installation script
# Put commands you need to run after the FSL archive has been installed
# in here

# Where is this script?
script_dir=$( cd $(dirname $0) ; pwd)

function syntax {
    echo "post_install.sh [-f <FSLDIR>] [-q]"
    echo "  -f <FSLDIR> Location of installed FSL, e.g. /usr/local/fsl"
    echo "                  if not provided looks for FSLDIR in environment"
    echo "  -q          Install silently"
}

script_list=`"${script_dir}/scripts.manifest"`

OPTIND=1
quiet=''
while getopts "h?qf:" opt; do
    case "${opt}" in
    h|\?)
        syntax
        exit 0
        ;;
    q)  quiet='-q'
        ;;
    f)  fsl_dir=${OPTARG}
        ;;
    esac
done

shift $((OPTIND-1))

[ "$1" = "--" ] && shift

for script in ${script_list}; do
    ${script_dir}/${script} -f "${fsl_dir}" ${quiet}
    if [ $? -ne 0 ]; then
        echo "Post install setup failed"
        exit $?
    fi
done

if [ -z "${quiet}" ]; then
    echo "Post install setup complete"
fi