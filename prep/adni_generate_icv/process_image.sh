#!/bin/bash

##########################################################
# PARSE ARGUMENTS
##########################################################

# Define help text for the script.
usage() {
	echo ""
	echo "Usage:	process_image.sh -i <input> -o <ouput>"
	echo ""
	echo "Description:"
	echo "    Apply bias correction & skull stripping.  Images are output in same folder as original image."
	echo "Arguments:"
	echo "    -i <path>    path to an input T1 image"
	echo "    -P <name>    Name to give preskullstripping output (NO EXTENSION!).  Default='t1_corrected'"
    echo "    -B <name>    Name to give brain image (NO EXTENSION!).  Default='brain'"
    echo "    -M <name>    Name to give brain mask image (NO EXTENSION!).  Default='brain_mask'"
    echo "    -T <path>    Path to TTP2 code.  Defaults to /home/tom.earnest/scripts/tom_t1_pipeline2"
	echo ""
	}

# defaults
PSS_NAME='t1_corrected'
BRAIN_NAME='brain'
BRAINMASK_NAME='brain_mask'
TTP2='/home/tom.earnest/scripts/tom_t1_pipeline2'

# read arguments
while getopts hi:P:B:M:T: arg
do
	case $arg in
	h)	usage
		exit 0;;
	i)	T1PATH=$(realpath ${OPTARG});;
	P)	PSS_NAME=${OPTARG};;
    B)  BRAIN_NAME=${OPTARG};;
    M)  BRAINMASK_NAME=${OPTARG};;
    T)  TTP2=${OPTARG};;
	?)	echo ""
		echo "Unknown arguments passed; exiting."
		echo ""
		usage;
		exit 1;;
	esac
done

PSS_SCRIPT="${TTP2}/preskullstripping.sh"
SS_SCRIPT="${TTP2}/skullstripping_deepmrseg.sh"

##########################################################
# CHECKS
##########################################################

if [[ -z $T1PATH ]]
then
    echo ""
    echo "ERROR: T1 image must be provided (-i)."
    exit 1
fi

if [[ ! -d $TTP2 ]]
then
    echo ""
    echo "ERROR: Cannot find '$TTP2'.  Download code from here: https://github.com/sotiraslab/tom_t1_pipeline2/tree/main"
    exit 1
fi

if [[ ! -f $PSS_SCRIPT || ! -f $SS_SCRIPT ]]
then
    echo ""
    echo "ERROR: Cannot find either preskullstripping.sh or skullstripping_deepmrseg.sh in the TTP2 folder."
    exit 1
fi

##########################################################
# DEPENDENCIES
##########################################################

# this could be more contditional on what tools are being used
# but keeping it more general for now

module load afni > /dev/null 2>&1
module load ants > /dev/null 2>&1
module load fsl > /dev/null 2>&1

if [[ -z $(command -v afni) ]]
then
	echo ""
	echo "ERROR: Cannot find required dependency, AFNI.  Exiting."
	echo ""
	exit 1
else
	echo  "- Loaded dependency AFNI, version: $(afni --version)"
fi

if [[ -z $(command -v ANTS) ]]
then
	echo ""
	echo "ERROR: Cannot find required dependency, ANTS.  Exiting."
	echo ""
	exit 1
else
	echo  "- Loaded dependency ANTS."
fi

if [[ -z $(command -v fsl) ]]
then
	echo ""
	echo "ERROR: Cannot find required dependency, FSL.  Exiting."
	echo ""
	exit 1
else
	echo  "- Loaded dependency FSL, version: `cat $FSLDIR/etc/fslversion`"
fi


##########################################################
# SETUP OUTPUTS
##########################################################

OUTFOLDER=$(dirname $T1PATH)
PSS_PATH=${OUTFOLDER}/${PSS_NAME}.nii.gz
BRAIN_PATH=${OUTFOLDER}/${BRAIN_NAME}.nii.gz
BRAINMASK_PATH=${OUTFOLDER}/${BRAINMASK_NAME}.nii.gz

##########################################################
# PRINT STATUS
##########################################################

echo ''
echo "-------------------"
echo 'Skull Stripping'
echo '-------------------'
echo ""

echo "Input Image: $T1PATH"
echo "Output w/ bias correction: $PSS_PATH"
echo "Output brain: $BRAIN_PATH"
echo "Output brain mask: $BRAINMASK_PATH"

##########################################################
# MAIN
##########################################################

PSS_COMMAND=(${PSS_SCRIPT} -i ${T1PATH} -o ${PSS_PATH})
SS_COMMAND=(${SS_SCRIPT} -i ${PSS_PATH} -m ${BRAINMASK_PATH} -b ${BRAIN_PATH})

# preskullstripping
echo ""
echo "Running bias correction and reorientation command:"
echo "     ${PSS_COMMAND[@]}"
"${PSS_COMMAND[@]}"

# skullstripping
echo ""
echo "Running Deep MR Segmentation command:"
echo "     ${SS_COMMAND[@]}"
"${SS_COMMAND[@]}"

echo ""