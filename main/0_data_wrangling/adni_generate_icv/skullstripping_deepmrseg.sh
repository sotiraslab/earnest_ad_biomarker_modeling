#!/bin/sh

# Script for running skull stripping using the CBICA DeepMRSeg

# https://github.com/CBICA/DeepMRSeg/tree/main/DeepMRSeg

##########################################################

##########################################################
# HELPER FUNCTIONS
##########################################################

# add to JSON
add_to_json_usage() {
	echo ""
	echo "Usage:"
	echo "add_to_json [-i <json_file>] -j <json_string> -o <output_file> [-l level] "
	echo ""
}


add_to_json() {
	arginput=""
	arglevel="."
	local OPTIND o a
	while getopts hi:j:o:l: arg
	do
		case $arg in
		h)	handle_fast_output_usage
			exit 0;;
		i) 	arginput=${OPTARG};;
		j) 	argjson=${OPTARG};;
		o)  argoutput=${OPTARG};;
		l)  arglevel=${OPTARG};;
		?)	echo >&2 ""
			echo >&2 "Unknown arguments passed; exiting."
			echo >&2 ""
			add_to_json_usage
			exit 1;;
		esac
	done
	shift $((OPTIND-1))
	
	if [[ ! -f $arginput ]]
	then
		arginput=""
		argnull='--null-input'
	fi
	
	tmpdir=$(mktemp -d)
	tmpjson=$tmpdir/temp.json
	jq \
		$argnull \
		"$arglevel += $argjson" \
		$arginput \
		> $tmpjson

	mv $tmpjson $argoutput
	rm -rf $tmpdir
}


##########################################################
# PARSE ARGUMENTS
##########################################################

# Define help text for the script.
usage() {
	echo "Usage:	skullstripping_deepmrseg.sh -i <t1path> -m <brainmask> [-b <brain>] "
	echo ""
	echo "Description:"
	echo "    Skullstrip a T1 image using DeepMRSeg, a deep-learning based approach."
	echo "    A T1 image is taken as an input (-i), and a skull stripped image (-b) "
	echo "    and / or brain mask (-m) are produced."
	echo ""
	echo "Required Arguments:"
	echo "    -i <path>    path to an input T1 image"
	echo "    -m <path>    where to output the binary brain mask"
	echo ""
	echo "Output Arguments:"
	echo "    -b <path>    where to output skullstripped T1 image"
	echo "    -C           For CHPC useage; loads the default softwares/libraries."
	echo ""
	}

# read arguments
CHPCLOAD=0
while getopts hi:b:m:C arg
do
	case $arg in
	h)	usage
		exit 0;;
	i)	t1path=${OPTARG};;
	b)	brain=${OPTARG};;
	m)  brainmask=${OPTARG};;
	C)  CHPCLOAD=1;;
	?)	echo ""
		echo "Unknown arguments passed; exiting."
		echo ""
		usage;
		exit 1;;
	esac
done

##########################################################
# CHECK INPUT
##########################################################

if [[ -z $t1path || -z $brainmask ]]
then
	echo ""
	echo "ERROR: An input image (-i) and an output path (-m) are both required.  Exiting."
	echo ""
	usage;
	exit 1
fi

if [[ ! -f $t1path ]]
then
	echo ""
	echo "ERROR: Cannot find input (-i) '$t1path'.  Exiting."
	echo ""
	exit 1
fi

##########################################################
# LOAD REQUIREMENTS
##########################################################

if [[ $CHPCLOAD == '1' ]]
then
	module load python > /dev/null 2>&1
	module load cuda/10.1 > /dev/null 2>&1
 	module load cudnn/8.1.1 > /dev/null 2>&1
	source activate DeepMRSeg > /dev/null 2>&1
fi

##########################################################
# CHECK DEPENDENCIES
##########################################################

SOFTWARE_NEEDED="deepmrseg_apply"

for tool in $SOFTWARE_NEEDED
do
	found=`command -v $tool`
	if [[ -z $found ]]
	then
		echo ""
		echo "Error; cannot find at least one required software dependency: $tool"
		echo "Exiting."
		echo ""
		exit 1
	fi
done

# must have pre-trained model!

if [[ ! -d $HOME/.deepmrseg/trained_models/dlicv/ ]]
then
	echo ""
	echo "Error; cannot find pre-trained 'dlicv' model for DeepMRSeg."
	echo "You can download with the command 'deepmrseg_downloadmodel --model dlicv'."
	echo "See here for more information: https://github.com/CBICA/DeepMRSeg"
	echo "Exiting."
	echo ""
	exit 1
fi

##########################################################
# PRINT STATUS
##########################################################

echo ""
echo "--------------------"
echo "Deep MR Segmentation"
echo "--------------------"
echo ""

echo "Input: $t1path"
echo "Output brain mask: $brainmask"
echo "Output brain: $brain"

# find original jsonpath
t1name=`basename $t1path`
jsonpath=`dirname $t1path`/`basename $t1name .nii.gz`.json
if [[ ! -f $jsonpath ]]
then
	echo "WARNING: No JSON sidecar file found for $t1path; one will be created."
fi

##########################################################
# MAIN
##########################################################

DEEPMRSEG_COMMAND=(deepmrseg_apply --task dlicv --inImg $t1path --outImg $brainmask)

echo ""
echo "Running Deep MR Segmentation command:"
echo "     ${DEEPMRSEG_COMMAND[@]}"

# run!
"${DEEPMRSEG_COMMAND[@]}"

# create skullstripped image and save if requested

if [[ ! -z $brain ]]
then
	echo ""
	echo "Creating skull stripped T1 image with DLICV mask."
	
	fslmaths $brainmask -mul $t1path $brain
	
	# augment JSON
	brain_json=`dirname "$brain"`/`basename "$brain" .nii.gz`.json
	json="{\"SkullStripped\":{\"Program\": \"DeepMRSeg\", \"Input\": \"$t1path\", \"Date\": \"$(date +%Y-%m-%d)\"}}"
	add_to_json -i $jsonpath -j "$json" -o "$brain_json" -l '.Preprocessing'
fi

#####################################
#####################################
#####################################
#####################################
