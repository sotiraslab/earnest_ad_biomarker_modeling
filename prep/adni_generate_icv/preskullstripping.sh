#!/bin/bash

# Script to prepare a T1 image for skullstripping.  Includes bias correction & reorientation.

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
	echo ""
	echo "Usage:	preskullstripping.sh -i <input> -o <ouput>"
	echo ""
	echo "Description:"
	echo "    Apply preprocessing steps to prepare a T1 image for skullstripping.  Current steps are: "
	echo "        - reorientation with AFNI 3dresample"
	echo "        - bias correction with N4"
	echo "Arguments:"
	echo "    -i <path>    path to an input T1 image"
	echo "    -o <path>    output path for image that has all the preprocessing steps applied"
	echo ""
	}

# read arguments
while getopts hi:o: arg
do
	case $arg in
	h)	usage
		exit 0;;
	i)	t1path=${OPTARG};;
	o)	output=${OPTARG};;
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

if [[ -z $t1path ]] || [[ -z $output ]]
then
	echo ""
	echo "Both input (-i) and output (-o) MUST be set for $0; exiting."
	echo ""
	usage;
	exit 1
fi

if [[ ! -f $t1path ]]
then
	echo ""
	echo "Input file '$t1path' does not exist; exiting."
	echo ""
	usage;
	exit 1
fi

##########################################################
# CHECK DEPENDENCIES
##########################################################

SOFTWARE_NEEDED="3dresample N4BiasFieldCorrection"

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

##########################################################
# PRINT STATUS
##########################################################

echo ''
echo "-------------------"
echo 'Pre Skull Stripping'
echo '-------------------'
echo ""

echo "Input Image: $t1path"
echo "Output Image: $output"

#find the original JSON path
t1name=`basename $t1path`
jsonpath=`dirname $t1path`/`basename $t1name .nii.gz`.json
if [[ ! -f $jsonpath ]]
then
	echo ""
	echo "WARNING: No JSON sidecar file found for $t1path; one will be created."
	echo ""
fi

##########################################################
# TEMPORARY WORKING DIRECTORY
##########################################################

echo ""
echo "Creating temporary working directory for running pre-skullstripping functions."

STARTINGDIR=`pwd`
WORKINGDIR=`mktemp -d /tmp/pre_skullstripping_XXXXXXXX`

echo ""
echo "Moving into temporary working directory at '$WORKINGDIR'."

cd $WORKINGDIR

# cleanup function
cleanup() {
	echo ""
	echo "RETURNING TO STARTING DIRECTORY '$STARTINGDIR'."
	cd $STARTINGDIR
	echo "REMOVING TEMPORARY WORKING DIRECTORY '$WORKINGDIR'."
	rm -rf $WORKINGDIR
	echo ""
	exit
}
trap cleanup EXIT
trap cleanup INT

##########################################################
# REORIENTATION
##########################################################

echo ''
echo 'REORIENTATION'
echo '-------------'
echo ''

og_orient=`3dinfo -orient $t1path` # add into JSON
reoriented=${WORKINGDIR}/reoriented.nii.gz
REORIENT_COMMAND=(3dresample -orient RPI -prefix $reoriented -input $t1path)

echo ""
echo "Original Image Orientation: $og_orient"
echo "Running reorientation to RPI with AFNI 3dresample."
echo "     ${REORIENT_COMMAND[@]}"

# run
"${REORIENT_COMMAND[@]}"

echo "Incorporating orientation information into JSON..."
reorient_json=`dirname "$reoriented"`/`basename "$reoriented" .nii.gz`.json
json="{\"Reoriented\": {\"Program\": \"3dresample\", \"Input\": \"$t1path\", \"Date\": \"$(date +%Y-%m-%d)\", \"OldOrientation\": \"$og_orient\", \"NewOrientation\":\"RPI\"}}"
add_to_json -i $jsonpath -j "$json" -o "$reorient_json" -l ".Preprocessing"
echo "Done."

echo ''
echo 'BIAS CORRECTION'
echo '---------------'
echo ''

BIAS_COMMAND=(N4BiasFieldCorrection -d 3 -i $reoriented -o $output)


echo "Appling ANTS N4 bias correction"
echo "     ${BIAS_COMMAND[@]}"

"${BIAS_COMMAND[@]}"

echo "Incorporating bias-correction into JSON..."
output_json=`dirname "$output"`/`basename "$output" .nii.gz`.json
json="{\"BiasCorrected\": {\"Program\": \"N4BiasFieldCorrection\", \"Input\": \"$t1path\", \"Date\": \"$(date +%Y-%m-%d)\"}}"
add_to_json -i $reorient_json -j "$json" -o "$output_json" -l ".Preprocessing"
echo "Done."

echo ""
