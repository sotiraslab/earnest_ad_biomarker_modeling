#!/bin/bash

# variables
INPUT_DIRECTORY='REPLACE ME'
TTP2='.'
SCRIPT='process_image.sh'

PATTERN='ADNI*.nii.gz'
FORCE='1'

# setup
module load afni
module load ants
module load fsl

# get image list
IMAGES=$(find ${INPUT_DIRECTORY} -name "$PATTERN")
IMAGE_COUNT=$(echo ${IMAGES} | wc -w)

# check
echo ""
echo "IMAGE COUNT: ${IMAGE_COUNT}; Proceed? [y/n]"
read ANS

if [[ ! $ANS == 'y' ]]
then
    echo ""
    echo "Exiting."
    echo ""
    exit 0
fi

# main loop
for IMAGE in $IMAGES
do
    NAME=$(basename ${IMAGE})
    DIRECTORY=$(dirname ${IMAGE})
    LOG=${DIRECTORY}/slurm.log
    SLURM_COMMAND=(sbatch -n 1 -N 1 --time='24:00:00' -o ${LOG} -J ${NAME} --mem=32GB)
    PROC_COMMAND=(${SCRIPT} -i ${IMAGE})
    FULL_COMMAND=("${SLURM_COMMAND[@]}" "${PROC_COMMAND[@]}")

    # report
    echo ""
    echo "Submitting command:"
    echo "    ${FULL_COMMAND[@]}"

    # run
    "${FULL_COMMAND[@]}"
    sleep 0.5
done

echo ""