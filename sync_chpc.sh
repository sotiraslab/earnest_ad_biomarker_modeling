#!/bin/bash

CHPC_PATH='/mnt/beegfs/scratch/tom.earnest/atn_cognition'

rsync -av 'data' $CHPC:${CHPC_PATH}
rsync -av --exclude="*.py" --exclude='*.R' --exclude='*.sh' $CHPC:${CHPC_PATH}/'main' '.'