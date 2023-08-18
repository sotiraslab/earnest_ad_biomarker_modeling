#!/bin/bash

# $CHPC is an enivornment variable with login address
LOGIN=$CHPC
CHPC_PATH='/mnt/beegfs/scratch/tom.earnest/atn_cognition'

UPLOADS='uploads.txt'
DOWNLOADS='downloads.txt'

# upload
rsync -avRr --include-from="${UPLOADS}" --exclude='*' . $LOGIN:${CHPC_PATH}/.

# download
rsync -avRr --files-from="${DOWNLOADS}" "$LOGIN:${CHPC_PATH}/." '.'
