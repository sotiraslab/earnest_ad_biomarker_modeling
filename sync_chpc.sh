
# $CHPC is an enivornment variable with login address
LOGIN=$CHPC2
CHPC_PATH='/mnt/beegfs/scratch/tom.earnest/atn_cognition'

# upload
# only the input data and processed tables are synced
rsync -av inputs $LOGIN:${CHPC_PATH}
rsync -av "outputs/maindata" "$LOGIN:${CHPC_PATH}/outputs"

# download
# outputs are downloaded
rsync -av "$LOGIN:${CHPC_PATH}/outputs" '.'