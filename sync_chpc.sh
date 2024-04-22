
# $CHPC is an enivornment variable with login address
LOGIN=$CHPC2
CHPC_PATH='/mnt/beegfs/scratch/tom.earnest/atn_cognition'

# reminder
echo ""
echo "!! REMINDER !!"
echo "    Make sure the CHPC repo is up to date (git pull)!"
echo ""

# upload
# only the input data and processed tables are synced
# rsync -av inputs $LOGIN:${CHPC_PATH}
rsync -av "main/inputs/" "$LOGIN:${CHPC_PATH}/main/inputs"
rsync -av "main/datasets/" "$LOGIN:${CHPC_PATH}/main/datasets"

# download
# outputs are downloaded
rsync -av "$LOGIN:${CHPC_PATH}/main/outputs" 'main'