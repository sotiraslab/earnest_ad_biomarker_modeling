
# $CHPC is an enivornment variable with login address
LOGIN=$CHPC2
CHPC_PATH='/home/tom.earnest/projects/atn_cognition'

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
