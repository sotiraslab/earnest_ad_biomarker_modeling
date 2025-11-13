
# $CHPC is an enivornment variable with login address
LOGIN=$CHPC
CHPC_PATH='/scratch/tom.earnest/earnest_ad_biomarker_modeling'

# reminder
echo ""
echo "!! REMINDER !!"
echo "    Make sure the CHPC repo is up to date (git pull)!"
echo ""

# download
# outputs are downloaded
rsync -av "$LOGIN:${CHPC_PATH}/main/outputs" 'main'
rsync -av "$LOGIN:${CHPC_PATH}/main/logs" 'main'
