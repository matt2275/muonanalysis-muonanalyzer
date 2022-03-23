ls -d $CMSSW_BASE/src/.git  $CMSSW_BASE/tmp
rm -rf $CMSSW_BASE/src/.git 
rm -rf $CMSSW_BASE/tmp
cd $CMSSW_BASE
cd ../
tar czvf $CI_PROJECT_DIR/build.tar.gz $CMSSW_VERSION
