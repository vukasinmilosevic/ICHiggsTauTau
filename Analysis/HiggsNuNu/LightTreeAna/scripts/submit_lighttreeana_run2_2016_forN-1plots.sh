#!/bin/sh
DOCERN=0
DOSUBMIT=1

DATE=170228

## Try and take the JOBWRAPPER and JOBSUBMIT commands
## from the environment if set, otherwise use these defaults
: ${JOBWRAPPER:="./scripts/generate_job.sh"}
: ${JOBSUBMIT:="eval"}

GRIDSETUP=1
if [ "$DOCERN" = "0" ]
  then
  JOBSCRIPT="./scripts/submit_ic_batch_job.sh"
else
  JOBSCRIPT="./scripts/submit_cern_batch_job.sh"
  GRIDSETUP=0
fi
export JOBSUBMIT=$JOBSCRIPT" "$JOBQUEUE


echo "Using job-wrapper: " $JOBWRAPPER
echo "Using job-submission: " $JOBSUBMIT

CONFIG=scripts/DefaultRun2Config_vetos_forN-1plots.cfg

QUEUEDIR=short #medium long
# no_METno_cut no_dijet_dphi_cut no_dijet_deta_cut no_alljetsmetno_mindphi_cut no_dijet_M_cut
JOBDIRPREFIX=jobs_run2ana_${DATE}_forN-1plots/no_dijet_M_cut
JOBDIR=$JOBDIRPREFIX/
OUTPUTPREFIX=output_run2ana_${DATE}_forN-1plots/no_dijet_M_cut
OUTPUTDIR=$OUTPUTPREFIX/

OUTPUTNAME="output.root"

echo "Config file: $CONFIG"
mkdir -p $JOBDIR
mkdir -p $OUTPUTDIR

if [ "$DOCERN" = "0" ]
  then
  if [ "$QUEUEDIR" = "medium" ]
    then
    JOBQUEUE="5:59:0"
  elif [ "$QUEUEDIR" = "long" ]
    then
    JOBQUEUE="47:59:0"
  else
    JOBQUEUE="2:59:0"
  fi
else
  if [ "$QUEUEDIR" = "medium" ]
    then
    JOBQUEUE="1nd"
  elif [ "$QUEUEDIR" = "long" ]
    then
    JOBQUEUE="2nd"
  else
    JOBQUEUE="1nh"
  fi
fi
export JOBSUBMIT=$JOBSCRIPT" "$JOBQUEUE
echo "Using job-submission: " $JOBSUBMIT

echo "JOB name = $JOB"
for syst in ""
do
  mkdir -p $JOBDIR$syst
  mkdir -p $OUTPUTDIR$syst
  for channels in enu munu taunu ee mumu qcd nunu
    do
    JOB=$channels
    #executable expect strings separated by "!"

    ## no_METno_cut
#     HISTSTRING=";E_{T,no-#mu}^{miss} (GeV);Events"
#     SHAPESTRING="metnomuons(50,0.,500.)"
#     if [ "$channels" = "ee" ]; then
#       HISTSTRING=";E_{T,no-el}^{miss} (GeV);Events"
#       SHAPESTRING="metnoelectrons(50,0.,500.)"
#     fi
#     if [ "$channels" = "enu" ]; then
#       HISTSTRING=";E_{T,no-el}^{miss} (GeV);Events"
#       SHAPESTRING="metnoelectrons(50,0.,500.)"
#     fi

    ## no_dijet_dphi_cut
#     HISTSTRING=";#Delta#phi_{jj};Events"
#     SHAPESTRING="dijet_dphi(20,0.,3.1416)"

    ## no_dijet_deta_cut
#     HISTSTRING=";#Delta#eta_{jj};Events"
#     SHAPESTRING="dijet_deta(45,0.,6.5)"

    ## no_alljetsmetnomu_mindphi_cut
#     HISTSTRING=";#Delta#phi(E_{T,no-#mu}^{miss},j);Events"
#     SHAPESTRING="alljetsmetnomu_mindphi(60,0.,3.1416)"
#     if [ "$channels" = "ee" ]; then
#       HISTSTRING=";#Delta#phi(E_{T,no-el}^{miss},j);Events"
#       SHAPESTRING="alljetsmetnoel_mindphi(60,0.,3.1416)"
#     fi
#     if [ "$channels" = "enu" ]; then
#       HISTSTRING=";#Delta#phi(E_{T,no-el}^{miss},j);Events"
#       SHAPESTRING="alljetsmetnoel_mindphi(60,0.,3.1416)"
#     fi

    ## no_dijet_M_cut
    HISTSTRING=";M_{jj} (GeV);Events"
    SHAPESTRING="dijet_M(22,600.,3500.)"

    echo "Making histograms: " $SHAPESTRING
    OUTPUTNAME="$channels.root"
    MINDPHICUT="alljetsmetnomu_mindphi\>=0.5"
    if [ "$channels" = "taunu" ]; then
	############MINDPHICUT="jetmetnomu_mindphi\>=1.0" #\&\&alljetsmetnomu_mindphi\<2.3"
	#MINDPHICUT="jetmetnomu_mindphi\>=1.0\&\&alljetsmetnomu_mindphi\<2.3"
      MINDPHICUT="alljetsmetnomu_mindphi\>=0.5"
    fi
    if [ "$channels" = "qcd" ]; then
      MINDPHICUT="alljetsmetnomu_mindphi\<0.5"
    fi
    if [ "$channels" = "ee" ]; then
      MINDPHICUT="alljetsmetnoel_mindphi\>=0.5"
    fi
    if [ "$channels" = "enu" ]; then
      MINDPHICUT="alljetsmetnoel_mindphi\>=0.5"
    fi
    if [ "$syst" = "" ]
      then
      $JOBWRAPPER "./bin/LTAnalysisRun2_2016 --cfg=$CONFIG --channel=$channels --histTitlePar='$HISTSTRING' --shapePar='$SHAPESTRING' -o $OUTPUTDIR$syst/$OUTPUTNAME --jetmetdphicut=$MINDPHICUT &> $JOBDIR$syst/$JOB.log" $JOBDIR$syst/$JOB.sh $GRIDSETUP
    else

      $JOBWRAPPER "./bin/LTAnalysisRun2_2016 --cfg=$CONFIG --channel=$channels --histTitlePar='$HISTSTRING' --shapePar='$SHAPESTRING' --syst=$syst -o $OUTPUTDIR$syst/$OUTPUTNAME --jetmetdphicut=$MINDPHICUT &> $JOBDIR$syst/$JOB.log" $JOBDIR$syst/$JOB.sh $GRIDSETUP

    fi
    if [ "$DOSUBMIT" = "1" ]; then 
      $JOBSUBMIT $JOBDIR$syst/$JOB.sh
    else 
      echo "$JOBSUBMIT $JOBDIR$syst/$JOB.sh"
    fi
  done
done
