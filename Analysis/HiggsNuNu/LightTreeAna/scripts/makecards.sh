#!/bin/bash

if [ "$#" -ne "2" ]; then
  echo "Usage: $0 <doSubmit> <do4params>"
  exit 0
fi


DOSUBMIT=$1
DO4PARAMS=$2
infolder=output_run2ana_170222/
outfolder=cards_run2ana_170222/
do_tau_veto_unc=false
wzqcd_syst=1.15
wzewk_syst=1.15
blind=true
#zvvstat=18
mkdir -p $outfolder

extraoptions="--do_ues=false" #--do_ggh=false --do_separate_qcdewk=false"

for channel in enu #munu taunu mumu ee qcd nunu
  do
  echo " ********************************"
  echo " *** Processing channel $channel"
  echo " ********************************"

  mkdir -p $outfolder/$channel

  echo "channel $channel"

  HistToIntegrate="jet1_pt"

  echo "channel $channel, HistToIntegrate: $HistToIntegrate"

  OUTNAME=$outfolder/$channel/vbfhinv_${channel}_13TeV.txt
  if (( "$DO4PARAMS" == "1" )); then
    OUTNAME=$outfolder/$channel/vbfhinv_${channel}_13TeV_4params.txt
  fi
  if (( "$DOSUBMIT" == "0" )); then
    echo "./bin/makeCountingCard -i $infolder --blind=$blind -o $OUTNAME -m 125 --channel $channel --do_latex true --do_datatop false --zvvstat 0 --qcdrate 0 --mcBkgOnly=true --do_run2=true --do_4params=$DO4PARAMS --histoToIntegrate=$HistToIntegrate --do_tau_veto_unc=$do_tau_veto_unc --wzqcd_syst=$wzqcd_syst --wzewk_syst=$wzewk_syst $extraoptions | tee $outfolder/$channel/card.log"
  else
    ./bin/makeCountingCard -i $infolder --blind=$blind -o $OUTNAME -m 125 --channel $channel --do_latex true --do_datatop false --zvvstat 0 --qcdrate 0 --mcBkgOnly=true --do_run2=true --do_4params=$DO4PARAMS --histoToIntegrate=$HistToIntegrate --do_tau_veto_unc=$do_tau_veto_unc --wzqcd_syst=$wzqcd_syst --wzewk_syst=$wzewk_syst $extraoptions | tee $outfolder/$channel/card_${channel}.log
  fi
done
