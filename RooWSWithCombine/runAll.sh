#!/bin/bash
### set this up in CMSSW combine setup###

###clean the area under your dir to have a clean run##
rm -f *.root *.txt *.pdf
### copy the input rootfiles and txt files to combine that were produced after running Florian's produceDataCards.py 
cp /home/ucl/cp3/$3/bamboodev/DatacardProducer/DataCard2_FH/$1/*.txt .
cp /home/ucl/cp3/$3/bamboodev/DatacardProducer/DataCard2_FH/$1/*.root .

#### make the histograms in rootworkspace format, have the fits done and update the root files ### 
python exampleHHWWgg.py

####### these root files are updated in combine datacards produced from florian's script #####  
source produceDataCards.sh $2

### manual combination of various datacards here ### 

combineCards.py bin1=Inv_mass_gghasOneL_DNN_1_HL.txt bin2=Inv_mass_gghasOneL_DNN_2_HL.txt bin3=Inv_mass_gghasOneL_DNN_3_HL.txt bin4=Inv_mass_gghasOneL_DNN_4_HL.txt  >& wwgg_oneL.txt
combineCards.py bin1=wwgg_oneL.txt bin2=Inv_mass_gghasTwoL_HL.txt  >& wwgg.txt
combineCards.py bin1=Mgg_c3_DNN_1_HL.txt bin2=Mgg_c3_DNN_2_HL.txt  >& ttgg_onetau.txt
combineCards.py bin1=ttgg_onetau.txt bin2=Mgg_c4_Zveto_HL.txt  >& ttgg.txt
combineCards.py bin1=Inv_mass_gghasOneL_DNN_1_HL.txt bin2=Inv_mass_gghasOneL_DNN_2_HL.txt bin3=Inv_mass_gghasOneL_DNN_3_HL.txt bin4=Inv_mass_gghasOneL_DNN_4_HL.txt  bin5=Inv_mass_gghasTwoL_HL.txt bin6=Mgg_c3_DNN_1_HL.txt bin7=Mgg_c3_DNN_2_HL.txt bin8=Mgg_c4_Zveto_HL.txt >& combined.txt



### creating workspace for various combinations
text2workspace.py wwgg_oneL.txt -o wwgg_oneL.root
text2workspace.py wwgg.txt -o wwgg.root
text2workspace.py ttgg_onetau.txt -o ttgg_onetau.root
text2workspace.py ttgg.txt -o ttgg.root
text2workspace.py combined.txt -o combined.root

echo "running Results from combination"
combine -M Significance -m 125 --rMin -10 --rMax 10 -t -1 --expectSignal 1 -d wwgg_oneL.root --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --cminDefaultMinimizerType Minuit2 --cminDefaultMinimizerStrategy 0 --cminDefaultMinimizerTolerance 0.1 --cminFallbackAlgo Minuit2,0:0.2 --cminFallbackAlgo Minuit2,0:0.4 > sig_wwgg_oneL.txt
combine -M Significance -m 125 --rMin -10 --rMax 10 -t -1 --expectSignal 1 -d wwgg.root --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --cminDefaultMinimizerType Minuit2 --cminDefaultMinimizerStrategy 0 --cminDefaultMinimizerTolerance 0.1 --cminFallbackAlgo Minuit2,0:0.2 --cminFallbackAlgo Minuit2,0:0.4 > sig_wwgg.txt
combine -M Significance -m 125 --rMin -10 --rMax 10 -t -1 --expectSignal 1 -d ttgg_onetau.root --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --cminDefaultMinimizerType Minuit2 --cminDefaultMinimizerStrategy 0 --cminDefaultMinimizerTolerance 0.1 --cminFallbackAlgo Minuit2,0:0.2 --cminFallbackAlgo Minuit2,0:0.4 > sig_ttgg_onetau.txt
combine -M Significance -m 125 --rMin -10 --rMax 10 -t -1 --expectSignal 1 -d ttgg.root --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --cminDefaultMinimizerType Minuit2 --cminDefaultMinimizerStrategy 0 --cminDefaultMinimizerTolerance 0.1 --cminFallbackAlgo Minuit2,0:0.2 --cminFallbackAlgo Minuit2,0:0.4 > sig_ttgg.txt
combine -M Significance -m 125 --rMin -10 --rMax 10 -t -1 --expectSignal 1 -d combined.root --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --cminDefaultMinimizerType Minuit2 --cminDefaultMinimizerStrategy 0 --cminDefaultMinimizerTolerance 0.1 --cminFallbackAlgo Minuit2,0:0.2 --cminFallbackAlgo Minuit2,0:0.4 > sig_combined.txt

./getSignificance_newsetup.sh

mkdir $2
mv *.root *.txt *.pdf $2/
cp *.sh *.py $2/
