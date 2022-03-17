#!/bin/bash

for channel in Inv_mass_gghasOneL_DNN_1_HL Inv_mass_gghasOneL_DNN_2_HL Inv_mass_gghasOneL_DNN_3_HL Inv_mass_gghasOneL_DNN_4_HL Inv_mass_gghasTwoL_HL Mgg_c3_DNN_1_HL Mgg_c3_DNN_2_HL Mgg_c4_Zveto_HL
do 	       
 sed -e 's/.*shapes.*/\
shapes data_obs      * roofit_'$channel'.root HHWWgg:data_obs_ \
shapes Continuum_Bkg * roofit_'$channel'.root HHWWgg:mc_ \
shapes GGHH          * roofit_'$channel'.root HHWWgg:mh_GGHH_ \
shapes GGH           * roofit_'$channel'.root HHWWgg:mh_GGH_ \
shapes VBFH          * roofit_'$channel'.root HHWWgg:mh_VBFH_ \
shapes VH            * roofit_'$channel'.root HHWWgg:mh_VH_ \
shapes tHq           * roofit_'$channel'.root HHWWgg:mh_tHq_ \
shapes ttH           * roofit_'$channel'.root HHWWgg:mh_ttH_/' ${channel}.txt > roofit_${channel}.txt

mv roofit_${channel}.txt ${channel}.txt
text2workspace.py ${channel}.txt -o roofit_${channel}_workspace.root
#cp roofit_${channel}_workspace.root /home/ucl/cp3/sjain/bamboodev/DatacardProducer/DataCard2_FH/$1/${channel}.root
#cp ${channel}.txt /home/ucl/cp3/sjain/bamboodev/DatacardProducer/DataCard2_FH/$1/
echo ${channel}
combine -M Significance -m 125 --rMin -10 --rMax 10 -t -1 --expectSignal 1 -d roofit_${channel}_workspace.root --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --cminDefaultMinimizerType Minuit2 --cminDefaultMinimizerStrategy 0 --cminDefaultMinimizerTolerance 0.1 --cminFallbackAlgo Minuit2,0:0.2 --cminFallbackAlgo Minuit2,0:0.4 > sig_${channel}.txt

done

#shapes Continuum_Bkg * roofit_Inv_mass_gghasOneL_DNN_2_HL.root HHWWgg:exp_gghasOneL_DNN_2_HL/' Inv_mass_gghasOneL_DNN_HL.txt
#shapes GGHH          * roofit_Inv_mass_gghasOneL_DNN_2_HL.root HHWWgg:mh_GGHH_gghasOneL_DNN_2_HL/' Inv_mass_gghasOneL_DNN_HL.txt 
#shapes GGH           * roofit_Inv_mass_gghasOneL_DNN_2_HL.root HHWWgg:mh_GGH_gghasOneL_DNN_2_HL \
#shapes VBFH          * roofit_Inv_mass_gghasOneL_DNN_2_HL.root HHWWgg:mh_VBFH_gghasOneL_DNN_2_HL \ 
#shapes VH            * roofit_Inv_mass_gghasOneL_DNN_2_HL.root HHWWgg:mh_VH_gghasOneL_DNN_2_HL \
#shapes tHq           * roofit_Inv_mass_gghasOneL_DNN_2_HL.root HHWWgg:mh_tHq_gghasOneL_DNN_2_HL \
#shapes ttH           * roofit_Inv_mass_gghasOneL_DNN_2_HL.root HHWWgg:mh_ttH_gghasOneL_DNN_2_HL/' Inv_mass_gghasOneL_DNN_HL.txt 


