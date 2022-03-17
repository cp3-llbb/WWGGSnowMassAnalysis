#!/bin/bash 
slash='\\\'
single='\'
prefit_stat=ResultswithNoFit_4GeV
prefit_exp=ResultswithFit_4GeV
prefit_exp1=RemovalResultsWithNoFit_4GeV
prefit=RemovalResultsWithFit_4GeV
echo "##############################################"
echo "WWGG SL categories"
echo "##############################################"
cat $prefit_stat/sig_Inv_mass_gghasOneL_DNN_1_HL.txt | grep "Significance:" | awk -v ORS=" " '{printf "Category 1 & %.4f & ", $2}'
cat $prefit_exp/sig_Inv_mass_gghasOneL_DNN_1_HL.txt| grep "Significance:" | awk -v ORS=" " '{printf "%.4f & ", $2}'
cat $prefit_exp1/sig_Inv_mass_gghasOneL_DNN_1_HL.txt| grep "Significance:" | awk -v ORS=" " '{printf "%.4f & ", $2}'
cat $prefit/sig_Inv_mass_gghasOneL_DNN_1_HL.txt | grep "Significance:" | awk -v slash="$slash" '{printf "%.4f %s \n", $2, slash}'
cat $prefit_stat/sig_Inv_mass_gghasOneL_DNN_2_HL.txt | grep "Significance:" | awk -v ORS=" " '{printf "Category 2 &%.4f & ", $2}'
cat $prefit_exp/sig_Inv_mass_gghasOneL_DNN_2_HL.txt | grep "Significance:" | awk -v ORS=" " '{printf "%.4f & ", $2}'
cat $prefit_exp1/sig_Inv_mass_gghasOneL_DNN_2_HL.txt | grep "Significance:" | awk -v ORS=" " '{printf "%.4f & ", $2}'
cat $prefit/sig_Inv_mass_gghasOneL_DNN_2_HL.txt | grep "Significance:" | awk -v slash="$slash" '{printf "%.4f %s \n", $2, slash}'
cat $prefit_stat/sig_Inv_mass_gghasOneL_DNN_3_HL.txt | grep "Significance:" | awk -v ORS=" " '{printf "Category 3 & %.4f & ", $2}'
cat $prefit_exp/sig_Inv_mass_gghasOneL_DNN_3_HL.txt | grep "Significance:" | awk -v ORS=" " '{printf "%.4f & ", $2}'
cat $prefit_exp1/sig_Inv_mass_gghasOneL_DNN_3_HL.txt | grep "Significance:" | awk -v ORS=" " '{printf "%.4f & ", $2}'
cat $prefit/sig_Inv_mass_gghasOneL_DNN_3_HL.txt | grep "Significance:" | awk -v slash="$slash" '{printf "%.4f %s \n", $2, slash}'
cat $prefit_stat/sig_Inv_mass_gghasOneL_DNN_4_HL.txt | grep "Significance:" | awk -v ORS=" " '{printf "Category 4 & %.4f & ", $2}'
cat $prefit_exp/sig_Inv_mass_gghasOneL_DNN_4_HL.txt | grep "Significance:" | awk -v ORS=" " '{printf "%.4f & ", $2}'
cat $prefit_exp1/sig_Inv_mass_gghasOneL_DNN_4_HL.txt | grep "Significance:" | awk -v ORS=" " '{printf "%.4f & ", $2}'
cat $prefit/sig_Inv_mass_gghasOneL_DNN_4_HL.txt | grep "Significance:" | awk -v slash="$slash" '{printf "%.4f %s \n", $2, slash}'
echo "\hline"
cat $prefit_stat/sig_wwgg_oneL.txt| grep "Significance:" | awk -v ORS=" " '{printf "One Leptonic & %.4f & ", $2}'
cat $prefit_exp/sig_wwgg_oneL.txt| grep "Significance:" | awk -v ORS=" " '{printf "%.4f & ", $2}'
cat $prefit_exp1/sig_wwgg_oneL.txt| grep "Significance:" | awk -v ORS=" " '{printf "%.4f & ", $2}'
cat $prefit/sig_wwgg_oneL.txt| grep "Significance:" | awk -v slash="$slash" '{printf "%.4f %s \n", $2, slash}'
echo "\hline"
echo "\hline"
cat $prefit_stat/sig_Inv_mass_gghasTwoL_HL.txt | grep "Significance:" | awk -v ORS=" " '{printf "Two Leptonic & %.4f & ", $2}'
cat $prefit_exp/sig_Inv_mass_gghasTwoL_HL.txt | grep "Significance:" | awk -v ORS=" " '{printf "%.4f & ", $2}'
cat $prefit_exp1/sig_Inv_mass_gghasTwoL_HL.txt | grep "Significance:" | awk -v ORS=" " '{printf "%.4f & ", $2}'
cat $prefit/sig_Inv_mass_gghasTwoL_HL.txt | grep "Significance:" | awk -v slash="$slash" '{printf "%.4f %s \n", $2, slash}'
cat $prefit_stat/sig_wwgg.txt | grep "Significance:" | awk -v ORS=" " '{printf "Combination & %.4f & ", $2}'
cat $prefit_exp/sig_wwgg.txt | grep "Significance:" | awk -v ORS=" " '{printf "%.4f & ", $2}'
cat $prefit_exp1/sig_wwgg.txt | grep "Significance:" | awk -v ORS=" " '{printf "%.4f & ", $2}'
cat $prefit/sig_wwgg.txt | grep "Significance:" | awk -v slash="$slash" '{printf "%.4f %s \n", $2, slash}'
echo "\hline"
echo "##############################################"
echo "tau categories"
echo "##############################################"
cat $prefit_stat/sig_Mgg_c3_DNN_1_HL.txt | grep "Significance:" | awk -v ORS=" " '{printf "Category 1 & %.4f & ", $2}'
cat $prefit_exp/sig_Mgg_c3_DNN_1_HL.txt | grep "Significance:" | awk -v ORS=" " '{printf "%.4f & ", $2}'
cat $prefit_exp1/sig_Mgg_c3_DNN_1_HL.txt | grep "Significance:" | awk -v ORS=" " '{printf "%.4f & ", $2}'
cat $prefit/sig_Mgg_c3_DNN_1_HL.txt | grep "Significance:" | awk -v slash="$slash" '{printf "%.4f %s \n", $2, slash}'
cat $prefit_stat/sig_Mgg_c3_DNN_2_HL.txt | grep "Significance:" | awk -v ORS=" " '{printf "Category 2 &%.4f & ", $2}'
cat $prefit_exp/sig_Mgg_c3_DNN_2_HL.txt | grep "Significance:" | awk -v ORS=" " '{printf "%.4f & ", $2}'
cat $prefit_exp1/sig_Mgg_c3_DNN_2_HL.txt | grep "Significance:" | awk -v ORS=" " '{printf "%.4f & ", $2}'
cat $prefit/sig_Mgg_c3_DNN_2_HL.txt | grep "Significance:" | awk -v slash="$slash" '{printf "%.4f %s \n", $2, slash}'
echo "\hline"
cat $prefit_stat/sig_ttgg_onetau.txt | grep "Significance:" | awk -v ORS=" " -v sslash="${single}" '{printf "1 $%stau$ & %.4f & ", sslash, $2}'
cat $prefit_exp/sig_ttgg_onetau.txt | grep "Significance:" | awk -v ORS=" " '{printf "%.4f & ", $2}'
cat $prefit_exp1/sig_ttgg_onetau.txt | grep "Significance:" | awk -v ORS=" " '{printf "%.4f & ", $2}'
cat $prefit/sig_ttgg_onetau.txt | grep "Significance:" | awk -v slash="$slash" '{printf "%.4f %s \n", $2, slash}'
echo "\hline"
#cat $prefit_stat/sig_ttgg/Mgg_c3_HL/log_combine.out | grep "Significance:" | awk -v ORS=" " -v sslash="${single}" '{printf "1$%stau$ + 0 lepton & %.4f & ", sslash,$2}'
#cat $prefit_exp/sig_ttgg/Mgg_c3_HL/log_combine.out | grep "Significance:" | awk -v ORS=" " '{printf "%.4f & ", $2}'
#cat $prefit/sig_ttgg/Mgg_c3_HL/log_combine.out | grep "Significance:" | awk -v slash="$slash" '{printf "%.4f %s \n", $2, slash}'
cat $prefit_stat/sig_Mgg_c4_Zveto_HL.txt | grep "Significance:" | awk -v ORS=" " -v sslash="${single}" '{printf "2 $%stau$s & %.4f & ", sslash, $2}'
cat $prefit_exp/sig_Mgg_c4_Zveto_HL.txt | grep "Significance:" | awk -v ORS=" " '{printf "%.4f & ", $2}'
cat $prefit_exp1/sig_Mgg_c4_Zveto_HL.txt | grep "Significance:" | awk -v ORS=" " '{printf "%.4f & ", $2}'
cat $prefit/sig_Mgg_c4_Zveto_HL.txt | grep "Significance:" | awk -v slash="$slash" '{printf "%.4f %s \n", $2, slash}'
echo "\hline"
cat $prefit_stat/sig_ttgg.txt| grep "Significance:" | awk -v ORS=" " '{printf "Combination & %.4f & ", $2}'
cat $prefit_exp/sig_ttgg.txt| grep "Significance:" | awk -v ORS=" " '{printf "%.4f & ", $2}'
cat $prefit_exp1/sig_ttgg.txt| grep "Significance:" | awk -v ORS=" " '{printf "%.4f & ", $2}'
cat $prefit/sig_ttgg.txt| grep "Significance:" | awk -v slash="$slash" '{printf "%.4f %s \n", $2, slash}'
echo "\hline"
echo "##############################################"
echo "Combination"
echo "##############################################"
cat $prefit_stat/sig_wwgg.txt | grep "Significance:" | awk -v ORS=" " -v sslash="${single}" '{printf "%swwgg & %.4f & ",sslash, $2}'
cat $prefit_exp/sig_wwgg.txt | grep "Significance:" | awk -v ORS=" " '{printf "%.4f & ", $2}'
cat $prefit_exp1/sig_wwgg.txt | grep "Significance:" | awk -v ORS=" " '{printf "%.4f & ", $2}'
cat $prefit/sig_wwgg.txt | grep "Significance:" | awk -v slash="$slash" '{printf "%.4f %s \n", $2, slash}'
cat $prefit_stat/sig_ttgg.txt| grep "Significance:" | awk -v ORS=" " -v sslash="${single}" '{printf "%sttgg & %.4f & ", sslash, $2}'
cat $prefit_exp/sig_ttgg.txt| grep "Significance:" | awk -v ORS=" " '{printf "%.4f & ", $2}'
cat $prefit_exp1/sig_ttgg.txt| grep "Significance:" | awk -v ORS=" " '{printf "%.4f & ", $2}'
cat $prefit/sig_ttgg.txt| grep "Significance:" | awk -v slash="$slash" '{printf "%.4f %s \n", $2, slash}'
cat $prefit_stat/sig_combined.txt  | grep "Significance:" | awk -v ORS=" " '{printf "Combination & %.4f & ", $2}'
cat $prefit_exp/sig_combined.txt  | grep "Significance:" | awk -v ORS=" " '{printf "%.4f & ", $2}'
cat $prefit_exp1/sig_combined.txt  | grep "Significance:" | awk -v ORS=" " '{printf "%.4f & ", $2}'
cat $prefit/sig_combined.txt  | grep "Significance:" | awk -v slash="$slash" '{printf "%.4f %s \n", $2, slash}'
