#!/bin/bash 
slash='\\\'
single='\'

prefix=$1

#for energy in 13TeV 14TeV
for energy in 14TeV
do 
    echo "energy:" $energy  
    for range in 180
    do 
	echo "range:" ${range}"GeV"
	#./testYieldsPerDir.sh ${prefix}_${energy}_1GeVBinning_allSys_newranges${range}
	
	for bin in 1 
	do 
	    echo "binning:" ${bin}"GeV"
	    for sys in all
	    do
		#echo range : $range, sys: $sys "in" $bin "GeV binning: "
		#cat ${dir_name}_${bin}GeVBinning_${sys}Sys_newranges${range}/limit_WWgg/combination_HL/log_combine.out | grep "Expected 50.0" | awk '{printf "%.3f \n" $5}'
		echo "OneL numbers:"
		cat sig_Inv_mass_gghasOneL_DNN_1_HL.txt | grep "Significance:" | awk -v syst="${sys}" '{printf "Cat1 with %s unc: %.4f \n ", syst, $2}'
		cat sig_Inv_mass_gghasOneL_DNN_2_HL.txt | grep "Significance:" | awk -v syst="${sys}" '{printf "Cat2 with %s unc: %.4f \n ", syst, $2}'
		cat sig_Inv_mass_gghasOneL_DNN_3_HL.txt | grep "Significance:" | awk -v syst="${sys}" '{printf "Cat3 with %s unc: %.4f \n ", syst, $2}'
		cat sig_Inv_mass_gghasOneL_DNN_4_HL.txt | grep "Significance:" | awk -v syst="${sys}" '{printf "Cat4 with %s unc: %.4f \n ", syst, $2}'
		cat sig_wwgg_oneL.txt | grep "Significance:" | awk -v syst="${sys}" '{printf "Significance with %s unc: %.4f \n ", syst, $2}'
		echo "TwoL numbers:"
		cat sig_Inv_mass_gghasTwoL_HL.txt | grep "Significance:" | awk -v syst="${sys}" '{printf "Significance with %s unc: %.4f \n ", syst, $2}'

		echo "WWgg combined numbers:"
		cat sig_wwgg.txt | grep "Significance:" | awk -v syst="${sys}" '{printf "Significance with %s unc: %.4f \n ", syst, $2}'

		echo "OneTau numbers:"
		cat sig_Mgg_c3_DNN_1_HL.txt | grep "Significance:" | awk -v syst="${sys}" '{printf "Cat1 with %s unc: %.4f \n ", syst, $2}'
		cat sig_Mgg_c3_DNN_2_HL.txt | grep "Significance:" | awk -v syst="${sys}" '{printf "Cat2 with %s unc: %.4f \n ", syst, $2}'
		cat sig_ttgg_onetau.txt | grep "Significance:" | awk -v syst="${sys}" '{printf "Significance with %s unc: %.4f \n ", syst, $2}'

		echo "2Tau numbers:"
		cat sig_Mgg_c4_Zveto_HL.txt | grep "Significance:" | awk -v syst="${sys}" '{printf "Significance with %s unc: %.4f \n ", syst, $2}'
		
		echo "ttgg combined numbers:"
		cat sig_ttgg.txt | grep "Significance:" | awk -v syst="${sys}" '{printf "Significance with %s unc: %.4f \n ", syst, $2}'

		echo "combined numbers:"
		cat sig_combined.txt | grep "Significance:" | awk -v syst="${sys}" '{printf "Significance with %s unc: %.4f \n ", syst, $2}'
	    done
	    printf "\n"
	done
	printf "\n"
    done
    printf "\n"
done
echo "==============="

