w
cd BambooDatacardProducer/
echo -e "Running: source inference/setup.sh"
source inference/setup.sh

echo -e "cd /afs/cern.ch/work/a/aguzel/private/CMSSW_10_2_13/src"
cd /afs/cern.ch/work/a/aguzel/private/CMSSW_10_2_13/src

echo -e "cmsenv"
cmsenv
cd -

path=/afs/cern.ch/work/a/aguzel/private/bamboodev/output-tautaugg-looseID-withPTmggRatio_split/

categories="c1_Zveto c2_Zveto c3 c4_Zveto"

for category in $categories
do
    echo -e "Running combine Significance for category: $category"
    mkdir -p ${path}combine_output_${category}/
    combine -M Significance ${path}forcombine_${category}/Mgg_HL-LHC.txt -m 125 --rMin -10 --rMax 10 -t -1 --expectSignal 1 > ${path}combine_output_${category}/${category}.txt
    mv higgsCombineTest.Significance.mH125.root ${path}combine_output_${category}/

    echo -e "Running combine  MultiDimFit for category: $category"
    combine -M MultiDimFit ${path}forcombine_${category}/Mgg_HL-LHC.txt -m 125 --rMin -10 --rMax 10 -t -1 --expectSignal 1 --algo grid --setParameterRanges r=-1,3 --saveNLL
    mv higgsCombineTest.MultiDimFit.mH125.root ${path}combine_output_${category}/higgsCombineTest.MultiDimFit.mH125.root
done

mkdir -p ${path}combine_output_combined/

combineCards.py bin1=${path}forcombine_c1_Zveto/Mgg_HL-LHC.txt bin2=${path}forcombine_c2_Zveto/Mgg_HL-LHC.txt bin3=${path}forcombine_c3/Mgg_HL-LHC.txt bin4=${path}forcombine_c4_Zveto/Mgg_HL-LHC.txt >& ${path}combine_output_combined/Mgg_HL-LHC_combined_datacard.txt

combine -M Significance ${path}combine_output_combined/Mgg_HL-LHC_combined_datacard.txt -m 125 --rMin -10 --rMax 10 -t -1 --expectSignal 1 > ${path}combine_output_combined/combined_significance.txt

mv higgsCombineTest.Significance.mH125.root ${path}combine_output_combined/

combine -M MultiDimFit ${path}combine_output_combined/Mgg_HL-LHC_combined_datacard.txt -m 125 --rMin -10 --rMax 10 -t -1 --expectSignal 1 --algo grid --setParameterRanges r=-1,3 --saveNLL

mv higgsCombineTest.MultiDimFit.mH125.root ${path}combine_output_combined/

