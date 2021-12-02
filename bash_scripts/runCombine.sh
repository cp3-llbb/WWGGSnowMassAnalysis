#!/bin/bash

cd BambooDatacardProducer/
echo -e "Running: source inference/setup.sh"
source inference/setup.sh

path=/afs/cern.ch/work/a/aguzel/private/bamboodev/output-tautaugg_withcategories6-2/

categories="hasTwoTausNoLept hasOneTauOneElec hasOneTauOneMuon hasOneTauNoLept"

for category in $categories
do
    echo -e "Running combine Significance for category: $category"
    combine -M Significance ${path}forcombine_${category}/Inv_massGG_HL-LHC.txt -m 125 --rMin -10 --rMax 10 -t 1 --expectSignal 1 > ${category}.txt
    mv ${category}.txt ${category}/${category}.txt
    mv higgsCombineTest.Significance.mH125.123456.root ${category}/higgsCombineTest.Significance.mH125.123456.root

    echo -e "Running combine  MultiDimFit for category: $category"
    combine -M MultiDimFit ${path}forcombine_${category}/Inv_massGG_HL-LHC.txt -m 125 --rMin -10 --rMax 10 -t 1 --expectSignal 1 --algo grid --setParameterRanges r=-1,3 --saveNLL
    mv higgsCombineTest.MultiDimFit.mH125.123456.root ${category}/higgsCombineTest.MultiDimFit.mH125.123456.root
done

