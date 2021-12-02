#!/bin/bash

echo -e "cd BambooDatacardProducer/"
cd BambooDatacardProducer/

echo -e "source inference/setup.sh"
source inference/setup.sh

unset PYTHONPATH

echo -e "source /cvmfs/sft.cern.ch/lcg/views/LCG_101/x86_64-centos7-gcc10-opt/setup.sh"
source /cvmfs/sft.cern.ch/lcg/views/LCG_101/x86_64-centos7-gcc10-opt/setup.sh

echo -e "python -m venv python3"
python -m venv python3

echo - e "source python3/bin/activate"
source python3/bin/activate

# pip install enlighten # if needed

path=/afs/cern.ch/work/a/aguzel/private/bamboodev/output-tautaugg_withcategories6-2/

categories="hasTwoTausNoLept hasOneTauOneElec hasOneTauOneMuon hasOneTauNoLept"

for category in $categories
do
    echo -e "python produceDataCards.py --yaml datacard_tautaugg.yml --pseudodata --custom category=${category}"
    python produceDataCards.py --yaml datacard_tautaugg.yml --pseudodata --custom category=${category} path=${path}
done
