#!/bin/bash

echo -e "cd BambooDatacardProducer/"
cd BambooDatacardProducer/

echo -e "source inference/setup.sh"
source inference/setup.sh

unset PYTHONPATH

echo -e "source /cvmfs/sft.cern.ch/lcg/views/LCG_101/x86_64-centos7-gcc10-opt/setup.sh"
source /cvmfs/sft.cern.ch/lcg/views/LCG_101/x86_64-centos7-gcc10-opt/setup.sh

echo -e "python -m venv python3"
python -m venv python3 # execute this at first trial and remove

echo - e "source python3/bin/activate"
source python3/bin/activate

pip install enlighten # execute this at first trial and remove

path=/afs/cern.ch/work/a/aguzel/private/bamboodev/output-tautaugg-looseID-withPTmggRatio_split/

categories="c1_Zveto c2_Zveto c3 c4_Zveto"

for category in $categories
do
    echo -e "python produceDataCards.py --yaml datacard_tautaugg.yml --pseudodata --custom category=${category}"
    python produceDataCards.py --yaml datacard_tautaugg.yml --pseudodata --custom category=${category} path=${path}
done
