#!/bin/bash 

sample=$1
echo $sample
xrdfs se01.indiacms.res.in ls /cms/store/group/Snowmass_2021_2022/DelphesNtuplizer/$sample | awk -v s=$sample '{print "root://se01.indiacms.res.in/"$1}' > $sample.txt
