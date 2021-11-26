#!/bin/bash 

sample=$1
echo $sample
xrdfs root://cmseos.fnal.gov/ ls /store/user/snowmass/Snowmass2021/DelphesNtuplizer/$sample | awk -v s=$sample '{print "root://cmseos.fnal.gov/"$1}' > $sample.txt
