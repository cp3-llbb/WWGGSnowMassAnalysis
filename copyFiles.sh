#!/bin/bash 

sample=$1
echo $sample
xrdfs root://cmseos.fnal.gov/ ls /store/user/snowmass/Snowmass2021/DelphesNtuplizer/$sample | awk -v s=$sample '{print "xrdcp root://cmseos.fnal.gov/"$1" /nfs/user/sdonerta/"$1}' > copy_$sample.txt
chmod 755 copy_$sample.txt
./copy_$sample.txt >& chk_$sample.txt &
