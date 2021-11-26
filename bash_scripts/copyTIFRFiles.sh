#!/bin/bash 

sample=$1
echo $sample
xrdfs se01.indiacms.res.in ls /cms/store/group/Snowmass_2021_2022/DelphesNtuplizer/$sample | awk -v s=$sample '{print "xrdcp root://se01.indiacms.res.in/"$1" /nfs/user/sdonerta/"$1}' > copy_$sample.txt
chmod 755 copy_$sample.txt
./copy_$sample.txt >& chk_$sample.txt &
