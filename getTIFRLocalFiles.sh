#!/bin/bash 

sample=$1
echo $sample
ls /nfs/user/sdonerta/cms/store/group/Snowmass_2021_2022/DelphesNtuplizer/$sample | awk -v s=$sample '{print "/nfs/user/sdonerta/cms/store/group/Snowmass_2021_2022/DelphesNtuplizer/'$sample'/"$1}' > $sample.txt
