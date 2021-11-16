#!/bin/bash 

sample=$1
echo $sample
ls /nfs/user/sdonerta/store/user/snowmass/Snowmass2021/DelphesNtuplizer/$sample | awk -v s=$sample '{print "/nfs/user/sdonerta/store/user/snowmass/Snowmass2021/DelphesNtuplizer/'$sample'/"$1}' > $sample.txt