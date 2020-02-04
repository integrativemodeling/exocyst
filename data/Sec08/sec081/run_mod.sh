#!/bin/bash

mdlName=$1

module load modeller
python model.py > model.log

rm *.V99* *.D00* *.sch *.rsr *.ini

bestModel=`cat model.log | grep B999 | sort -n -k 3 | awk '{print $1}' | head -n 1`
echo $bestModel

cp $bestModel $mdlName'.pdb'

#mkdir modeller
#mv *.seq *.ali *.py *B999* *.log modeller
