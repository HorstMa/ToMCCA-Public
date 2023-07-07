#!/bin/bash
#Give the JSON file with the setting's name (emit .json)
JSONName="JSONs/$1"
if [[ $JSONName != *".json" ]]; then
    JSONName=$JSONName".json"
fi
if [ -f $JSONName ]; then
    echo "Running JSON $JSONName"
else
    echo "JSON JSONs/$JSONName not found"
    exit
fi
Mult=$(eval jq .Mult $JSONName) #Multiplicity * 10 to get 1 digit precision, Mult dN/dEta
MultType=$(eval jq .MultType $JSONName) #0=Fixed, 1=DistributionFromParameterizazion, 2=DistributionFromFile, 3=Poissonian
NEvents=$(eval jq .NEvents $JSONName)


######## do not fiddle below here if you do not know what you are doing! ################



NoFoldersToday=$(find . -maxdepth 1 -type d | grep $(eval date +%b%d) | wc -l)
DateString=$(eval date +%b%d)
starttime=$(date)
FolderName="${DateString}_${NoFoldersToday}"
mkdir -p $FolderName
mkdir -p slurmFiles
mkdir -p temp
cp ToMCCA_Bash.cpp $FolderName
cp ToMCCA_func.h $FolderName
cp ToMCCA_Opts.sh $FolderName
cp $JSONName $FolderName
#cd $FolderName

#[[ $# -ge 1 ]] && echo "$(date): $@" >> RunPlan.txt
echo -n "Start: $FolderName $(date +%H:%M): ## $@ ## End: " >> RunPlan.txt

root -l -q 'ToMCCA_Bash.cpp( '${Mult}', '${MultType}', '${NEvents}')'
wait
mv Output.root $FolderName
echo "$(date +%H:%M)" >> RunPlan.txt
echo "all done"





