#!/bin/bash
#Check if root is installed and activated:
if ! [[ $(eval which root) ]]; then
    echo "ROOT not found. Please install the root framework!"
    return 1
fi


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

mode=$(eval jq .mode $JSONName) #mode: 1=Single, 2=serial Fit, 3=parallelFit, 4=parallelRun, 5=parallel Fit with power Law source
Mult=$(eval jq .Mult $JSONName) #Multiplicity * 10 to get 1 digit precision, Mult dN/dEta
MultType=$(eval jq .MultType $JSONName) #0=Fixed, 1=DistributionFromParameterizazion, 2=DistributionFromFile, 3=Poissonian
NEvents=$(eval jq .NEvents $JSONName)
outputNameNoQuote=ToMMCA
outputName='"'${outputNameNoQuote}'"'
HadronizationModeNoQuote=$(eval jq .HadronizationMode $JSONName) #QuarkRecombination
HadronizationMode=$(eval jq .HadronizationMode $JSONName)

######## do not fiddle below here if you do not know what you are doing! ################



NoFoldersToday=$(find ./Output/ -maxdepth 1 -type d | grep $(eval date +%b%d) | wc -l)
DateString=$(eval date +%b%d)
starttime=$(date)
FolderName="Output/${DateString}_${NoFoldersToday}"
mkdir -p $FolderName
mkdir -p slurmFiles
mkdir -p temp
cp ToMCCA_main.cpp $FolderName
cp -r include $FolderName
cp ToMCCA.sh $FolderName
cp $JSONName $FolderName
#cd $FolderName

#[[ $# -ge 1 ]] && echo "$(date): $@" >> RunPlan.txt
echo -n "Start: $FolderName $(date +%H:%M): ## $@ ## End: " >> RunPlan.txt

if [ $mode -eq 1 ]
then
    root -l -q 'ToMCCA_main.cpp++( '${mode}', '${Mult}', '${MultType}', '${NEvents}', '${HadronizationMode}')'
    wait
    mv Output.root $FolderName
    echo "Files can be found in $FolderName"
    echo "$(date +%H:%M)" >> RunPlan.txt
    echo "all done"
fi




