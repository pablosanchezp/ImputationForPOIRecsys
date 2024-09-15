#!/bin/bash
#Necessary to execute the script with ./ or bash

#Global variables for jvm
JAR=ImputationForPOIRecSys.jar
jvmMemory=-Xmx24G


javaCommand=java

#Directories
originalCities=OriginalCities
processedCities=ProcessedCities
extensionMap=_Mapping.txt

aggregateStrategyTime=LAST
aggregateStrategy=SUM

cities="MX_MexicoCity RU_Moscow CL_Santiago JP_Tokyo US_NewYork GB_London"

# This cities are for the original dataset of Foursquare
KCore=2
coordFile="POIS_Coords.txt"


prefix="_K"$KCore"_AgT"$aggregateStrategyTime"_AP"$aggregateStrategy"_T"
suffixnewTrain=$prefix"Train"
suffixnewTest=$prefix"Test"


unzip $processedCities/POIS_Mapping.zip
mv POIS_Mapping.txt $processedCities

unzip $processedCities/POIS_Coords.zip
mv POIS_Coords.txt $processedCities

for city in $cities
do
  fullPathOriginal="$originalCities"/"$city".txt


  # Generate the new file with timestamps
  fileNewIds=$processedCities/"$city"_Processed.txt
  if [ ! -f "$fileNewIds" ]; then
    $javaCommand $jvmMemory -jar $JAR -o generateNewCheckingFileWithTimeStamps -trf $fullPathOriginal -IMapping $processedCities/POIS"$extensionMap" -UMapping $processedCities/Users"$extensionMap" -newDataset $fileNewIds
  fi


  #KCore
  fileKCore=$processedCities/"$city"_Processed_KCore"$KCore".txt
  if [ ! -f "$fileKCore" ]; then
    $javaCommand $jvmMemory -jar $JAR -o Kcore -trf $fileNewIds -mru $KCore -mri $KCore -orf $fileKCore
  fi

  #Aggregate BEFORE splitting
  fileAggregate=$processedCities/"$city"_Processed_KCore"$KCore"_AggrT"$aggregateStrategyTime"_AggrP"$aggregateStrategy".dat
  if [ ! -f "$fileAggregate" ]; then
    $javaCommand $jvmMemory -jar $JAR -o AggregateWithWrapperTimeStamps -trf $fileKCore -wStrat $aggregateStrategy -wStratTime $aggregateStrategyTime -newDataset $fileAggregate
  fi

  #Generate TRAIN TEST SPLITS
  trainFile=$processedCities/"$city"_Processed_KCore"$KCore"_AggrT"$aggregateStrategyTime"_AggrP"$aggregateStrategy"_TempTrain.dat
  testFile=$processedCities/"$city"_Processed_KCore"$KCore"_AggrT"$aggregateStrategyTime"_AggrP"$aggregateStrategy"_TempTest.dat
  if [ ! -f "$trainFile" ]; then
    $javaCommand $jvmMemory -jar $JAR -o TemporalGlobalSplit -trf $fileAggregate $trainFile $testFile 0.8
  fi

  # Statistics Only in Train
  fileStatistics=$processedCities/"$city"_Processed_KCore"$KCore"_AggrT"$aggregateStrategyTime"_AggrP"$aggregateStrategy"_SatatisticsPercentageCheckinsByMidPoint.dat
  if [ ! -f "$fileStatistics" ]; then
    $javaCommand $jvmMemory -jar $JAR -o SatatisticsPercentageCheckinsByMidPoint -trf $trainFile -coordFile "$processedCities"/POIS_Coords.txt -orf $fileStatistics -listNumbers "1,5,10,15,20,2000"
  fi


  # Rename/move the files
  originaltrainFile=$processedCities/"$city"_Processed_KCore"$KCore"_AggrT"$aggregateStrategyTime"_AggrP"$aggregateStrategy"_TempTrain.dat
  originaltestFile=$processedCities/"$city"_Processed_KCore"$KCore"_AggrT"$aggregateStrategyTime"_AggrP"$aggregateStrategy"_TempTest.dat

  newtrainFile=$processedCities/"$city""$suffixnewTrain"".dat"
  newtestFile=$processedCities/"$city""$suffixnewTest"".dat"

  if [ ! -f "$newtrainFile" ]; then
    cp $originaltrainFile $newtrainFile
  fi

  if [ ! -f "$newtestFile" ]; then
    cp $originaltestFile $newtestFile
  fi

done # En cities
wait

# Generate ALL ImputationFiles

perUser=true
distanceMidpoint=10 #For imputing preferences, the maximum distance that the POI must be (obtained by the statistics file)


for city in $cities
do
  # For every city we generate new imputation train (random, pure distance, and recommenders IB and PopGeo)


  neighbours=10
  trainfile=$processedCities/"$city""$suffixnewTrain"".dat"
  testfile=$processedCities/"$city""$suffixnewTest"".dat"


  # Imputation for NOT reverse users (ordering the users from lower to higher number of checkins and impute them)
  for perIncr in 10 30
  do

    for reverse in false #true
    do

      # For random imputation
      outputImputationFile=$processedCities/"$city""$suffixnewTrain""_BINIMPU"$perIncr"_KmDist"$distanceMidpoint"_Rnd_AvgPrefUsersOrdRev"$reverse".dat"
      $javaCommand $jvmMemory -jar $JAR -o ImputationDataset -trf $trainfile -tsf $testfile -cIndex false -rr RandomRecommender -orf $outputImputationFile -perUser $perUser -perIncr $perIncr -coordFile "$processedCities"/POIS_Coords.txt -thr $distanceMidpoint -n $neighbours -reverse $reverse

      # For averageDistance
      outputImputationFile=$processedCities/"$city""$suffixnewTrain""_BINIMPU"$perIncr"_KmDist"$distanceMidpoint"_AvgDistFreq_AvgPrefUsersOrdRev"$reverse".dat"
      $javaCommand $jvmMemory -jar $JAR -o ImputationDataset -trf $trainfile -tsf $testfile -cIndex false -rr AverageDistanceUserGEO -rs "notUsed" -orf $outputImputationFile -perUser $perUser -perIncr $perIncr -coordFile "$processedCities"/POIS_Coords.txt -thr $distanceMidpoint -n $neighbours -scoreFreq FREQUENCY -reverse $reverse


      # For ITEM NEIGHBORHOOD dataset generation
      IBsimilarity=SJIS
      neighbours=90
      outputImputationFile=$processedCities/"$city""$suffixnewTrain""_BINIMPU"$perIncr"_KmDist"$distanceMidpoint"_IB_"$IBsimilarity"k"$neighbours"_AvgPrefUsersOrdRev"$reverse".dat"
      $javaCommand $jvmMemory -jar $JAR -o ImputationDataset -trf $trainfile -tsf $testfile -cIndex false -rr ItemNeighborhoodRecommender -rs $IBsimilarity -n $neighbours -orf $outputImputationFile -perUser $perUser -perIncr $perIncr -coordFile "$processedCities"/POIS_Coords.txt -thr $distanceMidpoint -reverse $reverse

      # For PopGeoNN generation
      UBsimilarity="SJUS"
      neighbours=100
      outputImputationFile=$processedCities/"$city""$suffixnewTrain""_BINIMPU"$perIncr"_KmDist"$distanceMidpoint"_PopGeoNN_"$UBsimilarity"k"$neighbours"_AvgPrefUsersOrdRev"$reverse".dat"
      $javaCommand $jvmMemory -jar $JAR -o ImputationDataset -trf $trainfile -tsf $testfile -cIndex false -rr PopGeoNN -rs $UBsimilarity -n $neighbours -orf $outputImputationFile -perUser $perUser -perIncr $perIncr -coordFile "$processedCities"/POIS_Coords.txt -thr $distanceMidpoint -reverse $reverse

    done
    wait

  done
  wait


  # Imputation all training Users by a percentage
  for perIncr2 in 10 30
  do
    #For random generation
    outputImputationFile=$processedCities/"$city""$suffixnewTrain""_BINIMPUTrUsers"$perIncr2"_KmDist"$distanceMidpoint"_Rnd.dat"
    $javaCommand $jvmMemory -jar $JAR -o ImputationDataset -trf $trainfile -tsf $testfile -cIndex false -rr RandomRecommender -orf $outputImputationFile -perUser $perUser -perIncr $perIncr2 -coordFile "$processedCities"/POIS_Coords.txt -thr $distanceMidpoint -n $neighbours -reverse true -imputeAllUsers true

    #For avgDis
    outputImputationFile=$processedCities/"$city""$suffixnewTrain""_BINIMPUTrUsers"$perIncr2"_KmDist"$distanceMidpoint"_AvgDistFreq.dat"
    $javaCommand $jvmMemory -jar $JAR -o ImputationDataset -trf $trainfile -tsf $testfile -cIndex false -rr AverageDistanceUserGEO -rs "notUsed" -orf $outputImputationFile -perUser $perUser -perIncr $perIncr2 -coordFile "$processedCities"/POIS_Coords.txt -thr $distanceMidpoint -n $neighbours -scoreFreq FREQUENCY -reverse true -imputeAllUsers true


    #For ITEM NEIGHBORHOOD dataset generation
    IBsimilarity=SJIS
    neighbours=90
    outputImputationFile=$processedCities/"$city""$suffixnewTrain""_BINIMPUTrUsers"$perIncr2"_KmDist"$distanceMidpoint"_IB_"$IBsimilarity"k"$neighbours".dat"
    $javaCommand $jvmMemory -jar $JAR -o ImputationDataset -trf $trainfile -tsf $testfile -cIndex false -rr ItemNeighborhoodRecommender -rs $IBsimilarity -n $neighbours -orf $outputImputationFile -perUser $perUser -perIncr $perIncr2 -coordFile "$processedCities"/POIS_Coords.txt -thr $distanceMidpoint -reverse true -imputeAllUsers true

    # For PopGeoNN generation
    UBsimilarity="SJUS"
    neighbours=100
    outputImputationFile=$processedCities/"$city""$suffixnewTrain""_BINIMPUTrUsers"$perIncr2"_KmDist"$distanceMidpoint"_PopGeoNN_"$UBsimilarity"k"$neighbours".dat"
    $javaCommand $jvmMemory -jar $JAR -o ImputationDataset -trf $trainfile -tsf $testfile -cIndex false -rr PopGeoNN -rs $UBsimilarity -n $neighbours -orf $outputImputationFile -perUser $perUser -perIncr $perIncr2 -coordFile "$processedCities"/POIS_Coords.txt -thr $distanceMidpoint -reverse true -imputeAllUsers true

  done
  wait

done
wait
