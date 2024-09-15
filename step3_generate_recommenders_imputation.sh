: '
  Script for generating the recommenders

'

# !/bin/bash

JAR=ImputationForPOIRecSys.jar
jvmMemory=-Xmx24G


javaCommand=java
originalCities=OriginalCities
processedCities=ProcessedCities
extensionMap=_Mapping.txt

aggregateStrategy=SUM
aggregateStrategyTime=LAST



# This cities are for the original dataset of Foursquare

KCore=2
coordFile="POIS_Coords.txt"


prefix="_K"$KCore"_AgT"$aggregateStrategyTime"_AP"$aggregateStrategy"_T"
suffixnewTrain=$prefix"Train"
suffixnewTest=$prefix"Test"


recPrefix=rec

# cities selected to perform the experiments
cities="MX_MexicoCity RU_Moscow CL_Santiago JP_Tokyo US_NewYork GB_London"



recFolder=RecommendationFolder
resultFolder=ResultFolder

allneighbours="10 20 30 40 50 60 70 80 90 100"
itemsRecommended=30

coordFile="$processedCities"/POIS_Coords.txt

allKFactorizerRankSys="10 50 100"
allLambdaFactorizerRankSys="0.1 1 10"
allAlphaFactorizerRankSys="0.1 1 10"


cutoffs="5,10,20"
cutoffsWrite="5-10-20"

mymedialitePath=MyMediaLite-3.11/bin
BPRFactors=$allKFactorizerRankSys
BPRBiasReg="0 0.5 1"
BPRLearnRate=0.05
BPRNumIter=50
BPRRegUs="0.0025 0.001 0.005 0.01 0.1"
BPRRegJ="0.00025 0.0001 0.0005 0.001 0.01"
extensionMyMediaLite=MyMedLt
geoBPRmaxDists="1 4"



# For rankGeoFM
cs_RANKGEOFM="1"
alphas_RANKGEOFM="0.1 0.2"
epsilons_RANKGEOFM="0.3"
kfactors_RANKGEOFM=$allKFactorizerRankSys

ns_RANKGEOFM="10 50 100 200"
iters_RANKGEOFM="50 120 200"
decays_RANKGEOFM="1"
isboldDriver_RANKGEOFM="true"
learnRates_RANKGEOFM="0.001"
maxRates_RANKGEOFM="0.001"


# For FMFMGM
alphas_FMFMGM="0.2 0.4"
thetas_FMFMGM="0.02 0.1"
distances_FMFMGM="15"
iters_FMFMGM="30"
kfactors_FMFMGM="50 100"
alphas2_FMFMGM="20 40"
betas_FMFMGM="0.2"
learningRates_FMFMGM="0.0001"
sigmoids_FMFMGM="false"


maxDiffTime=1814400
minDiffTime=60
minClosePrefBot=3

coordFile="$processedCities"/POIS_Coords.txt

#IRenMF parameters and data
#IRenMF parameters and data
# PLEASE, CHANGE THE PATH FOR YOUR MATLAB INSTALATION

fullpathMatlab="/home/pablosanchez/MATLABR2024aAll/bin/matlab"
#fullpathMatlab="/home/pablosanchez/MATLABR2024aAll/bin/matlab"

pathMatlab="$(dirname $fullpathMatlab)"

fullPath="$(pwd)"
pathIrenMF=$fullPath/IRenMFWithScore
pathIrenMFTesWithTrain=$fullPath/IRenMFTestWithTrain

pathDest=$fullPath

#Variables for IrenMF (paper of ExperimentalEvaluation. Most of these variables are updated in the configureFile)
lambda1=0.015
lambda2=0.015
GeoNN=10

allKIRENMF="100 50"
allAlphaIRENMF="0.4 0.6"
allLambda3IRENMF="1 0.1"
clustersIRENMF="5 50"
extensionCoords=_Coords.txt


# For IRenMF
function obtainConfigureFile() #The arguments are the k, the alpha and the lambda3
{
  configureFileSimple=""
  if [ "$1" == "100" ] && [ "$2" == "0.4" ] && [ "$3" == "1" ] && [ "$4" == "50" ] ; then
    configureFileSimple="configure00" #Configure 0
  fi

  if [ "$1" == "100" ] && [ "$2" == "0.4" ] && [ "$3" == "0.1" ] && [ "$4" == "50" ] ; then
    configureFileSimple="configure01" #Configure 1
  fi

  if [ "$1" == "100" ] && [ "$2" == "0.6" ] && [ "$3" == "1" ] && [ "$4" == "50" ] ; then
    configureFileSimple="configure02" #Configure 2
  fi

  if [ "$1" == "100" ] && [ "$2" == "0.6" ] && [ "$3" == "0.1" ] && [ "$4" == "50" ] ; then
    configureFileSimple="configure03" #Configure 3
  fi

  if [ "$1" == "50" ] && [ "$2" == "0.4" ] && [ "$3" == "1" ] && [ "$4" == "50" ] ; then
    configureFileSimple="configure04" #Configure 4
  fi

  if [ "$1" == "50" ] && [ "$2" == "0.4" ] && [ "$3" == "0.1" ] && [ "$4" == "50" ] ; then
    configureFileSimple="configure05" #Configure 5
  fi

  if [ "$1" == "50" ] && [ "$2" == "0.6" ] && [ "$3" == "1" ] && [ "$4" == "50" ] ; then
    configureFileSimple="configure06" #Configure 6
  fi

  if [ "$1" == "50" ] && [ "$2" == "0.6" ] && [ "$3" == "0.1" ] && [ "$4" == "50" ] ; then
    configureFileSimple="configure07" #Configure 7
  fi

  #Change of clusters
  if [ "$1" == "100" ] && [ "$2" == "0.4" ] && [ "$3" == "1" ] && [ "$4" == "5" ] ; then
    configureFileSimple="configure08" #Configure 0
  fi

  if [ "$1" == "100" ] && [ "$2" == "0.4" ] && [ "$3" == "0.1" ] && [ "$4" == "5" ] ; then
    configureFileSimple="configure09" #Configure 1
  fi

  if [ "$1" == "100" ] && [ "$2" == "0.6" ] && [ "$3" == "1" ] && [ "$4" == "5" ] ; then
    configureFileSimple="configure10" #Configure 2
  fi

  if [ "$1" == "100" ] && [ "$2" == "0.6" ] && [ "$3" == "0.1" ] && [ "$4" == "5" ] ; then
    configureFileSimple="configure11" #Configure 3
  fi

  if [ "$1" == "50" ] && [ "$2" == "0.4" ] && [ "$3" == "1" ] && [ "$4" == "5" ] ; then
    configureFileSimple="configure12" #Configure 4
  fi

  if [ "$1" == "50" ] && [ "$2" == "0.4" ] && [ "$3" == "0.1" ] && [ "$4" == "5" ] ; then
    configureFileSimple="configure13" #Configure 5
  fi

  if [ "$1" == "50" ] && [ "$2" == "0.6" ] && [ "$3" == "1" ] && [ "$4" == "5" ] ; then
    configureFileSimple="configure14" #Configure 6
  fi

  if [ "$1" == "50" ] && [ "$2" == "0.6" ] && [ "$3" == "0.1" ] && [ "$4" == "5" ] ; then
    configureFileSimple="configure15" #Configure 7
  fi


  echo "$configureFileSimple"
}


# For each city (single domain, each city individually)
for city in $cities
do
  # The test file is always the same

  testfile=$processedCities/"$city""$suffixnewTest"".dat"
  nonImputedTrainFile=$processedCities/"$city""$suffixnewTrain"".dat"
  completeFile=$processedCities/"$city"_Processed_KCore2.txt


  recommendationFolder=$recFolder"$city"
  mkdir -p $recommendationFolder


  # Find all the training files of the city EXCEPT THE ORIGINAL. That is, all files that matches "$city""$suffixnewTrain" but are not "$city""$suffixnewTrain"".dat"
  find $processedCities/ -name "$city""$suffixnewTrain"* -not -name "$city""$suffixnewTrain"".dat" | while read trainFile; do

  # Create spcific POI coords
   cityPOICoords="$processedCities"/"$city"POIS_Coords.txt
   $javaCommand $jvmMemory -jar $JAR -o SpecificPOICoords -trf $nonImputedTrainFile -coordFile $coordFile -orf $cityPOICoords

   cityCompletePOICoords="$processedCities"/"$city"CompletePOIS_Coords.txt
   $javaCommand $jvmMemory -jar $JAR -o SpecificPOICoords -trf $completeFile -coordFile $coordFile -orf $cityCompletePOICoords

   echo "Current trainFile is $trainFile"

   title=$(basename "$trainFile" .dat) #extension removed. We need to know in the recommender generation which is the trainFile used. Title is the variable to be used in the recommendation generation



   for ranksysRecommender in RandomRecommender PopularityRecommender
   do
     outputRecfile=$recommendationFolder/"$recPrefix"_"$title"_RSys_"$ranksysRecommender".txt
     $javaCommand $jvmMemory -jar $JAR -o ranksysOnlyComplete -trf $trainFile -trf2 $nonImputedTrainFile -tsf $testfile -cIndex true -rr "$ranksysRecommender" -rs "notUsed" -nI $itemsRecommended -n 20 -orf $outputRecfile
     $javaCommand $jvmMemory -jar $JAR -o CheckRecommendationsFile -trf $nonImputedTrainFile -tsf $testfile -rf $outputRecfile


   done # End ranksys recommenders
   wait



   # Only bests for every recommender
   UBsimilarity=""
   neighboursUB=""
   neighboursIB=""
   IBsimilarity=""

   #GeoBPRMF
   factor=""
   GEOBPRNumIter=""
   GEOBPRLearnRate=""
   bias_reg=""
   regU=""
   maxDist=""

   #HKV
   kFactor=""
   alphaValue=""
   lambdaValue=""

   #FMFMGM
   alpha_FMFMGM=""
   theta_FMFMGM=""
   distance_FMFMGM=""
   iter_FMFMGM=""
   kfactor_FMFMGM=""
   alpha2_FMFMGM=""
   beta_FMFMGM=""
   learningRate_FMFMGM=""
   sigmoid_FMFMGM=""

   #RankGeoFM
   c="1"
   kFactorizer=""
   kNeighbour=""
   decay=""
   isboldDriver=""
   numIter=""
   learnRate=""
   maxRate=""
   alphaF=""
   epsilon=""

   #PopGeoNN
   UBsimilarityPOP=""
   neighboursPOP=""

   #IRENMF


   #RankGeoTest

   if [[ $city == *CL_Santiago* ]]; then
     UBsimilarity="SJUS"
     neighboursUB="100"
     neighboursIB="90"
     IBsimilarity="SJIS"

     #GeoBPRMF
     factor="50"
     GEOBPRNumIter="50"
     GEOBPRLearnRate="0.05"
     bias_reg="1"
     regU="0.01"
     maxDist="4"

     #HKV
     kFactor="100"
     alphaValue="0.1"
     lambdaValue="0.1"

     #FMFMGM
     alpha_FMFMGM="0.2"
     theta_FMFMGM="0.02"
     distance_FMFMGM="15"
     iter_FMFMGM="30"
     kfactor_FMFMGM="100"
     alpha2_FMFMGM="20"
     beta_FMFMGM="0.2"
     learningRate_FMFMGM="0.0001"
     sigmoid_FMFMGM="false"

     #RankGeoFM
     kFactorizer="100"
     kNeighbour="100"
     decay="1"
     isboldDriver="true"
     numIter="200"
     learnRate="0.001"
     maxRate="0.001"
     alphaF="0.2"
     epsilon="0.3"

     #PopGeoNN
     UBsimilarityPOP="SJUS"
     neighboursPOP="100"

     #_BPRMF_nF
     BPRfactor="10"
     BPRNumIter="50"
     BPRLearnRate="0.05"
     BPRbias_reg="0"
     BPRregU="0.01"

     #IRENMF
     KIRENMF="100"
     alphaIRENMF="0.4"
     lambda3IRENMF="1"
     clusterIRENMF="50"


   fi

   if [[ $city == *JP_Tokyo* ]]; then
     UBsimilarity="SJUS"
     neighboursUB="100"
     neighboursIB="100"
     IBsimilarity="SJIS"

     #GeoBPRMF
     factor="100"
     GEOBPRNumIter="50"
     GEOBPRLearnRate="0.05"
     bias_reg="0"
     regU="0.001"
     maxDist="1"

     #HKV
     kFactor="10"
     alphaValue="0.1"
     lambdaValue="1"


     #FMFMGM
     alpha_FMFMGM="0.4"
     theta_FMFMGM="0.02"
     distance_FMFMGM="15"
     iter_FMFMGM="30"
     kfactor_FMFMGM="100"
     alpha2_FMFMGM="20"
     beta_FMFMGM="0.2"
     learningRate_FMFMGM="0.0001"
     sigmoid_FMFMGM="false"

     #RankGeoFM
     kFactorizer="100"
     kNeighbour="10"
     decay="1"
     isboldDriver="true"
     numIter="200"
     learnRate="0.001"
     maxRate="0.001"
     alphaF="0.2"
     epsilon="0.3"

     #PopGeoNN
     UBsimilarityPOP="SJUS"
     neighboursPOP="80"

     #_BPRMF_nF
     BPRfactor="50"
     BPRNumIter="50"
     BPRLearnRate="0.05"
     BPRbias_reg="0"
     BPRregU="0.1"

     #IRENMF
     KIRENMF="100"
     alphaIRENMF="0.4"
     lambda3IRENMF="1"
     clusterIRENMF="50"
   fi

   if [[ $city == *US_NewYork* ]]; then
     UBsimilarity="SJUS"
     neighboursUB="100"
     neighboursIB="90"
     IBsimilarity="SJIS"

     #GeoBPRMF
     factor="100"
     GEOBPRNumIter="50"
     GEOBPRLearnRate="0.05"
     bias_reg="0.5"
     regU="0.1"
     maxDist="4"

     #HKV
     kFactor="50"
     alphaValue="10"
     lambdaValue="0.1"

     #FMFMGM
     alpha_FMFMGM="0.4"
     theta_FMFMGM="0.1"
     distance_FMFMGM="15"
     iter_FMFMGM="30"
     kfactor_FMFMGM="100"
     alpha2_FMFMGM="20"
     beta_FMFMGM="0.2"
     learningRate_FMFMGM="0.0001"
     sigmoid_FMFMGM="false"

     #RankGeoFM
     kFactorizer="100"
     kNeighbour="100"
     decay="1"
     isboldDriver="true"
     numIter="200"
     learnRate="0.001"
     maxRate="0.001"
     alphaF="0.2"
     epsilon="0.3"


     #PopGeoNN
     UBsimilarityPOP="SJUS"
     neighboursPOP="100"

     #_BPRMF_nF
     BPRfactor="50"
     BPRNumIter="50"
     BPRLearnRate="0.05"
     BPRbias_reg="0"
     BPRregU="0.005"

     #IRENMF
     KIRENMF="100"
     alphaIRENMF="0.4"
     lambda3IRENMF="1"
     clusterIRENMF="5"


   fi

   if [[ $city == *GB_London* ]]; then
     UBsimilarity="SJUS"
     neighboursUB="100"
     neighboursIB="90"
     IBsimilarity="SJIS"

     #GeoBPRMF
     factor="50"
     GEOBPRNumIter="50"
     GEOBPRLearnRate="0.05"
     bias_reg="0.5"
     regU="0.01"
     maxDist="1"

     #HKV
     kFactor="100"
     alphaValue="10"
     lambdaValue="0.1"

     #FMFMGM
     alpha_FMFMGM="0.4"
     theta_FMFMGM="0.1"
     distance_FMFMGM="15"
     iter_FMFMGM="30"
     kfactor_FMFMGM="100"
     alpha2_FMFMGM="20"
     beta_FMFMGM="0.2"
     learningRate_FMFMGM="0.0001"
     sigmoid_FMFMGM="false"


     #RankGeoFM
     kFactorizer="100"
     kNeighbour="50"
     decay="1"
     isboldDriver="true"
     numIter="200"
     learnRate="0.001"
     maxRate="0.001"
     alphaF="0.2"
     epsilon="0.3"

     #PopGeoNN
     UBsimilarityPOP="SJUS"
     neighboursPOP="100"

     #_BPRMF_nF
     BPRfactor="50"
     BPRNumIter="50"
     BPRLearnRate="0.05"
     BPRbias_reg="0"
     BPRregU="0.1"

     KIRENMF="100"
     alphaIRENMF="0.4"
     lambda3IRENMF="1"
     clusterIRENMF="5"

   fi

   if [[ $city == *MX_MexicoCity* ]]; then
     UBsimilarity="SJUS"
     neighboursUB="100"
     neighboursIB="100"
     IBsimilarity="SJIS"

     #GeoBPRMF
     factor="50"
     GEOBPRNumIter="50"
     GEOBPRLearnRate="0.05"
     bias_reg="0.5"
     regU="0.01"
     maxDist="4"

     #HKV
     kFactor="100"
     alphaValue="1"
     lambdaValue="0.1"

     #FMFMGM
     alpha_FMFMGM="0.2"
     theta_FMFMGM="0.1"
     distance_FMFMGM="15"
     iter_FMFMGM="30"
     kfactor_FMFMGM="100"
     alpha2_FMFMGM="20"
     beta_FMFMGM="0.2"
     learningRate_FMFMGM="0.0001"
     sigmoid_FMFMGM="false"


     #RankGeoFM
     kFactorizer="100"
     kNeighbour="200"
     decay="1"
     isboldDriver="true"
     numIter="200"
     learnRate="0.001"
     maxRate="0.001"
     alphaF="0.1"
     epsilon="0.3"

     #PopGeoNN
     UBsimilarityPOP="SJUS"
     neighboursPOP="90"

     #_BPRMF_nF
     BPRfactor="50"
     BPRNumIter="50"
     BPRLearnRate="0.05"
     BPRbias_reg="0"
     BPRregU="0.1"

     KIRENMF="100"
     alphaIRENMF="0.4"
     lambda3IRENMF="1"
     clusterIRENMF="5"

   fi


   if [[ $city == *RU_Moscow* ]]; then
     UBsimilarity="SJUS"
     neighboursUB="100"
     neighboursIB="100"
     IBsimilarity="SJIS"

     #GeoBPRMF
     factor="100"
     GEOBPRNumIter="50"
     GEOBPRLearnRate="0.05"
     bias_reg="1"
     regU="0.005"
     maxDist="1"

     #HKV
     kFactor="100"
     alphaValue="1"
     lambdaValue="0.1"

     #FMFMGM
     alpha_FMFMGM="0.4"
     theta_FMFMGM="0.02"
     distance_FMFMGM="15"
     iter_FMFMGM="30"
     kfactor_FMFMGM="100"
     alpha2_FMFMGM="20"
     beta_FMFMGM="0.2"
     learningRate_FMFMGM="0.0001"
     sigmoid_FMFMGM="false"


     #RankGeoFM
     kFactorizer="100"
     kNeighbour="100"
     decay="1"
     isboldDriver="true"
     numIter="200"
     learnRate="0.001"
     maxRate="0.001"
     alphaF="0.1"
     epsilon="0.3"

     #PopGeoNN
     UBsimilarityPOP="SJUS"
     neighboursPOP="100"

     #_BPRMF_nF
     BPRfactor="100"
     BPRNumIter="50"
     BPRLearnRate="0.05"
     BPRbias_reg="1"
     BPRregU="0.01"

     KIRENMF="100"
     alphaIRENMF="0.4"
     lambda3IRENMF="1"
     clusterIRENMF="50"

   fi

   #_BPRMF_nF

   BPRregJ=$(echo "$BPRregU/10" | bc -l)

   outputRecfile=$recommendationFolder/"$recPrefix"_"$title"_"$extensionMyMediaLite"_BPRMF_nF"$BPRfactor"_nIt"$BPRNumIter"_LR"$BPRLearnRate"_BR"$BPRbias_reg"_RU"$BPRregU"_RI"$BPRregU"_RJ"$BPRregJ""Rep$repetition"NOTUSE.txt
   echo $outputRecfile

   realoutputRecfile=$recommendationFolder/"$recPrefix"_"$title"_"$extensionMyMediaLite"_BPRMF_nF"$BPRfactor"_nIt"$BPRNumIter"_LR"$BPRLearnRate"_BR"$BPRbias_reg"_RU"$BPRregU"_RI"$BPRregU"_RJ"$BPRregJ""Rep$repetition".txt
   if [ ! -f "$realoutputRecfile" ]; then

     outputRecfile2=$outputRecfile"Aux".txt
     # We need repeated items
     echo "./$mymedialitePath/item_recommendation --training-file=$trainFile --recommender=BPRMF --prediction-file=$outputRecfile2 --predict-items-number=200 --repeated-items --recommender-options="num_factors=$BPRfactor bias_reg=$BPRbias_reg reg_u=$BPRregU reg_i=$BPRregU reg_j=$BPRregJ learn_rate=$BPRLearnRate UniformUserSampling=false WithReplacement=false num_iter=$BPRNumIter""

     ./$mymedialitePath/item_recommendation --training-file=$trainFile --recommender=BPRMF --prediction-file=$outputRecfile2 --repeated-items --predict-items-number=200 --recommender-options="num_factors=$BPRfactor bias_reg=$BPRbias_reg reg_u=$BPRregU reg_i=$BPRregU reg_j=$BPRregJ learn_rate=$BPRLearnRate UniformUserSampling=false WithReplacement=false num_iter=$BPRNumIter"
     $javaCommand $jvmMemory -jar $JAR -o ParseMyMediaLite -trf $outputRecfile2 $testfile $outputRecfile
     rm $outputRecfile2
     echo "This should find repetitions"
     # $javaCommand $jvmMemory -jar $JAR -o CheckRecommendationsFile -trf $nonImputedTrainFile -tsf $testfile -rf $outputRecfile

     #Now, we filter
     $javaCommand $jvmMemory -jar $JAR -o filterRecFileByTrain -trf $nonImputedTrainFile -rf $outputRecfile -orf $realoutputRecfile -thr 100
     rm $outputRecfile

     #echo "This should NOT find repetitions"
     #$javaCommand $jvmMemory -jar $JAR -o CheckRecommendationsFile -trf $nonImputedTrainFile -tsf $testfile -rf $realoutputRecfile
   fi




   #RankGeoFM
   outputRecfile=$recommendationFolder/"$recPrefix"_"$title"_RSys_POI_"RankGeoFM""kFac"$kFactorizer"kNgh"$kNeighbour"dec"$decay"bDriv"$isboldDriver"It"$numIter"lR"$learnRate"mR"$maxRate"a"$alphaF"c"$c"eps"$epsilon"".txt
   $javaCommand $jvmMemory -jar $JAR -o ranksysOnlyComplete -trf $trainFile -trf2 $nonImputedTrainFile -tsf $testfile -cIndex false -rr "RankGeoFMRecommender" -coordFile $cityPOICoords -n $kNeighbour -nIFactorizer $numIter -kFactorizer $kFactorizer -aFactorizer $alphaF -epsilon $epsilon -c $c -svdDecay $decay -svdIsboldDriver $isboldDriver -svdLearnRate $learnRate -svdMaxLearnRate $maxRate -nI $itemsRecommended -orf $outputRecfile
   #$javaCommand $jvmMemory -jar $JAR -o CheckRecommendationsFile -trf $nonImputedTrainFile -tsf $testfile -rf $outputRecfile



    # UB

    outputRecfile=$recommendationFolder/"$recPrefix"_"$title"_RSys_UB_"$UBsimilarity"_k"$neighboursUB".txt
    $javaCommand $jvmMemory -jar $JAR -o ranksysOnlyComplete -trf $trainFile -trf2 $nonImputedTrainFile -tsf $testfile -cIndex false -rr UserNeighborhoodRecommender -rs $UBsimilarity -nI $itemsRecommended -n $neighboursUB -orf $outputRecfile
    #$javaCommand $jvmMemory -jar $JAR -o CheckRecommendationsFile -trf $nonImputedTrainFile -tsf $testfile -rf $outputRecfile

    # IB
    outputRecfile=$recommendationFolder/"$recPrefix"_"$title"_RSys_IB_"$IBsimilarity"_k"$neighboursIB".txt
    $javaCommand $jvmMemory -jar $JAR -o ranksysOnlyComplete -trf $trainFile -trf2 $nonImputedTrainFile -tsf $testfile -cIndex false -rr ItemNeighborhoodRecommender -rs $IBsimilarity -nI $itemsRecommended -n $neighboursIB -orf $outputRecfile
    #$javaCommand $jvmMemory -jar $JAR -o CheckRecommendationsFile -trf $nonImputedTrainFile -tsf $testfile -rf $outputRecfile





    #GeoBPRMF
    outputRecfile=$recommendationFolder/"$recPrefix"_"$title"_RSys_GeoBPRMF_nF"$factor"_nIt"$GEOBPRNumIter"_LR"$GEOBPRLearnRate"_BR"$bias_reg"_RU"$regU"_RI"$regU"_MaxD"$maxDist"_Rep"$repetition".txt
    $javaCommand $jvmMemory -jar $JAR -o ranksysOnlyComplete -trf $trainFile -trf2 $nonImputedTrainFile -tsf $testfile -cIndex false -rr GeoBPRMF -nI $itemsRecommended -n 20 -orf $outputRecfile -kFactorizer $factor -nIFactorizer $GEOBPRNumIter -svdLearnRate $GEOBPRLearnRate -svdRegBias $bias_reg -svdRegUser $regU -svdRegItem $regU -coordFile $cityPOICoords -maxDist $maxDist
    #$javaCommand $jvmMemory -jar $JAR -o CheckRecommendationsFile -trf $nonImputedTrainFile -tsf $testfile -rf $outputRecfile


    #HKV

    rankRecommenderNoSim=MFRecommenderHKV
    outputRecfile=$recommendationFolder/"$recPrefix"_"$title"_RSys_"$rankRecommenderNoSim"_kF"$kFactor"_aF"$alphaValue"_lF"$lambdaValue".txt
    #$javaCommand $jvmMemory -jar $JAR -o ranksysOnlyComplete -trf $trainFile -trf2 $nonImputedTrainFile -tsf $testfile -cIndex false -rr $rankRecommenderNoSim -rs "notUsed" -nI $itemsRecommended -n 20 -orf $outputRecfile -kFactorizer $kFactor -aFactorizer $alphaValue -lFactorizer $lambdaValue



    for poiRecommender in AverageDistanceUserGEO
    do
      outputRecfile=$recommendationFolder/"$recPrefix"_"$title"_RSys_POI_"$poiRecommender"SIMPLE.txt
      $javaCommand $jvmMemory -jar $JAR -o ranksysOnlyComplete -trf $trainFile -trf2 $nonImputedTrainFile -tsf $testfile -cIndex false -rr $poiRecommender -rs "notUsed" -nI $itemsRecommended -n 20 -orf $outputRecfile -coordFile $cityPOICoords
      #$javaCommand $jvmMemory -jar $JAR -o CheckRecommendationsFile -trf $nonImputedTrainFile -tsf $testfile -rf $outputRecfile
    done # End poiRecommender
    wait


    # Using the SUM

    for poiRecommender in AverageDistanceUserGEO
    do
      outputRecfile=$recommendationFolder/"$recPrefix"_"$title"_RSys_POI_"$poiRecommender"FREQUENCY.txt
      $javaCommand $jvmMemory -jar $JAR -o ranksysOnlyComplete -trf $trainFile -trf2 $nonImputedTrainFile -tsf $testfile -cIndex false -rr $poiRecommender -rs "notUsed" -nI $itemsRecommended -n 20 -orf $outputRecfile -coordFile $cityPOICoords -scoreFreq FREQUENCY
      #$javaCommand $jvmMemory -jar $JAR -o CheckRecommendationsFile -trf $nonImputedTrainFile -tsf $testfile -rf $outputRecfile
    done
    wait


    for poiRecommender in KDEstimatorRecommender
    do
      outputRecfile=$recommendationFolder/"$recPrefix"_"$title"_RSys_POI_"$poiRecommender".txt
      $javaCommand $jvmMemory -jar $JAR -o ranksysOnlyComplete -trf $trainFile -trf2 $nonImputedTrainFile -tsf $testfile -cIndex false -rr $poiRecommender -rs "notUsed" -nI $itemsRecommended -n 20 -orf $outputRecfile -coordFile $cityPOICoords
      #$javaCommand $jvmMemory -jar $JAR -o CheckRecommendationsFile -trf $nonImputedTrainFile -tsf $testfile -rf $outputRecfile

    done # End poiRecommender
    wait


    #FMFMGM
    outputRecfile=$recommendationFolder/"$recPrefix"_"$title"_RSys_POI_"FMFMGM"_a"$alpha_FMFMGM"_t"$theta_FMFMGM"_d"$distance_FMFMGM"_i"$iter_FMFMGM"_f"$kfactor_FMFMGM"_a2"$alpha2_FMFMGM"_b"$beta_FMFMGM"_lR"$learningRate_FMFMGM"_sig"$sigmoid_FMFMGM".txt
    $javaCommand $jvmMemory -jar $JAR -o ranksysOnlyComplete -trf $trainFile -trf2 $nonImputedTrainFile -tsf $testfile -cIndex false -rr "FMFMGM" -coordFile $cityPOICoords -aFactorizer $alpha_FMFMGM -thetaFactorizer $theta_FMFMGM -maxDist $distance_FMFMGM -nIFactorizer $iter_FMFMGM -kFactorizer $kfactor_FMFMGM -aFactorizer2 $alpha2_FMFMGM -svdBeta $beta_FMFMGM -svdLearnRate $learningRate_FMFMGM -useSigmoid $sigmoid_FMFMGM -n 20 -nI $itemsRecommended -orf $outputRecfile
    #$javaCommand $jvmMemory -jar $JAR -o CheckRecommendationsFile -trf $nonImputedTrainFile -tsf $testfile -rf $outputRecfile




    #PopGeoNN

    poiRecommender="PopGeoNN"
    outputRecfile=$recommendationFolder/"$recPrefix"_"$title"_RSys_POI_"$poiRecommender"_UBSim_"$UBsimilarityPOP"_k"$neighboursPOP".txt
    $javaCommand $jvmMemory -jar $JAR -o ranksysOnlyComplete -trf $trainFile -trf2 $nonImputedTrainFile -tsf $testfile -cIndex true -rr $poiRecommender -rs $UBsimilarityPOP -nI $itemsRecommended -n $neighboursPOP -orf $outputRecfile -coordFile $cityCompletePOICoords
    #$javaCommand $jvmMemory -jar $JAR -o CheckRecommendationsFile -trf $nonImputedTrainFile -tsf $testfile -rf $outputRecfile


    #IRENMF
    resultFileNNCityFile=$pathDest/"$processedCities"/POIS_""$city""$suffixnewTrain""_"$GeoNN"NN.txt
    poisCoordsOfCityFile=$pathDest/"$processedCities"/POIS_""$city""$suffixnewTrain""_"$extensionCoords"
    configureFileSimple=$(obtainConfigureFile "$KIRENMF" "$alphaIRENMF" "$lambda3IRENMF" "$clusterIRENMF")

    echo "Parameters K=$KIRENMF alpha=$alphaIRENMF lambda3=$lambda3IRENMF clusters=$clusterIRENMF"
    echo "File " $title " working with configure " $configureFileSimple

    configureFileToSelect="configure"_"K"$KIRENMF"_Alpha""$alphaIRENMF""_L1_"$lambda1"_L2_"$lambda2"_L3_"$lambda3IRENMF"_GeoNN"$GeoNN"_clusters"$clusterIRENMF
    outputRecfile=$pathDest/$recommendationFolder/"$recPrefix"_"$title"_IRENMF_"$configureFileToSelect".txt

    clusterFile=$pathDest/$recommendationFolder/$title"_"Clusters$clusterIRENMF"_"
    rm $clusterFile"_clusters".mat

    if [ ! -f $outputRecfile ]; then
        #ItemGroupPOI(configure_file, trainFile, poisCoordsFile,GeoNNFile,testingFile, candidatesFile,outFile,clusterFile)
        outputRecfile_IRENMF="$outputRecfile""Aux"
        $pathMatlab/./matlab -nodisplay -nodesktop -r "cd '$pathIrenMFTesWithTrain/'; ItemGroupPOI('$configureFileSimple', '$pathDest/$trainFile', '$poisCoordsOfCityFile', '$resultFileNNCityFile', '$pathDest/$testfile', '$resultFileNNCityFile', '$outputRecfile_IRENMF', '$clusterFile'); quit"

        #
        #$javaCommand $jvmMemory -jar $JAR -o CheckRecommendationsFile -trf $nonImputedTrainFile -tsf $testfile -rf $outputRecfile_IRENMF

        $javaCommand $jvmMemory -jar $JAR -o filterRecFileByTrain -trf $nonImputedTrainFile -rf $outputRecfile_IRENMF -orf $outputRecfile -thr 100
        rm $outputRecfile_IRENMF

        echo "This should NOT find repetitions"
        #$javaCommand $jvmMemory -jar $JAR -o CheckRecommendationsFile -trf $nonImputedTrainFile -tsf $testfile -rf $outputRecfile
    fi

    echo "Finished $configureFileSimple"




    if [[ $city == *GB_London* ]]; then

      for neighbours in $allneighbours
      do

        # ranksys with similarities UserBased
        for UBsimilarity in SJUS VCUS
        do
          outputRecfile=$recommendationFolder/"$recPrefix"_"$title"_RSys_UB_"$UBsimilarity"_k"$neighbours".txt
          $javaCommand $jvmMemory -jar $JAR -o ranksysOnlyComplete -trf $trainFile -trf2 $nonImputedTrainFile -tsf $testfile -cIndex false -rr UserNeighborhoodRecommender -rs $UBsimilarity -nI $itemsRecommended -n $neighbours -orf $outputRecfile
          $javaCommand $jvmMemory -jar $JAR -o CheckRecommendationsFile -trf $nonImputedTrainFile -tsf $testfile -rf $outputRecfile

        done # UB sim
        wait

        # Ranksys with similarities ItemBased
        for IBsimilarity in SJIS VCIS
        do
          outputRecfile=$recommendationFolder/"$recPrefix"_"$title"_RSys_IB_"$IBsimilarity"_k"$neighbours".txt
          $javaCommand $jvmMemory -jar $JAR -o ranksysOnlyComplete -trf $trainFile -trf2 $nonImputedTrainFile -tsf $testfile -cIndex false -rr ItemNeighborhoodRecommender -rs $IBsimilarity -nI $itemsRecommended -n $neighbours -orf $outputRecfile
          $javaCommand $jvmMemory -jar $JAR -o CheckRecommendationsFile -trf $nonImputedTrainFile -tsf $testfile -rf $outputRecfile
        done
        wait # End IB Sim

      done # End Neighbours
      wait


      if [[ "$title" == *"APSUM_TTrain_BINIMPUTrUsers10_KmDist10_PopGeoNN_SJUSk100"* ]] || [[ "$title" == *"APSUM_TTrain_BINIMPUTrUsers30_KmDist10_PopGeoNN_SJUSk100"* ]]; then


        for repetition in 1 #2 3 4
        do
          for BPRfactor in $BPRFactors
          do
            for BPRbias_reg in $BPRBiasReg
            do
              for BPRregU in $BPRRegUs #Regularization for items and users is the same
              do
                BPRregJ=$(echo "$BPRregU/10" | bc -l)

                outputRecfile=$recommendationFolder/"$recPrefix"_"$title"_"$extensionMyMediaLite"_BPRMF_nF"$BPRfactor"_nIt"$BPRNumIter"_LR"$BPRLearnRate"_BR"$BPRbias_reg"_RU"$BPRregU"_RI"$BPRregU"_RJ"$BPRregJ""Rep$repetition"NOTUSE.txt
                echo $outputRecfile

                realoutputRecfile=$recommendationFolder/"$recPrefix"_"$title"_"$extensionMyMediaLite"_BPRMF_nF"$BPRfactor"_nIt"$BPRNumIter"_LR"$BPRLearnRate"_BR"$BPRbias_reg"_RU"$BPRregU"_RI"$BPRregU"_RJ"$BPRregJ""Rep$repetition".txt
                if [ ! -f "$realoutputRecfile" ]; then

                  outputRecfile2=$outputRecfile"Aux".txt
                  # We need repeated items
                  echo "./$mymedialitePath/item_recommendation --training-file=$trainFile --recommender=BPRMF --prediction-file=$outputRecfile2 --predict-items-number=200 --repeated-items --recommender-options="num_factors=$BPRfactor bias_reg=$BPRbias_reg reg_u=$BPRregU reg_i=$BPRregU reg_j=$BPRregJ learn_rate=$BPRLearnRate UniformUserSampling=false WithReplacement=false num_iter=$BPRNumIter""

                  ./$mymedialitePath/item_recommendation --training-file=$trainFile --recommender=BPRMF --prediction-file=$outputRecfile2 --repeated-items --predict-items-number=200 --recommender-options="num_factors=$BPRfactor bias_reg=$BPRbias_reg reg_u=$BPRregU reg_i=$BPRregU reg_j=$BPRregJ learn_rate=$BPRLearnRate UniformUserSampling=false WithReplacement=false num_iter=$BPRNumIter"
                  $javaCommand $jvmMemory -jar $JAR -o ParseMyMediaLite -trf $outputRecfile2 $testfile $outputRecfile
                  rm $outputRecfile2
                  echo "This should find repetitions"
                  # $javaCommand $jvmMemory -jar $JAR -o CheckRecommendationsFile -trf $nonImputedTrainFile -tsf $testfile -rf $outputRecfile

                  #Now, we filter
                  $javaCommand $jvmMemory -jar $JAR -o filterRecFileByTrain -trf $nonImputedTrainFile -rf $outputRecfile -orf $realoutputRecfile -thr 100
                  rm $outputRecfile

                  echo "This should NOT find repetitions"
                  $javaCommand $jvmMemory -jar $JAR -o CheckRecommendationsFile -trf $nonImputedTrainFile -tsf $testfile -rf $realoutputRecfile
                fi

              done # Reg U
              wait
            done # Bias reg
            wait
          done # Factors
          wait
        done # Repetition
        wait

        for repetition in 1 #2 3 4 5
        do
          for factor in $BPRFactors
          do
            for bias_reg in $BPRBiasReg
            do
              for regU in $BPRRegUs #Regularization for items and users is the same
              do
                for maxDist in $geoBPRmaxDists
                do
                  outputRecfile=$recommendationFolder/"$recPrefix"_"$title"_RSys_GeoBPRMF_nF"$factor"_nIt"$GEOBPRNumIter"_LR"$GEOBPRLearnRate"_BR"$bias_reg"_RU"$regU"_RI"$regU"_MaxD"$maxDist"_Rep"$repetition".txt
                  $javaCommand $jvmMemory -jar $JAR -o ranksysOnlyComplete -trf $trainFile -trf2 $nonImputedTrainFile -tsf $testfile -cIndex false -rr GeoBPRMF -nI $itemsRecommended -n 20 -orf $outputRecfile -kFactorizer $factor -nIFactorizer $GEOBPRNumIter -svdLearnRate $GEOBPRLearnRate -svdRegBias $bias_reg -svdRegUser $regU -svdRegItem $regU -coordFile $cityPOICoords -maxDist $maxDist
                  $javaCommand $jvmMemory -jar $JAR -o CheckRecommendationsFile -trf $nonImputedTrainFile -tsf $testfile -rf $outputRecfile
                done
                wait
              done # Reg U
              wait
            done # Bias reg
            wait
          done # Factors
          wait
        done # Repetition
        wait

        for rankRecommenderNoSim in MFRecommenderHKV
        do
          for kFactor in $allKFactorizerRankSys
          do
            for lambdaValue in $allLambdaFactorizerRankSys
            do
              for alphaValue in $allAlphaFactorizerRankSys
              do
                # Neighbours is put to 20 because this recommender does not use it

                rankRecommenderNoSim=MFRecommenderHKV
                outputRecfile=$recommendationFolder/"$recPrefix"_"$title"_RSys_"$rankRecommenderNoSim"_kF"$kFactor"_aF"$alphaValue"_lF"$lambdaValue".txt
                $javaCommand $jvmMemory -jar $JAR -o ranksysOnlyComplete -trf $trainFile -trf2 $nonImputedTrainFile -tsf $testfile -cIndex false -rr $rankRecommenderNoSim -rs "notUsed" -nI $itemsRecommended -n 20 -orf $outputRecfile -kFactorizer $kFactor -aFactorizer $alphaValue -lFactorizer $lambdaValue



              done
              wait # End alpha values
            done
            wait # End lambda
          done
          wait # End KFactor
        done
        wait # End RankRecommender

        #FMFMGM


        for repetition in 1 #2 3 4 5
        do
          for alpha_FMFMGM in $alphas_FMFMGM
          do
            for theta_FMFMGM in $thetas_FMFMGM
            do
              for distance_FMFMGM in $distances_FMFMGM #
              do
                for iter_FMFMGM in $iters_FMFMGM
                do
                  for kfactor_FMFMGM in $kfactors_FMFMGM
                  do
                    for alpha2_FMFMGM in $alphas2_FMFMGM
                    do
                      for beta_FMFMGM in $betas_FMFMGM
                      do
                        for learningRate_FMFMGM in $learningRates_FMFMGM
                        do
                          for sigmoid_FMFMGM in $sigmoids_FMFMGM
                          do

                            outputRecfile=$recommendationFolder/"$recPrefix"_"$title"_RSys_POI_"FMFMGM"_a"$alpha_FMFMGM"_t"$theta_FMFMGM"_d"$distance_FMFMGM"_i"$iter_FMFMGM"_f"$kfactor_FMFMGM"_a2"$alpha2_FMFMGM"_b"$beta_FMFMGM"_lR"$learningRate_FMFMGM"_sig"$sigmoid_FMFMGM".txt
                            $javaCommand $jvmMemory -jar $JAR -o ranksysOnlyComplete -trf $trainFile -trf2 $nonImputedTrainFile -tsf $testfile -cIndex false -rr "FMFMGM" -coordFile $cityPOICoords -aFactorizer $alpha_FMFMGM -thetaFactorizer $theta_FMFMGM -maxDist $distance_FMFMGM -nIFactorizer $iter_FMFMGM -kFactorizer $kfactor_FMFMGM -aFactorizer2 $alpha2_FMFMGM -svdBeta $beta_FMFMGM -svdLearnRate $learningRate_FMFMGM -useSigmoid $sigmoid_FMFMGM -n 20 -nI $itemsRecommended -orf $outputRecfile
                            $javaCommand $jvmMemory -jar $JAR -o CheckRecommendationsFile -trf $nonImputedTrainFile -tsf $testfile -rf $outputRecfile

                          done # sigmoid
                          wait
                        done # learningRate
                        wait
                      done # beta
                      wait
                     done # alphas2
                     wait
                    done # kFactor
                    wait
                  done # iters
                  wait
                done # distance
                wait
              done # Theta
              wait
            done # alphas
            wait
          done # Repetition
          wait


  # RankGeoTest
    for c in $cs_RANKGEOFM
    do
      for alphaF in $alphas_RANKGEOFM
      do
        for epsilon in $epsilons_RANKGEOFM
        do
          for kFactorizer in $kfactors_RANKGEOFM
          do
            for kNeighbour in $ns_RANKGEOFM
            do
              for numIter in $iters_RANKGEOFM
              do

                for decay in $decays_RANKGEOFM
                do

                  for isboldDriver in $isboldDriver_RANKGEOFM
                  do

                    for learnRate in $learnRates_RANKGEOFM
                    do
                      for maxRate in $maxRates_RANKGEOFM
                      do

                        outputRecfile=$recommendationFolder/"$recPrefix"_"$title"_RSys_POI_"RankGeoFM""kFac"$kFactorizer"kNgh"$kNeighbour"dec"$decay"bDriv"$isboldDriver"It"$numIter"lR"$learnRate"mR"$maxRate"a"$alphaF"c"$c"eps"$epsilon"".txt
                        $javaCommand $jvmMemory -jar $JAR -o ranksysOnlyComplete -trf $trainFile -trf2 $nonImputedTrainFile -tsf $testfile -cIndex false -rr "RankGeoFMRecommender" -coordFile $cityPOICoords -n $kNeighbour -nIFactorizer $numIter -kFactorizer $kFactorizer -aFactorizer $alphaF -epsilon $epsilon -c $c -svdDecay $decay -svdIsboldDriver $isboldDriver -svdLearnRate $learnRate -svdMaxLearnRate $maxRate -nI $itemsRecommended -orf $outputRecfile
                        #$javaCommand $jvmMemory -jar $JAR -o CheckRecommendationsFile -trf $nonImputedTrainFile -tsf $testfile -rf $outputRecfile

                      done # ENd max learn
                      wait

                    done # Learn rate
                    wait
                  done # End bold
                  wait

                done # End decay
                wait
              done # End iter
              wait

            done # ENd neigh
            wait

          done # End kFactors
          wait
        done # Epsilon
        wait
      done # Alpha
      wait
    done # End c
    wait

        #Second baseline: popularity, knn and minimum distance
        for poiRecommender in PopGeoNN
        do
          for UBsimilarity in SJUS VCUS
          do
            for neighbours in $allneighbours
            do
              poiRecommender="PopGeoNN"
              outputRecfile=$recommendationFolder/"$recPrefix"_"$title"_RSys_POI_"$poiRecommender"_UBSim_"$UBsimilarity"_k"$neighbours".txt
              $javaCommand $jvmMemory -jar $JAR -o ranksysOnlyComplete -trf $trainFile -trf2 $nonImputedTrainFile -tsf $testfile -cIndex true -rr $poiRecommender -rs $UBsimilarity -nI $itemsRecommended -n $neighbours -orf $outputRecfile -coordFile $cityCompletePOICoords
              #$javaCommand $jvmMemory -jar $JAR -o CheckRecommendationsFile -trf $nonImputedTrainFile -tsf $testfile -rf $outputRecfile

            done # End neighbours
            wait
          done # End UB similarity
          wait
        done # End poiRecommender
        wait



        #IRENMF
        resultFileNNCityFile=$pathDest/"$processedCities"/POIS_""$city""$suffixnewTrain""_"$GeoNN"NN.txt
        poisCoordsOfCityFile=$pathDest/"$processedCities"/POIS_""$city""$suffixnewTrain""_"$extensionCoords"

        for KIRENMF in $allKIRENMF
        do
          for alphaIRENMF in $allAlphaIRENMF
          do
            for lambda3IRENMF in $allLambda3IRENMF
            do
              for clusterIRENMF in $clustersIRENMF
              do
                  configureFileSimple=$(obtainConfigureFile "$KIRENMF" "$alphaIRENMF" "$lambda3IRENMF" "$clusterIRENMF")

                  echo "Parameters K=$KIRENMF alpha=$alphaIRENMF lambda3=$lambda3IRENMF clusters=$clusterIRENMF"
                  echo "File " $title " working with configure " $configureFileSimple

                  configureFileToSelect="configure"_"K"$KIRENMF"_Alpha""$alphaIRENMF""_L1_"$lambda1"_L2_"$lambda2"_L3_"$lambda3IRENMF"_GeoNN"$GeoNN"_clusters"$clusterIRENMF
                  outputRecfile=$pathDest/$recommendationFolder/"$recPrefix"_"$title"_IRENMF_"$configureFileToSelect".txt

                  clusterFile=$pathDest/$recommendationFolder/$title"_"Clusters$clusterIRENMF"_"
                  rm $clusterFile"_clusters".mat

                  if [ ! -f $outputRecfile ]; then
                      outputRecfile_IRENMF="$outputRecfile""Aux"
                      #ItemGroupPOI(configure_file, trainFile, poisCoordsFile,GeoNNFile,testingFile, candidatesFile,outFile,clusterFile)
                      $pathMatlab/./matlab -nodisplay -nodesktop -r "cd '$pathIrenMFTesWithTrain/'; ItemGroupPOI('$configureFileSimple', '$pathDest/$trainFile', '$poisCoordsOfCityFile', '$resultFileNNCityFile', '$pathDest/$testfile', '$resultFileNNCityFile', '$outputRecfile_IRENMF', '$clusterFile'); quit"


                      #$javaCommand $jvmMemory -jar $JAR -o CheckRecommendationsFile -trf $nonImputedTrainFile -tsf $testfile -rf $outputRecfile_IRENMF

                      $javaCommand $jvmMemory -jar $JAR -o filterRecFileByTrain -trf $nonImputedTrainFile -rf $outputRecfile_IRENMF -orf $outputRecfile -thr 100
                      rm $outputRecfile_IRENMF

                      echo "This should NOT find repetitions"
                      #$javaCommand $jvmMemory -jar $JAR -o CheckRecommendationsFile -trf $nonImputedTrainFile -tsf $testfile -rf $outputRecfile

                      echo "Finished $configureFileSimple"
                  fi
              done #End cluster
              wait
          done
          wait #End lambda
        done
        wait #End alpha IrenMF
      done
      wait #End KIRENMF

      fi



    fi



  done # End find
  wait





done # End cities
wait







: '
Evaluation part
'
JAR=ImputationForPOIRecSys.jar
nonaccresultsPrefix=naev
evthreshold=1


for city in $cities
do


  # In every case
  testfile=$processedCities/"$city""$suffixnewTest"".dat"
  nonImputedTrainFile=$processedCities/"$city""$suffixnewTrain"".dat"


  resFolder=$resultFolder"$city"
  mkdir -p $resFolder

  recommendationFolder=$recFolder"$city"

  find $recommendationFolder/ -name "$recPrefix"*"$city"* | while read recFile; do


    recFileName=$(basename "$recFile" .txt) #extension removed

    #The trainFile must be the non imputed
    outputResultfile=$resFolder/"$nonaccresultsPrefix"_EvTh"$evthreshold"_"$recFileName"_C"$cutoffsWrite".txt
    $javaCommand $jvmMemory -jar $JAR -o ranksysNonAccuracyMetricsEvaluation -trf $nonImputedTrainFile -tsf $testfile -rf $recFile -thr $evthreshold -rc $cutoffs -orf $outputResultfile -onlyAcc false

  done # End Find
  wait

done # End cities
wait
