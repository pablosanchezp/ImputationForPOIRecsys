package es.uam.eps.ir.sr.utils;

import static org.ranksys.formats.parsing.Parsers.lp;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.Set;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.function.DoubleUnaryOperator;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

import org.jooq.lambda.tuple.Tuple2;
import org.jooq.lambda.tuple.Tuple3;
import org.jooq.lambda.tuple.Tuple4;
import org.ranksys.formats.feature.SimpleFeaturesReader;
import org.ranksys.formats.parsing.Parser;
import org.ranksys.lda.LDAModelEstimator;
import org.ranksys.lda.LDARecommender;

import com.google.common.util.concurrent.AtomicDouble;

import cc.mallet.topics.ParallelTopicModel;
import es.uam.eps.ir.antimetrics.recommenders.SkylineAntiRelevanceRecommender;
import es.uam.eps.ir.antimetrics.recommenders.SkylineRelevanceRecommender;
import es.uam.eps.ir.attrrec.datamodel.feature.UserFeatureData;
import es.uam.eps.ir.crossdomainPOI.datamodel.SimpleFastTemporalFeaturePreferenceData;
import es.uam.eps.ir.crossdomainPOI.datamodel.SimpleFastTemporalPreferenceData;
import es.uam.eps.ir.crossdomainPOI.datamodel.temporal.interfaces.FastTemporalPreferenceDataIF;
import es.uam.eps.ir.crossdomainPOI.datamodel.temporal.preferences.IdTimePref;
import es.uam.eps.ir.crossdomainPOI.datamodel.temporal.preferences.IdxTimePref;
import es.uam.eps.ir.crossdomainPOI.recommenders.AverageDistanceToUserGEO;
import es.uam.eps.ir.crossdomainPOI.recommenders.KDERecommender;
import es.uam.eps.ir.crossdomainPOI.recommenders.PopGeoNN;
import es.uam.eps.ir.crossdomainPOI.recommenders.SkylineTestOrder;
import es.uam.eps.ir.crossdomainPOI.recommenders.TrainRecommender;
import es.uam.eps.ir.crossdomainPOI.utils.KernelDensityEstimation;
import es.uam.eps.ir.crossdomainPOI.utils.KernelDensityEstimationCoordinates;
import es.uam.eps.ir.crossdomainPOI.utils.UsersMidPoints;
import es.uam.eps.ir.crossdomainPOI.utils.UsersMidPoints.SCORES_FREQUENCY;
import es.uam.eps.ir.ranksys.core.feature.FeatureData;
import es.uam.eps.ir.ranksys.core.feature.SimpleFeatureData;
import es.uam.eps.ir.ranksys.core.preference.PreferenceData;
import es.uam.eps.ir.ranksys.core.preference.SimplePreferenceData;
import es.uam.eps.ir.ranksys.fast.index.FastItemIndex;
import es.uam.eps.ir.ranksys.fast.index.FastUserIndex;
import es.uam.eps.ir.ranksys.fast.index.SimpleFastItemIndex;
import es.uam.eps.ir.ranksys.fast.index.SimpleFastUserIndex;
import es.uam.eps.ir.ranksys.fast.preference.FastPreferenceData;
import es.uam.eps.ir.ranksys.fast.preference.SimpleFastPreferenceData;
import es.uam.eps.ir.ranksys.metrics.rank.ExponentialDiscountModel;
import es.uam.eps.ir.ranksys.metrics.rank.LogarithmicDiscountModel;
import es.uam.eps.ir.ranksys.metrics.rank.NoDiscountModel;
import es.uam.eps.ir.ranksys.metrics.rank.RankingDiscountModel;
import es.uam.eps.ir.ranksys.metrics.rank.ReciprocalDiscountModel;
import es.uam.eps.ir.ranksys.metrics.rel.BackgroundBinaryRelevanceModel;
import es.uam.eps.ir.ranksys.metrics.rel.BinaryRelevanceModel;
import es.uam.eps.ir.ranksys.metrics.rel.NoRelevanceModel;
import es.uam.eps.ir.ranksys.metrics.rel.RelevanceModel;
import es.uam.eps.ir.ranksys.mf.Factorization;
import es.uam.eps.ir.ranksys.mf.als.HKVFactorizer;
import es.uam.eps.ir.ranksys.mf.als.PZTFactorizer;
import es.uam.eps.ir.ranksys.mf.plsa.PLSAFactorizer;
import es.uam.eps.ir.ranksys.mf.rec.MFRecommender;
import es.uam.eps.ir.ranksys.nn.item.ItemNeighborhoodRecommender;
import es.uam.eps.ir.ranksys.nn.item.neighborhood.CachedItemNeighborhood;
import es.uam.eps.ir.ranksys.nn.item.neighborhood.ItemNeighborhood;
import es.uam.eps.ir.ranksys.nn.item.neighborhood.TopKItemNeighborhood;
import es.uam.eps.ir.ranksys.nn.item.sim.SetCosineItemSimilarity;
import es.uam.eps.ir.ranksys.nn.item.sim.SetJaccardItemSimilarity;
import es.uam.eps.ir.ranksys.nn.item.sim.VectorCosineItemSimilarity;
import es.uam.eps.ir.ranksys.nn.item.sim.VectorJaccardItemSimilarity;
import es.uam.eps.ir.ranksys.nn.user.UserNeighborhoodRecommender;
import es.uam.eps.ir.ranksys.nn.user.neighborhood.TopKUserNeighborhood;
import es.uam.eps.ir.ranksys.nn.user.neighborhood.UserNeighborhood;
import es.uam.eps.ir.ranksys.nn.user.sim.SetCosineUserSimilarity;
import es.uam.eps.ir.ranksys.nn.user.sim.SetJaccardUserSimilarity;
import es.uam.eps.ir.ranksys.nn.user.sim.VectorCosineUserSimilarity;
import es.uam.eps.ir.ranksys.nn.user.sim.VectorJaccardUserSimilarity;
import es.uam.eps.ir.ranksys.rec.Recommender;
import es.uam.eps.ir.ranksys.rec.fast.basic.PopularityRecommender;
import es.uam.eps.ir.ranksys.rec.fast.basic.RandomRecommender;
import es.uam.eps.ir.seqawareev.rec.cb.BaselineCB;
import es.uam.eps.ir.seqawareev.rec.cb.UserItemProfileContentRecommender;
import es.uam.eps.ir.seqawareev.sim.AbstractTourCachedItemSimilarity;
import es.uam.eps.ir.seqawareev.sim.location.ItemLocationSimilarity;
import es.uam.eps.ir.seqawareev.sim.probability.ItemFeatureTransitionSimilarity;
import es.uam.eps.ir.seqawareev.sim.probability.ItemTransitionSimilarity;
import es.uam.eps.ir.seqawareev.smooth.JelineckMercer;
import es.uam.eps.ir.seqawareev.smooth.NoSmoothingCond;
import es.uam.eps.ir.seqawareev.smooth.NoSmoothingPrior;
import es.uam.eps.ir.seqawareev.smooth.SmoothingProbabilityIF;
import es.uam.eps.ir.seqawareev.structures.FeatureItemTransitionMatrix;
import es.uam.eps.ir.seqawareev.structures.ItemTransitionMatrix;
import es.uam.eps.ir.seqawareev.tour.recommenders.GenericItemSimilarityTourRecommender;
import es.uam.eps.ir.seqawareev.tour.recommenders.POI.RankGeoFMRecommender;
import es.uam.eps.ir.seqawareev.utils.CBBaselinesUtils.ITEMTRANSFORMATION;
import es.uam.eps.ir.seqawareev.utils.CBBaselinesUtils.USERTRANSFORMATION;
import es.uam.eps.ir.sr.data.POIProcessData;
import es.uam.eps.ir.sr.recommenders.POI.GeoBPRMF;
import es.uam.eps.ir.sr.recommenders.POI.RankGeoFMRecommenderEXSur;
import es.uam.eps.ir.sr.recommenders.POI.RankGeoFMRecommenderNotDense;
import es.uam.eps.ir.sr.recommenders.POI.FMFMGM.FMFMGM;

import es.uam.eps.ir.sr.recommenders.cf.mf.BPRMF;
import es.uam.eps.ir.sr.recommenders.skyline.SkylineRelevanceRecommenderUsersTrain;
import es.uam.eps.ir.sr.utils.comparators.PreferenceComparators;
import es.uam.eps.nets.rankfusion.combination.AnzCombiner;
import es.uam.eps.nets.rankfusion.combination.MNZCombiner;
import es.uam.eps.nets.rankfusion.combination.MaxCombiner;
import es.uam.eps.nets.rankfusion.combination.MedCombiner;
import es.uam.eps.nets.rankfusion.combination.MinCombiner;
import es.uam.eps.nets.rankfusion.combination.SumCombiner;
import es.uam.eps.nets.rankfusion.interfaces.IFCombiner;
import es.uam.eps.nets.rankfusion.interfaces.IFNormalizer;
import es.uam.eps.nets.rankfusion.normalization.DistributionNormalizer;
import es.uam.eps.nets.rankfusion.normalization.DummyNormalizer;
import es.uam.eps.nets.rankfusion.normalization.RankSimNormalizer;
import es.uam.eps.nets.rankfusion.normalization.StdNormalizer;
import es.uam.eps.nets.rankfusion.normalization.SumNormalizer;
import es.uam.eps.nets.rankfusion.normalization.TwoMUVNormalizer;
import es.uam.eps.nets.rankfusion.normalization.ZMUVNormalizer;
import es.uam.eps.nets.util.math.Function.Distribution;

/**
 * Final class that have useful methods in order to work in a more general way.
 *
 *
 * @author Pablo Sanchez (psperez@icai.comillas.edu)
 *
 */
public final class SequentialRecommendersUtils {
	
	//Threshold to allow differences between the expected value and the value obtained in the tests
	public static double thresholdError = 0.001;


	
	
	public static <F, V> void generateJaccardSimFileByUserFeatureData(FastUserIndex<Long> userIndex, UserFeatureData<Long, F, V> ufD, String outputFile) {
			Map<es.uam.eps.ir.attrrec.sim.SymmetricSimilarityBean<Long>, Tuple3<Double, Integer, Integer>> similarities = new HashMap<es.uam.eps.ir.attrrec.sim.SymmetricSimilarityBean<Long>, Tuple3<Double, Integer, Integer>>();

			userIndex.getAllUsers().forEach(u1 -> {
				Set<F> userFeatures1 = ufD.getUserFeatures(u1).map(t -> t.v1).collect(Collectors.toSet());
				userIndex.getAllUsers().filter(u2 -> !u1.equals(u2)).forEach(u2 -> {
					Set<F> userFeatures2 = ufD.getUserFeatures(u2).map(t -> t.v1).collect(Collectors.toSet());
					Set<F> union = new HashSet<>(userFeatures1);
					union.addAll(userFeatures2);
					Set<F> intersection = new HashSet<>(userFeatures1);
					intersection.retainAll(userFeatures2);
					if (intersection.size() != 0 && union.size() != 0) {
						double value = (double) intersection.size() / union.size();
						
						similarities.put(new es.uam.eps.ir.attrrec.sim.SymmetricSimilarityBean<Long>(u1, u2),
								new Tuple3<Double, Integer, Integer>(value, -1, -1));
					}
				});
			});
			es.uam.eps.ir.attrrec.sim.UserIndexReadingSimilarity<Long> userReadingSim = new es.uam.eps.ir.attrrec.sim.UserIndexReadingSimilarity<Long>(similarities, userIndex);
			userReadingSim.save(outputFile, 0.0);
			

		
	}
	
	
	
    public static <K extends Comparable<K>, V extends Comparable<V>> Map<K, V> sortByValue(Map<K, V> unsortMap, boolean reverse) {

        // 1. Convert Map to List of Map
        List<Map.Entry<K, V>> list = new LinkedList<Map.Entry<K, V>>(unsortMap.entrySet());

        // 2. Sort list with Collections.sort(), provide a custom Comparator
        //    Try switch the o1 o2 position for a different order
        Collections.sort(list, new Comparator<Map.Entry<K, V>>() {
            @Override
			public int compare(Map.Entry<K, V> o1,
                    Map.Entry<K, V> o2) {
                if (reverse) {
                    return (o2.getValue()).compareTo(o1.getValue());
                } else {
                    return (o1.getValue()).compareTo(o2.getValue());
                }
            }
        });

        // 3. Loop the sorted list and put it into a new insertion order Map LinkedHashMap
        Map<K, V> sortedMap = new LinkedHashMap<K, V>();
        for (Map.Entry<K, V> entry : list) {
            sortedMap.put(entry.getKey(), entry.getValue());
        }
        return sortedMap;
    }

    
    /**
     * Method to obtain the ranking discount model
     * @param model the model in a string
     * @param base only required for the exponential discount model
     * @return the RankSys Discount model
     */
    public static <U, I> RankingDiscountModel obtRankingDiscountModel(String model, double base) {
        switch (model) {
            case "ExponentialDiscountModel":
                return new ExponentialDiscountModel(base);
            case "LogarithmicDiscountModel":
                return new LogarithmicDiscountModel();
            case "NoDiscountModel":
                return new NoDiscountModel();
            case "ReciprocalDiscountModel":
                return new ReciprocalDiscountModel();
            default:
                return null;
        }
    }

    /**
     * Method to obtain the relevance model from RankSys
     * @param model String of the relevance model
     * @param prefData the preference data
     * @param threshold the threshold
     * @param background only required for background
     * @return the relevance model
     */
    public static <U, I> RelevanceModel<U, I> obtRelevanceModelRanksys(String model, PreferenceData<U, I> prefData, int threshold, double background) {
        switch (model) {
            case "BinaryRelevanceModel":
                return new BinaryRelevanceModel<>(false, prefData, threshold);
            case "NoRelevanceModel":
                return new NoRelevanceModel<U, I>();
            case "BackgroundBinaryRelevanceModel":
                return new BackgroundBinaryRelevanceModel<>(false, prefData, threshold, background);
            default:
                return null;
        }
    }

    /**
     * Obtain RankSys User Similarity
     *
     * @param data the preference data used to compute the similarities
     * @param similarity the similarity
     * @return a RankSys user similarity model
     */
    public static <U, I> es.uam.eps.ir.ranksys.nn.user.sim.UserSimilarity<U> obtRanksysUserSimilarity(FastPreferenceData<U, I> data, String similarity) {
        switch (similarity) {
            case "VectorCosineUserSimilarity":
            case "VCUS":
                // 0.5 to make it symmetrical.
                return new VectorCosineUserSimilarity<>(data, 0.5, true);
            case "VectorJaccardUserSimilarity":
            case "VJUS":
                return new VectorJaccardUserSimilarity<>(data, true);
            case "SetJaccardUserSimilarity":
            case "SJUS":
                return new SetJaccardUserSimilarity<>(data, true);
            case "SetCosineUserSimilarity":
            case "SCUS":
                return new SetCosineUserSimilarity<>(data, 0.5, true);
            default:
            	System.out.println("RankSys user similarity is null");
                return null;
        }
    }
    

    
    public static SmoothingProbabilityIF getSmoothingProbability(String smoothString, double lambda) {
    	if (smoothString == null) {
    		System.out.println("Smoothing probability is null");
			return null;
    	}
    	switch (smoothString) {
    		case "JelineckMercer":
    			return new JelineckMercer(lambda);
    		case "NoSmoothingPrior":
    			return new NoSmoothingPrior();
    		case "NoSmoothingCond":
    			return new NoSmoothingCond();
    		default:
    			System.out.println("Smoothing probability is null");
    			return null;
    	}
    		
    }
    
 
    

    
    
    

    

    

    
    public static <U, I, F> es.uam.eps.ir.ranksys.nn.item.sim.ItemSimilarity<I> obtRanksysItemSimilarity(
            FastPreferenceData<U, I> data, String similarity, FeatureData<I, F, Double> featureData){
    	return obtRanksysItemSimilarity(data, similarity, featureData, true);
    }
    /**
     * Obtain RankSys Item Similarity
     *
     * @param data the preference data used to compute the similarities
     * @param similarity the string identifier of the similarity
     * @return a RankSys item similarity model
     */
    public static <U, I, F> es.uam.eps.ir.ranksys.nn.item.sim.ItemSimilarity<I> obtRanksysItemSimilarity(
            FastPreferenceData<U, I> data, String similarity, FeatureData<I, F, Double> featureData, boolean dense) {
        switch (similarity) {
            case "VectorCosineItemSimilarity":
            case "VCIS":
                // 0.5 to make it symmetrical.
                return new VectorCosineItemSimilarity<>(data, 0.5, dense);
            case "VectorJaccardItemSimilarity":
            case "VJIS":
                return new VectorJaccardItemSimilarity<>(data, dense);
            case "SetJaccardItemSimilarity":
            case "SJIS":
                return new SetJaccardItemSimilarity<>(data, dense);
            case "SetCosineItemSimilarity":
            case "SCIS":
                return new SetCosineItemSimilarity<>(data, 0.5, dense);	
            default:
            	System.out.println(similarity + " not recognized");
                return null;
        }
    }

    public static Map<Long, Double> obtainMapOrdered(List<Tuple2<Double, List<Tuple2<Long, Double>>>> list) {
        Map<Long, Double> m = new HashMap<Long, Double>();
        for (Tuple2<Double, List<Tuple2<Long, Double>>> t : list) { //External list. First element of tuple is the neigh sim. Second element list of items rated by that neighbour 

            List<Tuple2<Long, Double>> items = t.v2;

            for (Tuple2<Long, Double> item : items) {
                double s = t.v1 * item.v2;
                if (m.get(item.v1) == null) {
                    //Rating of the neighbour multiplied by its similarity
                    m.put(item.v1, s);
                } else {
                    m.put(item.v1, s + m.get(item.v1));
                }
            }

        }

        return m;
    }
    
    public static <U, F, V> void analyzeAndSplitPerUserbyUserFeatures(UserFeatureData<U, F, V> ufd, String originalPerUserFile, String resultFile, Parser<U> up) {
    	//Value of feature -> metric -> user -> result
    	Map<String, Map<String, Map<U, Double>>> mappings = new HashMap<>();
    	
    	
    	//A file for each feature
    	Map<String, PrintStream> mapPrint = new HashMap<>();
    	ufd.getAllFeatures().forEach(feature -> {
        	String realResult = resultFile + feature.toString() + ".txt";
        	try {
				mapPrint.put(realResult, new PrintStream(realResult));
				mappings.put(feature.toString(), new HashMap<>());
			} catch (FileNotFoundException e) {
				e.printStackTrace();
			}
        });
    	//For each line, we obtain the feature of the user and send her to the correspondent file
    	//Format is metric, user and result
		try (Stream<String> stream = Files.lines(Paths.get(originalPerUserFile))) {
			stream.forEach(line -> {
				String [] data = line.split("\t");
				String metric = data[0];
				U user = up.parse(data[1]);
				Double result = Double.parseDouble(data[2]);
				
				ufd.getUserFeatures(user).forEach(featureUser -> {
					
					if (mappings.get(featureUser.toString()).get(metric) == null) {
						mappings.get(featureUser.toString()).put(metric, new HashMap<>());
					}
					mappings.get(featureUser.toString()).get(metric).put(user, result);
				});
			});
		} catch (IOException e) {
			e.printStackTrace();
		}
    	//Map read, now split in all files and add the all in each file
    	ufd.getAllFeatures().forEach(feature -> {
    		PrintStream pt = mapPrint.get(feature);
    		Map<String, Map<U, Double>> metricsAndResult = mappings.get(feature);
    		
    		for (String metric: metricsAndResult.keySet()) {
        		double acc = 0;
        		for (U user: metricsAndResult.get(metric).keySet()) {
        			double result = metricsAndResult.get(metric).get(user);
        			pt.println(metric + "\t" + user + "\t" + result);
        			acc += result;
        		}
        		acc /= metricsAndResult.get(metric).keySet().size();
    			pt.println(metric + "\t" + "all" + "\t" + acc);
    		}
    		pt.close();
    	});

		

    	
          
    }
    

    /***
     * Method to print the recommendations for all the recommenders implemented in the framework
     * @param user the user identifier
     * @param item the item identifier
     * @param valueItem the score of the recommendation
     * @param rank the rank
     * @return the rank
     */
    public static String formatRank(Object user, Object item, double valueItem, int rank) {
        return user.toString() + "\t" + item.toString() + "\t" + valueItem + "\t" + rank;
    }

    /**
     * For every item of u1 we will search in u2. If the item of U1 is in U2,
     * (not -1) and the total length (sum of length of indexes of U1 and U2) is
     * higher that the 2 indexes, then we update the values.
     *
     * @param timeStampOrderedU1 the list ordered by timestamp by the first user
     * (from old to new items)
     * @param timeStampOrderedU2 the list ordered by timestamp by the second
     * user (from old to new items)
     * @return a tuple with the last indexes of both users1 and 2
     */
    @Deprecated
    public static Tuple2<Integer, Integer> obtainBestIndexes(List<Long> timeStampOrderedU1, List<Long> timeStampOrderedU2) {
        int maxValue = Integer.MIN_VALUE;
        int indexU1 = 0;
        int indexU2 = 0;

        //Ordered from lower to higher timestamp, but we want to start in the back
        for (int i = timeStampOrderedU1.size() - 1; i >= 0; i--) {
            Long itemU1 = timeStampOrderedU1.get(i);
            boolean foundI2 = false;
            int indexu2 = -1;

            //Search in second user
            for (int j = timeStampOrderedU2.size() - 1; j >= 0; j--) {
                Long itemU2 = timeStampOrderedU2.get(j);

                //Found item of U1 to U2
                if (itemU1 == itemU2) {
                    foundI2 = true;
                    indexu2 = j;
                    break;
                }
                //Not possible to find better index
                if (maxValue >= i + j) {
                    break;
                }

            }

            if (foundI2 && maxValue < (i + indexu2)) {
                maxValue = i + indexu2;
                indexU1 = i;
                indexU2 = indexu2;
            }

            //Not possible to obtain better indexes
            if (maxValue >= i + timeStampOrderedU2.size() - 1) {
                break;
            }
        }

        //Users do not have items in common
        if (maxValue == Integer.MIN_VALUE) {
            indexU1 = -1;
            indexU2 = -1;
        }

        return new Tuple2<Integer, Integer>(indexU1, indexU2);

    }
    
    /***
     * Method to obtain the last indexes in a equivalent way that are obtained in the LCS algorithm
     * @param timeStampOrderedU1 the items consumed by the first user (ordered by timestamp, from the lower to the higher one)
     * @param timeStampOrderedU2 the items consumed by the second user (ordered by timestamp, from the lower to the higher one)
     * @return a tuple indicating the last indexes of both users
     */
    public static Tuple2<Integer, Integer> obtainBestIndexesLCSEquivalent(List<Long> timeStampOrderedU1,List<Long> timeStampOrderedU2){
		int indexU1 = -1;
		int indexU2 = -1;
		
		Set<Long> timeStampOrderedU1Set = new HashSet<>(timeStampOrderedU1);
		Set<Long> timeStampOrderedU2Set = new HashSet<>(timeStampOrderedU2);

		
		
		//Obtain the index of the first user
		for (int i = timeStampOrderedU1.size() - 1; i >= 0; i--){
			Long itemU1 = timeStampOrderedU1.get(i);
			
			if (timeStampOrderedU2Set.contains(itemU1)) {
				indexU1 = i;
				break;
			}
			
		}
		
		//Obtain the index of the second user
		for (int i = timeStampOrderedU2.size() - 1; i >= 0; i--){
			Long itemU2 = timeStampOrderedU2.get(i);			
			if (timeStampOrderedU1Set.contains(itemU2)) {
				indexU2 = i;
				break;
			}
			
		}
		return new Tuple2<Integer, Integer>(indexU1, indexU2);
		
		
	}
    

    
    /***
     * Method to obtain the RankSys recommender depending on the string. It requires the rest of the variables to build the model.
     * The parameters that are not used by our recommender can be null
     * 
     * @param rec the string denoting the recommender we want to instantiate
     * @param ranksysSimilarity similarity for RankSys
     * @param similarity2 the second similarity for RankSys
     * @param data the preference data to build the recommenderf
     * @param kNeighbours the k neighbours
     * @param userIndex the user index
     * @param itemIndex the item index
     * @param kFactorizer k factors for the MF recommender
     * @param alphaF the alpha value for the HKV recommender
     * @param lambdaF the lambda value for the HKV recommender
     * @param numInteractions the numInteractions value for a MF recommender
     * @param originalPOICoordsFile the POICoords file for POI recommenders
     * @param originalFeatureFile the Feature file to be exploited by CB recommender
     * @param poiCityFile the file containing all POIs of a city
     * @param citySimilarityFile 
     * @param reg1 regularization for some MF recommenders
     * @param reg2 second regularization for some MF recommenders
     * @param score if we are considering the scores as frequency or not. Only for POI recommender
     * @return a RankSys recommender
     */
    public static Recommender<Long, Long> obtRankSysRecommeder	(String rec, String additionalRecs, String ranksysSimilarity, String ranksysSimilarity2,
            FastPreferenceData<Long, Long> data, int kNeighbours, int kFactorizer, double alphaF, double lambdaF, int numInteractions,
            String originalPOICoordsFile, String poiCityFile, String citySimilarityFile, double reg1, double reg2, 
            SCORES_FREQUENCY score,  double regUser, double regItem, double learnRate, double maxRate, double regBias, double decay, boolean isboldDriver, 
            double regImpItem, String socialInfluenceFile, double eta, double beta, boolean normalize,
            USERTRANSFORMATION ut, ITEMTRANSFORMATION it, double epsilon, double c, FeatureData<Long, String, Double> featureData, 
            UserFeatureData<Long, String, Double> userFeatureData, String combiner, String normalizer, String weightsForHybrid, 
            String [] completeOrNot, String trainFileForHyrid, String testFileForHybrid, double maxDistance, String mode, double theta, 
            double alphaF2, boolean useSigmoid, String userFeaturesForRecommender, int maxPopItems, double threshold, boolean inverse, boolean negAnti) {
  
    	
    	
        switch (rec) {
        	case "RndRec":
            case "RandomRecommender": {
            	System.out.println("RandomRecommender");
                return new RandomRecommender<>(data, data);
            }
            case "PopRec":
            case "PopularityRecommender":{
            	System.out.println("PopularityRecommender");
                return new PopularityRecommender<>(data);
            }

            case "UBKnnRec":
            case "UserNeighborhoodRecommender": { // User based. Pure CF recommendation
                es.uam.eps.ir.ranksys.nn.user.sim.UserSimilarity<Long> simUNR = obtRanksysUserSimilarity(data, ranksysSimilarity);
                if (simUNR == null) {
                    return null;
                } else {
                	System.out.println("UserNeighborhoodRecommender");
                	System.out.println("kNeighs: "+ kNeighbours);
                	System.out.println("Sim: "+ ranksysSimilarity);
                }
                UserNeighborhood<Long> urneighborhood = new TopKUserNeighborhood<>(simUNR, kNeighbours);
                return new UserNeighborhoodRecommender<>(data, urneighborhood, 1);
            }
            

            
            //End anti-neighbours version
            ///

            case "IBKnnRec":
            case "ItemNeighborhoodRecommender": {// Item based. Pure CF recommendation
                es.uam.eps.ir.ranksys.nn.item.sim.ItemSimilarity<Long> simINR = obtRanksysItemSimilarity(data, ranksysSimilarity, featureData);
                if (simINR == null) {
                	return null;
                }else{
                	System.out.println("ItemNeighborhoodRecommender");
                	System.out.println("kNeighs: "+ kNeighbours);
                	System.out.println("Sim: "+ ranksysSimilarity);             
                }
                ItemNeighborhood<Long> neighborhood = new TopKItemNeighborhood<>(simINR, kNeighbours);
                neighborhood = new CachedItemNeighborhood<>(neighborhood);
                return new ItemNeighborhoodRecommender<>(data, neighborhood, 1);
            }

            case "TrainRec":
            case "TrainRecommender": {
            	System.out.println("TrainRecommender");
            	return new TrainRecommender<>(data, false);
            }
            case "TrainRecRev":
            case "TrainRecommenderReverse": {
            	System.out.println("TrainRecommenderReverse");
            	return new TrainRecommender<>(data, true);
            }
            case "MFRecHKV":
            case "MFRecommenderHKV": { // Matrix factorization
                int k = kFactorizer;
                double lambda = lambdaF;
                double alpha = alphaF;
                int numIter = numInteractions;
                System.out.println("MFRecommenderHKV");
                System.out.println("kFactors: " + k);
                System.out.println("lambda: " + lambda);
                System.out.println("alpha: " + alpha);
                System.out.println("numIter: " + numIter);


                DoubleUnaryOperator confidence = x -> 1 + alpha * x;
                Factorization<Long, Long> factorization = new HKVFactorizer<Long, Long>(lambda, confidence, numIter)
                        .factorize(k, data);
                return new MFRecommender<>(data, data, factorization);
            }
            case "MFRecPZT":
            case "MFRecommenderPZT": { // Probabilistic latent semantic analysis of
                // Hofmann 2004
                int k = kFactorizer;
                double lambda = lambdaF;
                double alpha = alphaF;
                int numIter = numInteractions;
                System.out.println("MFRecommenderPZT");
                System.out.println("kFactors: " + k);
                System.out.println("lambda: " + lambda);
                System.out.println("alpha: " + alpha);
                System.out.println("numIter: " + numIter);
                
                DoubleUnaryOperator confidence = x -> 1 + alpha * x;
                Factorization<Long, Long> factorization = new PZTFactorizer<Long, Long>(lambda, confidence, numIter)
                        .factorize(k, data);
                return new MFRecommender<>(data, data, factorization);
            }
            case "MFRecPLSA":
            case "MFRecommenderPLSA": { // Probabilistic latent semantic analysis of
                // Hofmann 2004
                int k = kFactorizer;
                int numIter = numInteractions;
                System.out.println("MFRecommenderPLSA");
                System.out.println("kFactors: " + k);
                System.out.println("numIter: " + numIter);
                
                Factorization<Long, Long> factorization = new PLSAFactorizer<Long, Long>(numIter).factorize(k, data);
                return new MFRecommender<>(data, data, factorization);
            }
            case "LDARec":
            case "LDARecommender": { // Probabilistic latent semantic analysis of
                // Hofmann 2004
                int k = kFactorizer;
                double alpha = alphaF;
                int numIter = numInteractions;
                int burninPeriod = 50;
                System.out.println("LDARecommender");
                System.out.println("kFactors: " + k);
                System.out.println("numIter: " + numIter);
                System.out.println("beta: " + beta);
                System.out.println("burninPeriod: " + burninPeriod);
                
                
                ParallelTopicModel topicModel = null;
                try {
                    topicModel = LDAModelEstimator.estimate(data, k, alpha, beta, numIter, burninPeriod);
                } catch (IOException e) {
                    return null;
                }
                return new LDARecommender<>(data, data, topicModel);
            }
                        
            case "BPRMF":{
            	int k = kFactorizer;
                int numIter = numInteractions;
                System.out.println("BPRMF");
            	System.out.println("kFactors: "+ k);
            	System.out.println("numIter: "+ numInteractions);
            	System.out.println("regUser: " + regUser);
            	System.out.println("regItem: " + regItem);
            	System.out.println("learnRate: " + learnRate);
            	System.out.println("regBias: " + regBias);
            	BPRMF<Long, Long> bprmf = new BPRMF<>(data, regUser, regItem, regBias, numIter, k, learnRate);
            	bprmf.fit();
            	return bprmf;
            }

            
            case "BaseCB":
            case "BaselineCB": {         		
            	return new BaselineCB<>(data, featureData);
            }

            case "UIProfCBRec":
            case "UserItemProfileContentRecommender": {	
            	return new UserItemProfileContentRecommender<>(data, featureData, normalize);
            	
            }
           
            
            
            /**
             * POI Recommendation algorithms
             */
            //Basic Geographical recommender
            case "AvgDis":
            case "AverageDistanceUserGEO": {
                System.out.println("AverageDistanceUserGEO");

                Map<Long, Tuple2<Double, Double>> mapCoordinates = POICoordinatesMap(originalPOICoordsFile, lp);
                
                UsersMidPoints<Long,Long> usermid = new UsersMidPoints<>(data, mapCoordinates, score);
                return new AverageDistanceToUserGEO<>(data, usermid);
            }
            case "KDERec":
            case "KDEstimatorRecommender" : {
                System.out.println("KDEstimatorRecommender");

                Map<Long, Tuple2<Double, Double>> mapCoordinates = POICoordinatesMap(originalPOICoordsFile, lp);
                Map<Long, List<Tuple2<Double, Double>>> userCoordinates = getCoordinatesUser(data, mapCoordinates);
                   
                KernelDensityEstimation<Long> kde = new KernelDensityEstimationCoordinates<>(userCoordinates);
                return new KDERecommender<>(data, kde, mapCoordinates);	            	
            }
            case "RankGeoFMRec":
            case "RankGeoFMRecommender": {
                System.out.println("RankGeoFMRecommender");
            	System.out.println("kFactors: " + kFactorizer);
            	System.out.println("kNeighbours " + kNeighbours);
            	System.out.println("decay: " + decay);
            	System.out.println("isboldDriver: " + isboldDriver);
            	System.out.println("numIter: " + numInteractions);
            	System.out.println("learnRate: " + learnRate);
            	System.out.println("maxRate: " + maxRate);
                System.out.println("alpha: " + alphaF);
                System.out.println("c: " + c);
                System.out.println("epsilon: " + epsilon);

                Map<Long, Tuple2<Double, Double>> mapCoordinates = SequentialRecommendersUtils.POICoordinatesMap(originalPOICoordsFile, lp);
				float [][] distanceMatrix = POIProcessData.distanceMatrix(data.getAllItems().collect(Collectors.toList()), data, mapCoordinates, false);
                System.out.println("Distance matrix size: " + data.getAllItems().count() + " x " + data.getAllItems().count());
				
        		System.out.println("--Version Reset GeoInfluenceMatrix--");

                
				//C and epsilon
                return new RankGeoFMRecommender<>(data, distanceMatrix, kFactorizer, kNeighbours, alphaF, c, epsilon, numInteractions, learnRate, maxRate, isboldDriver, decay);        	
            }
            case "RankGeoFMEXSurRec":
            case "RankGeoFMRecommenderEXSur": {
                System.out.println("RankGeoFMRecommenderEXSur");
            	System.out.println("kFactors: " + kFactorizer);
            	System.out.println("kNeighbours " + kNeighbours);
            	System.out.println("numIter: " + numInteractions);
            	System.out.println("learnRate (in the ExpSur model is alpha): " + learnRate);
                System.out.println("alpha (in the ExpSur model factorInf): " + alphaF);
                System.out.println("c: " + c);
                System.out.println("epsilon: " + epsilon);

                Map<Long, Tuple2<Double, Double>> mapCoordinates = SequentialRecommendersUtils.POICoordinatesMap(originalPOICoordsFile, lp);

/*(FastPreferenceData<U, I> prefData,
    			Map<Long, Tuple2<Double, Double>> mapCoordinates, int numberFactors, int knnPois, double factorInf, double C,
    			double epsilon, double alpha, int numberIterations)*/

                return new RankGeoFMRecommenderEXSur<>(data, mapCoordinates, kFactorizer, kNeighbours, alphaF, c, epsilon, learnRate , numInteractions);        	
            }
            
            
            case "RankGeoFMRecND":
            case "RankGeoFMRecommenderNotDense": {
                System.out.println("RankGeoFMRecommenderNotDense");
            	System.out.println("kFactors: " + kFactorizer);
            	System.out.println("kNeighbours " + kNeighbours);
            	System.out.println("decay: " + decay);
            	System.out.println("isboldDriver: " + isboldDriver);
            	System.out.println("numIter: " + numInteractions);
            	System.out.println("learnRate: " + learnRate);
            	System.out.println("maxRate: " + maxRate);
                System.out.println("alpha: " + alphaF);
                System.out.println("c: " + c);
                System.out.println("epsilon: " + epsilon);

                Map<Long, Tuple2<Double, Double>> mapCoordinates = SequentialRecommendersUtils.POICoordinatesMap(originalPOICoordsFile, lp);
				
        		System.out.println("--Version Reset GeoInfluenceMatrix--");

				//C and epsilon
                return new RankGeoFMRecommenderNotDense<>(data, mapCoordinates, kFactorizer, kNeighbours, alphaF, c, epsilon, numInteractions, learnRate, maxRate, isboldDriver, decay);        	
            }
            
            case "GeoBPRMF": {
            	int k = kFactorizer;
                int numIter = numInteractions;
                System.out.println("GeoBPRMF");
            	System.out.println("kFactors: "+ k);
            	System.out.println("numIter: "+ numInteractions);
            	System.out.println("regUser: " + regUser);
            	System.out.println("regItem: " + regItem);
            	System.out.println("learnRate: " + learnRate);
            	System.out.println("regBias: " + regBias);
            	System.out.println("maxDistance: " + maxDistance);
            	
                Map<Long, Tuple2<Double, Double>> mapCoordinates = SequentialRecommendersUtils.POICoordinatesMap(originalPOICoordsFile, lp);

            	GeoBPRMF<Long, Long> geoBPRMF = new GeoBPRMF<>(data, regUser, regItem, regBias, numIter, k, learnRate, mapCoordinates, maxDistance);
            	geoBPRMF.fit();
            	return geoBPRMF;
            	
            }
            

            
            case "FMFMGM": {
        		System.out.println("Alpha: " + alphaF);
        		System.out.println("Theta: " + theta);
        		System.out.println("Maximum distance" + maxDistance);
        		System.out.println("Iterations " + numInteractions);
        		System.out.println("Factors: " + kFactorizer);
        		System.out.println("Alpha2: " + alphaF2);
        		System.out.println("Beta: " + beta);
        		System.out.println("Lerarning rate: " + learnRate);
        		System.out.println("Sigmoid: " + useSigmoid);
        		Map<Long, Tuple2<Double, Double>> mapCoordinates = SequentialRecommendersUtils.POICoordinatesMap(originalPOICoordsFile, lp);

        		
        		FMFMGM<Long, Long> fmfMGMRecommender = new FMFMGM<>(data, mapCoordinates, alphaF, theta, maxDistance, numInteractions, kFactorizer, alphaF2, beta, learnRate, useSigmoid);
        		return fmfMGMRecommender;
            }
            
            
            

            
            //Algorithm that combines popularity, knn and geospatial
            case "PopGeoNN": {
                System.out.println("PopGeoNN");
                System.out.println("kNeighbours: " + kNeighbours);
                System.out.println("Sim: " + ranksysSimilarity);
                
                PopularityRecommender<Long, Long> popRec = new PopularityRecommender<Long, Long>(data);
                Map<Long, Tuple2<Double, Double>> mapCoordinates = POICoordinatesMap(originalPOICoordsFile, lp);
                UsersMidPoints<Long, Long> usermid = new UsersMidPoints<>(data, mapCoordinates, score);

                AverageDistanceToUserGEO<Long, Long> distanceRev = new AverageDistanceToUserGEO<Long, Long>(data, usermid);

                es.uam.eps.ir.ranksys.nn.user.sim.UserSimilarity<Long> simUNR = obtRanksysUserSimilarity(data, ranksysSimilarity);
                if (simUNR == null) {
                    return null;
                }
                UserNeighborhood<Long> urneighborhood = new TopKUserNeighborhood<>(simUNR, kNeighbours);
                UserNeighborhoodRecommender<Long, Long> ubRec = new UserNeighborhoodRecommender<>(data, urneighborhood, 1);

                
                return new PopGeoNN<>(data, popRec, distanceRev, ubRec);
            }


            default:
                System.out.println("Carefull. Recommender is null");
                return null;
        }

    }
    
    
   	public static <U> double getAverageUserTimeStamp(FastTemporalPreferenceDataIF<U, ?> data, U u) {
   		return data.getUidxTimePreferences(data.user2uidx(u)).mapToDouble(IdxTimePref::v3).average().getAsDouble();
   	}

   	public static <I> double getAverageItemTimeStamp(FastTemporalPreferenceDataIF<?, I> data, I i) {
        return data.getIidxTimePreferences(data.item2iidx(i)).mapToDouble(IdxTimePref::v3).average().getAsDouble();

   	}

    
    public static <U, I> Map<U, List<Tuple2<Double, Double>>> getCoordinatesUser(FastPreferenceData<U, I> data,
			Map<I, Tuple2<Double, Double>> mapCoordinates) {
    	Map<U, List<Tuple2<Double, Double>>> result = new HashMap<U, List<Tuple2<Double, Double>>>();
    	
    	data.getUsersWithPreferences().forEach(u -> {
    		result.put(u, new ArrayList<Tuple2<Double, Double>>());
    		List<Tuple2<Double, Double>> coordinatesUser = result.get(u);
    		data.getUserPreferences(u).forEach(pref ->{
    			coordinatesUser.add(mapCoordinates.get(pref.v1));
    		});
    	});
		return result;
	}


	/***
     * RankSys skyline recommenders
     * @param rec
     * @param testData
     * @param dataTemporalTest
     * @param userIndex
     * @param itemIndex
     * @param relevanceThreshold
     * @param antiRelevanceThreshold
     * @return
     */
    public static Recommender<Long, Long> obtRankSysSkylineRecommender(String rec, FastPreferenceData<Long, Long> trainData, FastPreferenceData<Long, Long> testData, FastTemporalPreferenceDataIF<Long, Long> dataTemporalTest, double relevanceThreshold, double antiRelevanceThreshold){
    	switch (rec){
    		case "skyRelRec":
    		case "SkyRelRec":
    		case "skylineRelevanceRecommender": 
    		case "SkylineRelevanceRecommender": {
            	System.out.println("SkylineRelevanceRecommender");
    			return new SkylineRelevanceRecommender<>(testData, relevanceThreshold);
    		}
    		case "SkylineRelevanceRecommenderUsersTrain": {
            	System.out.println("SkylineRelevanceRecommenderUsersTrain");
    			return new SkylineRelevanceRecommenderUsersTrain<>(trainData, testData, relevanceThreshold);
    		}
    		case "skyAntiRelRec":
    		case "SkyAntiRelRec":
    		case "skylineAntiRelevanceRecommender": 
    		case "SkylineAntiRelevanceRecommender": {
            	System.out.println("SkylineAntiRelevanceRecommender");
    			return new SkylineAntiRelevanceRecommender<>(testData, antiRelevanceThreshold);
    		}
    		case "skyTestOrder":
    		case "SkyTestOrder": 
    		case "skylineTestOrder": 
    		case "SkylineTestOrder": {
            	System.out.println("SkylineTestOrder");
            	return new SkylineTestOrder<>(dataTemporalTest, false);
            }
            case "skyTestOrderRev": 
            case "SkyTestOrderRev": 
            case "skylineTestOrderReverse": 
            case "SkylineTestOrderReverse": {
            	System.out.println("SkylineTestOrderReverse");
            	return new SkylineTestOrder<>(dataTemporalTest, true);

            }
    		default:
    			System.out.println("Recommender not recognized. Recommender is null");
    			return null;
    	}
    }
    
    
    
    




	public static es.uam.eps.ir.seqawareev.tour.abstracts.AbstractFastTourRecommender<Long, Long> obtRankSysTourRecommeder(String rec, FastTemporalPreferenceDataIF<Long, Long> prefDataTemporal, FastPreferenceData<Long, Long> prefData, String poiEstimatedTime, String fileCoordinatesItem, 
    		String FeatureFile, SmoothingProbabilityIF smooth, double lambda, int cached, int numberItemsCompute){
    	try {
	    	switch (rec) {

    		
	    		case "NNVenuesRecommender": {
	    			System.out.println("NNVenuesRecommender");
	    			/*
	    			 // Previous
	    			Map<Long, Tuple2<Double,Double>> mapCoordinates = SequentialRecommendersUtils.POICoordinatesMap(fileCoordinatesItem, lp);
	    			float [][] distanceMatrix = POIProcessData.distanceMatrix(prefData.getAllItems().collect(Collectors.toList()), prefData, mapCoordinates);
	    			NNVenuesRecommender<Long, Long> nnVenuesRec = new NNVenuesRecommender<Long, Long>(prefData, distanceMatrix); 
	    			return nnVenuesRec;
	    			*/
	    			Map<Long, Tuple2<Double, Double>> mapCoordinates = SequentialRecommendersUtils.POICoordinatesMap(fileCoordinatesItem, lp);
					float [][] distanceMatrix = POIProcessData.distanceMatrix(prefData.getAllItems().collect(Collectors.toList()), prefData, mapCoordinates, false);
					AbstractTourCachedItemSimilarity<Long> iSim = new ItemLocationSimilarity<>(prefData, distanceMatrix, cached, numberItemsCompute);
	    			return new GenericItemSimilarityTourRecommender<Long, Long>(prefData, iSim);
	    		} 

	    		case "ItemToItemMC" : {
	    			System.out.println("ItemToItemMC");
	    			/*
	    			// Previous
	    			ItemTransitionMatrix<Long, Long> im = new ItemTransitionMatrix<Long, Long>(prefDataTemporal);
	    			im.initialize();
	    			ItemToItemMC<Long, Long> itemToItem = new ItemToItemMC<>(prefDataTemporal, im, smooth);
	    			return itemToItem;
	    			*/
	    			ItemTransitionMatrix<Long, Long> im = new ItemTransitionMatrix<Long, Long>(prefDataTemporal);
	    			im.initialize();
	    			AbstractTourCachedItemSimilarity<Long> iSim = new ItemTransitionSimilarity<>(prefDataTemporal, im, smooth, cached, numberItemsCompute);
	    			return new GenericItemSimilarityTourRecommender<Long, Long>(prefDataTemporal, iSim);	    			
	    		}
	    		case "FeatureToFeature": {
	    			System.out.println("FeatureToFeature");
	    			/*
	    			// Previous
	    			FeatureData<Long, Long, Double> featureData = null;
					featureData = SimpleFeatureData.load(SimpleFeaturesReader.get().read(FeatureFile, lp, lp));
					FeatureItemTransitionMatrix<Long, Long, Long> fItemTrans = new FeatureItemTransitionMatrix<Long, Long, Long>(prefDataTemporal, featureData);
					fItemTrans.initialize();
					FeatureToFeatureMC<Long, Long, Long> fIt = new FeatureToFeatureMC<>(fItemTrans, prefDataTemporal, featureData, smooth);
					return fIt;
					*/
	    			FeatureData<Long, Long, Double> featureData = null;
					featureData = SimpleFeatureData.load(SimpleFeaturesReader.get().read(FeatureFile, lp, lp));
					FeatureItemTransitionMatrix<Long, Long, Long> fItemTrans = new FeatureItemTransitionMatrix<Long, Long, Long>(prefDataTemporal, featureData);
					fItemTrans.initialize();
					AbstractTourCachedItemSimilarity<Long> iSim = new ItemFeatureTransitionSimilarity<>(prefDataTemporal, fItemTrans, smooth, cached, numberItemsCompute);
	    			return new GenericItemSimilarityTourRecommender<Long, Long>(prefDataTemporal, iSim);	    			
	    		}
	    		
	    		default:
	    			System.out.println("Carefull. Tour Recommender is null (recommender not recognized)");
	                return null;
	    	}
    	} catch (Exception e) {
    		e.printStackTrace();
    	}
    	return null;
    }
    
    
    /***
     * Method to compute the midpoint from a list of coordinates
     * @param lstCoordinates
     * @return the midpoint of the coordinates
     */
    public static Tuple2<Double, Double> midPointCoordinates(List<Tuple2<Double, Double>> lstCoordinates){
    	double total = lstCoordinates.size();
    	double x = 0;
    	double y = 0;
    	double z = 0;
    	for (Tuple2<Double,Double> coordinate: lstCoordinates) {
    		double latitudeRad = Math.toRadians(coordinate.v1);
			double longitudeRad = Math.toRadians(coordinate.v2);
			x += Math.cos(latitudeRad) * Math.cos(longitudeRad);
			y += Math.cos(latitudeRad) * Math.sin(longitudeRad);
			z += Math.sin(latitudeRad);
    	}
		x /= total;
		y /= total;
		z /= total;
		
		//Now we have the average in the cartesian map, we must retrieve the original coordinates
		double longitude = Math.atan2(y,x);
		double hyp = Math.sqrt(x * x+ y * y);
		
		double lat2 = Math.atan2(z, hyp);
		
		return new Tuple2<Double,Double>(lat2 * (180.0/Math.PI), longitude * (180.0/Math.PI));
    }
   

    /**
     * *
     * Method to obtain a map of POI- coordinates
     *
     * @param originalPOICoordsFile the original file of POI coordinates
     * @param ip the item parser
     * @return the map
     */
    public static <I> Map<I, Tuple2<Double, Double>> POICoordinatesMap(String originalPOICoordsFile, Parser<I> ip) {
        try {
            Map<I, Tuple2<Double, Double>> result = new HashMap<I, Tuple2<Double, Double>>();
            Stream<String> stream = Files.lines(Paths.get(originalPOICoordsFile));
            stream.forEach(line -> {
                String[] split = line.split("\t");
                I item = ip.parse(split[0]);
                result.put(item, new Tuple2<Double, Double>(Double.parseDouble(split[1]), Double.parseDouble(split[2])));
            });
            stream.close();
            return result;
        } catch (IOException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
        return null;

    }
    
    
    
    
    
    /***
     * Method to store the estimated time 
     * @param originalPOIEstimatedTime the POI estimated time
     * @param ip 
     * @return the map
     */
    public static <I> Map<I, Double> POIEstimatedTime(String originalPOIEstimatedTime, Parser<I> ip){
    	try {
            Map<I, Double> result = new HashMap<I, Double>();
            Stream<String> stream = Files.lines(Paths.get(originalPOIEstimatedTime));
            stream.forEach(line -> {
                String[] split = line.split("\t");
                I item = ip.parse(split[0]);
                result.put(item, Double.parseDouble(split[1]));
            });
            stream.close();
            return result;
        } catch (IOException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
        return null;
    }
    
    /***
     * Method to read the friends of a user from a file
     * @param originalSocialInfluenceFile the path of the social influence file
     * @param up the user parser
     * @return
     */
    public static <U> Map<U, Set<U>> socialInfluenceMap(String originalSocialInfluenceFile, Parser<U> up){
    	try {
            Map<U, Set<U>> result = new HashMap<U, Set<U>>();
            Stream<String> stream = Files.lines(Paths.get(originalSocialInfluenceFile));
            stream.forEach(line -> {
                String[] split = line.split("\t");
                U user1 = up.parse(split[0]);
                U user2 = up.parse(split[1]);
                
                if (result.get(user1) == null) {
                	result.put(user1, new HashSet<>());
                }
                
                if (result.get(user2) == null) {
                	result.put(user2, new HashSet<>());
                }
                
                // I assume that this is symmetrical
                result.get(user1).add(user2);
                result.get(user2).add(user1);
                
            });
            stream.close();
            return result;
        } catch (IOException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
        return null;
    }
    
    

    public static <I, F> Map<I, Map<F, Double>> featureMap(String originalFeatureFile, Parser<I> ip, Parser<F> fp) {
        try {
            Map<I, Map<F, Double>> result = new HashMap<>();
            Stream<String> stream = Files.lines(Paths.get(originalFeatureFile));
            stream.forEach(line -> {
                // item feature weight
                String[] split = line.split("\t");
                I item = ip.parse(split[0]);
                F feat = fp.parse(split[1]);
                result.put(item, result.getOrDefault(item, new HashMap<>()));
                result.get(item).put(feat, Double.parseDouble(split[2]));
            });
            stream.close();
            return result;
        } catch (IOException e) {
            e.printStackTrace();
        }
        return null;
    }

    public static <I, C> Map<I, C> poiCityMap(String poiCityFile, Parser<I> ip, Parser<C> cp) {
        try {
            Map<I, C> result = new HashMap<>();
            Stream<String> stream = Files.lines(Paths.get(poiCityFile));
            stream.forEach(line -> {
                // item city
                String[] split = line.split("\t");
                I item = ip.parse(split[0]);
                C city = cp.parse(split[1]);
                result.put(item, city);
            });
            stream.close();
            return result;
        } catch (IOException e) {
            e.printStackTrace();
        }
        return null;
    }

    public static <C> Map<Tuple2<C, C>, Double> citySimilaritiesMap(String citySimilarityFile, Parser<C> cp) {
        try {
            Map<Tuple2<C, C>, Double> result = new HashMap<>();
            Stream<String> stream = Files.lines(Paths.get(citySimilarityFile));
            stream.forEach(line -> {
                // city1 city2 sim
                String[] split = line.split("\t");
                C city1 = cp.parse(split[0]);
                C city2 = cp.parse(split[1]);
                result.put(new Tuple2<>(city1, city2), Double.parseDouble(split[2]));
            });
            stream.close();
            return result;
        } catch (IOException e) {
            e.printStackTrace();
        }
        return null;
    }

    /**
     * *
     * Method that will read a file with n columns
     *
     * @param <I>
     * @param path the path of the file
     * @param catItem the category of the item (or the city) that the item must
     * match in order to be retrieved
     * @param parserItems the parser of the items
     * @param columnCategory the column category
     * @param columnItem the column of the item
     * @return the stream of the items
     */
    public static <I> Stream<I> obtainCandidatesFile(String path, String catItem, Parser<I> parserItems, int columnCategory, int columnItem) {
        try {
            List<I> result = new ArrayList<>();
            Stream<String> streamRead = Files.lines(Paths.get(path));
            streamRead.forEach(line -> {
                String[] data = line.split("\t");
                String cat = data[columnCategory];
                if (cat.equals(catItem)) {
                    result.add(parserItems.parse(data[columnItem]));
                }
            });
            streamRead.close();
            return result.stream();
        } catch (IOException e) {
            e.printStackTrace();
        }
        return null;
    }
    
    /***
     * 
     * @param itemIndex
     * @param featureData
     * @return
     */
    public static <I, F, V> Map<Integer, Map<F, Double>> featureDataToIndexMap(FastItemIndex<I> itemIndex, FeatureData<I, F, Double> featureData){
    	Map<Integer, Map<F, Double>> result = new HashMap<>();
    	featureData.getAllItems().forEach(item -> {
    		if (itemIndex.containsItem(item)) {
    			int iidx = itemIndex.item2iidx(item);
    			result.put(iidx, new HashMap<>());
    			Map<F, Double> featureItem = result.get(iidx);
    			featureData.getItemFeatures(item).forEach(f -> {
    				featureItem.put(f.v1, f.v2);
    			});
    		}
    		
    	});
    	return result;
    }

    /**
     * Method to obtain a k-core from a original datamodel
     *
     * @param <U>
     * @param <I>
     * @param original the original dataModel
     * @param minimumNumberRatingsUser the minimum number of ratings per user
     * @param minimumNumberRatingsItem the minimum number of ratings per item
     * @return
     */
    public static <U, I> SimpleFastTemporalPreferenceData<U, I> KCore(SimpleFastTemporalPreferenceData<U, I> original, int minimumNumberRatingsUser, int minimumNumberRatingsItem) {
        SimpleFastTemporalPreferenceData<U, I> result = original;

        while (!checkKCore(result, minimumNumberRatingsUser, minimumNumberRatingsItem)) {
            SimpleFastTemporalPreferenceData<U, I> aux = result;
            Set<U> usersFiltered = result.getAllUsers().filter(u -> aux.getUserPreferences(u).count() >= minimumNumberRatingsUser).collect(Collectors.toSet());
            Set<I> itemsFiltered = result.getAllItems().filter(i -> aux.getItemPreferences(i).count() >= minimumNumberRatingsItem).collect(Collectors.toSet());
            if (usersFiltered.isEmpty() || itemsFiltered.isEmpty()) {
                return null;
            }
            FastUserIndex<U> userIndex = SimpleFastUserIndex.load(usersFiltered.stream());
            FastItemIndex<I> itemIndex = SimpleFastItemIndex.load(itemsFiltered.stream());
            Stream<Tuple4<U, I, Double, Long>> prev = original.getAsTemporalTuples();
            Stream<Tuple4<U, I, Double, Long>> prevFiltered = prev.filter(t -> userIndex.containsUser(t.v1) && itemIndex.containsItem(t.v2));

            result = SimpleFastTemporalFeaturePreferenceData.loadTemporalFeature(prevFiltered, userIndex, itemIndex);
        }

        return result;
    }

    /**
     *
     * @param <U>
     * @param <I>
     * @param original
     * @param from
     * @param to
     * @return
     */
    public static <U, I> SimpleFastTemporalPreferenceData<U, I> temporalFilter(SimpleFastTemporalPreferenceData<U, I> original, Long from, Long to) {

        Stream<Tuple4<U, I, Double, Long>> prev = original.getAsTemporalTuples();
        Stream<Tuple4<U, I, Double, Long>> prevFiltered = prev.filter(t -> (t.v4 >= from) && (t.v4 < to));
        List<Tuple4<U, I, Double, Long>> lst = prevFiltered.collect(Collectors.toList());
        Set<U> usersFiltered = new HashSet<>();
        Set<I> itemsFiltered = new HashSet<>();
        lst.stream().forEach(t -> {
            usersFiltered.add(t.v1);
            itemsFiltered.add(t.v2);
        });
        FastUserIndex<U> userIndex = SimpleFastUserIndex.load(usersFiltered.stream());
        FastItemIndex<I> itemIndex = SimpleFastItemIndex.load(itemsFiltered.stream());

        SimpleFastTemporalPreferenceData<U, I> result = SimpleFastTemporalFeaturePreferenceData.loadTemporalFeature(lst.stream(), userIndex, itemIndex);

        return result;
    }

    public static <U, I> SimpleFastTemporalPreferenceData<U, I> removeDuplicates(SimpleFastTemporalPreferenceData<U, I> original) {
        List<Tuple4<U, I, Double, Long>> tuples = new ArrayList<>();
        original.getAllUsers().forEach(u -> {
            Set<I> items = new HashSet<>();
            Comparator<IdTimePref<I>> timeComparatorIdTimePref = (o1, o2) -> o1.v3 != o2.v3 ? Long.compare(o1.v3, o2.v3) : ((Comparable) o1.v1).compareTo(o2.v1);
            original.getUserPreferences(u).sorted(timeComparatorIdTimePref).forEach(p -> {
                if (!items.contains(p.v1)) {
                    tuples.add(new Tuple4<>(u, p.v1, p.v2, p.v3));
                    items.add(p.v1);
                }
            });
        });

        Set<U> usersFiltered = new HashSet<>();
        tuples.forEach(t -> usersFiltered.add(t.v1));
        Set<I> itemsFiltered = new HashSet<>();
        tuples.forEach(t -> itemsFiltered.add(t.v2));
        FastUserIndex<U> userIndex = SimpleFastUserIndex.load(usersFiltered.stream());
        FastItemIndex<I> itemIndex = SimpleFastItemIndex.load(itemsFiltered.stream());
        SimpleFastTemporalPreferenceData<U, I> model = SimpleFastTemporalFeaturePreferenceData.loadTemporalFeature(tuples.stream(), userIndex, itemIndex);
        return model;
    }

    public static <U, I> SimpleFastTemporalPreferenceData<U, I>[] splitData(SimpleFastTemporalPreferenceData<U, I> original, final Long[][] constraints, final Boolean[] allowedRepetitions, boolean filterTestByTrain) {
        SimpleFastTemporalPreferenceData<U, I>[] models = new SimpleFastTemporalPreferenceData[constraints.length];

        for (int i = 0; i < constraints.length; i++) {
            SimpleFastTemporalPreferenceData<U, I> split = temporalFilter(original, constraints[i][0], constraints[i][1]);
            if (!allowedRepetitions[i]) {
                split = removeDuplicates(split);
            }
            models[i] = split;
            // check if this model contains anything that is already in the first split (the test may or may not contain any rating from the train)
            if (i > 0 && filterTestByTrain) {
                Stream<Tuple4<U, I, Double, Long>> prev = split.getAsTemporalTuples();
                // only add the tuple if it is not present in model 0
                Stream<Tuple4<U, I, Double, Long>> prevFiltered = prev.filter(t -> (models[0].getPreferences(t.v1, t.v2) == null || models[0].getPreferences(t.v1, t.v2).count() == 0));
                List<Tuple4<U, I, Double, Long>> lst = prevFiltered.collect(Collectors.toList());
                Set<U> usersFiltered = new HashSet<>();
                Set<I> itemsFiltered = new HashSet<>();
                lst.stream().forEach(t -> {
                    usersFiltered.add(t.v1);
                    itemsFiltered.add(t.v2);
                });
                FastUserIndex<U> userIndex = SimpleFastUserIndex.load(usersFiltered.stream());
                FastItemIndex<I> itemIndex = SimpleFastItemIndex.load(itemsFiltered.stream());
                models[i] = SimpleFastTemporalFeaturePreferenceData.loadTemporalFeature(lst.stream(), userIndex, itemIndex);
            }
        }
        return models;
    }

    public static <U, I> SimpleFastTemporalPreferenceData<U, I>[] temporalPartition(SimpleFastTemporalPreferenceData<U, I> original, final double trainingPercentage) {
        SimpleFastTemporalPreferenceData<U, I>[] models = new SimpleFastTemporalPreferenceData[2];

        int nprefs = original.numPreferences();
        int nprefsTraining = (int) Math.floor(trainingPercentage * nprefs);
        System.out.println("Num preferences for train" + nprefsTraining);
        final List<Tuple4<U, I, Double, Long>> lstTraining = new ArrayList<>();
        final List<Tuple4<U, I, Double, Long>> lstTest = new ArrayList<>();

        AtomicInteger counter = new AtomicInteger(0);
        Stream<Tuple4<U, I, Double, Long>> prev = original.getAsTemporalTuples();
        // sort so first tuples are older than following ones
        Comparator<Tuple4<?, I, ?, Long>> timestampComparator = (o1, o2) -> o1.v4().compareTo(o2.v4());
        prev.sorted(timestampComparator).forEach(t -> {
            if (counter.get() <= nprefsTraining) {
                lstTraining.add(t);
            } else {
                lstTest.add(t);
            }
            counter.incrementAndGet();
        });
        System.out.println("Preferences train " + lstTraining.size());
        System.out.println("Preferences test " + lstTest.size());

        IntStream.range(0, 2).forEach(i -> {
            Set<U> usersFiltered = new HashSet<>();
            Set<I> itemsFiltered = new HashSet<>();
            List<Tuple4<U, I, Double, Long>> lst = (i == 0 ? lstTraining : lstTest);
            lst.stream().forEach(t -> {
                usersFiltered.add(t.v1);
                itemsFiltered.add(t.v2);
            });
            FastUserIndex<U> userIndex = SimpleFastUserIndex.load(usersFiltered.stream());
            FastItemIndex<I> itemIndex = SimpleFastItemIndex.load(itemsFiltered.stream());
            models[i] = SimpleFastTemporalFeaturePreferenceData.loadTemporalFeature(lst.stream(), userIndex, itemIndex);
        });
        return models;
    }

    public static <U, I> SimpleFastTemporalPreferenceData<U, I>[] splitPercentageData(SimpleFastTemporalPreferenceData<U, I> original, final double trainingPercentage, final Boolean[] allowedRepetitions) {
        SimpleFastTemporalPreferenceData<U, I>[] models = new SimpleFastTemporalPreferenceData[2];

        SimpleFastTemporalPreferenceData<U, I>[] split = temporalPartition(original, trainingPercentage);
        for (int i = 0; i < 2; i++) {
            if (!allowedRepetitions[i]) {
                split[i] = removeDuplicates(split[i]);
            }
            models[i] = split[i];
            // check if this model contains anything that is already in the first split
            if (i > 0) {
                Stream<Tuple4<U, I, Double, Long>> prev = split[i].getAsTemporalTuples();
                // only add the tuple if it is not present in model 0
                Stream<Tuple4<U, I, Double, Long>> prevFiltered = prev.filter(t -> (models[0].getPreferences(t.v1, t.v2) == null || models[0].getPreferences(t.v1, t.v2).count() == 0));
                List<Tuple4<U, I, Double, Long>> lst = prevFiltered.collect(Collectors.toList());
                Set<U> usersFiltered = new HashSet<>();
                Set<I> itemsFiltered = new HashSet<>();
                lst.stream().forEach(t -> {
                    usersFiltered.add(t.v1);
                    itemsFiltered.add(t.v2);
                });
                FastUserIndex<U> userIndex = SimpleFastUserIndex.load(usersFiltered.stream());
                FastItemIndex<I> itemIndex = SimpleFastItemIndex.load(itemsFiltered.stream());
                models[i] = SimpleFastTemporalFeaturePreferenceData.loadTemporalFeature(lst.stream(), userIndex, itemIndex);
            }
        }
        return models;
    }

    public static <U, I> SimpleFastTemporalPreferenceData<U, I>[] splitPercentageRemoveDuplicatesAndSplit(
            SimpleFastTemporalPreferenceData<U, I> fullDataRepetitions, double trainingPercentage, Boolean[] allowedRepetitions, int kCoreUsers, int kCoreItems) {
        SimpleFastTemporalPreferenceData<U, I> fullDataNoRepetitions = removeDuplicates(fullDataRepetitions);

        SimpleFastTemporalPreferenceData<U, I> fullDataNoRepetitionsKCore = KCore(fullDataNoRepetitions, kCoreUsers, kCoreItems);
        SimpleFastTemporalPreferenceData<U, I>[] models = new SimpleFastTemporalPreferenceData[2];

        System.out.println("Total preferences " + fullDataNoRepetitionsKCore.numPreferences());
        SimpleFastTemporalPreferenceData<U, I>[] split = temporalPartition(fullDataNoRepetitionsKCore, trainingPercentage);

        for (int i = 0; i < 2; i++) {
            models[i] = split[i];
            // check if this model contains anything that is already in the first split
            if (i > 0) {
                Stream<Tuple4<U, I, Double, Long>> prev = split[i].getAsTemporalTuples();
                System.out.println("Prev to filter " + split[i].getAsTemporalTuples().collect(Collectors.toList()).size());
                // only add the tuple if it is not present in model 0
                Stream<Tuple4<U, I, Double, Long>> prevFiltered = prev.filter(t -> (models[0].getPreferences(t.v1, t.v2) == null || models[0].getPreferences(t.v1, t.v2).count() == 0));
                List<Tuple4<U, I, Double, Long>> lst = prevFiltered.collect(Collectors.toList());
                System.out.println("after to filter " + lst.size());

                Set<U> usersFiltered = new HashSet<>();
                Set<I> itemsFiltered = new HashSet<>();
                lst.stream().forEach(t -> {
                    usersFiltered.add(t.v1);
                    itemsFiltered.add(t.v2);
                });
                FastUserIndex<U> userIndex = SimpleFastUserIndex.load(usersFiltered.stream());
                FastItemIndex<I> itemIndex = SimpleFastItemIndex.load(itemsFiltered.stream());
                models[i] = SimpleFastTemporalFeaturePreferenceData.loadTemporalFeature(lst.stream(), userIndex, itemIndex);
            }
        }
        return models;
    }
    

    


    /**
     * *
     * Method to check if all items and user of a specific dataModel have a
     * certain number of ratings
     *
     * @param dataModel the datamodel
     * @param minimumNumberRatingsUser the minimum number of ratings per user
     * @param minimumNumberRatingsItem the minimum number of ratings per item
     * @return boolean if it matches the condition, false otherwise
     */
    private static <U, I> boolean checkKCore(SimpleFastTemporalPreferenceData<U, I> dataModel, int minimumNumberRatingsUser, int minimumNumberRatingsItem) {
        Optional<U> opUser = dataModel.getAllUsers().filter(u -> dataModel.getUserPreferences(u).count() < minimumNumberRatingsUser).findAny();
        if (opUser.isPresent()) {
            return false;
        }

        Optional<I> opItem = dataModel.getAllItems().filter(i -> dataModel.getItemPreferences(i).count() < minimumNumberRatingsItem).findAny();
        if (opItem.isPresent()) {
            return false;
        }

        return true;
    }
    /***
     * Method to filter the preference data by providing the set of valid users and valid items 
     * @param original the original preferences
     * @param validUsers the set of valid users
     * @param validItems the set of valid items
     * @return the preference data filtered
     */
    public static <U, I> PreferenceData<U, I> filterPreferenceData(PreferenceData<U, I> original, Set<U> validUsers, Set<I> validItems) {
        final List<Tuple3<U, I, Double>> tuples = new ArrayList<>();
        original.getUsersWithPreferences().filter(u -> validUsers.contains(u)).forEach(u -> {
            if (validItems != null) {
                original.getUserPreferences(u).filter(t -> (validItems.contains(t.v1))).forEach(idPref -> {
                    tuples.add(new Tuple3<>(u, idPref.v1, idPref.v2));
                });
            } else {
                original.getUserPreferences(u).forEach(idPref -> {
                    tuples.add(new Tuple3<>(u, idPref.v1, idPref.v2));
                });
            }
        });
        System.out.println("Tuples original: " + original.numPreferences());
        System.out.println("Tuples original filtered: " + tuples.size());
        Stream<Tuple3<U, I, Double>> prev = tuples.stream();
        PreferenceData<U, I> result = SimplePreferenceData.load(prev);
        return result;
    }
    
    
    public static <U, I> List<PreferenceData<U,I>> splitTrainDataBasedOnRatings(final PreferenceData<U, I> originalTrain, double maxValue){
       
    	List<PreferenceData<U,I>> result = new ArrayList<>();


        	List<Tuple3<U, I, Double>> tuplesLess = new ArrayList<>();
        	List<Tuple3<U, I, Double>> tuplesMore = new ArrayList<>();

        	originalTrain.getUsersWithPreferences().forEach(
        			u -> originalTrain.getUserPreferences(u).forEach(pref -> {
        				double val = pref.v2;
        				if (val < maxValue) {
        					tuplesLess.add(new Tuple3<>(u, pref.v1, pref.v2));
        				} else {
        					tuplesMore.add(new Tuple3<>(u, pref.v1, pref.v2));
        				}
        					
        			}));
        	
        	PreferenceData<U, I> fstList = SimplePreferenceData.load(tuplesLess.stream());
        	result.add(fstList);
        	
        	PreferenceData<U, I> sndList = SimplePreferenceData.load(tuplesMore.stream());
        	result.add(sndList);
        
    	
    	return result;

    	
    }
    
    /***
     * Method to filter the test data based on the training interactions
     * @param originalTrain the original training set
     * @param originalTest the original test set
     * @param minUserRatings the minimum ratings of the user in the training set to consider in the test set
     * @param minItemRatings the minimum ratings of the item in the training set to consider in the test set
     * @return the preference data filtered
     */
    public static <U, I> PreferenceData<U, I> filterTestBasedOnTraining(final PreferenceData<U, I> originalTrain, final PreferenceData<U, I> originalTest,
            int minUserRatings, int minItemRatings) {
        // we can filter out from test: 
        //       users with less than n ratings in training
        //       items with less than m ratings in training
        PreferenceData<U, I> tmpTestData = null;
        if ((minUserRatings == 0) && (minItemRatings == 0)) {
            tmpTestData = originalTest;
        } else {
            Set<U> filteredUsers = originalTest.getUsersWithPreferences().filter(u -> originalTrain.numItems(u) >= minUserRatings).collect(Collectors.toSet());
            Set<I> filteredItems = originalTest.getItemsWithPreferences().filter(i -> originalTrain.numUsers(i) >= minItemRatings).collect(Collectors.toSet());
            System.out.println("filteredUsers (" + minUserRatings + ") " + filteredUsers.size());
            System.out.println("filteredItems (" + minItemRatings + ") " + filteredItems.size());
            System.out.println("Original Test: " + originalTest.numPreferences());
            tmpTestData = filterPreferenceData(originalTest, filteredUsers, filteredItems);
        }
        return tmpTestData;
    }
    
    /***
     * Method to filter a preference data using the users and item features.
     * 
     * The filtered feature data will contain only users or items matching the features send as argument
     * 
     * @param originalTrain the original training set
     * @param originalTest the original test set
     * @param minUserRatings the minimum ratings of the user in the training set to consider in the test set
     * @param minItemRatings the minimum ratings of the item in the training set to consider in the test set
     * @return the preference data filtered
     */
    public static <U, I, F> FastPreferenceData<U, I> filterPreferenceDataByUserItemsFeatures(final FastPreferenceData<U, I> originalTrain, 
    		UserFeatureData<U, F, Double> userFeatureData, FeatureData<I, F, Double> itemFeatureData, F userFeature, F itemFeature) {
        
    	FastPreferenceData<U, I> result = null;
        if (((userFeatureData == null) && (itemFeatureData == null)) || (userFeature == null) && (itemFeature == null)) {
        	result = originalTrain;
        } else {
        	//Filtering users
            Set<U> trainUsers = originalTrain.getUsersWithPreferences().collect(Collectors.toSet());
            Set<U> filteredUsers = trainUsers;
            if (userFeatureData != null) {
            	Set<U> usersHavingThatFeature = userFeatureData.getFeatureUsers(userFeature).map(t2 -> t2.v1).collect(Collectors.toSet());
            	filteredUsers.retainAll(usersHavingThatFeature);
            }
            //Filtering items
            Set<I> trainItems = originalTrain.getItemsWithPreferences().collect(Collectors.toSet());
            Set<I> filteredItems = trainItems;
            if (itemFeatureData != null) {
            	Set<I> itemsHavingThatFeature = itemFeatureData.getFeatureItems(itemFeature).map(t2 -> t2.v1).collect(Collectors.toSet());
            	filteredItems.retainAll(itemsHavingThatFeature);
            }
                        
            System.out.println("filteredUsers (" + userFeature + ") " + filteredUsers.size());
            System.out.println("filteredItems (" + itemFeature + ") " + filteredItems.size());
            System.out.println("Original Train: " + originalTrain.numPreferences());
            
            //Now, transform the normal preference data to FAST preference data
            List<U> usersList = filteredUsers.stream().sorted().collect(Collectors.toList());
			List<I> itemsList = filteredItems.stream().sorted().collect(Collectors.toList());
            
            FastUserIndex<U> usersIndex = SimpleFastUserIndex.load(usersList.stream());
			FastItemIndex<I> itemsIndex = SimpleFastItemIndex.load(itemsList.stream());

			List<Tuple3<U, I, Double>> lst = new ArrayList<>();
			
			originalTrain.getUsersWithPreferences().forEach(u -> {
				originalTrain.getUserPreferences(u).forEach(t2 -> {
					if (filteredUsers.contains(u) && filteredItems.contains(t2.v1)) {
						lst.add(new Tuple3<U, I, Double>(u, t2.v1, t2.v2));
					}
				});
			});            
			result = SimpleFastPreferenceData.load(lst.stream(), usersIndex, itemsIndex);

        }
        return result;
    }

    /**
     * *
     * Method to obtain a stream of tuples from a preference data
     *
     * @param prefData the preference data
     * 
     * @return a stream of tuples (user, item, rating)
     */
    public static <U, I> Stream<Tuple3<U, I, Double>> preferenceDataToTuples(PreferenceData<U, I> prefData) {
        List<Tuple3<U, I, Double>> lstPref = new ArrayList<>();
        prefData.getUsersWithPreferences().forEach(u -> {
            lstPref.addAll(prefData.getUserPreferences(u).map(pref -> new Tuple3<>(u, pref.v1, pref.v2)).collect(Collectors.toList()));
        });
        return lstPref.stream();
    }
    
    /***
     * Method to obtain the combiner for the Rank Fusion library
     * @param combiner string identifier of the combiner
     * @return the interface of the combiner
     */
    public static IFCombiner getCombiner(String combiner) {
        IFCombiner comb = null;
        if (combiner.compareTo("sumcomb") == 0) {
            comb = new SumCombiner();
        } else if (combiner.compareTo("mincomb") == 0) {
            comb = new MinCombiner();
        } else if (combiner.compareTo("mnzcomb") == 0) {
            comb = new MNZCombiner();
        } else if (combiner.compareTo("anzcomb") == 0) {
            comb = new AnzCombiner();
        } else if (combiner.compareTo("maxcomb") == 0) {
            comb = new MaxCombiner();
        } else if (combiner.compareTo("medcomb") == 0) {
            comb = new MedCombiner();
        } else if (combiner.compareTo("defaultcomb") == 0) {
            // Default combiner
            comb = new SumCombiner();
        } else {
            // Unrecognized combiner
        	System.out.println(combiner +" not recognized");
            comb = null;
        }
        return comb;
    }

    /**
     * Method to obtain the normalizer of the RankFusion algorithm
     * @param normalizer the string containing the normalizer
     * @param idealDist the ideal distribution
     * @param dataDist the data distribution
     * @return the Normalizer selected
     */
    public static IFNormalizer getNormalizer(String normalizer, Distribution idealDist, Distribution dataDist) {
        IFNormalizer norm = null;
        if (normalizer.compareTo("stdnorm") == 0) {
            norm = new StdNormalizer();
        } else if (normalizer.compareTo("zmuvnorm") == 0) {
            norm = new ZMUVNormalizer();
        } else if (normalizer.compareTo("ranksimnorm") == 0) {
            norm = new RankSimNormalizer();
        } else if (normalizer.compareTo("sumnorm") == 0) {
            norm = new SumNormalizer();
        } else if (normalizer.compareTo("2muvnorm") == 0) {
            norm = new TwoMUVNormalizer();
//            norm = new StdNormalizer();
        } else if (normalizer.compareTo("distrnorm") == 0) {
            Map<String, Distribution> mapDist = new HashMap<String, Distribution>();
            mapDist.put("", dataDist);
            norm = new DistributionNormalizer(idealDist, mapDist);
        } else if (normalizer.compareTo("defaultnorm") == 0) {
            // Default normalizer
            norm = new DummyNormalizer();
        } else {
            // Unrecognized normalizer
        	System.out.println(normalizer +" not recognized");
        	norm = null;
        }
        return norm;
    }
    
    /***
     * Method to obtain the minimum and the maximum ratings for a RankSysPreferenceData
     * @param data the preference data
     * @return a tuple containing the minimum and the maximum rating of the dataset
     */
    public static <U,I> Tuple2<Double, Double> getMinMaxDataRating(PreferenceData<U, I> data){
    	AtomicDouble minRating = new AtomicDouble(Double.MAX_VALUE);
    	AtomicDouble maxRating = new AtomicDouble(Double.MIN_VALUE);
    	data.getUsersWithPreferences().forEach(u -> {
    		double userMin = data.getUserPreferences(u).mapToDouble(timePref -> timePref.v2).min().getAsDouble();
    		double userMax = data.getUserPreferences(u).mapToDouble(timePref -> timePref.v2).max().getAsDouble();
    		minRating.set(userMin < minRating.get() ? userMin: minRating.get());
    		maxRating.set(userMax > maxRating.get() ? userMax: maxRating.get());
    	});
    	return new Tuple2<>(minRating.get(), maxRating.get());
    }
    
    
    /***
     * Method to generate a test set (for mymedialite) from a training set
     * The result file will contain for every user, the set of items that the user have not rated in the training ser
     * @param data
     * @param fileResult
     */
    public static <U, I> void generateTestFromTrain(FastPreferenceData<U,I> data, String fileResult){
    	try {
			PrintWriter writer = new PrintWriter(fileResult);
			data.getUsersWithPreferences().forEach(u -> {
				Set<I> usersItems = data.getUserPreferences(u).map(pref -> pref.v1).collect(Collectors.toSet());
				data.getItemsWithPreferences().forEach(i -> {
					if (!usersItems.contains(i)) {
						writer.println(u + "\t" + i +"\t" + 1);
					}	
				});
			});
			writer.close();
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
    	
    }

    /***
     * Method to order a recommendation raking
     * @param data the preference data
     * @param outputFile the output file
     * @param topN the topN to select
     */
	public static <U, I> void orderTopNRecommendationFile(FastPreferenceData<U, I> data, String outputFile, Integer topN) {
    	try {
			PrintWriter writer = new PrintWriter(outputFile);
			data.getUidxWithPreferences().forEach(uidx -> {
				data.getUidxPreferences(uidx).sorted(PreferenceComparators.preferenceComparatorIdxPref.reversed()).limit(topN).forEach(pref -> {
					writer.println(data.uidx2user(uidx) + "\t" + data.iidx2item(pref.v1) + "\t" + pref.v2);
				});
			});
			writer.close();
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
    
	
	

}
