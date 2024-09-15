
package es.uam.eps.ir.sr.mains;

import static org.ranksys.formats.parsing.Parsers.lp;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintStream;
import java.lang.reflect.InvocationTargetException;
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
import java.util.Map.Entry;
import java.util.Properties;
import java.util.Set;
import java.util.TreeMap;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import org.jooq.lambda.tuple.Tuple2;
import org.jooq.lambda.tuple.Tuple3;
import org.jooq.lambda.tuple.Tuple4;
import org.ranksys.formats.parsing.Parser;
import org.ranksys.formats.preference.SimpleRatingPreferencesReader;
import org.ranksys.formats.rec.RecommendationFormat;
import org.ranksys.formats.rec.SimpleRecommendationFormat;

import com.google.common.util.concurrent.AtomicDouble;

import es.uam.eps.ir.antimetrics.AntiNDCG;
import es.uam.eps.ir.antimetrics.AntiRecommendationMetricComplement;
import es.uam.eps.ir.antimetrics.AntiSuccess;
import es.uam.eps.ir.antimetrics.Success;
import es.uam.eps.ir.antimetrics.antirel.BinaryAntiRelevanceModel;
import es.uam.eps.ir.antimetrics.ratio.AntiRelevanceRatio;
import es.uam.eps.ir.antimetrics.ratio.BorderlineRatio;
import es.uam.eps.ir.antimetrics.ratio.UnknownItemsRatio;
import es.uam.eps.ir.attrrec.datamodel.feature.UserFeatureData;
import es.uam.eps.ir.attrrec.metrics.recommendation.averages.WeightedModelUser.UserMetricWeight;
import es.uam.eps.ir.crossdomainPOI.datamodel.SimpleFastTemporalFeaturePreferenceData;
import es.uam.eps.ir.crossdomainPOI.datamodel.temporal.interfaces.FastTemporalPreferenceDataIF;
import es.uam.eps.ir.crossdomainPOI.datamodel.temporal.utils.ExtendedPreferenceReader;
import es.uam.eps.ir.crossdomainPOI.utils.UsersMidPoints;
import es.uam.eps.ir.crossdomainPOI.utils.UsersMidPoints.SCORES_FREQUENCY;
import es.uam.eps.ir.ranksys.core.feature.FeatureData;
import es.uam.eps.ir.ranksys.core.preference.ConcatPreferenceData;
import es.uam.eps.ir.ranksys.core.preference.PreferenceData;
import es.uam.eps.ir.ranksys.core.preference.SimplePreferenceData;
import es.uam.eps.ir.ranksys.diversity.distance.metrics.EILD;
import es.uam.eps.ir.ranksys.diversity.intentaware.IntentModel;
import es.uam.eps.ir.ranksys.diversity.intentaware.metrics.AlphaNDCG;
import es.uam.eps.ir.ranksys.diversity.intentaware.metrics.ERRIA;
import es.uam.eps.ir.ranksys.fast.index.FastItemIndex;
import es.uam.eps.ir.ranksys.fast.index.FastUserIndex;
import es.uam.eps.ir.ranksys.fast.index.SimpleFastItemIndex;
import es.uam.eps.ir.ranksys.fast.index.SimpleFastUserIndex;
import es.uam.eps.ir.ranksys.fast.preference.FastPreferenceData;
import es.uam.eps.ir.ranksys.fast.preference.SimpleFastPreferenceData;
import es.uam.eps.ir.ranksys.metrics.RecommendationMetric;
import es.uam.eps.ir.ranksys.metrics.basic.AveragePrecision;
import es.uam.eps.ir.ranksys.metrics.basic.NDCG;
import es.uam.eps.ir.ranksys.metrics.basic.ReciprocalRank;
import es.uam.eps.ir.ranksys.metrics.rank.RankingDiscountModel;
import es.uam.eps.ir.ranksys.metrics.rel.BinaryRelevanceModel;
import es.uam.eps.ir.ranksys.metrics.rel.RelevanceModel;
import es.uam.eps.ir.ranksys.nn.item.sim.ItemSimilarity;
import es.uam.eps.ir.ranksys.novdiv.distance.ItemDistanceModel;
import es.uam.eps.ir.ranksys.novelty.longtail.FDItemNovelty;
import es.uam.eps.ir.ranksys.novelty.longtail.PCItemNovelty;
import es.uam.eps.ir.ranksys.novelty.longtail.metrics.EFD;
import es.uam.eps.ir.ranksys.novelty.longtail.metrics.EPC;
import es.uam.eps.ir.ranksys.novelty.temporal.ItemFreshness;
import es.uam.eps.ir.ranksys.novelty.temporal.ItemFreshness.FreshnessMetricNorm;
import es.uam.eps.ir.ranksys.novelty.temporal.ItemFreshness.MetricScheme;
import es.uam.eps.ir.ranksys.novelty.temporal.TimestampCalculator;
import es.uam.eps.ir.ranksys.novelty.temporal.metrics.GenericFreshness;
import es.uam.eps.ir.ranksys.novelty.unexp.PDItemNovelty;
import es.uam.eps.ir.ranksys.novelty.unexp.metrics.EPD;
import es.uam.eps.ir.seqawareev.metrics.DistanceEvaluationMetric;
import es.uam.eps.ir.seqawareev.metrics.FeatureLCSPrecision;
import es.uam.eps.ir.seqawareev.metrics.TestFeaturePrecision;
import es.uam.eps.ir.seqawareev.rerankers.LambdaRerankerRelevance.FILL_RERANKING_STRATEGY;
import es.uam.eps.ir.seqawareev.smooth.JelineckMercer;
import es.uam.eps.ir.seqawareev.smooth.NoSmoothingCond;
import es.uam.eps.ir.seqawareev.smooth.NoSmoothingPrior;
import es.uam.eps.ir.seqawareev.smooth.SmoothingProbabilityIF;
import es.uam.eps.ir.seqawareev.utils.CBBaselinesUtils.ITEMTRANSFORMATION;
import es.uam.eps.ir.seqawareev.utils.CBBaselinesUtils.USERTRANSFORMATION;
import es.uam.eps.ir.seqawareev.utils.UsersSessions.UserSession;
import es.uam.eps.ir.seqawareev.wrappers.SimpleFastTemporalFeaturePreferenceDataWrapper.RepetitionsStrategyPreference;
import es.uam.eps.ir.seqawareev.wrappers.SimpleFastTemporalFeaturePreferenceDataWrapper.RepetitionsStrategyTime;
import es.uam.eps.ir.sr.utils.PredicatesStrategies;
import es.uam.eps.ir.sr.utils.SequentialRecommendersUtils;
import es.uam.eps.ir.sr.utils.TimeStampUtils.TimestampStrategy;
import net.recommenders.rival.core.DataModel;
import net.recommenders.rival.evaluation.metric.EvaluationMetric;
import net.recommenders.rival.evaluation.metric.MultipleEvaluationMetricRunner;
import net.recommenders.rival.evaluation.metric.ranking.AbstractRankingMetric;
import net.recommenders.rival.evaluation.metric.ranking.MAP;
import net.recommenders.rival.evaluation.metric.ranking.PopularityStratifiedRecall;
import net.recommenders.rival.evaluation.metric.ranking.Precision;
import net.recommenders.rival.evaluation.metric.ranking.Recall;

/**
 * This class will have all methods that the experiment will use
 * 
 * @author Pablo Sanchez (psperez@icai.comillas.edu)
 *
 */
public class ExperimentUtils {
	
	private static int COLUMN_USER = 0;
	private static int COLUMN_ITEM = 0;


	public static void commonUsers(String file1, String path, String filesToProcess) throws IOException {
		String[] files = filesToProcess.split(",");
		for (String file : files) {
			commonUsers(file1, path + file);
		}
	}
	
	public static<U, F> boolean isValidUser(U user, Boolean computeFilter,  UserFeatureData<U, F, Double> userFeature, F feature) {
		if (userFeature == null || feature == null || !computeFilter) {
			return true;
		}
		
		return userFeature.getUserFeatures(user).map(tuples -> tuples.v1).anyMatch(f-> f.equals(feature));
	}
	
	
	/***
	 * Method to obtain the number of common users between 2 files
	 * @param file1 the first file
	 * @param file2 the second file
	 * @throws IOException
	 */
	public static void commonUsers(String file1, String file2) throws IOException {
		Set<String> users1 = new HashSet<String>();
		Set<String> users2 = new HashSet<String>();
		Set<String> users3 = new HashSet<String>();

		Stream<String> streamRatings1 = Files.lines(Paths.get(file1));
		streamRatings1.forEach(str -> {
			String[] data1 = str.split("\t");
			users1.add(data1[COLUMN_USER]);
		});
		streamRatings1.close();

		Stream<String> streamRatings2 = Files.lines(Paths.get(file2));
		streamRatings2.forEach(str -> {
			String[] data2 = str.split("\t");
			users2.add(data2[COLUMN_USER]);
		});
		streamRatings2.close();

		users3.addAll(users1);
		users3.addAll(users2);

		AtomicDouble intersection = new AtomicDouble(0);
		users1.forEach(user -> {
			if (users2.contains(user)) {
				intersection.addAndGet(1.0);
			}
		});

		System.out
				.println(file1 + " and " + file2 + " union: " + users3.size() + " intersection: " + intersection.get());
	}

	
	public static void generateLongsFeatureFile(String originalFeatureFile, String outputFeature) {
        File file = new File(outputFeature);
        Map<String, Integer> map = new HashMap<>();
		try (Stream<String> stream = Files.lines(Paths.get(originalFeatureFile))) {
			FileWriter fr = new FileWriter(file);
            
			AtomicInteger i = new AtomicInteger(1);
			stream.forEach(line -> {
				String [] data = line.split("\t"); 
				String feat = data[1];
				if (!map.containsKey(feat)) {
					map.put(feat, i.getAndIncrement());
				}
				try {
					fr.write(data[0] + "\t" + map.get(feat)  + "\t" + data[1]  + "\n");
				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
				
				
			});
			fr.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
	}
	
	public static<F> Map<F, Double> readFeatureDistribution(String featureDistributionFile, Parser<F> lp) {
		Map<F, Double> featureDistributionMap = new HashMap<F, Double>();
		
		try (Stream<String> stream = Files.lines(Paths.get(featureDistributionFile))) {

			stream.forEach(line -> {
				String [] data = line.split("\t"); //We create an array
				featureDistributionMap.put(lp.parse(data[0]), Double.parseDouble(data[1]));
			});

		} catch (IOException e) {
			e.printStackTrace();
		}

		
		
		return featureDistributionMap;
	}
	
	
	public static double estimateConfidence(String datapathTrain, double div) {
		BufferedReader br;
		try {
			br = new BufferedReader(new FileReader(datapathTrain));
			String line = br.readLine(); // First line
			Map<Integer, Integer> resMap = new HashMap<Integer, Integer>(); // Stores the id of the user and the number
			// of the item he has rated
			while (line != null) {
				if (line != null) {
					String[] data = line.split("\t");
					Integer id = Integer.parseInt(data[0]);
					if (resMap.get(id) == null) // New user
					{
						resMap.put(id, 1);
					} else {
						resMap.put(id, resMap.get(id) + 1);
					}
				}
				line = br.readLine();
			}
			br.close();
			// All computed. Sort by value
			Map<Integer, Integer> sorted = sortByComparator(resMap);
			// Now we compute the median
			int index;
			index = (int) (sorted.size() * div);

			int i = 0;
			for (Integer u : sorted.keySet()) {
				i++;
				if (i >= index) {
					return sorted.get(u);
				}
			}

		} catch (FileNotFoundException e) {
			System.out.println("Error loading " + datapathTrain + ". File not found");
		} catch (IOException e) {
			System.out.println("Error loading " + datapathTrain + ". IOException");
		}
		return 0;

	}

	public static Map<Integer, Integer> sortByComparator(Map<Integer, Integer> unsortMap) {
		List<Entry<Integer, Integer>> list = new LinkedList<Entry<Integer, Integer>>(unsortMap.entrySet());

		// Sorting the list based on values
		Collections.sort(list, new Comparator<Entry<Integer, Integer>>() {
			@Override
			public int compare(Entry<Integer, Integer> o1, Entry<Integer, Integer> o2) {
				return o1.getValue().compareTo(o2.getValue());

			}
		});

		// Maintaining insertion order with the help of LinkedList
		Map<Integer, Integer> sortedMap = new LinkedHashMap<Integer, Integer>();
		for (Entry<Integer, Integer> entry : list) {
			sortedMap.put(entry.getKey(), entry.getValue());
		}

		return sortedMap;
	}

	public static Map<String, Map<String, Double>> evaluateStrategy(final Properties properties,
			final DataModel<Long, Long> test, final DataModel<Long, Long> modelToEvaluate)
			throws ClassNotFoundException, IllegalAccessException, InstantiationException, InvocationTargetException,
			NoSuchMethodException {
		Map<String, Map<String, Double>> mapMetricResults = new HashMap<>();
		for (EvaluationMetric<Long> metric : MultipleEvaluationMetricRunner.instantiateEvaluationMetrics(properties,
				modelToEvaluate, test)) {
			metric.compute();
			Double all = metric.getValue();
			Map<String, Double> results = new HashMap<>();
			mapMetricResults.put(metric.toString(), results);
			results.put("all", all);
			Map<Long, Double> perUser = metric.getValuePerUser();
			for (Entry<Long, Double> e : perUser.entrySet()) {
				Long u = e.getKey();
				results.put(u.toString(), e.getValue());
			}
			// is this a ranking metric?
			if (metric instanceof AbstractRankingMetric) {
				@SuppressWarnings("unchecked")
				AbstractRankingMetric<Long, Long> rankingMetric = (AbstractRankingMetric<Long, Long>) metric;
				for (int n : rankingMetric.getCutoffs()) {
					all = rankingMetric.getValueAt(n);
					results = new HashMap<>();
					mapMetricResults.put(metric.toString() + "@" + n, results);
					results.put("all", all);
					perUser = rankingMetric.getValuePerUser();
					for (Long u : perUser.keySet()) {
						results.put(u.toString(), rankingMetric.getValueAt(u, n));
					}
				}
			}
		}
		return mapMetricResults;
	}

	public static Map<String, Map<String, Double>> evaluateStrategy(final Properties properties,
			final DataModel<String, Long> test, final DataModel<String, Long> modelToEvaluate,
			final Map<Long, Integer> observedItemRelevance) throws ClassNotFoundException, IllegalAccessException,
			InstantiationException, InvocationTargetException, NoSuchMethodException {
		double threshold = Double
				.parseDouble(properties.getProperty(MultipleEvaluationMetricRunner.RELEVANCE_THRESHOLD));
		Map<String, Map<String, Double>> mapMetricResults = new HashMap<>();
		List<EvaluationMetric<String>> metrics = new ArrayList<EvaluationMetric<String>>();
		metrics.add(new Recall<>(modelToEvaluate, test, threshold, new int[] { 10 }));
		metrics.add(new Precision<>(modelToEvaluate, test, threshold, new int[] { 10 }));
		metrics.add(new MAP<>(modelToEvaluate, test, threshold, new int[] { 10 }));
		for (EvaluationMetric<String> metric : metrics) {
			metric.compute();
			Double all = metric.getValue();
			Map<String, Double> results = new HashMap<>();
			mapMetricResults.put(metric.toString(), results);
			results.put("all", all);
			Map<String, Double> perUser = metric.getValuePerUser();
			for (Entry<String, Double> e : perUser.entrySet()) {
				String u = e.getKey();
				results.put(u.toString(), e.getValue());
			}
			// is this a ranking metric?
			if (metric instanceof AbstractRankingMetric) {
				@SuppressWarnings("unchecked")
				AbstractRankingMetric<String, Long> rankingMetric = (AbstractRankingMetric<String, Long>) metric;
				for (int n : rankingMetric.getCutoffs()) {
					all = rankingMetric.getValueAt(n);
					results = new HashMap<>();
					mapMetricResults.put(metric.toString() + "@" + n, results);
					results.put("all", all);
					perUser = rankingMetric.getValuePerUser();
					for (String u : perUser.keySet()) {
						results.put(u.toString(), rankingMetric.getValueAt(u, n));
					}
				}
			}
		}
		if (observedItemRelevance != null) {
			// add stratified-recall
			EvaluationMetric<String> metric = new PopularityStratifiedRecall<>(modelToEvaluate, test, threshold,
					new int[] { 10 }, 1.0, observedItemRelevance);
			metric.compute();
			Double all = metric.getValue();
			Map<String, Double> results = new HashMap<>();
			mapMetricResults.put(metric.toString(), results);
			results.put("all", all);
			Map<String, Double> perUser = metric.getValuePerUser();
			for (Entry<String, Double> e : perUser.entrySet()) {
				String u = e.getKey();
				results.put(u.toString(), e.getValue());
			}
			AbstractRankingMetric<String, Long> rankingMetric = (AbstractRankingMetric<String, Long>) metric;
			for (int n : rankingMetric.getCutoffs()) {
				all = rankingMetric.getValueAt(n);
				results = new HashMap<>();
				mapMetricResults.put(metric.toString() + "@" + n, results);
				results.put("all", all);
				perUser = rankingMetric.getValuePerUser();
				for (String u : perUser.keySet()) {
					results.put(u.toString(), rankingMetric.getValueAt(u, n));
				}
			}
		}
		return mapMetricResults;
	}

	/***
	 * Method to obtain the user and the items from the users and items files
	 * 
	 * @param fileTrain the train file
	 * @param fileTest  the test file
	 * @param up        the user parser
	 * @param ip        the item parser
	 * @return
	 */
	public static <U, I> Tuple2<List<U>, List<I>> getCompleteUserItems(String fileTrain, String fileTest, Parser<U> up,
			Parser<I> ip) {
		try {
			PreferenceData<U, I> trainData = SimplePreferenceData
					.load(SimpleRatingPreferencesReader.get().read(fileTrain, up, ip));
			PreferenceData<U, I> testData = SimplePreferenceData
					.load(SimpleRatingPreferencesReader.get().read(fileTest, up, ip));
			PreferenceData<U, I> totalData = new ConcatPreferenceData<>(trainData, testData);
			List<U> usersList = totalData.getAllUsers().collect(Collectors.toList());
			List<I> itemsList = totalData.getAllItems().collect(Collectors.toList());
			
			return new Tuple2<>(usersList.stream().sorted().collect(Collectors.toList()),
					itemsList.stream().sorted().collect(Collectors.toList()));
		} catch (IOException e) {
			e.printStackTrace();
			return null;
		}

	}

	/***
	 * Method to obtain a tuple of users and items from a file
	 * 
	 * @param file the file
	 * @param up   the user parser
	 * @param ip   the item parser
	 * @return
	 */
	public static <U, I> Tuple2<List<U>, List<I>> getUserItemsFromFile(String file, Parser<U> up, Parser<I> ip) {
		try {
			PreferenceData<U, I> data = SimplePreferenceData
					.load(SimpleRatingPreferencesReader.get().read(file, up, ip));
			List<U> usersList = data.getAllUsers().collect(Collectors.toList());
			List<I> itemsList = data.getAllItems().collect(Collectors.toList());

			System.out.println("Ordering by longs");

			return new Tuple2<>(usersList.stream().sorted().collect(Collectors.toList()),
					itemsList.stream().sorted().collect(Collectors.toList()));
		} catch (IOException e) {
			System.out.println("Exception catched");
			e.printStackTrace();
			return null;
		}

	}

	/***
	 * Method to retrieve the train test indexes from a file
	 * 
	 * @param completeIndexes flag indicating if we are going to concatenate the
	 *                        indexes of train and test or not
	 * @param trainFile       the train file
	 * @param testFile        the test file
	 * @param up              the user parser
	 * @param ip              the item parser
	 * @return
	 */
	public static <U, I> Tuple4<List<U>, List<I>, List<U>, List<I>> retrieveTrainTestIndexes(String trainFile, String testFile, boolean completeIndexes, Parser<U> up, Parser<I> ip) {
		List<U> usersTrain = null;
		List<I> itemsTrain = null;
		List<U> usersTest = null;
		List<I> itemsTest = null;
		if (completeIndexes) {
			Tuple2<List<U>, List<I>> userItems = ExperimentUtils.getCompleteUserItems(trainFile, testFile, up, ip);
			usersTrain = userItems.v1;
			itemsTrain = userItems.v2;
			usersTest = userItems.v1;
			itemsTest = userItems.v2;
		} else {
			Tuple2<List<U>, List<I>> userItemsTrain = ExperimentUtils.getUserItemsFromFile(trainFile, up, ip);
			usersTrain = userItemsTrain.v1;
			itemsTrain = userItemsTrain.v2;
			Tuple2<List<U>, List<I>> userItemsTest = ExperimentUtils.getUserItemsFromFile(testFile, up, ip);
			usersTest = userItemsTest.v1;
			itemsTest = userItemsTest.v2;

		}
		return new Tuple4<>(usersTrain, itemsTrain, usersTest, itemsTest);
	}

	/***
	 * Method to obtain a subset of the results ordered by the metric
	 * 
	 * @param resultFile            the result file obtained by the post eval
	 * @param indexesResMetricValue the indexes of the maps. The indexes of the
	 *                              result file name, the index of the metrics and
	 *                              the value
	 * @param filter                the set of strings that the resultFile must
	 *                              match
	 * @param recFamilies           the recommender families to return
	 * @param metricToOrder         the metric that will be used to order the
	 *                              results
	 */
	public static void parseResultsFile(String resultFile, int[] indexesResMetricValue, String[] filter,
			String[] recFamilies, String metricToOrder, String outputFile) {
		try {
			Stream<String> streamCitiesFile = Files.lines(Paths.get(resultFile));
			int indexRecFile = indexesResMetricValue[0];
			int indexMetric = indexesResMetricValue[1];
			int indexResult = indexesResMetricValue[2];
			Map<String, Map<String, Double>> mapMetricsResult = new HashMap<>();

			streamCitiesFile.forEach(line -> {
				String[] data = line.split("\t");
				String recFile = data[indexRecFile];
				String metric = data[indexMetric];
				double value = Double.parseDouble(data[indexResult]);
				if (matches(recFile, filter)) {
					if (mapMetricsResult.get(recFile) == null)
						mapMetricsResult.put(recFile, new HashMap<>());
					mapMetricsResult.get(recFile).put(metric, value);
				}
			});
			streamCitiesFile.close();
			// Map created, now we need to obtain the maximum by the recFamilies

			Map<String, Map<String, Double>> mapBestResults = new HashMap<>();

			// For each recommender in the family
			for (String family : recFamilies) {
				AtomicDouble bestValue = new AtomicDouble(Double.MIN_VALUE);
				String bestRec = " ";
				for (String recommender : mapMetricsResult.keySet()) {
					if (recommender.contains(family)
							&& bestValue.get() < mapMetricsResult.get(recommender).get(metricToOrder)) {
						bestValue.set(mapMetricsResult.get(recommender).get(metricToOrder));
						bestRec = recommender;
					}
				}
				if (bestRec != " ")
					mapBestResults.put(bestRec, mapMetricsResult.get(bestRec));
			}
			// Now, print all the recommenders
			PrintStream out = new PrintStream(outputFile);
			out.print("Recommender\t");
			for (String metric : mapBestResults.get(mapBestResults.keySet().stream().findFirst().get()).keySet()) {
				out.print(metric + "\t");
			}
			out.println();
			for (String recommender : mapBestResults.keySet()) {
				out.print(recommender + "\t");
				for (String metric : mapBestResults.get(recommender).keySet())
					out.print(mapBestResults.get(recommender).get(metric) + "\t");
				out.println();
			}
			out.close();

		} catch (IOException e) {
			System.out.println("Result file not found");
			e.printStackTrace();
		}

	}

	/****
	 * Method to obtain the PerUser recommendations
	 * 
	 * @param recFile            the recommender files
	 * @param perUserMetrics     the metric to compute the value per each user
	 * @param perUserEvaluations the map of metrics, user and result of the per user
	 *                           evaluation
	 */
	public static void evaluateRecommenderFile(final String recFile,
			final Map<String, RecommendationMetric<Long, Long>> perUserMetrics,
			final Map<String, Map<Long, Double>> perUserEvaluations) {
		RecommendationFormat<Long, Long> format = new SimpleRecommendationFormat<>(lp, lp);
		try {
			format.getReader(recFile).readAll().forEach(rec -> {
				Long u = rec.getUser();

				if ((perUserMetrics != null) && (perUserEvaluations != null)) {
					perUserMetrics.entrySet().forEach(e -> {
						perUserEvaluations.put(e.getKey(),
								perUserEvaluations.getOrDefault(e.getKey(), new HashMap<>()));
						perUserEvaluations.get(e.getKey()).put(u, e.getValue().evaluate(rec));
					});
				}

			});
		} catch (IOException ex) {

		}
	}

	private static boolean matches(String origin, String[] matchings) {
		for (String str : matchings) {
			if (!origin.contains(str))
				return false;
		}
		return true;
	}

	
	


	

	/***
	 * Method to obtain the type of user session (distance or time)
	 * @param userSession
	 * @return
	 */
	public static UserSession obtUserSession(String userSession) {
		for (UserSession m : UserSession.values()) {
			if (m.toString().equals(userSession)) {
				return m;
			}
		}
		return null;
	}



	

	/***
	 * Method to obtain the TimeStamp strategy
	 * 
	 * @param the time strategy used in svd read
	 * @return the time strategy for SVD
	 */
	public static TimestampStrategy obsSVDTimeStrategy(String svdTimeStrategy) {
		for (TimestampStrategy svdStrat : TimestampStrategy.values()) {
			if (svdStrat.toString().equals(svdTimeStrategy)) {
				return svdStrat;
			}
		}
		return TimestampStrategy.MAX_TIMESTAMP_TRAIN;
	}

	/**
	 * Method to obtain the content based strategy
	 * 
	 * @param stringUT the content based user transformation read
	 * @return
	 */
	public static USERTRANSFORMATION obtCBUserTransformation(String stringUT) {
		for (USERTRANSFORMATION ut : USERTRANSFORMATION.values()) {
			if (ut.toString().equals(stringUT)) {
				return ut;
			}
		}
		return null;
	}
	
	public static FILL_RERANKING_STRATEGY otFillRerankingStrategy(String stringUT) {
		for (FILL_RERANKING_STRATEGY strat : FILL_RERANKING_STRATEGY.values()) {
			if (strat.toString().equals(stringUT) || strat.getShortName().equals(stringUT)) {
				return strat;
			}
		}
		return null;
	}

	/**
	 * Method to obtain the recommendation strategy
	 * 
	 * @param stringRS the recomendation strategy read
	 * @return the enumeration value or null
	 */
	public static PredicatesStrategies.recommendationStrategy obtRecommendationStrategy(String stringRS) {
		for (PredicatesStrategies.recommendationStrategy ut : PredicatesStrategies.recommendationStrategy.values()) {
			if (ut.toString().equals(stringRS) || ut.getShortName().equals(stringRS)) {
				System.out.println("Selected " + ut.toString());
				return ut;
			}
		}
		return null;
	}
	


	/**
	 * Method to obtain the repetition strategy for the datamodels
	 * 
	 * @param obtScoreFreq the score Frequency read
	 * @return the enumeration value or null
	 */
	public static ITEMTRANSFORMATION obtCBItemTransformation(String stringiT) {
		for (ITEMTRANSFORMATION it : ITEMTRANSFORMATION.values()) {
			if (it.toString().equals(stringiT)) {
				return it;
			}
		}
		return null;
	}

	/**
	 * Method to obtain the repetition strategy for the datamodels (repetition
	 * strategy)
	 * 
	 * @param repStrategy the repetition strategy read
	 * @return the enumeration value or null
	 */
	public static RepetitionsStrategyPreference obtRepetitionsStrategy(String repStrategy) {
		for (RepetitionsStrategyPreference rs : RepetitionsStrategyPreference.values()) {
			if (rs.toString().equals(repStrategy)) {
				return rs;
			}
		}
		return null;
	}

	/**
	 * Method to obtain the repetition strategy for the datamodels (repetition
	 * strategy for the time)
	 * 
	 * @param repStrategyTime with the time read
	 * @return the enumeration value or null
	 */
	public static RepetitionsStrategyTime obtRepetitionsStrategyTime(String repStrategyTime) {
		for (RepetitionsStrategyTime rs : RepetitionsStrategyTime.values()) {
			if (rs.toString().equals(repStrategyTime)) {
				return rs;
			}
		}
		return null;
	}

	/**
	 * Method to obtain the enumeration of the scoreFreq
	 * 
	 * @param obtScoreFreq the score frequency read
	 * @return the enumeration value or null
	 */
	public static SCORES_FREQUENCY obtScoreFreq(String obtScoreFreq) {
		for (SCORES_FREQUENCY sf : SCORES_FREQUENCY.values()) {
			if (sf.toString().equals(obtScoreFreq)) {
				return sf;
			}
		}
		return null;
	}

	

	/**
	 * Method to obtain an array of name of files matching specific patterns
	 * 
	 * @param dirPath      the directory path
	 * @param filePatterns the file patterns
	 * @return the array of name of files mathing that pattern
	 */
	/*
	private static File[] obtainFiles(String dirPath, String filePatterns) {
		File dir = new File(dirPath);

		FilenameFilter fileFilter = new FilenameFilter() {
			@Override
			public boolean accept(File dir, String name) {
				String[] patterns = filePatterns.split(",");
				for (String pattern : patterns) {
					if (!name.contains(pattern))
						return false;
				}
				if (!name.contains("AVERAGE"))
					return true;
				else
					return false;
			}
		};

		File[] files = dir.listFiles(fileFilter);
		return files;
	}*/

	/***
	 * Method to obtain the average of different files mathing a pattern
	 * 
	 * @param dirPath     the directory of the files
	 * @param filePattern the file patterns (patterns that the name of the files
	 *                    have to match). If there is more than one, they are
	 *                    separated by commas
	 * @param resFile     the resultFile of the average of the results
	 */
	/*
	public static void averageResults(String dirPath, String filePattern, String resFile) {

		try {
			File[] files = obtainFiles(dirPath, filePattern);

			double divide = files.length;
			Map<String, Double> map = new LinkedHashMap<String, Double>();
			for (File file : files) {
				System.out.println(file.getPath());
				Stream<String> stream = Files.lines(file.toPath());
				stream.forEach(line -> {
					String[] fullData = line.split("\t");
					if (map.get(fullData[0]) == null)
						map.put(fullData[0], Double.parseDouble(fullData[1]));
					else
						map.put(fullData[0], map.get(fullData[0]) + Double.parseDouble(fullData[1]));

				});
				stream.close();
			}
			PrintStream resultFile = new PrintStream(resFile);
			for (String metric : map.keySet()) {
				resultFile.println(metric + "\t" + (map.get(metric) / divide));
			}
			resultFile.close();
			System.out.println("Writed on " + resFile);

		} catch (IOException e) {
			e.printStackTrace();
		}

	}*/

	/***
	 * Method to create a timeStamp calculator based on the files that we are
	 * working
	 * 
	 * @param trainFile           the train file
	 * @param itemsTimeAvgFile    the file of the average of the timestamps of the
	 *                            items
	 * @param itemsTimeFirstFile  the file containing the first timestamp of the
	 *                            items
	 * @param itemsTimeMedianFile the file of the median of the timestamps of the
	 *                            items
	 * @param itemsTimeLastFile   the file containing the last timestamps of the
	 *                            items
	 * @throws IOException
	 */
	public static TimestampCalculator<Long, Long> createTimeStampCalculator(String trainFile, String itemsTimeAvgFile,
			String itemsTimeFirstFile, String itemsTimeMedianFile, String itemsTimeLastFile) throws IOException {
		TimestampCalculator<Long, Long> its = new TimestampCalculator<>();
		if (itemsTimeAvgFile != null && itemsTimeFirstFile != null && itemsTimeMedianFile != null
				&& itemsTimeLastFile != null) {
			File fAvg = new File(itemsTimeAvgFile);
			File ffirst = new File(itemsTimeFirstFile);
			File fmed = new File(itemsTimeMedianFile);
			File flast = new File(itemsTimeLastFile);
			Map<MetricScheme, String> mapScheme = new HashMap<>();
			mapScheme.put(MetricScheme.AVG, itemsTimeAvgFile);
			mapScheme.put(MetricScheme.FIRST, itemsTimeFirstFile);
			mapScheme.put(MetricScheme.LAST, itemsTimeLastFile);
			mapScheme.put(MetricScheme.MEDIAN, itemsTimeMedianFile);
			if (!fAvg.exists() || !ffirst.exists() || !fmed.exists() || !flast.exists()) {
				System.out.println("COMPUTING all files for time-aware evaluation");
				its.computeInteractionTimes(trainFile, lp, "\t", 1, 3);
				its.save(mapScheme);
			} else {
				System.out.println("Reading all files for time-aware evaluation");
				its.load(mapScheme, lp);
			}
		}
		return its;
	}
	
	/***
	 * Method to get the user model weight
	 * @param userModelWeight
	 * @return the user model weight
	 */
	public static UserMetricWeight obtUserMetricWeight (String userModelWeight) {
		for (UserMetricWeight m : UserMetricWeight.values()) {
			if (m.toString().equals(userModelWeight) || m.getShortName().equals(userModelWeight) ) {
				return m;
			}
		}
		return null;
	}

	

	/***
	 * Method that will be in charge of add the specific recommendation metrics that
	 * we will use in our experiments
	 * 
	 * @param recMetricsAvgRelUsers the Map of recommendation metrics to add
	 * @param threshold             the relevance threshold
	 * @param cutoff                the cutOff of the recommendations
	 * @param trainData             the training data
	 * @param testData              the test data
	 * @param its                   the timeStampCalculator
	 * @param itemsTimeAvgFile      the file containing the avg of timestamps of the
	 *                              items
	 * @param itemsTimeFirstFile    the file containing the first timestamps of the
	 *                              items
	 * @param itemsTimeMedianFile   the file containing the median of timestamps of
	 *                              the items
	 * @param itemsTimeLastFile     the file containing the last timestamps of the
	 *                              items
	 * @param selectedRelevance     the selected relevance for the non accuracy
	 *                              metrics
	 * @param binRel                the binary relevance model (for accuracy
	 *                              metrics)
	 * @param discModel             the discount model
	 * @param featureData           the feature data for other non accuracy metrics
	 * @param dist                  the distance model
	 * @param intentModel           the intent model
	 * @param itemsReleaseDate      the release date of the items (if available)
	 * @param step                  the actual step of the recommendations
	 * @throws IOException
	 */
	public static void addMetrics(Map<String, RecommendationMetric<Long, Long>> recMetricsAvgRelUsers,
			Map<String, RecommendationMetric<Long, Long>> recMetricsAvgAntiRelUsers,
			Map<String, RecommendationMetric<Long, Long>> recMetricsAllRecUsers, 
			Map<String, RecommendationMetric<Long, Long>> recMetricsAllRecUsersForStd, int threshold,
			int thresholdAntiRelevance, int cutoff, PreferenceData<Long, Long> trainData, FastPreferenceData<Long, Long> trainDataFast,
			PreferenceData<Long, Long> testData, TimestampCalculator<Long, Long> its, String itemsTimeAvgFile,
			String itemsTimeFirstFile, String itemsTimeMedianFile, String itemsTimeLastFile,
			RelevanceModel<Long, Long> selectedRelevance, BinaryRelevanceModel<Long, Long> binRel,
			BinaryAntiRelevanceModel<Long, Long> antiBinRel, RankingDiscountModel discModel,
			FeatureData<Long, String, Double> featureData, ItemDistanceModel<Long> dist,
			IntentModel<Long, Long, String> intentModel, String itemsReleaseDate, String step, boolean antimetricsFlags,
			boolean computeOnlyAcc, FastTemporalPreferenceDataIF<Long, Long> dataTestTemporal,
			boolean computeLCSEvaluation, boolean computeRatios, boolean computeDistance, 
			String pathCoordinates, boolean computeNonExact, List<ItemSimilarity<Long>> itemSimilarities, List<String> itemSims, Boolean computePopEvMetrics, Integer bins,
			Map<Long, Double> prefsItem, Map<Long, Double> usersItem, Map<Long, Double> popByRatingItems,
			Map<String, Double> featuresCounterItems, Map<String, Double> featuresCounterPrefs, Map<String, Double> featuresCounterPopByRatingsItems) throws IOException {

		//Classic IR metrics
		recMetricsAvgRelUsers.put("Precision@" + cutoff + "_" + threshold,
				new es.uam.eps.ir.ranksys.metrics.basic.Precision<>(cutoff, binRel));
		recMetricsAvgRelUsers.put("MAP@" + cutoff + "_" + threshold, new AveragePrecision<>(cutoff, binRel));
		recMetricsAvgRelUsers.put("Recall@" + cutoff + "_" + threshold,
				new es.uam.eps.ir.ranksys.metrics.basic.Recall<>(cutoff, binRel));
		recMetricsAvgRelUsers.put("MRR@" + cutoff + "_" + threshold, new ReciprocalRank<>(cutoff, binRel));
		recMetricsAvgRelUsers.put("NDCG@" + cutoff + "_" + threshold,
				new NDCG<>(cutoff, new NDCG.NDCGRelevanceModel<>(false, testData, threshold)));
		recMetricsAvgRelUsers.put("Success@" + cutoff + "_" + threshold, new Success<>(cutoff, binRel));

		
		// Ratios 
		if (computeRatios) {
			recMetricsAllRecUsers.put("BorderlineRatio@" + cutoff + "_" + threshold + "_" + thresholdAntiRelevance,
					new BorderlineRatio<>(cutoff, testData, threshold, thresholdAntiRelevance));
			recMetricsAllRecUsers.put("UnknownItemsRatio@" + cutoff + "_" + threshold + "_" + thresholdAntiRelevance,
					new UnknownItemsRatio<>(cutoff, testData));
			recMetricsAllRecUsers.put("AntiRelevanceRatio@" + cutoff + "_" + thresholdAntiRelevance,
					new AntiRelevanceRatio<>(cutoff, antiBinRel));
			recMetricsAllRecUsers.put("RelevanceRatio@" + cutoff + "_" + threshold,
					new es.uam.eps.ir.ranksys.metrics.basic.Precision<>(cutoff, binRel));
		}
		//Novelty and diversity and popularity metric
		if (!computeOnlyAcc) {
			recMetricsAllRecUsers.put("epc@" + cutoff,
					new EPC<>(cutoff, new PCItemNovelty<>(trainData), selectedRelevance, discModel));
			recMetricsAllRecUsers.put("efd@" + cutoff,
					new EFD<>(cutoff, new FDItemNovelty<>(trainData), selectedRelevance, discModel));	
			
		}
		

			
		//Computing distance between items (if we have coordinates)
		if (computeDistance) {
			Map<Long, Tuple2<Double, Double>> mapCoordinates = SequentialRecommendersUtils.POICoordinatesMap(pathCoordinates, lp);

			
			recMetricsAllRecUsers.put("DistanceEvMetric@" + cutoff, new DistanceEvaluationMetric<>(cutoff, mapCoordinates));
			recMetricsAllRecUsersForStd.put("DistanceEvMetric@" + cutoff, new DistanceEvaluationMetric<>(cutoff, mapCoordinates));
			
			
			
		}
		
		
		// If we have files of first, last median and average of the timestamps of each
		// item (temporal novelty)
		if (itemsTimeAvgFile != null && itemsTimeFirstFile != null && itemsTimeMedianFile != null
				&& itemsTimeLastFile != null) {
			recMetricsAllRecUsers.put("adi@" + cutoff,
					new GenericFreshness<>(cutoff,
							new ItemFreshness<>(trainData, its, FreshnessMetricNorm.NO_NORM, MetricScheme.AVG),
							selectedRelevance, discModel));
			recMetricsAllRecUsers.put("fdi@" + cutoff,
					new GenericFreshness<>(cutoff,
							new ItemFreshness<>(trainData, its, FreshnessMetricNorm.NO_NORM, MetricScheme.FIRST),
							selectedRelevance, discModel));
			recMetricsAllRecUsers.put("mdi@" + cutoff,
					new GenericFreshness<>(cutoff,
							new ItemFreshness<>(trainData, its, FreshnessMetricNorm.NO_NORM, MetricScheme.MEDIAN),
							selectedRelevance, discModel));
			recMetricsAllRecUsers.put("ldi@" + cutoff,
					new GenericFreshness<>(cutoff,
							new ItemFreshness<>(trainData, its, FreshnessMetricNorm.NO_NORM, MetricScheme.LAST),
							selectedRelevance, discModel));

			// Freshness metrics: ADI, FDI, MDI and LDI (MinMaxNormalization)
			recMetricsAllRecUsers.put("adiMinMaxNorm@" + cutoff,
					new GenericFreshness<>(cutoff,
							new ItemFreshness<>(trainData, its, FreshnessMetricNorm.MINMAXNORM, MetricScheme.AVG),
							selectedRelevance, discModel));
			recMetricsAllRecUsers.put("fdMinMaxNorm@" + cutoff,
					new GenericFreshness<>(cutoff,
							new ItemFreshness<>(trainData, its, FreshnessMetricNorm.MINMAXNORM, MetricScheme.FIRST),
							selectedRelevance, discModel));
			recMetricsAllRecUsers.put("mdiMinMaxNorm@" + cutoff,
					new GenericFreshness<>(cutoff,
							new ItemFreshness<>(trainData, its, FreshnessMetricNorm.MINMAXNORM, MetricScheme.MEDIAN),
							selectedRelevance, discModel));
			recMetricsAllRecUsers.put("ldiMinMaxNorm@" + cutoff,
					new GenericFreshness<>(cutoff,
							new ItemFreshness<>(trainData, its, FreshnessMetricNorm.MINMAXNORM, MetricScheme.LAST),
							selectedRelevance, discModel));
		}
		// If we have the item release dates of the items
		if (itemsReleaseDate != null) {
			Map<MetricScheme, String> mapScheme = new HashMap<>();
			mapScheme.put(MetricScheme.FIRST_RELEASE, itemsReleaseDate);
			its.load(mapScheme, lp);
			recMetricsAllRecUsers.put("fdiRelease@" + cutoff, new GenericFreshness<>(cutoff,
					new ItemFreshness<>(trainData, its, FreshnessMetricNorm.NO_NORM, MetricScheme.FIRST_RELEASE),
					selectedRelevance, discModel));
			recMetricsAllRecUsers.put("fdiReleaseMinMaxNorm@" + cutoff, new GenericFreshness<>(cutoff,
					new ItemFreshness<>(trainData, its, FreshnessMetricNorm.MINMAXNORM, MetricScheme.FIRST_RELEASE),
					selectedRelevance, discModel));
		}
		// If we have features (novelty and diversity with features)
		if (!computeOnlyAcc && featureData != null) {
			recMetricsAllRecUsers.put("eild@" + cutoff, new EILD<>(cutoff, dist, selectedRelevance, discModel));
			recMetricsAllRecUsers.put("epd@" + cutoff,
					new EPD<>(cutoff, new PDItemNovelty<>(false, trainData, dist), selectedRelevance, discModel));
			recMetricsAllRecUsers.put("err-ia@" + cutoff,
					new ERRIA<>(cutoff, intentModel, new ERRIA.ERRRelevanceModel<>(false, testData, threshold)));
			recMetricsAllRecUsers.put("a-ndcg@" + cutoff, new AlphaNDCG<>(cutoff, 0.5, featureData, binRel));
			recMetricsAllRecUsers.put("TFP@" + cutoff, new TestFeaturePrecision<>(cutoff, featureData, testData, binRel));
		}
		//Anti relevance metrics
		if (antimetricsFlags) {
			recMetricsAvgAntiRelUsers.put("AntiSuccess@" + cutoff + "_" + thresholdAntiRelevance,
					new AntiRecommendationMetricComplement<>(new AntiSuccess<>(cutoff, antiBinRel)));
			recMetricsAvgAntiRelUsers.put("AntiPrecision@" + cutoff + "_" + thresholdAntiRelevance,
					new AntiRecommendationMetricComplement<>(
							new es.uam.eps.ir.ranksys.metrics.basic.Precision<>(cutoff, antiBinRel)));
			recMetricsAvgAntiRelUsers.put("AntiMAP@" + cutoff + "_" + thresholdAntiRelevance,
					new AntiRecommendationMetricComplement<>(new AveragePrecision<>(cutoff, antiBinRel)));
			recMetricsAvgAntiRelUsers.put("AntiRecall@" + cutoff + "_" + thresholdAntiRelevance,
					new AntiRecommendationMetricComplement<>(
							new es.uam.eps.ir.ranksys.metrics.basic.Recall<>(cutoff, antiBinRel)));
			recMetricsAvgAntiRelUsers.put("AntiMRR@" + cutoff + "_" + thresholdAntiRelevance,
					new AntiRecommendationMetricComplement<>(new ReciprocalRank<>(cutoff, antiBinRel)));
			recMetricsAvgAntiRelUsers.put("AntiNDCG@" + cutoff + "_" + thresholdAntiRelevance,
					new AntiRecommendationMetricComplement<>(new AntiNDCG<>(cutoff,
							new AntiNDCG.NDCGAntiRelevanceModel<>(false, testData, thresholdAntiRelevance))));
		}
		
		
		if (computeLCSEvaluation && !computeOnlyAcc) {
			recMetricsAvgRelUsers.put("FeatureLCSP@" + cutoff + "_" + threshold, new FeatureLCSPrecision<>(cutoff, featureData, dataTestTemporal, binRel));
		}

	}
	
	public static TreeMap<Double, Double> defaultSimsForNonExactTreeMap(){
		TreeMap<Double, Double> result = new TreeMap<>();
		result.put(0.0, 0.0);
		result.put(0.5, 0.25);
		result.put(0.75, 0.5);
		result.put(1.0, 0.75);
		return result;
		
	}

	/***
	 * Method to obtain the type of smoothing probability
	 * @param smooth the smooth string format
	 * @param lambda the lambda value
	 * @return
	 */
	public static SmoothingProbabilityIF obtSmoothingProbability(String smooth, double lambda) {
		if (smooth == null) {
			return null;
		}
		switch (smooth) {
			case "NoSmoothingPrior":
				return new NoSmoothingPrior();
			case "NoSmoothingCond":
				return new NoSmoothingCond();
			case "JelineckMercer":
				return new JelineckMercer(lambda);
			default:
				return null;
		}
	}
	
	
	
	/***
	 * Method to create a SimpleFastTemporalFeaturePreferenceData
	 * @param trainFile the train file
	 * @param testFile the test file
	 * @param completeOrNot if we will also use the test indexes to build the data
	 * @return
	 */
	public static FastTemporalPreferenceDataIF<Long, Long> loadTrainFastTemporalFeaturePreferenceData(String trainFile, String testFile, boolean completeOrNot, boolean useTrainTest){
		try {
			Tuple4<List<Long>, List<Long>, List<Long>, List<Long>> indexes = ExperimentUtils
					.retrieveTrainTestIndexes(trainFile, testFile, completeOrNot, lp, lp);
			List<Long> usersTrain = indexes.v1;
			List<Long> itemsTrain = indexes.v2;
			
			List<Long> usersTest = indexes.v3;
			List<Long> itemsTest = indexes.v4;

			//Retrieving train
			if (useTrainTest) {
				FastUserIndex<Long> userIndexTrain = SimpleFastUserIndex.load(usersTrain.stream());
				FastItemIndex<Long> itemIndexTrain = SimpleFastItemIndex.load(itemsTrain.stream());
			
				Stream<Tuple4<Long, Long, Double, Long>> tuplesStreamTrain = ExtendedPreferenceReader.get()
						.readPreferencesTimeStamp(trainFile, lp, lp);
				
				SimpleFastTemporalFeaturePreferenceData<Long, Long, Long, Double> ranksysTrainTemporal = SimpleFastTemporalFeaturePreferenceData
						.loadTemporalFeature(tuplesStreamTrain, userIndexTrain, itemIndexTrain);
				return ranksysTrainTemporal;
			} else {
				//Retrieving test
				FastUserIndex<Long> userIndexTest = SimpleFastUserIndex.load(usersTest.stream());
				FastItemIndex<Long> itemIndexTest = SimpleFastItemIndex.load(itemsTest.stream());
			
				Stream<Tuple4<Long, Long, Double, Long>> tuplesStreamTest = ExtendedPreferenceReader.get()
						.readPreferencesTimeStamp(testFile, lp, lp);
				
				SimpleFastTemporalFeaturePreferenceData<Long, Long, Long, Double> ranksysTestTemporal = SimpleFastTemporalFeaturePreferenceData
						.loadTemporalFeature(tuplesStreamTest, userIndexTest, itemIndexTest);
				return ranksysTestTemporal;
			}
			
		} catch (IOException e) {
			System.out.println("Error while creating SimpleFastTemporalFeaturePreferenceData");
			e.printStackTrace();
		}
		return null;
	}
	
	
	public static <U,I> FastPreferenceData<U, I> filterFastPreferenceData(FastPreferenceData<U, I> original, Set<U> userCandidates, Set<I> itemCandidates){
		FastUserIndex<U> userIndexFiltered = original;
		FastItemIndex<I> itemIndexFiltered = original;
		
		Set<U> finalUsers = new HashSet<>();
		Set<I> finalItems = new HashSet<>();
		
		if (userCandidates != null) {
			userIndexFiltered = SimpleFastUserIndex.load(original.getAllUsers().filter(u -> userCandidates.contains(u)));
		}
		
		if (itemCandidates != null) {
			itemIndexFiltered = SimpleFastItemIndex.load(original.getAllItems().filter(i -> itemCandidates.contains(i)));
		}
		
		FastUserIndex<U> userIndexFilteredAux = userIndexFiltered;
		FastItemIndex<I> itemIndexFilteredAux = itemIndexFiltered;
		List<Tuple3<U, I, Double>> preferences = new ArrayList<>();
		original.getUsersWithPreferences().forEach(u -> {
			if (userIndexFilteredAux.containsUser(u)) {
				original.getUserPreferences(u).forEach(pref -> {
					if (itemIndexFilteredAux.containsItem(pref.v1)) {
						preferences.add(new Tuple3<>(u, pref.v1, pref.v2));
						finalUsers.add(u);
						finalItems.add(pref.v1);
					}
				});
			}
		});
		
		// final filtering for preferenceData
		FastUserIndex<U> finalUsersIndx = SimpleFastUserIndex.load(finalUsers.stream());
		FastItemIndex<I> finalItemsIndx = SimpleFastItemIndex.load(finalItems.stream());
		
		return SimpleFastPreferenceData.load(preferences.stream(), finalUsersIndx, finalItemsIndx);
	}
	
	/***
	 * Method to create a SimpleFastTemporalFeaturePreferenceData
	 * @param trainFile the train file
	 * @param testFile the test file
	 * @param completeOrNot if we will also use the test indexes to build the data
	 * @return
	 */
	public static FastPreferenceData<Long, Long> loadTrainFastPreferenceData(String trainFile, String testFile, boolean completeOrNot, boolean useTrainTest){
		try {
			Tuple4<List<Long>, List<Long>, List<Long>, List<Long>> indexes = ExperimentUtils
					.retrieveTrainTestIndexes(trainFile, testFile, completeOrNot, lp, lp);
			List<Long> usersTrain = indexes.v1;
			List<Long> itemsTrain = indexes.v2;
			List<Long> usersTest = indexes.v3;
			List<Long> itemsTest = indexes.v4;
			
			//Retrieving train
			if (useTrainTest) {
				FastUserIndex<Long> userIndexTrain = SimpleFastUserIndex.load(usersTrain.stream());
				FastItemIndex<Long> itemIndexTrain = SimpleFastItemIndex.load(itemsTrain.stream());
				
				return SimpleFastPreferenceData.load(SimpleRatingPreferencesReader.get().read(trainFile, lp, lp), userIndexTrain, itemIndexTrain);
			} else {
				//Retrieving test
				FastUserIndex<Long> userIndexTest = SimpleFastUserIndex.load(usersTest.stream());
				FastItemIndex<Long> itemIndexTest = SimpleFastItemIndex.load(itemsTest.stream());

				return SimpleFastPreferenceData.load(SimpleRatingPreferencesReader.get().read(testFile, lp, lp), userIndexTest, itemIndexTest);
			}
		} catch (IOException e) {
			System.out.println("Error while creating SimpleFastTemporalFeaturePreferenceData");
			e.printStackTrace();
		}
		return null;
	}
	
	/***
	 * Method for checking that the reranked file and the original recommender file contain the same number of items for every user and that the first item in
	 * both list match and that all the items in the reranked file are in the original file
	 * 
	 * @param originalRecommender the original recommender
	 * @param rerankerRecommender the reranked recommender
	 */
	public static void checkRerankCorrectionFirstItemAndSameCandidates(String originalRecommender, String rerankerRecommender) {
		RecommendationFormat<Long, Long> format = new SimpleRecommendationFormat<>(lp, lp);
		try {
			Map<Long, List<Long>> mapUserItemsRatedOriginal = new HashMap<>();
			format.getReader(originalRecommender).readAll().forEach(rec -> {
				List<Long> itemsRatedOriginal = rec.getItems().stream().map(t -> t.v1).collect(Collectors.toList());
				mapUserItemsRatedOriginal.put(rec.getUser(), itemsRatedOriginal);
			});
			
			AtomicBoolean firstItemMatch = new AtomicBoolean(true);
			AtomicBoolean allCandidates = new AtomicBoolean(true);
			AtomicBoolean sameNumber = new AtomicBoolean(true);

			format.getReader(rerankerRecommender).readAll().forEach(rec -> {
				List<Long> originalCandidates = mapUserItemsRatedOriginal.get(rec.getUser());
				if (!originalCandidates.get(0).equals(rec.getItems().get(0).v1)) {
					firstItemMatch.set(false);
				}
				
				if (originalCandidates.size() != rec.getItems().size()) {
					sameNumber.set(false);
				}
				
				rec.getItems().stream().forEach(tuple -> {
					if (!originalCandidates.contains(tuple.v1)) {
						allCandidates.set(false);
					}
				});
				
				
			});
			if (firstItemMatch.get() && allCandidates.get() && sameNumber.get()) {
				System.out.println("All items of each user in the reranked recommender are in the original recommendation for that user");
				System.out.println("The first item recommended in the original recommender is also the first item in the reranked recommender");
				System.out.println("The number of items retrieved for each user in both recommenders in the same");

			}
			
			if (!firstItemMatch.get()) {
				System.out.println("For some users, the first item in the reranked recommended is not the same as the first in the original recommender");
			}
			
			if (!allCandidates.get()) {
				System.out.println("For some users, some of the reranked items are not in the original file");
			}
			
			if (!sameNumber.get()) {
				System.out.println("For some users, the number of reranked items and the number of recommended items is not the same");
			}
			
			
		} catch (Exception e) {
			e.printStackTrace();
			
		}
			
		
	}

	

	
	
	

}
