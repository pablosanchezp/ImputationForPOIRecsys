package es.uam.eps.ir.sr.data;

import static org.ranksys.formats.parsing.Parsers.lp;
import static org.ranksys.formats.parsing.Parsers.sp;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintStream;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Calendar;
import java.util.Collections;
import java.util.Comparator;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;
import java.util.TimeZone;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.AtomicLong;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import org.apache.commons.lang.time.DateUtils;
import org.jooq.lambda.tuple.Tuple2;
import org.jooq.lambda.tuple.Tuple3;
import org.ranksys.core.util.tuples.Tuple2od;
import org.ranksys.formats.feature.SimpleFeaturesReader;
import org.ranksys.formats.parsing.Parser;
import org.ranksys.formats.rec.RecommendationFormat;
import org.ranksys.formats.rec.SimpleRecommendationFormat;

import es.uam.eps.ir.attrrec.datamodel.feature.SimpleUserFeatureData;
import es.uam.eps.ir.attrrec.datamodel.feature.SimpleUserFeaturesReader;
import es.uam.eps.ir.attrrec.datamodel.feature.UserFeatureData;
import es.uam.eps.ir.crossdomainPOI.datamodel.FastTemporalPreferenceData;
import es.uam.eps.ir.crossdomainPOI.datamodel.temporal.interfaces.FastTemporalPreferenceDataIF;
import es.uam.eps.ir.crossdomainPOI.datamodel.temporal.preferences.IdTimePref;
import es.uam.eps.ir.crossdomainPOI.datamodel.temporal.preferences.IdxTimePref;

import com.google.common.util.concurrent.AtomicDouble;

import es.uam.eps.ir.ranksys.core.Recommendation;
import es.uam.eps.ir.ranksys.core.feature.FeatureData;
import es.uam.eps.ir.ranksys.core.feature.SimpleFeatureData;
import es.uam.eps.ir.ranksys.core.preference.PreferenceData;
import es.uam.eps.ir.ranksys.fast.index.FastUserIndex;
import es.uam.eps.ir.ranksys.fast.preference.FastPreferenceData;
import es.uam.eps.ir.ranksys.nn.item.sim.ItemSimilarity;
import es.uam.eps.ir.ranksys.novdiv.distance.ItemDistanceModel;
import es.uam.eps.ir.seqawareev.utils.UsersSessions;
import es.uam.eps.ir.seqawareev.utils.UsersSessions.UserSession;
import es.uam.eps.ir.sr.data.ProcessData.Preference;
import es.uam.eps.ir.sr.mains.ExperimentUtils;
import es.uam.eps.ir.sr.utils.TimeStampUtils;
import es.uam.eps.ir.sr.utils.comparators.PreferenceComparators;

/***
 * Static class that contains methods to process data. Methods to work with
 * externalFiles:
 * 
 * -DatasetReductionRatings: performs a k-core of a dataset (i.e only consider
 * items that have been rated at least by k users and considers only users that
 * have rated at least k ratings. There cannot be users or items with less than
 * that ratings)
 * 
 * -DatasetTransformation: method to transform a datasets to new ids (both users
 * and items)
 * 
 * -Stats: method to compute the statistics of the dataset (i.e. number of
 * items, number of ratings, averages etc).
 * 
 * -TimeSplitPerUser: method to compute a dataset split per user (all users
 * having at least k ratings will have other items in test)
 * 
 * -DatasetTemporalGlobalSplit: method to compute a global split of the dataset
 * (a percentage of newest ratings will go to test)
 * 
 * -obtainBestRecommenders: method to obtain a summary file of the best
 * recommenders. It also works with the user coverage non accuracy metric.
 *
 * -parseMyMediaLite: method to parse myMedialite file of recommendations and
 * prints an output that may be interpreted by rival
 * 
 * @author Pablo Sanchez (psperez@icai.comillas.edu)
 *
 */
public final class ProcessData {
	
	/****
	 * Method to obtain the number of Preferences of the items in the preference data
	 * @param prefData
	 * @return
	 */
	public static <U, I, F, V> Map<F, Double> featureCountsByPreferences(FeatureData<I, F, V> featureData, PreferenceData<U, I> prefData) {
		Map<F, Double> counterMap = new HashMap<>();
		
		featureData.getFeaturesWithItems().forEach(f -> {
			AtomicInteger counterFeature = new AtomicInteger();
			featureData.getFeatureItems(f).forEach(t -> {
				counterFeature.addAndGet((int) prefData.getItemPreferences(t.v1).count());
			});
			counterMap.put(f, (double) counterFeature.get());
		});
		
		return counterMap;
	}
	
	
	/****
	 * Method to obtain the number of Preferences of the items in the preference data
	 * @param prefData
	 * @return
	 */
	public static <U, I, F, V> Map<F, Double> featureCountsByItems(FeatureData<I, F, V> featureData, PreferenceData<U, I> prefData) {
		Map<F, Double> counterMap = new HashMap<>();
		
		featureData.getFeaturesWithItems().forEach(f -> {
			Set<I> itemsFeature = new HashSet<>();
			featureData.getFeatureItems(f).forEach(t -> {
				if (prefData.containsItem(t.v1)){
					itemsFeature.add(t.v1);
				}
			});
			counterMap.put(f, (double) itemsFeature.size());
		});
		
		return counterMap;
	}
	
	/****
	 * Method to obtain the preferences of a feature by adding all the ratings of the items
	 * @param prefData
	 * @return
	 */
	public static <U, I, F, V> Map<F, Double> featureCountsPopByRatings(FeatureData<I, F, V> featureData, PreferenceData<U, I> prefData) {
		Map<F, Double> counterMap = new HashMap<>();
		
		featureData.getFeaturesWithItems().forEach(f -> {
			Set<I> itemsFeature = new HashSet<>();
			featureData.getFeatureItems(f).forEach(t -> {
				if (prefData.containsItem(t.v1)){
					itemsFeature.add(t.v1);
				}
			});
			
			//Now, for every item, we add all the preferences of the feature (by ratings)
			Double totalPop = 0.0;
			for (I item: itemsFeature) {
				totalPop += prefData.getItemPreferences(item).mapToDouble(pref -> pref.v2).sum();
			}
			
			counterMap.put(f, totalPop);
		});
		
		return counterMap;
	}
	
	
	/****
	 * Method to obtain the number of Preferences of the items in the preference data
	 * @param prefData
	 * @return
	 */
	public static <U, I> Map<I, Double> numberPreferencesItems(PreferenceData<U, I> prefData){
		Map<I, Double> counterMap = new HashMap<>();
		prefData.getItemsWithPreferences().forEach(item -> {
			counterMap.put(item, (double) prefData.getItemPreferences(item).count());
		});
		return counterMap;
	}
	
	/***
	 * Method to obtain the number of users that have rated that item
	 * @param prefData
	 * @return
	 */
	public static <U, I> Map<I, Double> numberUsersRatedItems(PreferenceData<U, I> prefData){
		Map<I, Double> counterMap = new HashMap<>();
		prefData.getItemsWithPreferences().forEach(item -> {
			counterMap.put(item, (double) getNumberUsersWithPreferencesInItem(prefData, item));
		});
		return counterMap;
	}
	
	/***
	 * Method to obtain the sum of preferences of the items by rating
	 * @param prefData
	 * @return
	 */
	public static <U, I> Map<I, Double> popByRatingItems(PreferenceData<U, I> prefData){
		Map<I, Double> counterMap = new HashMap<>();
		prefData.getItemsWithPreferences().forEach(item -> {
			counterMap.put(item, prefData.getItemPreferences(item).mapToDouble(pref -> pref.v2).sum());
		});
		return counterMap;
	}
	
	
	/***
	 * Method to obtain the number of users that have a preference for item I
	 * @param prefData
	 * @param item
	 * @return
	 */
	public static <U, I> int getNumberUsersWithPreferencesInItem(PreferenceData<U, I> prefData, I item) {
		Set<U> differentUsers = new HashSet<>();
		prefData.getItemPreferences(item).forEach(pref -> {
			differentUsers.add(pref.v1);
		});
		return differentUsers.size();
	}
	
	
	
	/***
	 * Method to obtain the sum of all preferences of a preference data
	 * @param prefData
	 * @param item
	 * @return
	 */
	public static <U, I> double getTotalSumPreferences(PreferenceData<U, I> prefData) {
		AtomicDouble sumAllPrefs = new AtomicDouble(0);
		
		prefData.getUsersWithPreferences().forEach(u -> {
			double sumPrefsUser = prefData.getUserPreferences(u).mapToDouble(pref -> pref.v2).sum();
			sumAllPrefs.addAndGet(sumPrefsUser);
		});
		
		return sumAllPrefs.get();
	}
	
	
	
    public static <U, I> void popularityAnalisisAt(FastPreferenceData<U, I> preferenceDataTrain, String recommendedFile, String outputFile, int cutOff, Parser<U> parserUser, Parser<I> parserItem) {
		
    	
    	//Will store the number of times that every item have been recommended
    	Map<I, Integer> numberUserRecommendedItem = new HashMap<>();
    	Set<U> totalUsers = new HashSet<>();
    	
    	RecommendationFormat<U, I> format = new SimpleRecommendationFormat<>(parserUser, parserItem);

    	
    	try {
			format.getReader(recommendedFile).readAll()
					.forEach(rec -> {
						totalUsers.add(rec.getUser());
						rec.getItems().stream().limit(cutOff).forEach(tuple -> {
							I item = tuple.v1;
							if (!numberUserRecommendedItem.containsKey(item)) {
								numberUserRecommendedItem.put(item, 0);
							}
							numberUserRecommendedItem.put(item, numberUserRecommendedItem.get(item) + 1);
						});
					});
			
			double counterPopByUsers = 0;
			double counterPopByPrefs = 0;
			PrintStream out = new PrintStream(new File(outputFile));
			
			out.println("ItemId" + "\t" + "NºUsersRecommendedItem" + "\t" + "NºUsersRatedItem" + "\t" + "NºPreferencesItem");

			
	    	Map<I, Double> itemCounterPopByPrefs = new HashMap<>();
	    	Map<I, Double> itemCounterPopByUsers = new HashMap<>();


			for (I item: numberUserRecommendedItem.keySet()) {
				
				//Total preferences of the item
				Double preferencesItem = (double) preferenceDataTrain.getItemPreferences(item).count();
				
				
			
				//Real number of users rated that item (number of different users)
				Double numUsersItem = (double) ProcessData.getNumberUsersWithPreferencesInItem(preferenceDataTrain, item);

				//To accumulate
				double popByPrefs = preferencesItem * numberUserRecommendedItem.get(item);
				counterPopByPrefs += popByPrefs;
				itemCounterPopByPrefs.put(item, popByPrefs);
				
				double popByUsers = numUsersItem * numberUserRecommendedItem.get(item);
				counterPopByUsers += popByUsers;
				itemCounterPopByUsers.put(item, popByUsers);
				
				
				out.println(item + "\t" + numberUserRecommendedItem.get(item) + "\t" + numUsersItem + "\t" + preferencesItem);
			}
			double stdPopByPrefs = 0;
			double stdPopByUsers = 0;
			
			double avgPopByPrefs = counterPopByPrefs / numberUserRecommendedItem.keySet().size();
			double avgPopByUsers = counterPopByUsers / numberUserRecommendedItem.keySet().size();

			for (I item: numberUserRecommendedItem.keySet()) {
				stdPopByPrefs += Math.pow(itemCounterPopByPrefs.get(item) - avgPopByPrefs, 2);
				stdPopByUsers += Math.pow(itemCounterPopByUsers.get(item) - avgPopByUsers, 2);
			}
			
			stdPopByPrefs = Math.sqrt(stdPopByPrefs /(numberUserRecommendedItem.keySet().size() - 1));
			stdPopByUsers = Math.sqrt(stdPopByUsers /(numberUserRecommendedItem.keySet().size() - 1));



			
			
			out.println("all" + "\t" + "NºUsersRecommended" + "\t" + "2º*3ºColumn" + "\t" + "std(2º*3ºColumn)" + "\t" + "2º*4ºColumn" + "\t" + "std(2º*4ºColumn)");
			out.println("all" + "\t" + totalUsers.size() + "\t" + counterPopByUsers + "\t" + stdPopByUsers + "\t" + counterPopByPrefs + "\t" + stdPopByPrefs);
			out.close();
			
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
    }

	
	/**
	 * Method that will obtain a list of items ordered by their aggregated similarity 
	 * @param itemSims the list of similarities
	 * @param maxItems the maximum items to retreive
	 * @return
	 */
	public static <I extends Comparable<I>> List<I> majoritySims(List<ItemSimilarity<I>> itemSims, I sourceItem, int maxItems){
		Map<I, Double> counterMap = new HashMap<>();
		
		for (ItemSimilarity<I> iSim: itemSims) {
			iSim.getAllItems().forEach(item -> {
				Double sim = iSim.similarity(item, sourceItem);
				if (counterMap.get(item) == null) {
					counterMap.put(item, 0.0);
				}
				counterMap.put(item , counterMap.get(item) + sim); 
			});
		}
		
		List<Tuple2od<I>> lst = counterMap.keySet().stream().map(i -> new Tuple2od<>(i , counterMap.get(i))).collect(Collectors.toList());
		
		Collections.sort(lst, PreferenceComparators.recommendationComparatorTuple2od());
		Collections.reverse(lst);
		
		List<I> result = lst.stream().map(t -> t.v1).limit(maxItems).collect(Collectors.toList());
		return result;
	}
	
	
	
	
	public static void metricPair(String testFile, String perUserFile, String [] metrics, String outputFile){
		try {
			PrintStream outStream = new PrintStream(outputFile);
			Stream<String> stream = Files.lines(Paths.get(testFile));
			
			Map<String, Integer> userPreferences = new HashMap<>();
			Map<String, Map<String, Double>> mapUserMetrics = new HashMap<>();
			Set<String> metricsSet = new LinkedHashSet<>();
			for (String metric: metrics) {
				metricsSet.add(metric);
			}
			
			//Store the number of preferences of the users
			stream.forEach(line -> {
				String [] data = line.split("\t");
				String user = data[0];
				if (!userPreferences.containsKey(user)) {
					userPreferences.put(user, 0);
				}
				userPreferences.put(user, userPreferences.get(user) + 1);		
			});
			stream.close();
			
			
			//Store the metrics that we will use for each user
			Stream<String> stream2 = Files.lines(Paths.get(perUserFile));
			stream2.forEach(line -> {
				String data [] = line.split("\t");
				String metric = data[0];
				String user = data[1];
				Double value = Double.parseDouble(data[2]);
				
				// only store the metrics the we will use
				if (metricsSet.contains(metric)) {
					if (!mapUserMetrics.containsKey(user)) {
						mapUserMetrics.put(user, new HashMap<>());
					}
					mapUserMetrics.get(user).put(metric, value);
				}
			});
			stream2.close();
			//header
			outStream.print("User\tTestLength");
			for (String metric : metricsSet) {
				outStream.print("\t" + metric);
			}
			outStream.println();
			
			//each line
			for (String user: mapUserMetrics.keySet()) {
				outStream.print(user + "\t" + userPreferences.get(user));
				for (String metric : metricsSet) {
					outStream.print("\t" + mapUserMetrics.get(user).get(metric));
				}
				outStream.println();
			}
			
			

			outStream.close();
			
			
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	/***
	 * Method to "clean" the feature of the items changing them tolower and deleting spaces
	 * @param fileFeatureItem
	 * @param outputFeatureItem
	 */
	public static void parseFeatureItem(String fileFeatureItem, String outputFeatureItem) {

		try {
			PrintStream outStream = new PrintStream(outputFeatureItem);
			Stream<String> stream = Files.lines(Paths.get(fileFeatureItem));
			stream.forEach(line -> {
				String data [] = line.split("\t");
				if (data.length  == 2) {
					String itemId = data [0];
					String feature = data[1];
					feature = feature.toLowerCase();
					feature = feature.replace(" ", "_");
					outStream.println(itemId + "\t" + feature);
				}
			});
			outStream.close();
			stream.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		
	}
	
	
	/***
	 * Method to group the user age into an interval of ages
	 * @param featureUserFile the feature user file
	 * @param outputUserFile the output user file with the new ages
	 * @param agesOrdered the set of ages
	 */
	public static void groupUserAges(String featureUserFile, String outputUserFile, TreeSet<Integer> agesOrdered) {
		Stream<String> stream;
		try {
			PrintStream outStream = new PrintStream(outputUserFile);
			stream = Files.lines(Paths.get(featureUserFile));
			stream.forEach(line -> {	
				String data [] = line.split("\t");
				String user = data[0];
				try {
					Integer age = Integer.parseInt(data[1]);
					Integer newAge = agesOrdered.floor(age);
					if (newAge != null) {
						outStream.println(user + "\t" + newAge);
					}
				} catch (NumberFormatException exception) {
					
				}

			});
			stream.close();
			outStream.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		
	}
	
	
	/**
	 * Method to parse a preference file and change the item ids to a new mapping (new idsItems)
	 * 
	 * @param preferenceFile the original preference file
	 * @param outputItemMapping the output item mapping
	 * @param outputResultFile the result file with new items Ids
	 */
	public static void itemNewIds(String preferenceFile, String outputItemMapping, String outputResultFile) {
		Map<String, Integer> itemMapping = new HashMap<>();
		AtomicInteger itemCounter = new AtomicInteger(1);
		try {
			PrintStream outputFile = new PrintStream(outputResultFile);
			PrintStream outputItemParsedFile = new PrintStream(outputItemMapping);

			Stream<String> stream;
			stream = Files.lines(Paths.get(preferenceFile));
			stream.forEach(line -> {
				String [] data = line.split("\t");
				String user = data[0];
				String item = data[1];
				String rating = data[2];

				if (!itemMapping.containsKey(item)) {
					itemMapping.put(item, itemCounter.getAndIncrement());
				}
				
				outputFile.println(user + "\t" + itemMapping.get(item) + "\t" + rating);
				
			});
			for (String item: itemMapping.keySet()) {
				outputItemParsedFile.println(item + "\t" + itemMapping.get(item));
			}
			
			stream.close();
			outputFile.close();
			outputItemParsedFile.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}
	
	/***
	 * Method to filter a preference file by set of specific feature of the users.
	 * The result preference file will not have any preference of the users having any feature indicated in the set of invalidInfoFeatures
	 *  
	 * @param trainFile the original preference file
	 * @param userInfoFile the user info file
	 * @param invalidInfoFeatures the set of invalid features
	 * @param outputResultFile the result preference file
	 */
	public static void filterPreferenceFileByUserFeature(String preferenceFile, String userInfoFile, Set<String> invalidInfoFeatures, String outputResultFile) {
		Map<String, String> userFeature = new HashMap<>();
		Stream<String> stream;
		try {
			stream = Files.lines(Paths.get(userInfoFile));
			stream.forEach(line -> {
				String [] data = line.split("\t");
				String user = data[0];
				String feature = data[1];
				userFeature.put(user, feature);
			});
			stream.close();
			PrintStream ptr = new PrintStream(outputResultFile);
			stream = Files.lines(Paths.get(preferenceFile));
			stream.forEach(line -> {
				String [] data = line.split("\t");
				String user = data[0];
				String feature = userFeature.get(user);
				if (feature != null && !invalidInfoFeatures.contains(feature)) {
					ptr.println(line);
				}
			});
			stream.close();
			ptr.close();
			
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		
	}
	
	public static void statsUserFeature(String recFile, String userInfoFile) {
		Map<String, String> userFeature = new HashMap<>();
		Stream<String> stream;
		Set<String> featuresUser = new HashSet<>();
		Map<String, Integer> featuresUsersCount = new HashMap<>();
		try {
			stream = Files.lines(Paths.get(userInfoFile));
			stream.forEach(line -> {
				String [] data = line.split("\t");
				String user = data[0];
				String feature = data[1];
				userFeature.put(user, feature);
				featuresUser.add(feature);
			});
			stream.close();
			
			int none = 0;
			for (String feature: featuresUser) {
				featuresUsersCount.put(feature, 0);
			}
			
			Set<String> distinctUsers = new HashSet<>();
			stream = Files.lines(Paths.get(recFile));
			stream.forEach(line -> {
				String [] data = line.split("\t");
				String user = data[0];
				distinctUsers.add(user);
			});
			
			for (String user: distinctUsers) {
				String feat = userFeature.get(user);
				if (feat == null) {
					none++;
				} else {
					featuresUsersCount.put(feat, featuresUsersCount.get(feat) + 1);
				}				
			}
			for (String feature: featuresUser) {
				System.out.println("Users with feature " + feature + " " + featuresUsersCount.get(feature));
			}
			System.out.println("Users none feature " + none);


			
			
			
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
			
	}

	
	
	
	/***
	 * Method to filter a preference file by rating. The result preference file will have
	 * ratings higher than the scoreValue 
	 * @param inputPreferenceFile the input preference file
	 * @param scoreValue the minimum value to consider
	 * @param outputFile the result output file
	 */
	public static void filterPreferenceFileByRating(String inputPreferenceFile, int scoreValue, String outputFile) {
		try {
			Stream<String> stream;
			PrintStream outPutResult = new PrintStream(outputFile);

			stream = Files.lines(Paths.get(inputPreferenceFile));
			stream.forEach(line -> {
				String [] data = line.split("\t");
				double score = Double.parseDouble(data[2]);
				if (score > scoreValue) {
					outPutResult.println(line);
				}
				
			});
			outPutResult.close();
			stream.close();
			
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	/***
	 * Method to generate a file of users marking them as cold start users by the number of preferences
	 * @param prefData the preference data
	 * @param outputFile the output file (format userid \t COLD\NOCOLD)
	 * @param maxPreferencesForCold the maximum number of preferences to be considered as a cold start user
	 */
	public static<U, I> void markColdStart(PreferenceData<U, I> prefData, String outputFile, int maxPreferencesForCold) {
    	try {
			PrintStream outPutResult = new PrintStream(outputFile);
			prefData.getUsersWithPreferences().forEach(u -> {
				if (prefData.getUserPreferences(u).count() <= maxPreferencesForCold) {
					outPutResult.println(u + "\t" + "COLD");
				} else {
					outPutResult.println(u + "\t" + "NOCOLD");
				}
				
			});
			outPutResult.close();
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	public static void markUsersUFValuesQuartiles(String userFeatureFile, String outputResult, int column) {

		try {
			PrintStream outPutResult = new PrintStream(outputResult);

			Stream<String> stream = Files.lines(Paths.get(userFeatureFile));
			Map<String, Integer> valuesUser = new HashMap<String, Integer>();
			List<Number> listQuartiles = new ArrayList<Number>();
			
			stream.forEach(line -> {
				String [] data = line.split("\t");
				String idUser = data[0];
				Integer valueUser = Integer.parseInt(data[column]);
				valuesUser.put(idUser, valueUser);
				listQuartiles.add(valueUser);
			});
			
			List<Number> quartiles = TimeStampUtils.getQuartiles(listQuartiles);
			for (String user: valuesUser.keySet()) {
				Integer valueUser = valuesUser.get(user);
				
				if (valueUser <= quartiles.get(0).intValue()) {
					outPutResult.println(user + "\t" + "Q1");
				} else {
					if (valueUser <= quartiles.get(1).intValue()) {
						outPutResult.println(user + "\t" + "Q2");
					} else {
						if (valueUser <= quartiles.get(2).intValue()) {
							outPutResult.println(user + "\t" + "Q3");
						} else {
							outPutResult.println(user + "\t" + "Q4");
						}
					}
					
				}
			}
			
			
			
			stream.close();
			outPutResult.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		
	}
	
	public static<U, I> void markUsersByNumberPreferencesQuartiles(PreferenceData<U, I> prefData, String outputFile) {
		try {
			PrintStream outPutResult = new PrintStream(outputFile);
			List<Number> lstPreferences = new ArrayList<>();
			prefData.getUsersWithPreferences().forEach(u -> {
				lstPreferences.add(prefData.getUserPreferences(u).count());
			});
			
			List<Number> quartiles = TimeStampUtils.getQuartiles(lstPreferences);
			System.out.println("Quartiles: " + quartiles);
			prefData.getUsersWithPreferences().forEach(u -> {
				long prefsUser = prefData.getUserPreferences(u).count();
				//First
				if (prefsUser <= quartiles.get(0).longValue()) {
					outPutResult.println(u + "\t" + "Q1");
				} else {
					if (prefsUser <= quartiles.get(1).longValue()) {
						outPutResult.println(u + "\t" + "Q2");
					} else {
						if (prefsUser <= quartiles.get(2).longValue()) {
							outPutResult.println(u + "\t" + "Q3");
						} else {
							outPutResult.println(u + "\t" + "Q4");
						}
					}
					
				}
				
				
			});
			outPutResult.close();
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	
	
	/***
	 * Method to cumpute statistics using a file of user information
	 * @param inputPreferenceFile the input preference file
	 * @param userInfo the user info File
	 * @param outputStats the output file for statistics
	 */
	public static void userInfoFileStats(String inputPreferenceFile, String userInfo, String outputStats) {
		//User -> feature
		Map<String, String> mapUser = new HashMap<>();
		//User -> number of preferences
		Map<String, Integer> userCounter = new HashMap<>();
		//Feature -> counterUsers
		Map<String, Integer> mapFeatureCounterUserInfoFile = new HashMap<>();

		
		AtomicInteger numberPreferences = new AtomicInteger();
		
		
		try {
			Stream<String> stream;
			//First, analyze the preferences of the users
			stream = Files.lines(Paths.get(inputPreferenceFile));
			stream.forEach(line -> {
				String data [] = line.split("\t");
				String user = data[0];
				
				if (userCounter.get(user) == null) {
					userCounter.put(user, 0);
				}
				userCounter.put(user, userCounter.get(user) + 1);
				numberPreferences.incrementAndGet();
			});
			stream.close();
			
			Stream<String> stream2;
			stream2 = Files.lines(Paths.get(userInfo));
			//We will take into account only the users that are in the training set
			stream2.forEach(line -> {
				String data [] = line.split("\t");
				String user = data[0];
				String feature = data[1];
				
				//Only for the users appearing in the file
				if (userCounter.containsKey(user)) {
					if (mapFeatureCounterUserInfoFile.get(feature) == null) {
						mapFeatureCounterUserInfoFile.put(feature, 0);
					}
					mapUser.put(user, feature);
					
					mapFeatureCounterUserInfoFile.put(feature, mapFeatureCounterUserInfoFile.get(feature) + 1);
				}
			});
			stream2.close();
			PrintStream out = new PrintStream(outputStats);
			out.println("User features and their percentages");
			for (String feature: mapFeatureCounterUserInfoFile.keySet()) {
				out.println(mapFeatureCounterUserInfoFile.get(feature) + " users have feature " + feature + " from a total of " + userCounter.keySet().size());
			}
			out.println("Preferences and user features");
			out.println("Number of TOTAL preferences " + numberPreferences.get());
			for (String feature: mapFeatureCounterUserInfoFile.keySet()) {
				out.println(feature + " feature");
				int preferencesThatFeature = 0;
				for (String user: mapUser.keySet()) {
					if (mapUser.get(user).equals(feature)) {
						preferencesThatFeature+= userCounter.get(user);
					}
				}
				out.println(preferencesThatFeature + " preferences associated with user feature " + feature);
			}
			out.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		

		
	}
	
	
	/***
	 * Method to group the preferences of lastfm by a specific column. The lastfm dataset contains
	 * repeated listetings, so we can group all the data by an specific column to create a more suitable
	 * data for recommenders
	 * @param inputFile
	 * @param columnGroup
	 * @param outputFile
	 */
	public static void lastFmGroupByColumn(String inputFile, int columnGroup, String outputFile) {
		Map<String, Map<String, Integer>> mapping = new HashMap<>();
		
    	PrintStream outPutResult;
		try {
    		Stream<String> stream = Files.lines(Paths.get(inputFile));
    		stream.forEach(line -> {
    			String data [] = line.split("\t");
    			String user = data[0];
    			String groupingItem = data[columnGroup];
    			if (data.length != 6 || groupingItem.equals("") || groupingItem.equals(" ")) {
    				System.out.println(data);
    			}
    			
    			if (mapping.get(user) == null) {
    				mapping.put(user, new HashMap<>());
    			}
    			
    			
    			if (mapping.get(user).get(groupingItem) == null) {
    				mapping.get(user).put(groupingItem, 0);
    			}
    			int c = mapping.get(user).get(groupingItem);
				mapping.get(user).put(groupingItem, c + 1);
    		});
    		stream.close();
			outPutResult = new PrintStream(outputFile);
			for (String user: mapping.keySet()) {
				for (Map.Entry<String, Integer> entry : mapping.get(user).entrySet()) {
					outPutResult.println(user + "\t" + entry.getKey() + "\t" + entry.getValue());
				}
			}
	    	outPutResult.close();
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	
	/***
     * Method for printing the session of the users of a dataset
     * @param us the user sessions 
     * @param pathResult the result file
     * @param maxDiff the maximum difference allowed (distance or time)
     * @param userSession if the session is for distance or time
     */
    public static <U, I extends Comparable<I>> void printSessions(UsersSessions<U, I> us, String pathResult, double maxDiff, double maxDiffFstLst, UserSession userSession){
    	PrintStream outPutResult;
		try {
			outPutResult = new PrintStream(pathResult);
    		AtomicInteger sessionId = new AtomicInteger(1);
    		
    		//Print the sessions ordered by user Id
	    	us.getData().getUsersWithPreferences().sorted().forEach(user -> {
	    		List<List<IdTimePref<I>>> listUsession = us.getUserSession(user, maxDiff, maxDiffFstLst, userSession);
	    		for (List<IdTimePref<I>> lst : listUsession) {
		    		for (IdTimePref<I> pref: lst) {
		    			outPutResult.println(user + "\t" + pref.v1 + "\t" + pref.v2 + "\t" + pref.v3 + "\t" + sessionId.get());
		    		}
		    		sessionId.incrementAndGet();
	    		}
	    	});
	    	outPutResult.close();
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
    }
	
    
    /***
     * Method to obtain a user feature file with the new Ids of the users
     * 
     * @param checkinFile the checkin file with the new ids
     * @param userMapping the user mapping in the format oldid \t new id
     * @param userFeatureFile the original feature file with the old user ids
     * @param outputUserFeatureFile the output feature id, containing the new ids
     * of the users and the rest of the features
     */
    public static void filterUserInfoFile(String checkinFile, String userMapping, String userFeatureFile, String outputUserFeatureFile) {
    	Map<String, Long> oldNewUserMapping = new HashMap<>();
    	POIProcessData.readMap(userMapping, oldNewUserMapping);
    	
    	Set<Long> usersCheckinFile = new HashSet<>();
    	Set<String> usersProcessed = new HashSet<>(); //avoid repetitions
    	
    	try {
        	PrintStream outPutResult = new PrintStream(outputUserFeatureFile);
    		Stream<String> stream = Files.lines(Paths.get(checkinFile));
    		stream.forEach(line -> {
    			Long newIdUser= Long.parseLong(line.split("\t")[0]);
    			usersCheckinFile.add(newIdUser);
    		});
    		stream.close();
    		
    		// Now we have the users of our checkin file stored, we now see how many of them appear in the user file    		
    		Stream<String> stream2 = Files.lines(Paths.get(userFeatureFile));
    		stream2.forEach(line -> {
    			String data [] = line.split("\t");
    			String oldIdUser = data[0];
    			
    			Long userNewId = oldNewUserMapping.get(oldIdUser);
    			if (userNewId == null || !usersCheckinFile.contains(userNewId)) {
    				//System.out.println("User " + data[0] + " is not in our dataset");
    			} else {
    				if (!usersProcessed.contains(oldIdUser)) {
	    				usersProcessed.add(oldIdUser);
	    				outPutResult.print(userNewId);
	    				for (int i = 1; i < data.length; i++) {
	    					outPutResult.print("\t" + data[i]);
	    				}
	    				outPutResult.println();
    				}
    				
    			}
    		});
    		
    		stream2.close();
    		outPutResult.close();
    	} catch (IOException e) {
			e.printStackTrace();
		}
    }
    
    public static void filterCategoryInfoFile(String checkinFile, String completePoiFeatureFile, String outputPOIFeatureFile) {

    	Set<Long> poisProcessed = new HashSet<>(); //avoid repetitions
    	
		FeatureData<Long, String, Double> featureData = null;

    	
    	try {
    		featureData = SimpleFeatureData.load(SimpleFeaturesReader.get().read(completePoiFeatureFile, lp, sp));
        	PrintStream outPutResult = new PrintStream(outputPOIFeatureFile);
    		Stream<String> stream = Files.lines(Paths.get(checkinFile));
    		stream.forEach(line -> {
    			Long item = Long.parseLong(line.split("\t")[1]);
    			poisProcessed.add(item);
    		});
    		stream.close();
    		
    		for (Long item : poisProcessed) {
    			featureData.getItemFeatures(item).forEach(feat -> {
    				outPutResult.println(item + "\t" + feat.v1 + "\t" + feat.v2);
    			});
    		}
    		
    		
    		
    		outPutResult.close();
    	} catch (IOException e) {
			e.printStackTrace();
		}
    }
    
    
    
    
    /***
     * Method that will receive am user feature file (containing more than one feature) 
     * and will obtain a set of user features files
     * @param inputFeatureFile the input user file
     * @param outputFeatureFiles 
     */
    public static void splitFeatureUserFile(String inputFeatureFile, String [] outputFeatureFiles) {
		try {
			List<PrintStream> ptrs = new ArrayList<>();			
			Stream<String> stream = Files.lines(Paths.get(inputFeatureFile));
			for (String destPath: outputFeatureFiles) {
				ptrs.add(new PrintStream(destPath));
			}
			stream.forEach(line -> {
				String data [] = line.split("\t");
				
				//data [0] is the user
				for (int i = 0; i < outputFeatureFiles.length; i++) {
					ptrs.get(i).println(data[0] + "\t" + data[i + 1]);
				}
				
			});
			stream.close();
			for (PrintStream ptr: ptrs) {
				ptr.close();
			}
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

    	
    }
    
    
    /***
     * Method to parse the feature file of movielens to a format that can be used in ranksys
     * @param featureFile the feature file of movielens
     * @param outputResultFile the output file
     * @param outputfeatureMap the feature map 
     */
    public static void parseFeaturesMovielens(String moviesInfoFiles, String outputResultFile, String outputfeatureMap) {
    	AtomicInteger counter = new AtomicInteger(0);
    	Map<String, Integer> featuresIndex = new HashMap<>();
    	try {
    		PrintStream outputResultFilePS = new PrintStream(outputResultFile);
    		PrintStream outputfeatureMapPS = new PrintStream(outputfeatureMap);
    		
    		Stream<String> stream = Files.lines(Paths.get(moviesInfoFiles));
    		stream.forEach(line -> {
    			String data [] = line.split("::");
    			String genres [] = data[2].split("\\|");
    			
    			for (String genre: genres) {
    				if (featuresIndex.get(genre) == null) {
    					featuresIndex.put(genre, counter.incrementAndGet());
    				}
    				outputResultFilePS.println(data[0] + "\t" + featuresIndex.get(genre));
    			}
    		});
    		outputResultFilePS.close();
    		stream.close();
    		
    		for (Map.Entry<String, Integer> entry : featuresIndex.entrySet()) {
    			outputfeatureMapPS.println(entry.getKey() + "\t" + entry.getValue());
    		}    		
    		outputfeatureMapPS.close();
    		
    	} catch (IOException e) {
			e.printStackTrace();
		}
    }
    
    /**
     * Method to filter an original file by the users id
     * @param originalFile the original file
     * @param userFile the user file
     * @param outputResFile the output result file
     */
    public static void filterFileByUserFile(String originalFile, String userFile, String outputResFile) {
    	Set<Long> usersCandidates = new HashSet<>();

    	try {
    		PrintStream outPutResult = new PrintStream(outputResFile);
    		Stream<String> stream = Files.lines(Paths.get(userFile));
    		stream.forEach(line -> {
    			Long idUser = Long.parseLong(line.split("\t")[0]);
    			usersCandidates.add(idUser);
    		});
    		stream.close();
    		// we have all the candidates, now we filter the original file
    		
    		Stream<String> stream2 = Files.lines(Paths.get(originalFile));
    		stream2.forEach(line -> {
    			Long idUser = Long.parseLong(line.split("\t")[0]);
    			if (usersCandidates.contains(idUser)) {
    				outPutResult.println(line);
    			}
    		});
    		stream2.close();
    		outPutResult.close();
    	} catch (IOException e) {
			e.printStackTrace();
		}

    }
    
    
    
	/***
	 * This method will print the items appearing in both featureFile and trainFile in the outputFeatureFile
	 * It is basically a filter of a big featureFile
	 * @param trainFile the trainFile
	 * @param featureFile the featureFile
	 * @param outputFeatureFile the outputFeatureFile
	 */
    public static void specificFeatureFile(String trainFile, String featureFile, String itemMappingFile, String outputFeatureFile) {
    	int columnItemFeat = 0;
    	int columnItemTrain = 1;
    	int columnFeat = 1;
    	Map<String, Long> oldNewItemMapping = new HashMap<>();
    	if (itemMappingFile != null) {
    		POIProcessData.readMap(itemMappingFile, oldNewItemMapping);
    	}
    	
    	Set<String> itemsTrainNoMapping = new HashSet<>();
    	Set<Long> itemsTrainMapping = new HashSet<>();

    	
    	try {
        	PrintStream out = new PrintStream(outputFeatureFile);
    		Stream<String> stream = Files.lines(Paths.get(trainFile));
    		stream.forEach(line -> {
    			String [] data = line.split("\t");
    			String idItem = data[columnItemTrain];
    			// If the mapping is null, then the ids of the items do not need to be parsed
    			if (itemMappingFile == null) {
    				itemsTrainNoMapping.add(idItem);
    			} else {
    				itemsTrainMapping.add(Long.parseLong(idItem));
    			}
    		});
    		stream.close();
    		Stream<String> stream2 = Files.lines(Paths.get(featureFile));
    		stream2.forEach(line -> {
    			String [] data = line.split("\t");
    			//The item ids of the featuee file are the old ones
    			String idItem = data[columnItemFeat];
    			if (itemMappingFile == null) {
	    			if (itemsTrainNoMapping.contains(idItem)) {
	    				out.println(data[columnItemFeat] + "\t" + data[columnFeat]);
	    			}
    			} else {
    				Long idItemParsed = oldNewItemMapping.get(idItem);
    				if (itemsTrainMapping.contains(idItemParsed)) {
	    				out.println(idItemParsed + "\t" + data[columnFeat]);
	    			}
    			}
    		});
    		stream2.close();

    		
    		
    		out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
    	
    	
    }
    
	

	/***
	 * Method to perform a simple random split
	 * 
	 * @param source          the source file
	 * @param destTrain       the destination train
	 * @param destTest        the destination test
	 * @param percentajeTrain the percentage of ratings that will go to the train
	 *                        and the test set
	 */
	public static void randomSplit(String source, String destTrain, String destTest, double percentageTrain, long seed) {
		try {
			PrintStream resultFile1 = new PrintStream(destTrain);
			PrintStream resultFile2 = new PrintStream(destTest);

			Stream<String> stream = Files.lines(Paths.get(source));
			Random n = new Random();
			n.setSeed(seed);
			stream.forEach(line -> {
				if (n.nextFloat() > percentageTrain) {
					resultFile2.println(line);
				}
				else {
					resultFile1.println(line);
				}
			});
			stream.close();
			resultFile1.close();
			resultFile2.close();

		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	/***
	 * Method to reduce the data of a dataset (making a k core)
	 * 
	 * @param srcpath              the source path of dataset
	 * @param dstPath              the destination route of the dataset
	 * @param numberOfRatingsUser  ratings necessary of a user to be maintained in
	 *                             the destination dataset
	 * @param numberOfRatingsItems ratings necessary of an item to be maintained in
	 *                             the destination dataset
	 */
	public static void DatasetReductionRatings(String srcpath, String dstPath, int numberOfRatingsUser,
			int numberOfRatingsItems, int numberOfUniqueRatingsUser, int numberOfUniqueRatingsItems) {
		boolean implicit = true;
		// Assume structure of the file is: UserdId(long) ItemID(long)
		// rating(double) timestamp(long)

		// Map of users -> userID and list of ItemID-rating-timestamp
		Map<String, List<Tuple3<String, Double, String>>> users = new TreeMap<String, List<Tuple3<String, Double, String>>>();
		// Map of items -> itemID and list of UserID-rating-timestamp
		Map<String, List<Tuple3<String, Double, String>>> items = new TreeMap<String, List<Tuple3<String, Double, String>>>();

		String characterSplit = "\t";

		// Parameters to configure
		int columnUser = 0;
		int columnItem = 1;
		int columnRating = 2;
		int columnTimeStamp = 2;
		// switch to true

		PrintStream writer = null;

		try (Stream<String> stream = Files.lines(Paths.get(srcpath))) {
			stream.forEach(line -> {
				String[] data = line.split(characterSplit);
				String idUser = data[columnUser];
				String idItem = data[columnItem];
				
				double rating = 0;

				if (!implicit) {
					rating = Double.parseDouble((data[columnRating]));
				}
				else {
					rating = 1.0;
				}
				String timestamp;
				
				if (data.length > 3) {
					timestamp = data[columnTimeStamp] + "\t" + data[columnTimeStamp + 1];
				} else {
					timestamp = "1";
				}
				
				
				List<Tuple3<String, Double, String>> lstu = users.get(idUser);
				if (lstu == null) { // New user
					lstu = new ArrayList<>();
					// Add item to list
					users.put(idUser, lstu);
				}
				lstu.add(new Tuple3<String, Double, String>(idItem, rating, timestamp));

				List<Tuple3<String, Double, String>> lsti = items.get(idItem);
				if (lsti == null) { // New item
					lsti = new ArrayList<>();
					// Add item to list
					items.put(idItem, lsti);
				}
				lsti.add(new Tuple3<String, Double, String>(idUser, rating, timestamp));
			});
			stream.close();

			while (true) {
				if (checkStateCore(users, items, numberOfRatingsUser, numberOfRatingsItems, numberOfUniqueRatingsUser, numberOfUniqueRatingsItems)) {
					break;
				}
				updateMaps(users, items, numberOfRatingsUser, numberOfRatingsItems, numberOfUniqueRatingsUser, numberOfUniqueRatingsItems);
			}
			writer = new PrintStream(dstPath);

			// Writing the new data to the file
			for (String user : users.keySet()) {
				List<Tuple3<String, Double, String>> lst = users.get(user);
				for (Tuple3<String, Double, String> t : lst) {
					writer.println(user + "\t" + t.v1 + "\t" + t.v3);
				}
			}
			writer.close();

		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}

	/***
	 * Method to remove the users and items that have not a certain number of
	 * ratings
	 * 
	 * @param users              the preferences of the users
	 * @param items              the preference of the items
	 * @param minimumRatingUsers minimum ratings that the users may have
	 * @param minimumRatingItems minimum ratings that the items may have
	 */
	private static void updateMaps(Map<String, List<Tuple3<String, Double, String>>> users,
			Map<String, List<Tuple3<String, Double, String>>> items, int minimumRatingUsers, int minimumRatingItems, 
			int numberOfUniqueRatingsUser, int numberOfUniqueRatingsItems) {
		Set<String> usersWithLess = new TreeSet<>(
				users.entrySet().stream().filter(t -> t.getValue().size() < minimumRatingUsers).map(t -> t.getKey())
						.collect(Collectors.toSet()));
		Set<String> itemsWithLess = new TreeSet<>(
				items.entrySet().stream().filter(t -> t.getValue().size() < minimumRatingItems).map(t -> t.getKey())
						.collect(Collectors.toSet()));
		
		Set<String> usersWithLessUnique = new TreeSet<>(
				users.entrySet().stream().filter(t -> users.get(t.getKey()).stream().map(t3 -> t3.v1).collect(Collectors.toSet()).size() < numberOfUniqueRatingsUser).map(t -> t.getKey())
						.collect(Collectors.toSet()));
		
		Set<String> itemsWithLessUnique = new TreeSet<>(
				items.entrySet().stream().filter(t -> items.get(t.getKey()).stream().map(t3 -> t3.v1).collect(Collectors.toSet()).size() < numberOfUniqueRatingsItems).map(t -> t.getKey())
						.collect(Collectors.toSet()));
		
		
		usersWithLess.addAll(usersWithLessUnique);
		itemsWithLess.addAll(itemsWithLessUnique);

		// Now for the items we must remove the users that have rated that item
		// AND the items that have less than that number of ratings
		Set<String> totalItems = new TreeSet<>(items.keySet());
		for (String item : totalItems) {
			List<Tuple3<String, Double, String>> lst = items.get(item);
			items.put(item, new ArrayList<>(
					lst.stream().filter(t -> !usersWithLess.contains(t.v1)).collect(Collectors.toList())));
		}

		// Now remove the items from the users
		Set<String> totalUsers = new TreeSet<>(users.keySet());
		for (String user : totalUsers) {
			List<Tuple3<String, Double, String>> lst = users.get(user);
			users.put(user, new ArrayList<>(
					lst.stream().filter(t -> !itemsWithLess.contains(t.v1)).collect(Collectors.toList())));
		}

		// Now remove all items and all users that have less than the ratings
		for (String user : usersWithLess) {
			users.remove(user);
		}

		for (String item : itemsWithLess) {
			items.remove(item);
		}

		//
		totalUsers = new TreeSet<>(users.keySet());
		for (String user : totalUsers) {
			List<Tuple3<String, Double, String>> lst = users.get(user);
			users.put(user,
					new ArrayList<>(lst.stream().filter(t -> (items.get(t.v1) != null)).collect(Collectors.toList())));
		}

		//
		totalItems = new TreeSet<>(items.keySet());
		for (String item : totalItems) {
			List<Tuple3<String, Double, String>> lst = items.get(item);
			items.put(item,
					new ArrayList<>(lst.stream().filter(t -> (users.get(t.v1) != null)).collect(Collectors.toList())));
		}

		System.out.println("Iteration: " + users.size() + " " + items.size());

	}

	/***
	 * Method to check if we are satisfying the k-core property
	 * 
	 * @param users              the preference of the users
	 * @param items              the preference of the items
	 * @param minimumRatingUsers minimum ratings that the users may have
	 * @param minimumRatingItems minimum ratings that the items may have
	 * @return true if the k-core is satisfied, false if not
	 */
	private static boolean checkStateCore(Map<String, List<Tuple3<String, Double, String>>> users,
			Map<String, List<Tuple3<String, Double, String>>> items, int minimumRatingUsers, int minimumRatingItems, int numberOfUniqueRatingsUser, int numberOfUniqueRatingsItem) {
		for (String user : users.keySet()) {
			if (users.get(user).size() < minimumRatingUsers) {
				return false;
			}
			
			Set<String> uniqueRatingsUser = new HashSet<>(users.get(user).stream().map(t -> t.v1).collect(Collectors.toSet()));
			if (uniqueRatingsUser.size() < numberOfUniqueRatingsUser) {
				return false;
			}
		}

		for (String item : items.keySet()) {
			if (items.get(item).size() < minimumRatingItems) {
				return false;
			}
			
			Set<String> uniqueRatingsItem = new HashSet<>(items.get(item).stream().map(t -> t.v1).collect(Collectors.toSet()));
			if (uniqueRatingsItem.size() < numberOfUniqueRatingsItem) {
				return false;
			}
		}
		return true;
	}

	/***
	 * Method to transform an original dataset into another one substituting the
	 * user and item identifiers to long values
	 * 
	 * @param srcPath          the source path
	 * @param dstPath          the destination path
	 * @param characterSplit   the character used to split the data
	 * @param removeDuplicates boolean to remove duplicates
	 */
	public static void DatasetTransformation(String srcPath, String dstPath, String characterSplit,
			boolean removeDuplicates) {
		Map<String, Long> element1s = new TreeMap<String, Long>();
		Map<String, Long> element2s = new TreeMap<String, Long>();

		// Id user, id item rating and timestamp
		// The data to store is a map of users, with a map of items and a list
		// associated with the user and the item
		Map<Long, Map<Long, List<Tuple2<Double, Long>>>> storationUsers = new TreeMap<Long, Map<Long, List<Tuple2<Double, Long>>>>();

		// Variables to customize
		boolean ignoreFirstLine = false; // If we want to ignore first line, switch to true
		int columnFirstElement = 0;
		int columnSecondElement = 1;
		int columnRating = 2;
		int columnTimeStamp = 3;

		BufferedWriter writer = null;
		AtomicLong totalElements1 = new AtomicLong(1L);
		AtomicLong totalElements2 = new AtomicLong(1L);

		Stream<String> stream = null;
		try {
			if (ignoreFirstLine) {
				stream = Files.lines(Paths.get(srcPath)).skip(1); // skipping the first line
			} else {
				stream = Files.lines(Paths.get(srcPath));
			}

			stream.forEach(line -> {
				String[] data = line.split(characterSplit);
				// Transforming the new ids
				if (element1s.get(data[columnFirstElement]) == null) {
					element1s.put(data[columnFirstElement], totalElements1.get());
					totalElements1.incrementAndGet();
				}
				if (element2s.get(data[columnSecondElement]) == null) {
					element2s.put(data[columnSecondElement], totalElements2.get());
					totalElements2.incrementAndGet();
				}

				Long idUserTransformed = element1s.get(data[columnFirstElement]);
				Long idItemTransformed = element2s.get(data[columnSecondElement]);
				double rating = Double.parseDouble((data[columnRating]));
				Long timestamp = Long.parseLong(data[columnTimeStamp]);
				// With the new ids, we will put the the ratings to be printed in the dst file

				if (storationUsers.get(idUserTransformed) == null) { // New user
					Map<Long, List<Tuple2<Double, Long>>> map = new TreeMap<Long, List<Tuple2<Double, Long>>>();
					List<Tuple2<Double, Long>> lst = new ArrayList<Tuple2<Double, Long>>();
					lst.add(new Tuple2<Double, Long>(rating, timestamp));
					map.put(idItemTransformed, lst);
					// Add item to list
					storationUsers.put(idUserTransformed, map);
				} else { // User exits
					Map<Long, List<Tuple2<Double, Long>>> map = storationUsers.get(idUserTransformed);
					List<Tuple2<Double, Long>> lst = map.get(idItemTransformed);
					if (lst == null) { // New item for the user
						lst = new ArrayList<Tuple2<Double, Long>>();
						lst.add(new Tuple2<Double, Long>(rating, timestamp));
						map.put(idItemTransformed, lst);
					} else {
						if (removeDuplicates) {
							if (lst.get(0).v2 < timestamp) { // We obtain the newest timestamp
								lst.clear();
								lst.add(new Tuple2<Double, Long>(rating, timestamp));
							}
						} else {
							lst.add(new Tuple2<Double, Long>(rating, timestamp));
						}
					}
				}
			});
			writer = new BufferedWriter(new FileWriter(dstPath));
			// print the new ratings
			for (Long user : storationUsers.keySet()) {

				for (Long item : storationUsers.get(user).keySet()) {
					for (Tuple2<Double, Long> t : storationUsers.get(user).get(item)) {
						writer.write(user + "\t" + item + "\t" + t.v1 + "\t" + t.v2 + "\n");
					}

				}
			}

			writer.close();
		} catch (Exception e) {
			// TODO: handle exception
		}

	}
	
	public static void DatasetTransformationTime(String srcPath, String dstPath, String characterSplit,
			boolean removeDuplicates) {
		Map<String, Long> element1s = new TreeMap<String, Long>();
		Map<String, Long> element2s = new TreeMap<String, Long>();

		// Id user, id item rating and timestamp
		// The data to store is a map of users, with a map of items and a list
		// associated with the user and the item
		Map<Long, Map<Long, List<Tuple2<Double, Long>>>> storationUsers = new TreeMap<Long, Map<Long, List<Tuple2<Double, Long>>>>();

		// Variables to customize
		boolean ignoreFirstLine = true; // If we want to ignore first line, switch to true
		int columnFirstElement = 0;
		int columnSecondElement = 1;
		int columnRating = 2;
		int columnTimeStamp = 3;

		BufferedWriter writer = null;
		AtomicLong totalElements1 = new AtomicLong(1L);
		AtomicLong totalElements2 = new AtomicLong(1L);
		SimpleDateFormat dateFormat = new SimpleDateFormat("yyyy-MM-dd");


		Stream<String> stream = null;
		try {
			if (ignoreFirstLine) {
				stream = Files.lines(Paths.get(srcPath)).skip(1); // skipping the first line
			} else {
				stream = Files.lines(Paths.get(srcPath));
			}

			stream.forEach(line -> {
				String[] data = line.split(characterSplit);
					// Transforming the new ids
					if (element1s.get(data[columnFirstElement]) == null) {
						element1s.put(data[columnFirstElement], totalElements1.get());
						totalElements1.incrementAndGet();
					}
					if (element2s.get(data[columnSecondElement]) == null) {
						element2s.put(data[columnSecondElement], totalElements2.get());
						totalElements2.incrementAndGet();
					}
	
					Long idUserTransformed = element1s.get(data[columnFirstElement]);
					Long idItemTransformed = element2s.get(data[columnSecondElement]);
					double rating = Double.parseDouble((data[columnRating]));
					String time = data[columnTimeStamp];
					Long timestamp = null;
					try {
						timestamp = dateFormat.parse(time).getTime() / 1000;
					} catch (ParseException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
					// With the new ids, we will put the the ratings to be printed in the dst file
	
					if (storationUsers.get(idUserTransformed) == null) { // New user
						Map<Long, List<Tuple2<Double, Long>>> map = new TreeMap<Long, List<Tuple2<Double, Long>>>();
						List<Tuple2<Double, Long>> lst = new ArrayList<Tuple2<Double, Long>>();
						lst.add(new Tuple2<Double, Long>(rating, timestamp));
						map.put(idItemTransformed, lst);
						// Add item to list
						storationUsers.put(idUserTransformed, map);
					} else { // User exits
						Map<Long, List<Tuple2<Double, Long>>> map = storationUsers.get(idUserTransformed);
						List<Tuple2<Double, Long>> lst = map.get(idItemTransformed);
						if (lst == null) { // New item for the user
							lst = new ArrayList<Tuple2<Double, Long>>();
							lst.add(new Tuple2<Double, Long>(rating, timestamp));
							map.put(idItemTransformed, lst);
						} else {
							if (removeDuplicates) {
								if (lst.get(0).v2 < timestamp) { // We obtain the newest timestamp
									lst.clear();
									lst.add(new Tuple2<Double, Long>(rating, timestamp));
								}
							} else {
								lst.add(new Tuple2<Double, Long>(rating, timestamp));
							}
						}
					} 

			});
			writer = new BufferedWriter(new FileWriter(dstPath));
			// print the new ratings
			for (Long user : storationUsers.keySet()) {

				for (Long item : storationUsers.get(user).keySet()) {
					for (Tuple2<Double, Long> t : storationUsers.get(user).get(item)) {
						writer.write(user + "\t" + item + "\t" + t.v1 + "\t" + t.v2 + "\n");
					}

				}
			}

			writer.close();
		} catch (Exception e ) {
			System.out.println(e);
		}

	}
	
	/***
	 * Method to transform an original dataset into another one substituting the
	 * user and item identifiers to long values
	 * 
	 * @param srcPath          the source path
	 * @param dstPath          the destination path
	 * @param characterSplit   the character used to split the data
	 * @param removeDuplicates boolean to remove duplicates
	 */
	public static void DatasetTransformationRetailRocket(String srcPath, String dstPath, String characterSplit,
			boolean removeDuplicates) {
		Map<String, Long> element1s = new TreeMap<String, Long>();
		Map<String, Long> element2s = new TreeMap<String, Long>();

		// Id user, id item rating and timestamp
		// The data to store is a map of users, with a map of items and a list
		// associated with the user and the item
		Map<Long, Map<Long, List<Tuple2<Double, Long>>>> storationUsers = new TreeMap<Long, Map<Long, List<Tuple2<Double, Long>>>>();

		// Variables to customize
		boolean ignoreFirstLine = true; // If we want to ignore first line, switch to true
		int columnFirstElement = 1;
		int columnSecondElement = 3;
		int columnTimeStamp = 0; //Timestamp in this case is in milliseconds

		BufferedWriter writer = null;
		AtomicLong totalElements1 = new AtomicLong(1L);
		AtomicLong totalElements2 = new AtomicLong(1L);

		Stream<String> stream = null;
		try {
			if (ignoreFirstLine) {
				stream = Files.lines(Paths.get(srcPath)).skip(1); // skipping the first line
			} else {
				stream = Files.lines(Paths.get(srcPath));
			}

			stream.forEach(line -> {
				String[] data = line.split(",");
				// Transforming the new ids
				if (element1s.get(data[columnFirstElement]) == null) {
					element1s.put(data[columnFirstElement], totalElements1.get());
					totalElements1.incrementAndGet();
				}
				if (element2s.get(data[columnSecondElement]) == null) {
					element2s.put(data[columnSecondElement], totalElements2.get());
					totalElements2.incrementAndGet();
				}

				Long idUserTransformed = element1s.get(data[columnFirstElement]);
				Long idItemTransformed = element2s.get(data[columnSecondElement]);
				double rating = 1.0;
				Double timestampD = Math.floor(Double.parseDouble(data[columnTimeStamp]) / 1000);
				Long timestamp = timestampD.longValue();
				// With the new ids, we will put the the ratings to be printed in the dst file

				if (storationUsers.get(idUserTransformed) == null) { // New user
					Map<Long, List<Tuple2<Double, Long>>> map = new TreeMap<Long, List<Tuple2<Double, Long>>>();
					List<Tuple2<Double, Long>> lst = new ArrayList<Tuple2<Double, Long>>();
					lst.add(new Tuple2<Double, Long>(rating, timestamp));
					map.put(idItemTransformed, lst);
					// Add item to list
					storationUsers.put(idUserTransformed, map);
				} else { // User exits
					Map<Long, List<Tuple2<Double, Long>>> map = storationUsers.get(idUserTransformed);
					List<Tuple2<Double, Long>> lst = map.get(idItemTransformed);
					if (lst == null) { // New item for the user
						lst = new ArrayList<Tuple2<Double, Long>>();
						lst.add(new Tuple2<Double, Long>(rating, timestamp));
						map.put(idItemTransformed, lst);
					} else {
						if (removeDuplicates) {
							if (lst.get(0).v2 < timestamp) { // We obtain the newest timestamp
								lst.clear();
								lst.add(new Tuple2<Double, Long>(rating, timestamp));
							}
						} else {
							lst.add(new Tuple2<Double, Long>(rating, timestamp));
						}
					}
				}
			});
			writer = new BufferedWriter(new FileWriter(dstPath));
			// print the new ratings
			for (Long user : storationUsers.keySet()) {

				for (Long item : storationUsers.get(user).keySet()) {
					for (Tuple2<Double, Long> t : storationUsers.get(user).get(item)) {
						writer.write(user + "\t" + item + "\t" + t.v1 + "\t" + t.v2 + "\n");
					}

				}
			}

			writer.close();
		} catch (Exception e) {
			e.printStackTrace();
		}

	}
	/***
	 * Method to transform an original dataset into another one substituting the
	 * user and item identifiers to long values. 
	 * This transformation is for lastfm 1kUsers
	 * 
	 * @param srcPath          the source path
	 * @param dstPath          the destination path
	 * @param characterSplit   the character used to split the data
	 * @param removeDuplicates boolean to remove duplicates
	 */
	public static void DatasetTransformLastfm1k(String srcPath, String dstPath) {
		Map<String, Long> element1s = new TreeMap<String, Long>();
		Map<String, Long> element2s = new TreeMap<String, Long>();

		// Id user, id item rating and timestamp
		// The data to store is a map of users, with a map of items and a list
		// associated with the user and the item
		Map<Long, Map<Long, List<Tuple2<Double, Long>>>> storationUsers = new TreeMap<Long, Map<Long, List<Tuple2<Double, Long>>>>();

		// Variables to customize
		boolean ignoreFirstLine = false; // If we want to ignore first line, switch to true
		int columnFirstElement = 0;
		int columnSecondElement = 2;
		int columnSecondElement2 = 3;
		int columnTimeStamp = 1;

		BufferedWriter writer = null;
		AtomicLong totalElements1 = new AtomicLong(1L);
		AtomicLong totalElements2 = new AtomicLong(1L);

		Stream<String> stream = null;
		try {
			if (ignoreFirstLine) {
				stream = Files.lines(Paths.get(srcPath)).skip(1); // skipping the first line
			} else {
				stream = Files.lines(Paths.get(srcPath));
			}
			writer = new BufferedWriter(new FileWriter(dstPath));
			BufferedWriter writer2 = writer;

			stream.forEach(line -> {
				String[] data = line.split("\t");
				// Transforming the new ids
				if (element1s.get(data[columnFirstElement]) == null) {
					element1s.put(data[columnFirstElement], totalElements1.get());
					totalElements1.incrementAndGet();
				}
				String prevItemIDItemName = data[columnSecondElement] + data[columnSecondElement2];
				if (element2s.get(prevItemIDItemName) == null) {
					element2s.put(prevItemIDItemName, totalElements2.get());
					totalElements2.incrementAndGet();
				}

				Long idUserTransformed = element1s.get(data[columnFirstElement]);
				Long idItemTransformed = element2s.get(prevItemIDItemName);
				String dateRead = data[columnTimeStamp];
				SimpleDateFormat simpleDateFormat = new SimpleDateFormat("yyyy-MM-dd'T'HH:mm:ss'Z'");
				//simpleDateFormat.setTimeZone(TimeZone.getTimeZone("GMT"));
				//simpleDateFormat.setTimeZone(TimeZone.getTimeZone("UTC")); //UTC is the same as GMT
				Date date = null;
				try {
					date = simpleDateFormat.parse(dateRead);
				} catch (ParseException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
				
				long timestamp = date.getTime() / 1000;
				// With the new ids, we will put the the ratings to be printed in the dst file

				try {
					writer2.write(idUserTransformed + "\t" + idItemTransformed + "\t" + 1.0 + "\t" + timestamp + "\n");
				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			});


			writer.close();
		} catch (Exception e) {
			e.printStackTrace();

		}

	}
	
	public static void analysisOfTest(String testFile, String trainFile, String userFeatureFile, String dstPath) {
		try {
			UserFeatureData<Long, String, Double> ufD = SimpleUserFeatureData.load(SimpleUserFeaturesReader.get().read(userFeatureFile, lp, sp));
			Map<Long, Map<String, Integer>> itemsRatedByTypeOfUsersTest = null;
			Map<Long, Set<Long>> itemAndUsersRatedTest = null;
			Map<Long, Integer> itemRatingsTest = null;

		
			Tuple3<Map<Long, Map<String, Integer>>, Map<Long, Set<Long>>, Map<Long, Integer>> t = obtainMapsUserFeatureFromFile(ufD, testFile);
			itemsRatedByTypeOfUsersTest = t.v1;
			itemAndUsersRatedTest = t.v2;
			itemRatingsTest = t.v3;

			Map<Long, Map<String, Integer>> itemsRatedByTypeOfUsersTrain = null;
			Map<Long, Set<Long>> itemAndUsersRatedTrain = null;
			Map<Long, Integer> itemRatingsTrain = null;

			
			Tuple3<Map<Long, Map<String, Integer>>, Map<Long, Set<Long>>, Map<Long, Integer>> t2 = obtainMapsUserFeatureFromFile(ufD, trainFile);
			itemsRatedByTypeOfUsersTrain = t2.v1;
			itemAndUsersRatedTrain = t2.v2;
			itemRatingsTrain = t2.v3;
			
			
			/*
			Map<Long, Map<String, Integer>> itemsRatedByTypeOfUsersRec = null;
			Map<Long, Set<Long>> itemAndUsersRatedRec = null;
			
			
			Tuple2<Map<Long, Map<String, Integer>>, Map<Long, Set<Long>>> t3 = obtainMapsUserFeatureFromFile(ufD, recFile);
			itemsRatedByTypeOfUsersRec = t3.v1;
			itemAndUsersRatedRec = t3.v2;
			*/
			
			BufferedWriter writer = null;
			BufferedWriter writer2 = null;
			try {
				writer = new BufferedWriter(new FileWriter(dstPath));
				writer2 = writer;
			
				//Write header
				writer2.write("TestItem" + "\t" + "TotalUsersRatedInTrain" + "\t");
				List<String> userFeatures = ufD.getAllFeatures().collect(Collectors.toList());
				for (String userFeat: userFeatures) {
					writer2.write("N-UsersTestRated-" + userFeat + "\t");
					writer2.write("N-UsersTrainRated-" + userFeat + "\t");
					//writer2.write("N-UsersRecommended-" + userFeat + "\t");
				}
				writer2.write("\n");
				//Now print the stats for every item in the TEST set
				for (Long idItem: itemAndUsersRatedTest.keySet()) {
					writer2.write(idItem + "\t");
					if (itemAndUsersRatedTrain.get(idItem) != null) {
						writer2.write(itemAndUsersRatedTrain.get(idItem).size() + "\t");
					} else {
						writer2.write(0 + "\t");
					}
					
					
					for (String userFeat: userFeatures) {
						writer2.write(itemsRatedByTypeOfUsersTest.get(idItem).get(userFeat) + "\t");
						
						if (itemsRatedByTypeOfUsersTrain.get(idItem) == null || itemsRatedByTypeOfUsersTrain.get(idItem).get(userFeat) == null) {
							writer2.write(0 + "\t");
						}
						else {
							writer2.write(itemsRatedByTypeOfUsersTrain.get(idItem).get(userFeat) + "\t");
						}
						/*
						if (itemsRatedByTypeOfUsersRec.get(idItem) == null || itemsRatedByTypeOfUsersRec.get(idItem).get(userFeat) == null) {
							writer2.write(0 + "\t");
						}
						else {
							writer2.write(itemsRatedByTypeOfUsersRec.get(idItem).get(userFeat) + "\t");
						}
						*/
					}
					writer2.write("\n");
				}
				
				double avgCheckinsTest = 0.0;
				double avgCheckinsTrain = 0.0;
				for (Long idItem : itemRatingsTest.keySet()) {
					avgCheckinsTest+=itemRatingsTest.get(idItem);
				}
				for (Long idItem : itemRatingsTrain.keySet()) {
					avgCheckinsTrain+=itemRatingsTrain.get(idItem);
				}
				System.out.println("Total number of items in train: " + itemRatingsTrain.size());
				System.out.println("NºItems/Checkins in train: " + avgCheckinsTrain/itemRatingsTrain.size());
				System.out.println("Total number of items in test: " + itemRatingsTest.size());
				System.out.println("NºItems/Checkins in test: " + avgCheckinsTest/itemRatingsTest.size());

				
				
				//Now, I will find the percentage of items rated by each user feature
				for (String userFeat: userFeatures) {
					int callUserFeatTrain = 0;
					int conlyUserFeatTrain = 0;

					
					int callUserFeatTest = 0;
					double avgRatingsItemFeat = 0.0;
					double avgRatingsItemOnlyFeat = 0.0;
					

					
					//For each item in Train, we see if we have that feature with a value higher than 0
					for (Long idItem: itemsRatedByTypeOfUsersTrain.keySet()) {
						if (itemsRatedByTypeOfUsersTrain.get(idItem).get(userFeat) > 0) {
							callUserFeatTrain++;
							
							//Now, we filter out the rest of features (features that are not userFeat)
							List<String> userFeaturesNotThis = userFeatures.stream().filter(f -> !f.equals(userFeat)).collect(Collectors.toList());
							boolean found = false;
							for (String userFeat2: userFeaturesNotThis) {
								if (itemsRatedByTypeOfUsersTrain.get(idItem).get(userFeat2) > 0){
									found = true;
									break;
								}
							}
							if (!found) {
								conlyUserFeatTrain++;
							}
						}
						
					}
					System.out.println("Total number of items in train rated by " + userFeat + " " + callUserFeatTrain);
					System.out.println("Total number of items in train rated by ONLY " + userFeat + " " + conlyUserFeatTrain);

					//For each item from test, we check if it has been rated in train by any user of this feature
					for (Long idItem: itemsRatedByTypeOfUsersTest.keySet()) {

						if (itemsRatedByTypeOfUsersTrain.get(idItem) != null && itemsRatedByTypeOfUsersTrain.get(idItem).get(userFeat) > 0) {
							callUserFeatTest++;
							avgRatingsItemFeat += itemAndUsersRatedTrain.get(idItem).size();
							
							//Filer out the users that do not have that feature (we only want to obtain the interactions of the user of that feature)
							Set<Long> onlyCheckinsUF = itemAndUsersRatedTrain.get(idItem).stream().
									filter(idUser -> ufD.getUserFeatures(idUser).
											map(tuple-> tuple.v1).anyMatch(feature -> feature.equals(userFeat))).
									collect(Collectors.toSet());
							
							avgRatingsItemOnlyFeat += onlyCheckinsUF.size();
							
						}
					}
					
					System.out.println("Total number of test items appearing in train rated by: " + userFeat + " " + callUserFeatTest);
					System.out.println("Avg ratings of test items appearing in train rated by: " + userFeat + " : " + avgRatingsItemFeat / callUserFeatTest);
					System.out.println("Avg ratings of test items appearing in train (of only users with the same feature) rated by: "  + userFeat + ": " + avgRatingsItemOnlyFeat / callUserFeatTest);


					
					
				}
				

			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			
			
			writer.close();
		
		} catch (IOException e) {
			e.printStackTrace();
		}		
	}
	
	private static Tuple3< Map<Long, Map<String, Integer>>, Map<Long, Set<Long>>, Map<Long, Integer>> obtainMapsUserFeatureFromFile(UserFeatureData<Long, String, Double> ufD, String file){
		int indexUser = 0;
		int indexItem = 1;
		
		//For each item, we store the number of users of each feature that have rated it in test
		
		//Id Item -> String feature -> number of users having that feature that rated the item
		Map<Long, Map<String, Integer>> itemsRatedByTypeOfUsers = new HashMap<>();	
		
		//Id Item -> Set of users that have rated
		Map<Long, Set<Long>> itemAndUsersRated = new HashMap<>();
		
		//Id Item -> Number of interactions associated with that item
		Map<Long, Integer> itemRatings = new HashMap<>();
		
		Stream<String> stream;
		try {
			stream = Files.lines(Paths.get(file));
			stream.forEach(line -> {
				String data [] = line.split("\t");
				Long idUser = Long.parseLong(data[indexUser]);
				Long idItem = Long.parseLong(data[indexItem]);
				
				//For each item we store the maps of users 
				if (!itemsRatedByTypeOfUsers.containsKey(idItem)) {
					itemsRatedByTypeOfUsers.put(idItem, new HashMap<>());
					itemRatings.put(idItem, 0);
					ufD.getAllFeatures().forEach(uf -> {
						itemsRatedByTypeOfUsers.get(idItem).put(uf, 0);
					});
				}
				itemRatings.put(idItem, itemRatings.get(idItem) + 1);
				
				if (!itemAndUsersRated.containsKey(idItem)) {
					itemAndUsersRated.put(idItem, new HashSet<Long>());
				}
				
				//Now, store the info of the test set
				if (!itemAndUsersRated.get(idItem).contains(idUser)) {
					itemAndUsersRated.get(idItem).add(idUser);
					ufD.getUserFeatures(idUser).forEach(uf -> {
						itemsRatedByTypeOfUsers.get(idItem).put(uf.v1, itemsRatedByTypeOfUsers.get(idItem).get(uf.v1) + 1);
					});
					
				}				
			});
			stream.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return new Tuple3<>(itemsRatedByTypeOfUsers, itemAndUsersRated, itemRatings);
		
	}
	
	

	/***
	 * Method to compute the statistics of a dataset
	 * 
	 * @param srcPath the source dataset
	 * @param srcDest the destination dataset
	 */
	public static void Stats(String srcPath, String srcDest) {
		boolean implicit = true;
		// Parameters to configure
		String characterSplit = "\t"; // Substitute this character to the others to split
		int indexUser = 0;
		int indexItem = 1;
		int indexRating = 2;
		//

		// Users. ID of element 1 and number of items rated
		Map<String, Integer> element1s = new TreeMap<String, Integer>();
		// Items. ID of element 2 and number of items rated
		Map<String, Integer> element2s = new TreeMap<String, Integer>();
		Map<String, Map<String, List<Double>>> preferences = new TreeMap<String, Map<String, List<Double>>>();

		AtomicLong totalRatings = new AtomicLong(0);
		AtomicLong totalUsers = new AtomicLong(0);
		AtomicLong totalItems = new AtomicLong(0);
		AtomicDouble avgRatings = new AtomicDouble(0);
		double density = 0;
		try (Stream<String> stream = Files.lines(Paths.get(srcPath))) {

			stream.forEach(line -> {

				String[] data = line.split(characterSplit);
				if (element1s.get(data[indexUser]) == null) {
					element1s.put(data[indexUser], 1);
					totalUsers.incrementAndGet();
				} else
					element1s.put(data[indexUser], element1s.get(data[indexUser]) + 1);

				if (element2s.get(data[indexItem]) == null) {
					element2s.put(data[indexItem], 1);
					totalItems.incrementAndGet();
				} else
					element2s.put(data[indexItem], element2s.get(data[indexItem]) + 1);

				if (preferences.get(data[indexUser]) == null) { // New user
					Map<String, List<Double>> prefUser = new TreeMap<String, List<Double>>();
					List<Double> lstPrefRepeated = new ArrayList<Double>();
					if (!implicit)
						lstPrefRepeated.add(Double.parseDouble(data[indexRating]));
					else
						lstPrefRepeated.add(1.0);
					prefUser.put(data[indexItem], lstPrefRepeated);
					preferences.put(data[indexUser], prefUser);
				} else {
					if (preferences.get(data[indexUser]).get(data[indexItem]) == null) { // New item for user
						List<Double> lstPrefRepeated = new ArrayList<Double>();
						if (!implicit)
							lstPrefRepeated.add(Double.parseDouble(data[indexRating]));
						else
							lstPrefRepeated.add(1.0);
						preferences.get(data[indexUser]).put(data[indexItem], lstPrefRepeated);
					} else {
						List<Double> lstPrefRepeated = preferences.get(data[indexUser]).get(data[indexItem]);
						if (!implicit)
							lstPrefRepeated.add(Double.parseDouble(data[indexRating]));
						else
							lstPrefRepeated.add(1.0);
					}
				}
				if (!implicit)
					avgRatings.addAndGet(Double.parseDouble(data[indexRating]));
				else
					avgRatings.addAndGet(1.0);

				totalRatings.incrementAndGet();

				if (totalRatings.get() % 10000 == 0)
					System.out.println("10000 processed");
			});
			density = (double) totalRatings.get() / (double) (totalUsers.get() * totalItems.get());
			BufferedWriter writer = null;
			writer = new BufferedWriter(new FileWriter(srcDest));

			writer.write("Analyzing \n");
			writer.write("Average of ratings " + avgRatings.get() / totalRatings.get() + "\n");
			double onlyones = 0;
			double onlytwos = 0;
			double onlythrees = 0;
			double onlyfours = 0;
			double onlyfives = 0;
			double morethanFive = 0;

			for (String user : element1s.keySet()) {
				Integer setUser = element1s.get(user);
				if (setUser == 1)
					onlyones++;
				else if (setUser == 2)
					onlytwos++;
				else if (setUser == 3)
					onlythrees++;
				else if (setUser == 4)
					onlyfours++;
				else if (setUser == 5)
					onlyfives++;
				if (setUser > 5)
					morethanFive++;
			}

			double onlyonesI = 0;
			double onlytwosI = 0;
			double onlythreesI = 0;
			double onlyfoursI = 0;
			double onlyfivesI = 0;
			double morethanFiveI = 0;
			for (String item : element2s.keySet()) {
				Integer setItem = element2s.get(item);
				if (setItem == 1)
					onlyonesI++;
				else if (setItem == 2)
					onlytwosI++;
				else if (setItem == 3)
					onlythreesI++;
				else if (setItem == 4)
					onlyfoursI++;
				else if (setItem == 5)
					onlyfivesI++;
				if (setItem > 5)
					morethanFiveI++;
			}

			writer.write("Total users: " + totalUsers + "\n");
			writer.write("Total items:  " + totalItems + "\n");
			writer.write("Total ratings: " + totalRatings + "\n");
			writer.write(
					"Average ratings per user: " + ((double) totalRatings.get() / (double) totalUsers.get()) + "\n");
			writer.write(
					"Average ratings per item: " + ((double) totalRatings.get() / (double) totalItems.get()) + "\n");
			writer.write("Users with 1 rating: " + onlyones + "\n");
			writer.write("Users with 2 rating: " + onlytwos + "\n");
			writer.write("Users with 3 rating: " + onlythrees + "\n");
			writer.write("Users with 4 rating: " + onlyfours + "\n");
			writer.write("Users with 5 rating: " + onlyfives + "\n");
			writer.write("Users with more than 5 rating: " + morethanFive + "\n");

			writer.write("Items with 1 rating: " + onlyonesI + "\n");
			writer.write("Items with 2 rating: " + onlytwosI + "\n");
			writer.write("Items with 3 rating: " + onlythreesI + "\n");
			writer.write("Items with 4 rating: " + onlyfoursI + "\n");
			writer.write("Items with 5 rating: " + onlyfivesI + "\n");
			writer.write("Items with more than 5 rating: " + morethanFiveI + "\n");

			writer.write("Repeated ratings (One user rating the same item more than one time) +\n");
			int counter = 0;
			for (String user : preferences.keySet()) {
				for (String item : preferences.get(user).keySet()) {
					if (preferences.get(user).get(item).size() > 1)
						counter += preferences.get(user).get(item).size() - 1;
				}
			}
			writer.write("Total repetitions: " + counter + "\n");
			writer.write("Density :" + density);
			writer.close();

		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	public static void StatsPerMonth(String srcPath, String srcDest) {
		boolean implicit = false;
		// Parameters to configure
		String characterSplit = "\t"; // Substitute this character to the others to split
		int indexUser = 0;
		int indexItem = 1;
		int indexRating = 2;
		int indexTimeStamp = 3;
		
		//
		Map<Long, Set<Integer>> usersInTimeStamp = new TreeMap<>();
		Map<Long, Set<Integer>> itemInTimeStamp = new TreeMap<>();
		Map<Long, Integer>  ratingsTimeStamp = new TreeMap<>();

		AtomicLong minTimeStamp = new AtomicLong(Long.MAX_VALUE);
		AtomicLong maxTimeStamp = new AtomicLong(Long.MIN_VALUE);


		//First pass, compute min and max timestamp
		try (Stream<String> stream = Files.lines(Paths.get(srcPath))) {			
			stream.forEach(line -> {
				long readTimestamp = Long.parseLong(line.split(characterSplit)[indexTimeStamp]);
				if (minTimeStamp.get() > readTimestamp) {
					minTimeStamp.set(readTimestamp);
				}
				
				if (maxTimeStamp.get() < readTimestamp) {
					maxTimeStamp.set(readTimestamp);
				}
				
			});
			stream.close();
			//We have the maximum and the minimum timestamp
			long starting = minTimeStamp.get() * 1000;
			long starting2 = starting;

			Calendar calendar = Calendar.getInstance();
			Date minimum = new Date(starting); //In java dates are in ms
			calendar.setTime(minimum);
			calendar.setTimeZone(TimeZone.getTimeZone("GMT"));
			
			//First second o start new day
			calendar.set(Calendar.DAY_OF_MONTH, 1);
			calendar.set(Calendar.HOUR_OF_DAY, 0);
			calendar.set(Calendar.MINUTE, 0);
			calendar.set(Calendar.SECOND, 1);


			starting = calendar.getTimeInMillis();
			
			long ending = maxTimeStamp.get() * 1000;
			
			
			while (starting < ending) {
				usersInTimeStamp.put(starting, new HashSet<>());
				itemInTimeStamp.put(starting, new HashSet<>());
				ratingsTimeStamp.put(starting, 0);
				
				Date actual = new Date(starting);
				actual = DateUtils.addMonths(actual, 1);
				starting = actual.getTime();
			}
			
			//Now, second pass and add the complete information of the specified maps.
			
			Stream<String> stream2 = Files.lines(Paths.get(srcPath));
			stream2.forEach(line -> {
				String [] fields = line.split(characterSplit);
				int user = Integer.parseInt(fields[indexUser]);
				int item = Integer.parseInt(fields[indexItem]);
				Long time = Long.parseLong(fields[indexTimeStamp]) * 1000;
				
				long prev = starting2;
				for (Long possibleTime: usersInTimeStamp.keySet()) {
					if (time < possibleTime) {
						break;
					}
					prev = possibleTime;
				}
				usersInTimeStamp.get(prev).add(user);
				itemInTimeStamp.get(prev).add(item);
				ratingsTimeStamp.put(prev, ratingsTimeStamp.get(prev) + 1 );
			});
			stream2.close();
			
			//Now, we print
			
			BufferedWriter writer = null;
			writer = new BufferedWriter(new FileWriter(srcDest));
			
			List<Long> orderedTimeStamps = new ArrayList<>(usersInTimeStamp.keySet());
			for (int i = 0; i < orderedTimeStamps.size() - 1; i++) {
				long possibleTime = orderedTimeStamps.get(i);
				long possibleTime2 = orderedTimeStamps.get(i+1);

				writer.write("Stats between: \n");
				writer.write("TimeStamp: " + possibleTime + " with date " + new Date(possibleTime) + "\n");
				writer.write("TimeStamp: " + possibleTime2 + " with date " + new Date(possibleTime2) + "\n");

				writer.write("Number of users: " + usersInTimeStamp.get(possibleTime).size() + "\n");
				writer.write("Number of items: " + itemInTimeStamp.get(possibleTime).size() + "\n");
				writer.write("Number of ratings: " + ratingsTimeStamp.get(possibleTime) + "\n");

				
				
			}
			long last = orderedTimeStamps.get(orderedTimeStamps.size() - 1);
			writer.write("Stats between: \n");
			writer.write("TimeStamp: " + last + " with date " + new Date(last) + "\n");
			writer.write("TimeStamp: " + ending + " with date " + new Date(ending) + "\n");

			writer.write("Number of users: " + usersInTimeStamp.get(last).size() + "\n");
			writer.write("Number of items: " + itemInTimeStamp.get(last).size() + "\n");
			writer.write("Number of ratings: " + ratingsTimeStamp.get(last) + "\n");
			
			
			writer.close();

			
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	

			
	}

	public static void Stats(String srcPathTrain, String srcPathTest, String srcDest) {
		boolean implicit = true;
		// Parameters to configure
		String characterSplit = "\t"; // Substitute this character to the others to split
		int indexUser = 0;
		int indexItem = 1;
		int indexRating = 2;
		//

		// Users. ID of element 1 and number of items rated
		Map<String, Integer> element1s = new TreeMap<String, Integer>();
		// Items. ID of element 2 and number of items rated
		Map<String, Integer> element2s = new TreeMap<String, Integer>();
		Map<String, Map<String, List<Double>>> preferences = new TreeMap<String, Map<String, List<Double>>>();

		AtomicLong totalRatings = new AtomicLong(0);
		AtomicLong totalUsers = new AtomicLong(0);
		AtomicLong totalItems = new AtomicLong(0);
		AtomicDouble avgRatings = new AtomicDouble(0);

		try {
			Stream<String> stream = Files.lines(Paths.get(srcPathTrain));
			Stream<String> stream2 = Files.lines(Paths.get(srcPathTest));
			// For every line
			Stream.concat(stream, stream2).forEach(line -> {
				if (line != null) {
					String[] data = line.split(characterSplit);
					if (element1s.get(data[indexUser]) == null) {
						element1s.put(data[indexUser], 1);
						totalUsers.incrementAndGet();
					} else {
						element1s.put(data[indexUser], element1s.get(data[indexUser]) + 1);
					}
					if (element2s.get(data[indexItem]) == null) {
						element2s.put(data[indexItem], 1);
						totalItems.incrementAndGet();
					} else {
						element2s.put(data[indexItem], element2s.get(data[indexItem]) + 1);
					}
					if (preferences.get(data[indexUser]) == null) { // New user
						Map<String, List<Double>> prefUser = new TreeMap<String, List<Double>>();
						List<Double> lstPrefRepeated = new ArrayList<Double>();
						if (!implicit) {
							lstPrefRepeated.add(Double.parseDouble(data[indexRating]));
						}
						else {
							lstPrefRepeated.add(1.0);
						}
						prefUser.put(data[indexItem], lstPrefRepeated);
						preferences.put(data[indexUser], prefUser);
					} else {
						if (preferences.get(data[indexUser]).get(data[indexItem]) == null) { // New item for user
							List<Double> lstPrefRepeated = new ArrayList<Double>();
							if (!implicit) {
								lstPrefRepeated.add(Double.parseDouble(data[indexRating]));
							}
							else {
								lstPrefRepeated.add(1.0);
							}
							preferences.get(data[indexUser]).put(data[indexItem], lstPrefRepeated);
						} else {
							List<Double> lstPrefRepeated = preferences.get(data[indexUser]).get(data[indexItem]);
							if (!implicit) {
								lstPrefRepeated.add(Double.parseDouble(data[indexRating]));
							}
							else {
								lstPrefRepeated.add(1.0);
							}
						}
					}
					avgRatings.set(avgRatings.get() + Double.parseDouble(data[indexRating]));

				}
				totalRatings.incrementAndGet();

				if (totalRatings.get() % 10000 == 0)
					System.out.println("10000 processed");

			});
			BufferedWriter writer = null;
			writer = new BufferedWriter(new FileWriter(srcDest));

			writer.write("Analyzing \n");
			writer.write("Average of ratings " + avgRatings.get() / totalRatings.get() + "\n");
			double onlyones = 0;
			double onlytwos = 0;
			double onlythrees = 0;
			double onlyfours = 0;
			double onlyfives = 0;
			double moreThan15U = 0;
			double morethanFive = 0;

			for (String user : element1s.keySet()) {
				Integer setUser = element1s.get(user);
				if (setUser == 1) {
					onlyones++;
				}
				else if (setUser == 2) {
					onlytwos++;
				}
				else if (setUser == 3) {
					onlythrees++;
				}
				else if (setUser == 4) {
					onlyfours++;
				}
				else if (setUser == 5) {
					onlyfives++;
				}
				else if (setUser >= 15) {
					moreThan15U++;
				}
				if (setUser > 5) {
					morethanFive++;
				}
			}

			double onlyonesI = 0;
			double onlytwosI = 0;
			double onlythreesI = 0;
			double onlyfoursI = 0;
			double onlyfivesI = 0;
			double morethatTenI = 0;
			double morethanFiveI = 0;
			for (String item : element2s.keySet()) {
				Integer setItem = element2s.get(item);
				if (setItem == 1) {
					onlyonesI++;
				}
				else if (setItem == 2) {
					onlytwosI++;
				}
				else if (setItem == 3) {
					onlythreesI++;
				}
				else if (setItem == 4) {
					onlyfoursI++;
				}
				else if (setItem == 5) {
					onlyfivesI++;
				}
				else if (setItem >= 10) {
					morethatTenI++;
				}
				if (setItem > 5) {
					morethanFiveI++;
				}
			}

			writer.write("Total users: " + totalUsers + "\n");
			writer.write("Total items:  " + totalItems + "\n");
			writer.write("Total ratings: " + totalRatings + "\n");
			writer.write("Average ratings per user: " + (totalRatings.get() / totalUsers.get()) + "\n");
			writer.write("Average ratings per item: " + (totalRatings.get() / totalItems.get()) + "\n");
			writer.write("Users with 1 rating: " + onlyones + "\n");
			writer.write("Users with 2 rating: " + onlytwos + "\n");
			writer.write("Users with 3 rating: " + onlythrees + "\n");
			writer.write("Users with 4 rating: " + onlyfours + "\n");
			writer.write("Users with 5 rating: " + onlyfives + "\n");
			writer.write("Users with 5 rating: " + onlyfives + "\n");
			writer.write("Users with more than 5 rating: " + morethanFive + "\n");
			writer.write("Users with more than 15 rating: " + moreThan15U + "\n");

			writer.write("Items with 1 rating: " + onlyonesI + "\n");
			writer.write("Items with 2 rating: " + onlytwosI + "\n");
			writer.write("Items with 3 rating: " + onlythreesI + "\n");
			writer.write("Items with 4 rating: " + onlyfoursI + "\n");
			writer.write("Items with 5 rating: " + onlyfivesI + "\n");
			writer.write("Items with more than 5 rating: " + morethanFiveI + "\n");
			writer.write("Items with more than 10 rating: " + morethatTenI + "\n");

			writer.write("Repeated ratings (One user rating the same item more than one time) +\n");
			int counter = 0;
			for (String user : preferences.keySet()) {
				for (String item : preferences.get(user).keySet()) {
					if (preferences.get(user).get(item).size() > 1) {
						counter += preferences.get(user).get(item).size() - 1;
					}
				}
			}
			writer.write("Total repetitions: " + counter + "\n");

			writer.close();

		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	/**
	 * It assumes that the file has at least: userID, itemID, preference, timestamp,
	 * session
	 * 
	 * @param srcPathTrain
	 * @param resultFile
	 */
	public static void StatsSession(String srcPathTrain, String pathUserLocals, String resultFile, boolean ignoreFstLine) {

		// Map -> Id of user -> see if it is local or tourist
		Map<String, String> usersLocalOrTourist = new HashMap<>();

		// Map -> Id user -> set of sessions
		Map<String, Set<String>> totalUsersNumberSessions = new HashMap<>();

		// Map -> Id user -> Id session -> numberItems in that session
		Map<String, Map<String, Integer>> sessionsNumberItems = new HashMap<>();

		//Map -> Id user -> Id session -> List of timestamps
		Map<String, Map<String, List<Long>>> sessionsTimeStamps = new HashMap<>();
		
		Set<String> items = new HashSet<>();
		
		AtomicInteger totalPreferences = new AtomicInteger(0);
		try {
			PrintStream output = new PrintStream(resultFile);

			if (pathUserLocals != null) {

				Stream<String> stream = Files.lines(Paths.get(pathUserLocals));
				if (ignoreFstLine) {
					stream.skip(1).forEach(line -> {
						String data[] = line.split("\t");
						String userId = data[0];
						String touristOrNot = data[1];
						usersLocalOrTourist.put(userId, touristOrNot);
					});
					stream.close();
				} else {
					stream.forEach(line -> {
						String data[] = line.split("\t");
						String userId = data[0];
						String touristOrNot = data[1];
						usersLocalOrTourist.put(userId, touristOrNot);
					});
					stream.close();
				}
				 
			}

			Stream<String> stream = Files.lines(Paths.get(srcPathTrain));
			stream.forEach(line -> {
				String data[] = line.split("\t");
				String userId = data[0];
				String itemId = data[1];
				Double score = Double.parseDouble(data[2]);
				Long time = Long.parseLong(data[3]);
				String session = data[4];
				
				items.add(itemId);
				totalPreferences.incrementAndGet();

				// Map of number of sessions
				if (totalUsersNumberSessions.get(userId) == null) {
					totalUsersNumberSessions.put(userId, new HashSet<>());
				}
				totalUsersNumberSessions.get(userId).add(session);
				

				// Map of users, sessions and items per session. First, see if the user exists
				if (sessionsNumberItems.get(userId) == null) {
					sessionsNumberItems.put(userId, new HashMap<>());
					sessionsTimeStamps.put(userId, new HashMap<>());
				}

				// See if the session exist
				if (sessionsNumberItems.get(userId).get(session) == null) {
					sessionsNumberItems.get(userId).put(session, 1);
					sessionsTimeStamps.get(userId).put(session, new ArrayList<>());
					sessionsTimeStamps.get(userId).get(session).add(time);
				} else {
					sessionsNumberItems.get(userId).put(session, sessionsNumberItems.get(userId).get(session) + 1);
					sessionsTimeStamps.get(userId).get(session).add(time);
				}

			});
			stream.close();

			// Some statistics
			int totalUsers = totalUsersNumberSessions.keySet().size();
			int totalSessions = totalUsersNumberSessions.values().stream().mapToInt(i -> i.size()).sum();
			int totalLengthOfSessions = 0;
			double totalHours = 0;
			List<Number> hours = new ArrayList<>();
			for (String user : sessionsNumberItems.keySet()) {
				totalLengthOfSessions += sessionsNumberItems.get(user).values().stream().mapToInt(i -> i.intValue()).sum();
				
				//Compute the total hours
				for (String session : sessionsTimeStamps.get(user).keySet()) {
					double h = TimeStampUtils.getHoursListTime(sessionsTimeStamps.get(user).get(session));
					totalHours += h;
					hours.add(h);
				}
			}

			output.println("Total users: " + totalUsers);
			output.println("Total items: " + items.size());
			output.println("Total preferences: " + totalPreferences.get());
			output.println("Preferences/Users: " + (double) totalPreferences.get() / (double) totalUsers);
			output.println("Preferences/Items: " + (double) totalPreferences.get() / (double) items.size());
			output.println("Median Hours-Sessions: " + (double) TimeStampUtils.getMedian(hours));
			output.println("\n\n");
			
			output.println("Number of sessions: " + totalSessions);
			output.println("TotalSessions/Users: " + (double) totalSessions / (double) totalUsers);
			output.println("TotalHoursSessions/Sessions: " + totalHours / totalSessions);
			output.println("TotalInterations/Sessions: " + (double) totalLengthOfSessions / (double) totalSessions);
			output.println("TotalIterations/Users: " + (double) totalLengthOfSessions / (double) totalUsers);
			output.println("\n\n");

			if (pathUserLocals != null) {
				// Show same stats but with locals and tourists
				int totalSessionsLocals = 0;
				int totalSessionsTourists = 0;
				int totalSessionsBots = 0;
				
				int totalLocals = 0;
				int totalTourists = 0;
				int totalBots = 0;
				
				int totalLengthOfSessionsLocals = 0;
				int totalLengthOfSessionsTourists = 0;
				int totalLengthOfSessionsBots = 0;
				
				double totalHoursLocals = 0;
				double totalHoursTourists = 0;
				double totalHoursBots = 0;
				
				for (String user : sessionsNumberItems.keySet()) {
					if (usersLocalOrTourist.get(user).equals("L")) {
						totalLengthOfSessionsLocals += sessionsNumberItems.get(user).values().stream()
								.mapToInt(i -> i.intValue()).sum();
					} 
					if (usersLocalOrTourist.get(user).equals("T")){
						totalLengthOfSessionsTourists += sessionsNumberItems.get(user).values().stream()
								.mapToInt(i -> i.intValue()).sum();
					}
					
					if (usersLocalOrTourist.get(user).equals("B")){
						totalLengthOfSessionsBots += sessionsNumberItems.get(user).values().stream()
								.mapToInt(i -> i.intValue()).sum();
					}
					
				}
				List<Number> hoursLocals = new ArrayList<>();
				List<Number> hoursTourist = new ArrayList<>();
				List<Number> hoursBots = new ArrayList<>();
				
				for (String user : totalUsersNumberSessions.keySet()) {
					if (usersLocalOrTourist.get(user).equals("L")) {
						totalSessionsLocals += totalUsersNumberSessions.get(user).size();
						totalLocals ++;
						for (String session: sessionsTimeStamps.get(user).keySet()) {
							double h = TimeStampUtils.getHoursListTime(sessionsTimeStamps.get(user).get(session));
							totalHoursLocals += h;
							hoursLocals.add(h);
						}
						
					} 
					if (usersLocalOrTourist.get(user).equals("T")) {
						totalSessionsTourists += totalUsersNumberSessions.get(user).size();
						totalTourists ++;
						for (String session: sessionsTimeStamps.get(user).keySet()) {
							double h = TimeStampUtils.getHoursListTime(sessionsTimeStamps.get(user).get(session));
							totalHoursTourists += h;
							hoursTourist.add(h);
						}
						
					}
					
					if (usersLocalOrTourist.get(user).equals("B")) {
						totalSessionsBots += totalUsersNumberSessions.get(user).size();
						totalBots ++;
						for (String session: sessionsTimeStamps.get(user).keySet()) {
							double h = TimeStampUtils.getHoursListTime(sessionsTimeStamps.get(user).get(session));
							totalHoursBots += h;
							hoursBots.add(h);
						}
						
					}
					
					
				}
				output.println("Total locals: " + totalLocals);
				output.println("Total tourists: " + totalTourists);
				output.println("Total bots: " + totalBots);
				
				output.println("Number of sessions locals: " + totalSessionsLocals);
				output.println("Number of sessions tourists: " + totalSessionsTourists);
				output.println("Number of sessions bots: " + totalSessionsBots);
				
				output.println("SessionsLocals/locals: " + (double) totalSessionsLocals / (double) totalLocals);
				output.println("SessionsTourists/tourists: " + (double) totalSessionsTourists / (double) totalTourists);
				output.println("SessionsBots/Bots: " + (double) totalSessionsBots / (double) totalBots);
				
				output.println("HoursLocals/SessionsLocals: " + totalHoursLocals / totalSessionsLocals);
				output.println("HoursTourists/SessionTurists: " + totalHoursTourists / totalSessionsTourists);
				output.println("HoursBots/SessionBots: " + totalHoursBots / totalSessionsBots);
				
				output.println("Median Hours-Sessions-Locals: " + (double) TimeStampUtils.getMedian(hoursLocals));
				output.println("Median Hours-Sessions-Tourist: " + (double) TimeStampUtils.getMedian(hoursTourist));
				output.println("Median Hours-Sessions-Bots: " + (double) TimeStampUtils.getMedian(hoursBots));
				
				output.println("POIsLocals/SessionLocals: " + (double) totalLengthOfSessionsLocals / (double) totalSessionsLocals);
				output.println("POIsTourists/SessionTourists: " + (double) totalLengthOfSessionsTourists / (double) totalSessionsTourists);
				output.println("POIsBots/SessionBots: " + (double) totalLengthOfSessionsBots / (double) totalSessionsBots);

			}
			output.close();
			
			

		} catch (Exception e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}

	}
	
	/***
	 * Method to obtain a result file for tourists and locals
	 * @param dataTrain
	 * @param outputLocalTourist
	 * @param maxDiffBetweenFirstAndLastDate
	 */
	public static void obtainLocalsAndTouristsByTime(String dataTrain, String outputLocalTourist, long maxDiffBetweenFirstAndLastDate, long minDiffforBots, int minPrefsToBeCOnsideredBot) {
		
		//userid -> list of timestamps
		Map<String, List<Long>> usersTimestamps = new HashMap<>();
		
		//userid -> set of sessions
		Map<String, Set<String>> usersSessions = new HashMap<>();
		
		//userid -> List of POIS (can be repeated)
		Map<String, List<String>> usersPois = new HashMap<>();
		Stream<String> stream;
		try {
			PrintStream out = new PrintStream(outputLocalTourist);
			stream = Files.lines(Paths.get(dataTrain));
			stream.forEach(line -> {
				String data[] = line.split("\t");
				String userId = data[0];
				String itemId = data[1];
				Double score = Double.parseDouble(data[2]);
				Long time = Long.parseLong(data[3]);
				String session = "NA"; 
				if (data.length > 4) {
					session = data[4];
				}

				if (usersTimestamps.get(userId) == null) {
					usersTimestamps.put(userId, new ArrayList<>());
					usersSessions.put(userId, new HashSet<>());
					usersPois.put(userId, new ArrayList<>());
				}
				
				usersTimestamps.get(userId).add(time);
				usersSessions.get(userId).add(session);
				usersPois.get(userId).add(itemId);
			});
			stream.close();
			//All timestamps are stored -> sort and see if it is local (L) or Tourist (T)
			
			out.println("UserId" + "\t" + "L-T-B" + "\t" + "FirstTimestamp" + "\t" + "LastTimestamp" + "\t"+ "HoursDiffFstLstTime" + "\t" + "DaysDiff" + "\t" + "PoisVisited" + "\t" + "UniquePoisVisited" + "\t" + "NumberSessions" + "\t" + "PoisPerSession" + "\t" + "UniquePoisPerSession");
			for (String user: usersTimestamps.keySet()) {
				Collections.sort(usersTimestamps.get(user));
				
				Long fst = usersTimestamps.get(user).get(0);
				Long lst = usersTimestamps.get(user).get(usersTimestamps.get(user).size() - 1);
				double hoursDiff = (lst - fst) / 3600.0;
				double daysDiff = hoursDiff / 24.0;
				String localTouristBot = "";
				Set<String> uniquePoisUser = new HashSet<>(usersPois.get(user));
				boolean isbot = isBot(usersTimestamps.get(user), minDiffforBots, minPrefsToBeCOnsideredBot);
				if (isbot) {
					localTouristBot = "B";
				} else {
				
					if ((lst - fst) > maxDiffBetweenFirstAndLastDate) {
						localTouristBot = "L";
					} else {
						localTouristBot = "T";
					}
				}
				out.println(user + "\t" + localTouristBot + "\t" + fst + "\t" + lst + "\t" + hoursDiff + "\t" + daysDiff + "\t" + usersPois.get(user).size() + 
						"\t" + uniquePoisUser.size() + "\t" + usersSessions.get(user).size() + "\t" + (double)usersPois.get(user).size() / (double) usersSessions.get(user).size() + "\t" + (double)uniquePoisUser.size() / (double) usersSessions.get(user).size() );
			}
			out.close();
			
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	/***
	 * Method to obtain a result file for tourists and locals
	 * @param dataTrain
	 * @param outputLocalTourist
	 * @param maxDiffBetweenFirstAndLastDate
	 */
	public static void obtainLocalsAndTouristsByAverageSessionLenghts(String dataTrain, String outputLocalTourist, long minDiffforBots, int minPrefsToBeConsideredBot, double averageSessionLocalTourist) {
		
		//userid -> list of timestamps
		Map<String, List<Long>> usersTimestamps = new HashMap<>();
		
		//userid -> number of session -> length of session
		Map<String, Map<String, Integer>> usersSessions = new HashMap<>();
		
		//userid -> List of POIS (can be repeated)
		Map<String, List<String>> usersPois = new HashMap<>();
		Stream<String> stream;
		try {
			PrintStream out = new PrintStream(outputLocalTourist);
			stream = Files.lines(Paths.get(dataTrain));
			stream.forEach(line -> {
				String data[] = line.split("\t");
				String userId = data[0];
				String itemId = data[1];
				Double score = Double.parseDouble(data[2]);
				Long time = Long.parseLong(data[3]);
				String session = "NA"; 
				if (data.length > 4) {
					session = data[4];
				}

				if (usersTimestamps.get(userId) == null) {
					usersTimestamps.put(userId, new ArrayList<>());
					usersSessions.put(userId, new HashMap<>());
					usersPois.put(userId, new ArrayList<>());
				}
				
				usersTimestamps.get(userId).add(time);
				
				if (usersSessions.get(userId).get(session) == null) {
					usersSessions.get(userId).put(session, 0);
				}
				usersSessions.get(userId).put(session, usersSessions.get(userId).get(session) + 1);
				
				
				usersPois.get(userId).add(itemId);
			});
			stream.close();
			//All timestamps are stored -> sort and see if it is local (L) or Tourist (T)
			
			out.println("UserId" + "\t" + "L-T-B" + "\t" + "FirstTimestamp" + "\t" + "LastTimestamp" + "\t" + "HoursDiffFstLstTime" + "\t" + "PoisVisited" + "\t" + "UniquePoisVisited" + "\t" + "NumberSessions" + "\t" + "PoisPerSession" + "\t" + "UniquePoisPerSession");
			for (String user: usersTimestamps.keySet()) {
				Collections.sort(usersTimestamps.get(user));
				
				Long fst = usersTimestamps.get(user).get(0);
				Long lst = usersTimestamps.get(user).get(usersTimestamps.get(user).size() - 1);
				double hoursDiff = (lst - fst) / 3600.0;
				String localTouristBot = "";
				Set<String> uniquePoisUser = new HashSet<>(usersPois.get(user));
				boolean isbot = isBot(usersTimestamps.get(user), minDiffforBots, minPrefsToBeConsideredBot);
				
				Map<String, Integer> lengthSessionsUser = usersSessions.get(user);
				double averageSessionLengthUser = 0;
				
				for (String sessionID: lengthSessionsUser.keySet()) {
					averageSessionLengthUser += lengthSessionsUser.get(sessionID);
				}
				averageSessionLengthUser /= (double) lengthSessionsUser.keySet().size();
				
				
				if (isbot) {
					localTouristBot = "B";
				} else {
				
					if ( averageSessionLocalTourist < averageSessionLengthUser ) {
						localTouristBot = "L";
					} else {
						localTouristBot = "T";
					}
				}
				out.println(user + "\t" + localTouristBot + "\t" + fst + "\t" + lst + "\t" + hoursDiff + "\t" + usersPois.get(user).size() + 
						"\t" + uniquePoisUser.size() + "\t" + usersSessions.get(user).size() + "\t" + (double)usersPois.get(user).size() / (double) usersSessions.get(user).size() + "\t" + (double)uniquePoisUser.size() / (double) usersSessions.get(user).size() );
			}
			out.close();
			
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	/***
	 * Method to obtain a result file for tourists and locals by the categories they visit
	 * @param dataTrain
	 * @param outputLocalTourist
	 * @param maxNumberOfDifferentCategories
	 */
	public static void obtainLocalsAndTouristsByCategories(String dataTrain, String outputLocalTourist, long maxNumberVisitedCategories, String itemFeatureFile) {	

		
		//For each user we store the features of the Items she consumes
		Map<String, Set<String>> usersFeatures = new HashMap<>();


		Stream<String> stream;
		try {
			PrintStream out = new PrintStream(outputLocalTourist);
			FeatureData<String, String, Double> featureData = SimpleFeatureData.load(SimpleFeaturesReader.get().read(itemFeatureFile, sp, sp));

			stream = Files.lines(Paths.get(dataTrain));
			stream.forEach(line -> {
				String data[] = line.split("\t");
				String userId = data[0];
				String itemId = data[1];

				if (usersFeatures.get(userId) == null) {
					usersFeatures.put(userId, new HashSet<>());
				}
				Set<String> features = featureData.getItemFeatures(itemId).map(t -> t.v1).collect(Collectors.toSet());
				if (features != null) {
					usersFeatures.get(userId).addAll(features);
				}
			});
			stream.close();
			
			out.println("UserId" + "\t" + "L-T-B" + "\t" + "Number of different Categories visited");
			for (String user: usersFeatures.keySet()) {
					String localTouristBot="";
					int numberCat = usersFeatures.get(user).size();
				
					if (maxNumberVisitedCategories > numberCat) {
						localTouristBot = "L";
					} else {
						localTouristBot = "T";
					}
				
				out.println(user + "\t" + localTouristBot + "\t" + numberCat);
			}
			out.close();
			
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		
		
	}
	
	
	public static <U extends Comparable<U>, I extends Comparable<I>> void filterOutBotsTimestampUsers(FastTemporalPreferenceDataIF<U, I> data, String outputFile, long minDiffforBots, int minPrefsToBeConsideredBot) {
		PrintStream out;
		try {
			out = new PrintStream(outputFile);
			data.getUsersWithPreferences().forEach(user -> {
				List<Long> userTimestamps = data.getUserPreferences(user).map(timepref -> timepref.v3).collect(Collectors.toList());
				Collections.sort(userTimestamps);				
				if (!isBot(userTimestamps, minDiffforBots, minPrefsToBeConsideredBot)) {
					data.getUserPreferences(user).sorted(PreferenceComparators.timeComparatorIdTimePref()).forEach(pref -> {
						out.println(user + "\t" + pref.v1 + "\t" + pref.v2 + "\t" + pref.v3);
					});
				}
			});
			out.close();
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	
		
	}
	
	/***
	 * Method to detect if a user is a bot by a list of timestamps
	 * @param timestamps the lsit of timestamps
	 * @param minDiffforBots the difference between 2 timestamps to be considered as bot
	 * @param minPrefsToBeConsideredBot the number of times that the difference need to be lower so it is considered as bot
	 * @return
	 */
	public static boolean isBot(List<Long> timestamps, long minDiffforBots, int minPrefsToBeConsideredBot) {
		Long actualTime = timestamps.get(0);
		int count = 0;
		for (int i = 1; i < timestamps.size(); i++) {
			Long newTime = timestamps.get(i);
			Long diff = Math.abs(actualTime - newTime);
			if (diff <= minDiffforBots) {
				count ++;
			}
			actualTime = newTime;
		}
		return count >= minPrefsToBeConsideredBot;
	}
	
	
	

	/***
	 * Method to make a global split (oldest ratings to train, newest to test)
	 * 
	 * @param srcPath      the source path
	 * @param dstPathTrain the destination of the train subset
	 * @param dstPathTest  the destination of the test subset
	 * @param trainPercent the percentage of ratings that the train subset must have
	 */
	public static void DatasetTemporalGlobalSplit(String srcPath, String dstPathTrain, String dstPathTest,
			double trainPercent) {

		List<Preference> allRatings = new ArrayList<>();
		String characterSplit = "\t";

		// Parameters to configure
		int columnUser = 0;
		int columnItem = 1;
		int columnRating = 2;
		int columnTimeStamp = 3;
		boolean ignoreFirstLine = false; // If we want to ignore first line, switch to true

		Stream<String> stream = null;
		try {
			if (ignoreFirstLine) {
				stream = Files.lines(Paths.get(srcPath)).skip(1);
			} else {
				stream = Files.lines(Paths.get(srcPath));
			}
			stream.forEach(line -> {

				String[] data = line.split(characterSplit);
				Long idUser = Long.parseLong(data[columnUser]);
				Long idItem = Long.parseLong(data[columnItem]);
				Double rating = Double.parseDouble((data[columnRating]));
				Long timeStamp = Long.parseLong(data[columnTimeStamp]);
				Preference p = new Preference(idUser, idItem, rating, timeStamp);
				allRatings.add(p);

			});
			stream.close();

			PrintStream writer = null;
			PrintStream writer2 = null;

			writer = new PrintStream(dstPathTrain);
			writer2 = new PrintStream(dstPathTest);
			// Now temporal split
			Collections.sort(allRatings);
			int indexLimit = (int) (trainPercent * allRatings.size());

			for (int i = 0; i < allRatings.size(); i++) {
				Preference p = allRatings.get(i);
				if (i < indexLimit)
					writer.println(
							p.getIdUser() + "\t" + p.getIdItem() + "\t" + p.getRating() + "\t" + p.getTimeStamp());
				else
					writer2.println(
							p.getIdUser() + "\t" + p.getIdItem() + "\t" + p.getRating() + "\t" + p.getTimeStamp());
			}

			writer2.close();
			writer.close();

		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}
	
	public static void DatasetTemporalGlobalSplitByTimeStamps(String srcPath, String dstPathTrain, String dstPathTest, 
			long startingTrainingTimestamp, long startingTestTimestamp, long endingTestTimestamp) {

		String characterSplit = "\t";

		// Parameters to configure
		int columnUser = 0;
		int columnItem = 1;
		int columnRating = 2;
		int columnTimeStamp = 3;
		boolean ignoreFirstLine = false; // If we want to ignore first line, switch to true

		Stream<String> stream = null;
		try {
			if (ignoreFirstLine) {
				stream = Files.lines(Paths.get(srcPath)).skip(1);
			} else {
				stream = Files.lines(Paths.get(srcPath));
			}
			
			PrintStream writer = null;
			PrintStream writer2 = null;
			
			writer = new PrintStream(dstPathTrain);
			writer2 = new PrintStream(dstPathTest);
			
			PrintStream writerTraining = writer;
			PrintStream writerTest = writer2;

			
			stream.forEach(line -> {

				String[] data = line.split(characterSplit);
				Long idUser = Long.parseLong(data[columnUser]);
				Long idItem = Long.parseLong(data[columnItem]);
				Double rating = Double.parseDouble((data[columnRating]));
				Long timeStamp = Long.parseLong(data[columnTimeStamp]);
				
				//if the timestamp is between  startingTrainingTimestamp and startingTestTimestamp, then train
				
				if (timeStamp >= startingTrainingTimestamp && timeStamp < startingTestTimestamp) {
					writerTraining.println(idUser + "\t" +idItem + "\t" + rating + "\t" + timeStamp);
				} else if (timeStamp >= startingTestTimestamp && timeStamp < endingTestTimestamp) {
					writerTest.println(idUser + "\t" +idItem + "\t" + rating + "\t" + timeStamp);
				}
				
				

			});
			stream.close();

			writer2.close();
			writer.close();

			
			

			

		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}
	
	public static void DatasetOrderRatingsPerUser(String srcPath, String dstPathDest) {

		List<Preference> allRatings = new ArrayList<>();
		String characterSplit = "\t";

		// Parameters to configure
		int columnUser = 0;
		int columnItem = 1;
		int columnRating = 2;
		int columnTimeStamp = 3;
		boolean ignoreFirstLine = false; // If we want to ignore first line, switch to true

		Stream<String> stream = null;
		try {
			if (ignoreFirstLine) {
				stream = Files.lines(Paths.get(srcPath)).skip(1);
			} else {
				stream = Files.lines(Paths.get(srcPath));
			}
			stream.forEach(line -> {

				String[] data = line.split(characterSplit);
				Long idUser = Long.parseLong(data[columnUser]);
				Long idItem = Long.parseLong(data[columnItem]);
				Double rating = Double.parseDouble((data[columnRating]));
				Long timeStamp = Long.parseLong(data[columnTimeStamp]);
				Preference p = new Preference(idUser, idItem, rating, timeStamp);
				allRatings.add(p);

			});
			stream.close();

			PrintStream writer = null;

			writer = new PrintStream(dstPathDest);
			// Now temporal split
			Collections.sort(allRatings, new PreferenceByUserComparator());

			for (int i = 0; i < allRatings.size(); i++) {
				Preference p = allRatings.get(i);
					writer.println(
							p.getIdUser() + "\t" + p.getIdItem() + "\t" + p.getRating() + "\t" + p.getTimeStamp());
			}

			writer.close();

		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}
	
	

	/***
	 * Method to obtain splits user-dependent. For each user, if the number of
	 * ratings in train are lower or equal to the total preferences of the user
	 * minus the number stipulated to the test set, then the last n ratings will be
	 * put to the train subset and the rest to the train subset. If the user does
	 * not match the previous condition will be put only to the train. In order to
	 * be in the test set, the user must have at least numberOfTrain+numberOfTest
	 * preferences
	 * 
	 * @param srcpath            source route of the path.
	 * @param dstPathTrain       destination of the train subset
	 * @param dstPathTest        destination of the test subset
	 * @param numberOfTest       the number of test preferences obligatory
	 * @param minimumNumberTrain the minimum number of train preferences that the
	 *                           user must have
	 */
	public static void TimeSplitPerUser(String srcpath, String dstPathTrain, String dstPathTest, int numberOfTest,
			int minimumNumberTrain) {
		// Assume structure of the file is: UserdId(long) ItemID(long) rating(double)
		// timestamp(long)

		Map<Long, ArrayList<Tuple3<Long, Double, Long>>> element1s = new TreeMap<Long, ArrayList<Tuple3<Long, Double, Long>>>();
		String characterSplit = "\t";

		// Parameters to configure
		int columnUser = 0;
		int columnItem = 1;
		int columnRating = 2;
		int columnTimeStamp = 3;

		BufferedWriter writer = null;
		BufferedWriter writer2 = null;

		Stream<String> stream = null;
		try {
			stream = Files.lines(Paths.get(srcpath));
			writer = new BufferedWriter(new FileWriter(dstPathTrain));
			writer2 = new BufferedWriter(new FileWriter(dstPathTest));

			stream.forEach(line -> {
				String[] data = line.split(characterSplit);
				Long idUser = Long.parseLong(data[columnUser]);
				if (element1s.get(idUser) == null) { // New user
					ArrayList<Tuple3<Long, Double, Long>> lst = new ArrayList<Tuple3<Long, Double, Long>>();
					Long idItem = Long.parseLong(data[columnItem]);
					double rating = Double.parseDouble((data[columnRating]));
					Long timestamp = Long.parseLong(data[columnTimeStamp]);
					lst.add(new Tuple3<Long, Double, Long>(idItem, rating, timestamp));
					// Add item to list
					element1s.put(idUser, lst);
				} else {
					ArrayList<Tuple3<Long, Double, Long>> lst = element1s.get(idUser);
					Long idItem = Long.parseLong(data[columnItem]);
					double rating = Double.parseDouble((data[columnRating]));
					Long timestamp = Long.parseLong(data[columnTimeStamp]);
					// Add item to list
					lst.add(new Tuple3<Long, Double, Long>(idItem, rating, timestamp));
				}
			});
			// Now temporal split

			for (Long user : element1s.keySet()) {
				ArrayList<Tuple3<Long, Double, Long>> lst = element1s.get(user);
				// Order data
				Collections.sort(lst, new TimeStampComparator());

				// Number of elements in test must be lower than the total size - numberOfTest.
				// If not, then go to train
				if (minimumNumberTrain <= lst.size() - numberOfTest) {
					// Store test items
					for (int i = 0; i < numberOfTest; i++) {
						Tuple3<Long, Double, Long> t = lst.remove(lst.size() - 1);
						writer2.write(user + "\t" + t.v1 + "\t" + t.v2 + "\t" + t.v3 + "\n");
					}
				}
				// Rest of the list, to train
				for (Tuple3<Long, Double, Long> t : lst)
					writer.write(user + "\t" + t.v1 + "\t" + t.v2 + "\t" + t.v3 + "\n");

			}
			writer2.close();
			writer.close();

		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}
	
	public static void TimeSplitPerUserPercentage(String srcpath, String dstPathTrain, String dstPathTest, float percentageTrain) {
		
		Map<Long, ArrayList<Tuple3<Long, Double, Long>>> element1s = new TreeMap<Long, ArrayList<Tuple3<Long, Double, Long>>>();
		String characterSplit = "\t";

		// Parameters to configure
		int columnUser = 0;
		int columnItem = 1;
		int columnRating = 2;
		int columnTimeStamp = 3;

		BufferedWriter writer = null;
		BufferedWriter writer2 = null;

		Stream<String> stream = null;
		try {
			stream = Files.lines(Paths.get(srcpath));
			writer = new BufferedWriter(new FileWriter(dstPathTrain));
			writer2 = new BufferedWriter(new FileWriter(dstPathTest));

			stream.forEach(line -> {
				String[] data = line.split(characterSplit);
				Long idUser = Long.parseLong(data[columnUser]);
				if (element1s.get(idUser) == null) { // New user
					ArrayList<Tuple3<Long, Double, Long>> lst = new ArrayList<Tuple3<Long, Double, Long>>();
					Long idItem = Long.parseLong(data[columnItem]);
					double rating = Double.parseDouble((data[columnRating]));
					Long timestamp = (long) Double.parseDouble(data[columnTimeStamp]);
					lst.add(new Tuple3<Long, Double, Long>(idItem, rating, timestamp));
					// Add item to list
					element1s.put(idUser, lst);
				} else {
					ArrayList<Tuple3<Long, Double, Long>> lst = element1s.get(idUser);
					Long idItem = Long.parseLong(data[columnItem]);
					double rating = Double.parseDouble((data[columnRating]));
					Long timestamp = (long) Double.parseDouble(data[columnTimeStamp]);
					// Add item to list
					lst.add(new Tuple3<Long, Double, Long>(idItem, rating, timestamp));
				}
			});
			// Now temporal split

			for (Long user : element1s.keySet()) {
				ArrayList<Tuple3<Long, Double, Long>> lst = element1s.get(user);
				// Order data
				Collections.sort(lst, new TimeStampComparator());

				int total = lst.size();
				int numberTrain = Math.round(total * percentageTrain);
				int numberTest = total - numberTrain;
				

					// Store test items
					for (int i = 0; i < numberTest; i++) {
						Tuple3<Long, Double, Long> t = lst.remove(lst.size() - 1);
						writer2.write(user + "\t" + t.v1 + "\t" + t.v2 + "\t" + t.v3 + "\n");
					}
				
				// Rest of the list, to train
				for (Tuple3<Long, Double, Long> t : lst)
					writer.write(user + "\t" + t.v1 + "\t" + t.v2 + "\t" + t.v3 + "\n");

			}
			writer2.close();
			writer.close();

		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
	}
	
	
	private static boolean stringContainsAtLeastOne(String totalString, Set<String> matchings) {
		for (String s: matchings) {
			if (totalString.contains(s)) {
				return true;
			}
		}
		return false;
	}
	
	private static boolean stringContainsAll(String totalString, Set<String> matchings) {
		for (String s: matchings) {
			if (!totalString.contains(s)) {
				return false;
			}
		}
		return true;
	}

	/**
	 * Method to obtain the best recommenders of a different string matchings
	 * For the file containing the results of all recommenders, it will selected for each string matching recommender, the one
	 * that maximizes the metric stated.
	 * 
	 * @param sourceResultFile source file containing result of all the recommenders
	 * @param destinationFile		the destination file to write the best recommenders
	 * @param selectedMetric		the metric to maximize
	 * @param stringMatchingRecommenders the different recommenders to find the best
	 * @param matchingStrings set of strings tht all recommenders MUST have in order to be processed
	 * @param recommenderColumn     the index of the recommender file
	 * @param resultsColumn   the index of the result file
	 * @param metricColumn   the index of the metric
	 */
	public static void obtainBestRecommenders(String sourceResultFile, String destFile, String selectedMetric, String[] stringMatchingRecommenders, 
			String [] matchingStrings, int metricColumn, int resultsColumn, int recommenderColumn) {
		
		Set<String> possibleMatchingRecommenders = new HashSet<String>(Arrays.asList(stringMatchingRecommenders));
		Set<String> requiredMatchingStrings = new HashSet<String>(Arrays.asList(matchingStrings));

		
		// Map of recommenders -> metric, result
		Map<String, Map<String, Double>> results = new TreeMap<String, Map<String, Double>>();
		//First, read the whole file and store the recommenders that contains any string in the array stringMatchingRecommenders
		try (Stream<String> stream = Files.lines(Paths.get(sourceResultFile))) {
			stream.forEach(line -> {
				String data [] = line.split("\t");
				String recRead = data [recommenderColumn];
				String metricRead = data [metricColumn];
				Double metricValue = Double.parseDouble(data[resultsColumn]);
				
				//It is a valid recommender
				if (stringContainsAtLeastOne(recRead, possibleMatchingRecommenders) && stringContainsAll(recRead, requiredMatchingStrings)) {
					if (results.get(recRead) == null) {
						results.put(recRead, new TreeMap<>());
					}
					results.get(recRead).put(metricRead, metricValue);

				}
				
			});
			//Now we have the valid recommenders
		
			PrintStream writer = null;
			writer = new PrintStream(destFile);

			// Now, show the best recommenders matching
			for (String recommendersMathing : possibleMatchingRecommenders) {
				//Obtain the recommenders matching
				List<String> recToProcess = results.keySet().stream().filter(pref -> pref.contains(recommendersMathing))
						.collect(Collectors.toList());
				
				//There can be ties, so we will retrieve all the recommenders matching
				List<String> selectedBestRec = new ArrayList<>();
				double maxVal = Double.MIN_VALUE;
				for (String rec : recToProcess) {
					double val = results.get(rec).get(selectedMetric);
					if (maxVal < val) {
						selectedBestRec.clear();
					}
					if (maxVal <= val) {
						maxVal = val;
						selectedBestRec.add(rec);
					}
				}
				
				//Print
				for (String rec: selectedBestRec) {
					for (String m : results.get(rec).keySet()) {
						writer.println(rec + "\t" + m + "\t" + results.get(rec).get(m));
					}
				}
				
			}
			writer.close();

		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}

	/***
	 * Method to parse MyMedialLte recommendations
	 * 
	 * @param sourceFile the source file
	 * @param destFile   the new recommedner file
	 */
	public static void parseMyMediaLite(String sourceFile, FastUserIndex<Long> userIndexTest, String destFile) {
		BufferedReader br;
		PrintStream writer = null;
		try {
			br = new BufferedReader(new FileReader(sourceFile));
			String line = br.readLine();
			writer = new PrintStream(destFile);

			while (line != null) {
				if (line != null) {
					String[] fullData = line.split("\t");
					Long user = Long.parseLong(fullData[0]);
					line = fullData[1].replace("[", "");
					line = line.replace("]", "");
					String[] data = line.split(",");
					if (userIndexTest.containsUser(user)) {
						for (String itemAndPref : data) {
							Long item = Long.parseLong(itemAndPref.split(":")[0]);
							Double preference = Double.parseDouble(itemAndPref.split(":")[1]);
							writer.println(user + "\t" + item + "\t" + preference);
						}
					}

					line = br.readLine();
				}
			}
			br.close();
			writer.close();

		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	/**
	 * Method to filter out the users and items that does not have a minimum number
	 * of ratings It is NOT a K-core
	 * 
	 * @param sourceFile      the original (full dataset)
	 * @param destFile        the destination file
	 * @param minimumRatUsers the minimum ratings per user
	 * @param minimumRatItems the minimum ratings per item
	 */
	public static void usersItemsMinimumRatings(String sourceFile, String destFile, int minimumRatUsers,
			int minimumRatItems) {
		// Assume structure of the file is: UserdId(String) ItemID(String)
		// rating(double) timestamp(long)

		BufferedReader br;
		// Map of users String(user), Integer (number of ratings per User)
		Map<String, Integer> userMap = new HashMap<String, Integer>();

		// Map of users String(item), Integer (number of ratings per Item)
		Map<String, Integer> itemMap = new HashMap<String, Integer>();

		String characterSplit = "\t";

		// Parameters to configure
		int columnUser = 0;
		int columnItem = 1;

		// switch to true

		PrintStream writer = null;

		try {

			// First pass of the file
			br = new BufferedReader(new FileReader(sourceFile));

			String line = br.readLine();
			System.out.println("First pass");
			while (line != null) {
				if (line != null) {
					String[] data = line.split(characterSplit);
					String idUser = data[columnUser];
					String idItem = data[columnItem];

					// If the user does exist, we increment in 1 the number of ratings
					if (userMap.containsKey(idUser))
						userMap.put(idUser, userMap.get(idUser) + 1);
					else
						userMap.put(idUser, 1);

					if (itemMap.containsKey(idItem))
						itemMap.put(idItem, itemMap.get(idItem) + 1);
					else
						itemMap.put(idItem, 1);
				}
				line = br.readLine();
			}
			br.close();

			System.out.println("Second pass and writing in " + destFile);
			br = new BufferedReader(new FileReader(sourceFile));
			writer = new PrintStream(destFile);

			line = br.readLine();
			while (line != null) {
				if (line != null) {
					String[] data = line.split(characterSplit);
					String idUser = data[columnUser];
					String idItem = data[columnItem];

					if (userMap.get(idUser) >= minimumRatUsers && itemMap.get(idItem) >= minimumRatItems) {
						writer.println(line);
					}

				}
				line = br.readLine();
			}
			br.close();
			writer.close();

		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}

	public static void implicitToExplicitBins(String trainFile, String outputFile, String charactersSplit,
			int numberBins) {

		int columnUser = 0;
		int columnItem = 1;
		int columnRatingImplicits = 2;
		boolean ignoreFirstLine = true;
		Stream<String> stream = null;
		Map<String, List<Tuple2<String, Double>>> userPreferences = new HashMap<>();

		try {
			if (ignoreFirstLine) {
				stream = Files.lines(Paths.get(trainFile)).skip(1);
			}
			else {
				stream = Files.lines(Paths.get(trainFile));
			}
			stream.forEach(line -> {
				String[] data = line.split(charactersSplit);
				String user = data[columnUser];
				if (!userPreferences.containsKey(user)) {
					userPreferences.put(user, new ArrayList<>());
				}
				userPreferences.get(user)
						.add(new Tuple2<>(data[columnItem], Double.parseDouble(data[columnRatingImplicits])));
			});
			stream.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

		// Now, for every user, we order the list of tuples
		try {
			final PrintStream writer = new PrintStream(outputFile);
			for (String user : userPreferences.keySet()) {
				List<Tuple2<String, Double>> preferences = userPreferences.get(user);
				// If we have less than number of bins, we discard that user
				if (preferences.size() < numberBins)
					continue;
				Collections.sort(preferences, new WeightComparator());
				// All preferences sorted from the lower values to the higher

				double percent = 1.0 / numberBins;
				double percents = percent;

				double score = 1;

				for (int i = 0; i < preferences.size(); i++) {
					if (i >= percents * preferences.size()) {
						percents += percent;
						score++;
					}
					writer.println(user + "\t" + preferences.get(i).v1 + "\t" + score);
				}

			}

			writer.close();
			stream.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	/***
	 * Method to read a session file in the format we want
	 * @param originalFile the original file
	 * @return a map by userID, sessionId and time preferences
	 */
	public static Map<Long, Map<Long, List<IdTimePref<Long>>>> readSessionFormat(String originalFile){
		Map<Long, Map<Long, List<IdTimePref<Long>>>> mapSessions = new TreeMap<>();
		try {
			Stream<String> stream = Files.lines(Paths.get(originalFile));
			stream.forEach(line -> {
				String [] data = line.split("\t");
				Long user = Long.parseLong(data[0]);
				Long itemId = Long.parseLong(data[1]);
				Double pref = Double.parseDouble(data[2]);
				Long timestamp = Long.parseLong(data[3]);
				Long sessionId = Long.parseLong(data[4]);
				
				if (mapSessions.get(user) == null) {
					mapSessions.put(user, new TreeMap<>());
				}
				
				if (mapSessions.get(user).get(sessionId) == null) {
					mapSessions.get(user).put(sessionId, new ArrayList<>());
				}
				
				mapSessions.get(user).get(sessionId).add(new IdTimePref<>(itemId, pref, timestamp));
				
				
			});
			stream.close();
			return mapSessions;
		} catch (IOException e) {
			e.printStackTrace();
		}
		return null;

	}
	
	/****
	 * Method to filter the sessions of the user by taking the most n recent sessions
	 * @param mapSessions
	 * @param L the number of the lst N sessions
	 * @return the same map filtered
	 */
	public static Tuple2<Map<Long, Map<Long, List<IdTimePref<Long>>>>, Map<Long, Long>> filterMostLRecentSessions(Map<Long, Map<Long, List<IdTimePref<Long>>>> mapSessions, int L){
		Map<Long, Map<Long, List<IdTimePref<Long>>>> mapSessionsResult = new TreeMap<>();
		Map<Long, Long> userLastSession = new HashMap<>();
		
		//For each user
		for (Long user: mapSessions.keySet()) {
			Map<Long, List<IdTimePref<Long>>> sessionsUser = mapSessions.get(user);
			
			List<Tuple2<Long, Long>> idSessionFstTimeStamp = new ArrayList<Tuple2<Long,Long>>();
			for (Long idSession: sessionsUser.keySet()) {
				long lowestSession = Long.MAX_VALUE;
				List<IdTimePref<Long>> sessionUser = sessionsUser.get(idSession);
				for (IdTimePref<Long> preference: sessionUser) {
					if (preference.v3 < lowestSession)
						lowestSession = preference.v3;
				}
				idSessionFstTimeStamp.add(new Tuple2<Long, Long>(idSession, lowestSession));
			}
			Collections.sort(idSessionFstTimeStamp, PreferenceComparators.comparatorTuple2Second());
			mapSessionsResult.put(user, new TreeMap<>());
			
			userLastSession.put(user, idSessionFstTimeStamp.get(idSessionFstTimeStamp.size() - 1).v1);
			
			if (idSessionFstTimeStamp.size() <= L) {
				mapSessionsResult.get(user).putAll(sessionsUser);
			} else {
				//Only take the n most recent
				List<Tuple2<Long, Long>> sublist = idSessionFstTimeStamp.subList(idSessionFstTimeStamp.size() - L, idSessionFstTimeStamp.size());
				for (Tuple2<Long, Long> t: sublist) {
					mapSessionsResult.get(user).put(t.v1, sessionsUser.get(t.v1));
				}
			}
			
			
		}
		return new Tuple2<>(mapSessionsResult, userLastSession);
	}
	
	
	
	public static Map<Integer, Map<Long, List<Integer>>> readSessionFormatIdx(String originalFile, FastTemporalPreferenceDataIF<Long, Long> prefData){
		Map<Long, Map<Long, List<IdTimePref<Long>>>> mapSessions = readSessionFormat(originalFile);
	
		Map<Integer, Map<Long, List<Integer>>> result = new HashMap<>();
		for (Long user: mapSessions.keySet()) {
			int uidx = prefData.user2uidx(user);
			if (!result.containsKey(uidx)) {
				result.put(uidx, new HashMap<>());
			}
			for (Long sessionUser: mapSessions.get(user).keySet()) {
				List<Integer> itemsRated = mapSessions.get(user).get(sessionUser).stream().mapToInt(t -> prefData.item2iidx(t.v1)).boxed().collect(Collectors.toList());
				result.get(uidx).put(sessionUser, itemsRated);
			}
				
			
			
		}
		
		
		return result;

	}
	
	
	
	
	
	
	
	public static void fixSessionSplit(String originaFile, String trainFile, String testFile, int minNumberOfSessionsTest, int minimumNumberPoistoTest) {
		Map<Long, Map<Long, List<IdTimePref<Long>>>> mapSessions = readSessionFormat(originaFile);
		//Assume the original file is in this way
		try {
			PrintStream writerTrain = new PrintStream(trainFile);
			PrintStream writerTest = new PrintStream(testFile);
			PrintStream actual = null;
			for (Long userid: mapSessions.keySet()) {
				actual = writerTrain;
				int nSessions = mapSessions.get(userid).keySet().size();
				if (nSessions < minNumberOfSessionsTest) {
					for (Long sessionId: mapSessions.get(userid).keySet()) {
						for (IdTimePref<Long> t: mapSessions.get(userid).get(sessionId)) {
							actual.println(userid + "\t" + t.v1 + "\t" + t.v2 + "\t" + t.v3 + "\t" + sessionId);
						}
					}
				} else {
					int count = 1;
					for (Long sessionId: mapSessions.get(userid).keySet()) {
						if (count == nSessions && mapSessions.get(userid).get(sessionId).size() >= minimumNumberPoistoTest ) {
							actual = writerTest;
						}
							
						for (IdTimePref<Long> t: mapSessions.get(userid).get(sessionId)) {
							actual.println(userid + "\t" + t.v1 + "\t" + t.v2 + "\t" + t.v3 + "\t" + sessionId);	
						}
						count ++;
					}
				}
				
			}

			writerTest.close();
			writerTrain.close();

		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	/***
	 * Method to test the session file construction. It will check that:
	 * -The user sessions are ordered (the session with lower id have the lower timestamps)
	 * -The timestamps of the user are ordered
	 * @param originalFile 
	 */
	public static void checkSessionFiles(String originalFile) {
		Map<Long, Map<Long, List<IdTimePref<Long>>>> mapSessions = readSessionFormat(originalFile);
	
		Boolean correct = true;
		for (Long user: mapSessions.keySet()) {
			List<Long> sessions = new ArrayList<Long>(mapSessions.get(user).keySet());
			Collections.sort(sessions);
			Long maxTimestampUserPrevSess = Long.MIN_VALUE;
			for (Long sessionId: sessions) {
				Long minTimestampUserSess = mapSessions.get(user).get(sessionId).get(0).v3();
				Long maxTimestampUserSess = mapSessions.get(user).get(sessionId).get(mapSessions.get(user).get(sessionId).size() - 1).v3();
				if (minTimestampUserSess > maxTimestampUserSess) {
					System.out.println("User " + user + " session " + sessionId + " is not ordered");
					correct = false;
				}

				if (maxTimestampUserPrevSess > minTimestampUserSess) {
					System.out.println("User " + user + " session " + sessionId + " have a minimum timestamp higher than the highest timestamp of the previous session");
					correct = false;
				}
				
				maxTimestampUserPrevSess = maxTimestampUserSess;
				
			}
		}
		
		if (correct) {
			System.out.println("Session file is correct");
		}

		
	}
	
	/***
	 * Method to test the session file construction. It will check that:
	 * -The user sessions are ordered (the session with lower id have the lower timestamps)
	 * -The timestamps of the user are ordered
	 * @param originalFile 
	 */
	public static void checkSessionFilesTrainTest(String trainFile, String testFile) {
		Map<Long, Map<Long, List<IdTimePref<Long>>>> mapSessions = readSessionFormat(trainFile);
		Map<Long, Map<Long, List<IdTimePref<Long>>>> mapSessionsTest = readSessionFormat(testFile);

	
		Boolean correct = true;
		for (Long userTest: mapSessionsTest.keySet()) {
			List<Long> timestampsUserTest = new ArrayList<Long>();
			for (Long session: mapSessionsTest.get(userTest).keySet()) {
				timestampsUserTest.addAll(mapSessionsTest.get(userTest).get(session).stream().mapToLong(t -> t.v3).boxed().collect(Collectors.toList()));
			}
			
			Collections.sort(timestampsUserTest);
			
			long lowest = timestampsUserTest.get(0);
			
			List<Long> timestampsUserTrain = new ArrayList<Long>();
			for (Long session: mapSessions.get(userTest).keySet()) {
				timestampsUserTrain.addAll(mapSessions.get(userTest).get(session).stream().mapToLong(t -> t.v3).boxed().collect(Collectors.toList()));
			}
			
			Collections.sort(timestampsUserTrain);
			
			for (Long time: timestampsUserTrain) {
				if (time > lowest) {
					correct = false;
					System.out.println("User: " + userTest + " has a timestmap in train higher than in test");
				}
			}

			
			
			
		}
		
		if (correct) {
			System.out.println("For  each user, the timestamps in test are higher than the ones in train");
		}

		
	}
	
	
	
	
	public static void analyzeTour(String sessionFile, String pathItemOldIdNewID, String dstFile, int minNumberPois, double minHours, double maxHours) {
		Map<String, String> newOld = new HashMap<>();
		
		//Id user -> Id session -> List tuples ids, time 
		Map<String, Map<String, List<Tuple2<String, Long>>>> sessionsTimeStamps = new HashMap<>();

		try {
			Stream<String> stream = Files.lines(Paths.get(pathItemOldIdNewID));
			stream.forEach(line -> {
				String [] data = line.split("\t");
				newOld.put(data[1], data[0]);
			});
			stream.close();
			Stream<String> streamTr = Files.lines(Paths.get(sessionFile));
			streamTr.forEach(line -> {
				String data [] = line.split("\t");
				String userId = data[0];
				String itemId = data[1];
				Long time = Long.parseLong(data[3]);
				String session = data[4];
				
				// Map of users, sessions and items per session. First, see if the user exists
				if (sessionsTimeStamps.get(userId) == null) {
					sessionsTimeStamps.put(userId, new HashMap<>());
				}

				// See if the session exist
				if (sessionsTimeStamps.get(userId).get(session) == null) {
					sessionsTimeStamps.get(userId).put(session, new ArrayList<>());
					sessionsTimeStamps.get(userId).get(session).add(new Tuple2<>(itemId, time));
				} else {
					sessionsTimeStamps.get(userId).get(session).add(new Tuple2<>(itemId, time));
				}
			});
			streamTr.close();
			PrintStream out = new PrintStream(dstFile);
			for (String user: sessionsTimeStamps.keySet()) {
				for (String session: sessionsTimeStamps.get(user).keySet()) {
					if (sessionsTimeStamps.get(user).get(session).size() >= minNumberPois) {
						Collections.sort(sessionsTimeStamps.get(user).get(session), new TimeStampComparatorTuple2());
						Long fstTime = sessionsTimeStamps.get(user).get(session).get(0).v2;
						Long lstTime = sessionsTimeStamps.get(user).get(session).get(sessionsTimeStamps.get(user).get(session).size() - 1).v2;
						double hours = lstTime - fstTime;
						hours /= 3600.0;
						if (hours >= minHours && hours <= maxHours) {
							for (Tuple2<String, Long> pref : sessionsTimeStamps.get(user).get(session)) {
								out.println(user + "\t" + pref.v1 + "\t" + newOld.get(pref.v1) + "\t" + pref.v2() + "\t" + session);
							}
						}						
					}
				}
				
			}		
			out.close();

		} catch (IOException e) {
			e.printStackTrace();
		}

	}
	
	public static void checkOrdering(String trainFile) {
		Stream<String> stream = null;
		try {
			stream = Files.lines(Paths.get(trainFile));
			Map<Long, List<Long>> userTimestamps = new HashMap<Long, List<Long>>(); 
			Set<Long> prevUsers = new HashSet<Long>();
			AtomicLong prevUser = new AtomicLong(-1);
			stream.forEach(line -> {
				String data [] = line.split("\t");
				Long user = Long.parseLong(data[0]);
				Long time = Long.parseLong(data[3]);
				
				if (prevUser.get() != user) {
					prevUser.set(user);
					if (prevUsers.contains(user)) {
						System.out.println("User has previusly be here");
					}
					prevUsers.add(user);
				}
				
				if (!userTimestamps.containsKey(user)) {
					userTimestamps.put(user, new ArrayList<>());
				}
				userTimestamps.get(user).add(time);
				
				
			});
			
			//Now, check
			boolean correct = true;
			for (Long user: userTimestamps.keySet()) {
				List<Long> timestamps = userTimestamps.get(user);
				if (timestamps.size() >= 2) {
					for (int i = 0; i < timestamps.size() - 1 ; i++) {
						
						if (timestamps.get(i) > timestamps.get(i + 1)) {
							correct = false;
							System.out.println("User " + user + " ratings are not ordered");
						}
						
					}
				}
				
			}
			if (correct) {
				System.out.println("For every user, all timestamps in the training set are ordered (from lower to higher)");
			}
			

			
		} catch (IOException e) {
			e.printStackTrace();
		}
		
	}
	
	public static void checkNoItemRecommendedMoreThanOnceForUser(String trainFile) {
		Stream<String> stream = null;
		try {
			stream = Files.lines(Paths.get(trainFile));
			Map<Long, Set<Long>> userItems = new HashMap<Long, Set<Long>>(); 
			
			stream.forEach(line -> {
				String data [] = line.split("\t");
				Long user = Long.parseLong(data[0]);
				Long item = Long.parseLong(data[1]);
				
				if (!userItems.containsKey(user)) {
					userItems.put(user, new HashSet<>());
				}
				
				if (userItems.get(user).contains(item)) {
					System.out.println("User " + user + " Item recommended more than once " + item);
				}
				userItems.get(user).add(item);
				
			});
			


			
		} catch (IOException e) {
			e.printStackTrace();
		}
		
	}

	
	
	
	/***
	 * The way to obtain the explicit data is Round(w * 4 / Max_u) + 1
	 * 
	 * @param trainFile       the original train file
	 * @param outputFile      the output file (parsed to explicit)
	 * @param charactersSplit the character to be used as splitting
	 */
	public static void implicitToExplicitMaxRatings(String trainFile, String outputFile, String charactersSplit) {
		Map<String, Integer> userMaximums = new HashMap<String, Integer>();

		int columnUser = 0;
		int columnItem = 1;
		int columnRatingImplicits = 2;
		boolean ignoreFirstLine = true;
		Stream<String> stream = null;
		try {
			if (ignoreFirstLine) {
				stream = Files.lines(Paths.get(trainFile)).skip(1);
			}
			else {
				stream = Files.lines(Paths.get(trainFile));
			}
			stream.forEach(line -> {
				String[] data = line.split(charactersSplit);
				int numberList = Integer.parseInt(data[columnRatingImplicits]);
				if (userMaximums.get(data[columnUser]) == null) {
					userMaximums.put(data[columnUser], numberList);
				}
				if (userMaximums.get(data[columnUser]) < numberList) {
					userMaximums.put(data[columnUser], numberList);
				}
			});
			stream.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

		// Now, reading again and print in the output file the round dividing by the
		// largest number of listenings
		try {
			if (ignoreFirstLine) {
				stream = Files.lines(Paths.get(trainFile)).skip(1);
			}
			else {
				stream = Files.lines(Paths.get(trainFile));
			}
			final PrintStream writer = new PrintStream(outputFile);
			stream.forEach(line -> {
				String[] data = line.split(charactersSplit);
				// The new value is the actual number of listenings * 4 / maximum number of
				// listenings. As this value is between 0 and 4, we add a 1
				long newValue = Math.round(
						(Double.parseDouble(data[columnRatingImplicits]) * 4.0) / userMaximums.get(data[columnUser]))
						+ 1;
				writer.println(data[columnUser] + "\t" + data[columnItem] + "\t" + newValue);
			});
			writer.close();
			stream.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

	}
	
	
	 /***
     * 
     * @param originalFile
     * @param outputFile
     */
    public static void filterSessions(String originalFile, String outputFile, int minNumberLengthSession, int minNumberSecondsSession, int minimumPreferencesUser, long minDiffforBots, int minPrefsToBeConsideredBot) {
    	//Original file format = user id, item id, preference, timestamp, sessionId
    	//We wont assume that the items in the sessions are ordered
    	
    	//Id user -> sessionId -> Tuple ItemId, preference, timestamp 
    	Map<Long, Map<Long, List<IdTimePref<Long>>>> map =  readSessionFormat(originalFile);
    	try {			
			PrintStream out = new PrintStream(outputFile);			
			for (Long userid: map.keySet()) {
				if (!isBot(getTimestampUser(map.get(userid)), minDiffforBots, minPrefsToBeConsideredBot)) {
					int prefsUser = countNumberPreferencesSessionFormat(map.get(userid));
					if (prefsUser > minimumPreferencesUser) {
						for (Long sessionId: map.get(userid).keySet()) {
							List<IdTimePref<Long>> tuple3lst = map.get(userid).get(sessionId);
							Collections.sort(tuple3lst, PreferenceComparators.timeComparatorIdTimePref);
							double diff = Math.abs(tuple3lst.get(0).v3 - tuple3lst.get(tuple3lst.size() - 1).v3);
							if (diff >= minNumberSecondsSession && tuple3lst.size() >= minNumberLengthSession) {
								for (IdTimePref<Long> t: tuple3lst) {
									out.println(userid + "\t" + t.v1 + "\t" + t.v2 + "\t" + t.v3 + "\t" + sessionId);
								}
							} else {
								System.out.println("User: "+ userid + " with session: " + sessionId + " discarded because the sessions length is "
										+ "lower than " + minNumberSecondsSession + "or because the session length is lower than " + minNumberLengthSession);
							}
							
						}
					}
				}
			}
			out.close();
			
			
			
    	} catch (Exception e) {
			e.printStackTrace();

		}
			
    }
    
    private static List<Long> getTimestampUser(Map<Long, List<IdTimePref<Long>>> userSessions){
    	List<Long> userTimestampsOrdered = new ArrayList<>();
    	for (Long sessionId: userSessions.keySet()) {
    		userTimestampsOrdered.addAll(userSessions.get(sessionId).stream().map(pref -> pref.v3).collect(Collectors.toList()));
    	}
    	Collections.sort(userTimestampsOrdered);
    	return userTimestampsOrdered;
    }
    
    /***
     * Method to parse a test data from ranksys selecting just the first item in the test set
     * @param ranksysTestData testData from ranksys 
     * @param outputFile the output test file containing just the first item in the test set
     */
    public static <U, I> void generateTestOnlyFirstItem(FastTemporalPreferenceDataIF<U, I> ranksysTestData, String outputFile) {
    	try {
			PrintStream out = new PrintStream(outputFile);
	    	ranksysTestData.getUsersWithPreferences().forEach(u -> {
	    		IdxTimePref pref = ranksysTestData.getUidxTimePreferences(ranksysTestData.user2uidx(u)).sorted(PreferenceComparators.timeComparatorIdxTimePref).findFirst().get();
	    		out.println(u + "\t" + ranksysTestData.iidx2item(pref.v1) + "\t" + pref.v2 + "\t" + pref.v3);
	    	});
	    	out.close();
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

    	
    }
	
	private static int countNumberPreferencesSessionFormat(Map<Long, List<IdTimePref<Long>>> prefsUser) {
		int acc = 0;
		for (Long sessionId: prefsUser.keySet()) {
			acc += prefsUser.get(sessionId).size();
		}
		return acc;
	}

	static class Preference implements Comparable<Preference> {

		Long idUser;
		Long idItem;
		Double rating;
		Long timeStamp;

		public Preference(Long idUser, Long idItem, Double rating, Long timeStamp) {
			super();
			this.idUser = idUser;
			this.idItem = idItem;
			this.rating = rating;
			this.timeStamp = timeStamp;
		}

		public Long getIdUser() {
			return idUser;
		}

		public void setIdUser(Long idUser) {
			this.idUser = idUser;
		}

		public Long getIdItem() {
			return idItem;
		}

		public void setIdItem(Long idItem) {
			this.idItem = idItem;
		}

		public Double getRating() {
			return rating;
		}

		public void setRating(Double rating) {
			this.rating = rating;
		}

		public Long getTimeStamp() {
			return timeStamp;
		}

		public void setTimeStamp(Long timeStamp) {
			this.timeStamp = timeStamp;
		}

		@Override
		public int compareTo(Preference arg0) {
			return this.getTimeStamp().compareTo(arg0.getTimeStamp());
		}

	}
	
	public static void filterRecByTraining(FastTemporalPreferenceDataIF<Long, Long> ranksysTrainTemporal, String recommendedFile, int limit, String outputFile) {
		
		RecommendationFormat<Long, Long> format = new SimpleRecommendationFormat<>(lp, lp);
		
		try {
			PrintStream out = new PrintStream(outputFile);

			format.getReader(recommendedFile).readAll().forEach(rec -> {
				Long u = rec.getUser();
				
				Set<Long> itemsRatedByUser = ranksysTrainTemporal.getUserPreferences(u).map(t -> t.v1).collect(Collectors.toSet());
				
				List<Tuple2od<Long>> recommendation = rec.getItems();
				
				int c = 0;
				for (Tuple2od<Long> t: recommendation) {
					if (!itemsRatedByUser.contains(t.v1)) {
						out.println(u + "\t" + t.v1 + "\t" + t.v2);
					}
					c++;
					if  (c >= limit) {
						break;
					}
				}
				
			});
			out.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	public static <U, I> void filterTestByTraining(FastTemporalPreferenceDataIF<U, I> ranksysTrainTemporal,
			FastTemporalPreferenceDataIF<U, I> ranksysTestDataTemporal, boolean removeDuplicates, String outputFile) {
		try {
			Map<U, Set<I>> alreadyRatedItemsTest = new HashMap<>();
			PrintStream out = new PrintStream(outputFile);
			
			ranksysTestDataTemporal.getUsersWithPreferences().forEach(u -> {
				alreadyRatedItemsTest.put(u, new HashSet<I>());
				Set<I> itemsRatedByUserTraining = new HashSet<I>();
				
				if (ranksysTrainTemporal.containsUser(u)) {
					itemsRatedByUserTraining = ranksysTrainTemporal.getUserPreferences(u).map(t -> t.v1).collect(Collectors.toSet());
				}
				
				Set<I> itemsRatedByUserTrainingf = itemsRatedByUserTraining;
				
				ranksysTestDataTemporal.getUserPreferences(u).forEach(t -> {
					I itemTest = t.v1;
					
					if (!itemsRatedByUserTrainingf.contains(itemTest)) {
						//Item in test not in training, so print
						
						Set<I> processedItemsTest = alreadyRatedItemsTest.get(u); 
						if (!processedItemsTest.contains(itemTest) || !removeDuplicates) {
							out.println(u + "\t" + t.v1 + "\t" + t.v2 + "\t" + t.v3);
							processedItemsTest.add(t.v1);
						}
						
					}
				});
				
				
			});
			
			out.close();
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		
	}






	



}

class PreferenceByUserComparator implements Comparator<Preference>{

	@Override
	public int compare(Preference o1, Preference o2) {
		int a = o1.idUser.compareTo(o2.idUser);
		if (a == 0) {
			return o1.timeStamp.compareTo(o2.timeStamp);
		}
		return a;
	}
	
}

class WeightComparator implements Comparator<Tuple2<String, Double>> {

	@Override
	public int compare(Tuple2<String, Double> o1, Tuple2<String, Double> o2) {
		int a = o1.v2.compareTo(o2.v2);
		if (a == 0) {
			return o1.v1.compareTo(o2.v1);
		}
		return a;
	}

}

class TimeStampComparatorTuple2 implements Comparator<Tuple2<String, Long>> {
	@Override
	public int compare(Tuple2<String, Long> o1, Tuple2<String, Long> o2) {
		if (o1.v2 < o2.v2)
			return -1;
		 else
			return 1;
	}

}

class TimeStampComparator implements Comparator<Tuple3<Long, Double, Long>> {
	@Override
	public int compare(Tuple3<Long, Double, Long> o1, Tuple3<Long, Double, Long> o2) {
		if (o1.v3 < o2.v3)
			return -1;
		// Deterministic
		else if (o1.v3.equals(o2.v3)){
			if (o1.v1 < o2.v1)
				return -1;
			else if (o1.v1 > o2.v1)
				return 1;
			else
				return 0;
		} else
			return 1;
	}

}
