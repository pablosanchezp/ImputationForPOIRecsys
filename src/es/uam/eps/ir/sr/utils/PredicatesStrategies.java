package es.uam.eps.ir.sr.utils;

import static org.ranksys.formats.parsing.Parsers.lp;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintStream;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.function.Predicate;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import org.jooq.lambda.tuple.Tuple2;
import org.ranksys.core.util.tuples.Tuple2od;
import org.ranksys.formats.rec.RecommendationFormat;
import org.ranksys.formats.rec.SimpleRecommendationFormat;

import es.uam.eps.ir.crossdomainPOI.datamodel.temporal.interfaces.FastTemporalPreferenceDataIF;
import es.uam.eps.ir.crossdomainPOI.datamodel.temporal.interfaces.TemporalPreferenceDataIF;
import es.uam.eps.ir.crossdomainPOI.datamodel.temporal.preferences.IdxTimePref;
import es.uam.eps.ir.crossdomainPOI.utils.UsersMidPoints;
import es.uam.eps.ir.ranksys.core.Recommendation;
import es.uam.eps.ir.ranksys.core.preference.PreferenceData;
import es.uam.eps.ir.ranksys.fast.preference.FastPreferenceData;
import es.uam.eps.ir.ranksys.rec.Recommender;
import es.uam.eps.ir.seqawareev.tour.interfaces.TourRecommenderIF;
import es.uam.eps.ir.sr.data.POIProcessData;
import es.uam.eps.ir.sr.utils.comparators.PreferenceComparators;

/****
 * File that will contains all the predicates that may work as recommendation
 * strategies of the system.
 * 
 * @author Pablo Sanchez (psperez@icai.comillas.edu)
 *
 */
public final class PredicatesStrategies {

	/**
	 * Enumeration of the predicates: -NO_CONDITION: no conditions of the items
	 * recommended (can be rated by the user in train)
	 * 
	 * -TRAIN_ITEMS: every item recommended by the user must have a rating in train
	 * data and must have not been rated in the train subset by the user
	 * 
	 * -ALL_ITEMS: every item in the system is considered or the recommendation, but
	 * must not have been rated by the user in train
	 * 
	 * -ONLYITEMS_TEST: only make recommendations for the items that are in the user
	 * test
	 * 
	 * -ITEMS_IN_TRAIN: every item recommended by the user must appear in the
	 * training set
	 * 
	 * @author Pablo Sanchez (pablo.sanchezp@uam.es)
	 *
	 */
	public enum recommendationStrategy {
		NO_CONDITION("N_C"), TRAIN_ITEMS("T_I"), ALL_ITEMS("A_I"), ONLYITEMS_TEST("OI_T"), ITEMS_IN_TRAIN("I_I_T"), MIDPOINT_USER_TRAIN_ITEMS("M_U"), CLOSE_LAST_ITEM_TRAIN_ITEMS("C_L_I");
		
		private String shortName;
		
		private recommendationStrategy(String name) {
			this.shortName = name;
		}
		
		public String getShortName() {
			return shortName;
		}
	};
	
	/***
	 * Selection of the recommendation predicate
	 * 
	 * @param strat     the strategy of the recommendation
	 * @param user      the user
	 * @param trainData the train data
	 * @param testData  the test data
	 * @return a predicate for recommendation
	 */
	public static <U, I> Predicate<I> selectGeographicalPredicate(PreferenceData<U, I> trainData, recommendationStrategy strat, U user, Map<I, Tuple2<Double, Double>> coordinatesPOI, UsersMidPoints<U, I> usermid, I item, I prevItem, double limitKm) {
		switch (strat) {
		case MIDPOINT_USER_TRAIN_ITEMS:
			return isInRangeMidPointUser(coordinatesPOI, usermid, user, item, limitKm).and(trainItems(user, trainData));
		case CLOSE_LAST_ITEM_TRAIN_ITEMS:
			return isCloseLastItem(coordinatesPOI, item, prevItem, limitKm).and(trainItems(user, trainData));
		default:
			return p -> true;

		}

	}
	
	
	
	/***
	 * Predicate to test if the recommended venue is withtin the range of the user midpoint
	 * @param <U>
	 * @param <I>
	 * @param coordinatesPOI
	 * @param usermid
	 * @param user
	 * @param item
	 * @param limitKm
	 * @return
	 */
	public static <U, I> Predicate<I> isInRangeMidPointUser(Map<I, Tuple2<Double, Double>> coordinatesPOI, UsersMidPoints<U, I> usermid, U user, I item, double limitKm){
		Tuple2<Double, Double> coordPOI = coordinatesPOI.get(item);
		Tuple2<Double, Double> coordUser = usermid.getUserAvgCoordinates(user);
		return p -> POIProcessData.haversine(coordPOI.v1, coordPOI.v2, coordUser.v1, coordUser.v2, false) <= limitKm;
	}
	
	/***
	 * 
	 * @param <U>
	 * @param <I>
	 * @param coordinatesPOI
	 * @param usermid
	 * @param user
	 * @param item
	 * @param limitKm
	 * @return
	 */
	public static <U, I> Predicate<I> isCloseLastItem(Map<I, Tuple2<Double, Double>> coordinatesPOI, I nextItem, I lastRecItem, double limitKm){
		//For the first item recommended then we allow the recommendation, as the lastRecItem is null
		if (lastRecItem == null) {
			return p -> true;
		}
		Tuple2<Double, Double> coordPOI1 = coordinatesPOI.get(nextItem);
		Tuple2<Double, Double> coordPOI2 = coordinatesPOI.get(lastRecItem);
		return p -> POIProcessData.haversine(coordPOI1.v1, coordPOI1.v2, coordPOI2.v1, coordPOI2.v2, false) <= limitKm;
	}
	
	

	/**
	 * Predicate used in RankSys to make recommendations of items not rated by the
	 * user in train
	 *
	 * @param user      the user
	 * @param trainData the ranksysDatamodel of train
	 * @return a predicate
	 */
	public static <U, I> Predicate<I> isNotInTrain(U user, PreferenceData<U, I> trainData) {
		return p -> trainData.getUserPreferences(user).noneMatch(itemPref -> itemPref.v1.equals(p));
	}
	
	
	public static <U, I> Predicate<I> isNotInTrain(U user, TemporalPreferenceDataIF<U, I> trainData) {
		return p -> trainData.getUserPreferences(user).noneMatch(itemPref -> itemPref.v1.equals(p));
	}


	/**
	 * Predicate so see if an item is a candidate
	 * 
	 * @param candidates
	 * @return a predicate
	 */
	public static <I> Predicate<I> isCandidate(Set<I> candidates) {
		return p -> candidates.contains(p);
	}

	/**
	 * Predicate used in RankSys to make recommendations of items that appears in
	 * the test set
	 *
	 * @param user     the user
	 * @param testData the ranksysDatamodel of test
	 * @return a predicate
	 */
	public static <U, I> Predicate<I> isInTest(U user, PreferenceData<U, I> testData) {
		return p -> testData.getUserPreferences(user).anyMatch(itemPref -> itemPref.v1.equals(p));
	}
	
	public static <U, I> Predicate<I> isInTest(U user, TemporalPreferenceDataIF<U, I> testData) {
		return p -> testData.getUserPreferences(user).anyMatch(itemPref -> itemPref.v1.equals(p));
	}


	/***
	 * Predicate that will return true if the item has been rated in the training
	 * set (by any user)
	 * 
	 * @param trainData the trainData
	 * @return a predicate
	 */
	public static <U, I> Predicate<I> itemWithRatingInTrain(PreferenceData<U, I> trainData) {
		return p -> trainData.getItemPreferences(p).findFirst().isPresent();
	}
	public static <U, I> Predicate<I> itemWithRatingInTrain(TemporalPreferenceDataIF<U, I> trainData) {
		return p -> trainData.getItemPreferences(p).findFirst().isPresent();
	}
	
	/***
	 * Train items strategy of the predicate. -Only items with at least one rating
	 * in train should be recommended. -Only items that the user have not rated in
	 * train (new items for the user).
	 * 
	 * @param user      the user
	 * @param trainData the trainData
	 * @return a predicate
	 */
	public static <U, I> Predicate<I> trainItems(U user, PreferenceData<U, I> trainData) {
		return isNotInTrain(user, trainData).and(itemWithRatingInTrain(trainData));
	}
	public static <U, I> Predicate<I> trainItems(U user, TemporalPreferenceDataIF<U, I> trainData) {
		return isNotInTrain(user, trainData).and(itemWithRatingInTrain(trainData));
	}

	/***
	 * No condition predicate (every item is considered as a predicate)
	 * 
	 * @return a predicate
	 */
	public static <I> Predicate<I> noCondition() {
		return p -> true;
	}

	/**
	 * * Method that receiving an user, the train data and a stream of
	 * candidateItems, will return a stream of items that are not in the train of
	 * the user but that have at least 1 preference in the train subset and is in
	 * the list of candidate items
	 *
	 * @param user          the user to whom we are making recommendations
	 * @param trainData     the train data
	 * @param candidateItem the candidate items
	 * @return a stream of candidates items
	 */
	public static <U, I> Stream<I> generateCandidatesWithRatingsAndNotRatedByUserInTrain(U user,
			PreferenceData<U, I> trainData, Set<I> candidateItems) {
		return trainData.getItemsWithPreferences()
				.filter(i -> !trainData.getUserPreferences(user).map(pref -> pref.v1).anyMatch(i2 -> i2.equals(i)))
				.filter(i -> candidateItems.contains(i));

	}

	/**
	 * * Method to store a ranking file of a RankSys recommender
	 *
	 * @param ranksysTrainData     the train data
	 * @param ranksystestData      the test data
	 * @param rankSysrec           the recommender
	 * @param outputFile           the output file
	 * @param numberItemsRecommend the number of items to recommend
	 */
	public static <U, I> void ranksysWriteRanking(PreferenceData<U, I> ranksysTrainData, PreferenceData<U, I> ranksysTrainData2,
			PreferenceData<U, I> ranksystestData, Recommender<U, I> rankSysrec, String outputFile,
			int numberItemsRecommend, recommendationStrategy strat) {
		//Train data is the data used to train the recommenders
		//TrainData2 will be the one used for the candidates (normally trainData 1 and 2 will be the same), but for example, for data imputation
		//TrainData2 will be the original TrainFile
		
		PreferenceData<U, I> ranksysTrainDataUsed = ranksysTrainData2 != null ? ranksysTrainData2: ranksysTrainData;
		PrintStream out;
		try {
			out = new PrintStream(outputFile);
			Stream<U> targetUsers = ranksystestData.getUsersWithPreferences();
			System.out.println("Users in test data " + ranksystestData.numUsersWithPreferences());
			System.out.println("Users in test data appearing in train set "
					+ ranksystestData.getUsersWithPreferences().filter(u -> ranksysTrainDataUsed.containsUser(u)).count());
			if (strat.equals(recommendationStrategy.ONLYITEMS_TEST)) {
				numberItemsRecommend = Integer.MAX_VALUE;
			}
			final int numItemsRec = numberItemsRecommend;
			targetUsers.forEach(user -> {
				if (ranksystestData.getUserPreferences(user).count() != 0L && ranksysTrainDataUsed.containsUser(user)) {
					int rank = 1;
					/*
					System.out.println("USER" + user);
					System.out.println("numItemsRec" + numItemsRec);
					System.out.println("ranksysTrainDataUsed" + ranksysTrainDataUsed);
					System.out.println("ranksystestData" + ranksystestData);*/
					Recommendation<U, I> rec = rankSysrec.getRecommendation(user, numItemsRec,
							selectRecommendationPredicate(strat, user, ranksysTrainDataUsed, ranksystestData));
					for (Tuple2od<I> tup : rec.getItems()) { // Items recommended
						// We should see if the items are in train
						out.println(SequentialRecommendersUtils.formatRank(user, tup.v1, tup.v2, rank));
						rank++;
					}
				}
			});
			out.close();
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}
	
	/***
	 * Method to write the recommendation of a tour recommender in ranksys. For every recommendation, we send the first item rated
	 * in the test set, to the recommendation method in the recommender
	 * @param rankSysTrainData
	 * @param rankSysTestData
	 * @param rankSysrec
	 * @param outputFile
	 * @param numberItemsRecommend
	 * @param strat
	 */
	public static <U, I> void ranksysWriteRankingTourRecommender(FastTemporalPreferenceDataIF<U, I> rankSysTrainData, FastTemporalPreferenceDataIF<U, I> rankSysTestData, TourRecommenderIF<U, I> rankSysrec, String outputFile,
			int numberItemsRecommend, recommendationStrategy strat) {
		PrintStream out;
		try {
			out = new PrintStream(outputFile);
			Stream<U> targetUsers = rankSysTestData.getUsersWithPreferences();
			System.out.println("Users in test data " + rankSysTestData.numUsersWithPreferences());
			System.out.println("Users in test data appearing in train set "
					+ rankSysTestData.getUsersWithPreferences().filter(u -> rankSysTrainData.containsUser(u)).count());
			
			if (strat.equals(recommendationStrategy.ONLYITEMS_TEST)) {
				numberItemsRecommend = Integer.MAX_VALUE;
			}
			
			final int numItemsRec = numberItemsRecommend;
			targetUsers.forEach(user -> {
				if (rankSysTestData.getUserPreferences(user).count() != 0L && rankSysTrainData.containsUser(user)) {
					int rank = 1;
					//Get the first element from test
					IdxTimePref pref = rankSysTestData.getUidxTimePreferences(rankSysTestData.user2uidx(user)).sorted(PreferenceComparators.timeComparatorIdxTimePref).findFirst().get();
					Recommendation<U, I> rec = rankSysrec.getRecommendation(user, numItemsRec, selectRecommendationPredicateTour(strat, user, rankSysTrainData, rankSysTestData), pref.v3, rankSysTestData.iidx2item(pref.v1));
					
					for (Tuple2od<I> tup : rec.getItems()) { // Items recommended
						// We should see if the items are in traingetRecommendation
						out.println(SequentialRecommendersUtils.formatRank(user, tup.v1, tup.v2, rank));
						rank++;
					}
				}
			});
			out.close();

		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
	}
	/***
	 * Method to generate another recommendation file so that the new recommendation file have the first item of the test set as the first item recommended 
	 * 
	 * @param rankSysTestData the test data
	 * @param originRecommendation the original recommendation file
	 * @param destNewRecommendation the destination of the new recommendation file
	 */
	public static void parseRecWithFirstItemInTest(FastTemporalPreferenceDataIF<Long, Long> rankSysTestData,
			String originRecommendation, String destNewRecommendation) {
		try {
			RecommendationFormat<Long, Long> format = new SimpleRecommendationFormat<>(lp, lp);
			PrintStream out = new PrintStream(destNewRecommendation);
			format.getReader(originRecommendation).readAll().forEach(rec -> {
				Long u = rec.getUser();
				Long fstItemTest = rankSysTestData.getUserPreferences(u)
						.sorted(PreferenceComparators.timeComparatorIdTimePref).map(pref -> pref.v1).findFirst().get();
				Set<Long> addedItems = new HashSet<>();
				List<Tuple2od<Long>> lstRec = rec.getItems();

				int totalCount = lstRec.size();
				addedItems.add(fstItemTest);
				out.println(u + "\t" + fstItemTest + "\t" + totalCount);
				totalCount--;
				for (int i = 0; i < lstRec.size(); i++) {
					Long itemRec = lstRec.get(i).v1;
					if (!addedItems.contains(itemRec)) {
						addedItems.add(itemRec);
						int c = totalCount - i;
						out.println(u + "\t" + itemRec + "\t" + c);
					}
				}

			});
			out.close();
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}


	/**
	 * Method to create a recommender file in which every item recommended must be
	 * in the set of candidates
	 * 
	 * @param ranksysTrainData     the train data
	 * @param ranksystestData      the test data
	 * @param rankSysrec           the recommended
	 * @param outputFile           the output file
	 * @param numberItemsRecommend the number of items recommended
	 * @param candidates           the set of candidates
	 * @param strat                the strategy to apply
	 */
	public static <U, I> void ranksysWriteRankingFromCandidatesList(FastPreferenceData<U, I> ranksysTrainData,
			FastPreferenceData<U, I> ranksystestData, Recommender<U, I> rankSysrec, String outputFile,
			int numberItemsRecommend, Set<I> candidates, recommendationStrategy strat) {
		PrintStream out;
		try {
			out = new PrintStream(outputFile);
			Stream<U> targetUsers = ranksystestData.getUsersWithPreferences();
			System.out.println("Users in test data " + ranksystestData.numUsersWithPreferences());
			targetUsers.forEach(user -> {
				if (ranksystestData.getUserPreferences(user).count() != 0L && ranksysTrainData.containsUser(user)) {
					int rank = 1;
					Recommendation<U, I> rec = rankSysrec.getRecommendation(user,
							selectRecommendationPredicate(strat, user, ranksysTrainData, ranksystestData, candidates));
					if (rec == null || rec.getItems() == null) {
						System.out.println("Recommendation for user " + user + " is null");
					} else {
						for (Tuple2od<I> tup : rec.getItems()) { // Items recommended
							if (rank > numberItemsRecommend) {
								break;
							}
							out.println(SequentialRecommendersUtils.formatRank(user, tup.v1, tup.v2, rank));
							rank++;
						}
					}
				}
			});
			out.close();
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}

	/***
	 * Selection of the recommendation predicate
	 * 
	 * @param strat     the strategy of the recommendation
	 * @param user      the user
	 * @param trainData the train data
	 * @param testData  the test data
	 * @return a predicate for recommendation
	 */
	public static <U, I> Predicate<I> selectRecommendationPredicate(recommendationStrategy strat, U user,
			PreferenceData<U, I> trainData, PreferenceData<U, I> testData) {
		switch (strat) {
		case ALL_ITEMS:
			return isNotInTrain(user, trainData);
		case NO_CONDITION:
			return noCondition();
		case ONLYITEMS_TEST:
			return isInTest(user, testData);
		case ITEMS_IN_TRAIN:
			return itemWithRatingInTrain(trainData);
		case TRAIN_ITEMS:
		default:
			return trainItems(user, trainData);
		}

	}
	
	public static <U, I> Predicate<I> selectRecommendationPredicateTour(recommendationStrategy strat, U user,
			TemporalPreferenceDataIF<U, I> trainData, TemporalPreferenceDataIF<U, I> testData) {
		switch (strat) {
		case ALL_ITEMS:
			return isNotInTrain(user, trainData);
		case NO_CONDITION:
			return noCondition();
		case ONLYITEMS_TEST:
			return isInTest(user, testData);
		case ITEMS_IN_TRAIN:
			return itemWithRatingInTrain(trainData);
		case TRAIN_ITEMS:
		default:
			return trainItems(user, trainData);
		}

	}
	
	
	

	/**
	 * Selection of recommendation predicate combining the set of candidate files
	 * 
	 * @param strat     the strategy of the recommendation
	 * @param user      the user
	 * @param trainData the trainData
	 * @param testData  the testData
	 * @return the predicate
	 */
	public static <U, I> Predicate<I> selectRecommendationPredicate(recommendationStrategy strat, U user,
			PreferenceData<U, I> trainData, PreferenceData<U, I> testData, Set<I> candidates) {
		switch (strat) {
		case ALL_ITEMS:
			return isNotInTrain(user, trainData).and(isCandidate(candidates));
		case NO_CONDITION:
			return isCandidate(candidates);
		case ONLYITEMS_TEST:
			return isInTest(user, testData).and(isCandidate(candidates));
		case ITEMS_IN_TRAIN:
			return itemWithRatingInTrain(trainData);
		case TRAIN_ITEMS:
		default:
			return trainItems(user, trainData).and(isCandidate(candidates));
		}

	}

	public static <U, I> List<I> generateCandidateList(recommendationStrategy strat, U user, PreferenceData<U, I> trainData, PreferenceData<U, I> testData) {
		switch (strat) {
		case ALL_ITEMS:
			return trainData.getAllItems().filter(isNotInTrain(user, trainData)).collect(Collectors.toList());
		case NO_CONDITION:
			return trainData.getAllItems().collect(Collectors.toList());
		case ONLYITEMS_TEST:
			return testData.getAllItems().filter(isInTest(user, testData)).collect(Collectors.toList());
		case ITEMS_IN_TRAIN:
			return trainData.getAllItems().filter(itemWithRatingInTrain(trainData)).collect(Collectors.toList());
		case TRAIN_ITEMS:
		default:
			return trainData.getAllItems().filter(trainItems(user, trainData)).collect(Collectors.toList());
		}
	}
	
	

}
