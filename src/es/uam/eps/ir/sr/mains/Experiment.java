package es.uam.eps.ir.sr.mains;

import static org.ranksys.formats.parsing.Parsers.lp;
import static org.ranksys.formats.parsing.Parsers.sp;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.jooq.lambda.tuple.Tuple2;
import org.jooq.lambda.tuple.Tuple4;
import org.ranksys.formats.feature.SimpleFeaturesReader;
import org.ranksys.formats.preference.SimpleRatingPreferencesReader;
import org.ranksys.formats.rec.RecommendationFormat;
import org.ranksys.formats.rec.SimpleRecommendationFormat;

import es.uam.eps.ir.antimetrics.antirel.BinaryAntiRelevanceModel;
import es.uam.eps.ir.attrrec.datamodel.feature.SimpleUserFeatureData;
import es.uam.eps.ir.attrrec.datamodel.feature.SimpleUserFeaturesReader;
import es.uam.eps.ir.attrrec.datamodel.feature.UserFeatureData;
import es.uam.eps.ir.attrrec.metrics.recommendation.averages.WeightedAverageRecommendationMetricIgnoreNoRelevantUsersAndNaNs;
import es.uam.eps.ir.attrrec.metrics.recommendation.averages.WeightedModelUser;
import es.uam.eps.ir.attrrec.metrics.recommendation.averages.WeightedModelUser.UserMetricWeight;
import es.uam.eps.ir.attrrec.sim.ItemIndexReadingSimilarity;
import es.uam.eps.ir.attrrec.sim.ItemSim;
import es.uam.eps.ir.crossdomainPOI.datamodel.temporal.interfaces.FastTemporalPreferenceDataIF;
import es.uam.eps.ir.crossdomainPOI.metrics.system.RealAggregateDiversity;
import es.uam.eps.ir.crossdomainPOI.metrics.system.UserCoverage;
import es.uam.eps.ir.crossdomainPOI.utils.UsersMidPoints;
import es.uam.eps.ir.crossdomainPOI.utils.UsersMidPoints.SCORES_FREQUENCY;
import es.uam.eps.ir.ranksys.core.feature.FeatureData;
import es.uam.eps.ir.ranksys.core.feature.SimpleFeatureData;
import es.uam.eps.ir.ranksys.core.preference.ConcatPreferenceData;
import es.uam.eps.ir.ranksys.core.preference.PreferenceData;
import es.uam.eps.ir.ranksys.core.preference.SimplePreferenceData;
import es.uam.eps.ir.ranksys.diversity.intentaware.FeatureIntentModel;
import es.uam.eps.ir.ranksys.diversity.intentaware.IntentModel;
import es.uam.eps.ir.ranksys.diversity.sales.metrics.AggregateDiversityMetric;
import es.uam.eps.ir.ranksys.diversity.sales.metrics.GiniIndex;
import es.uam.eps.ir.ranksys.fast.index.FastUserIndex;
import es.uam.eps.ir.ranksys.fast.index.SimpleFastUserIndex;
import es.uam.eps.ir.ranksys.fast.preference.FastPreferenceData;
import es.uam.eps.ir.ranksys.metrics.RecommendationMetric;
import es.uam.eps.ir.ranksys.metrics.SystemMetric;
import es.uam.eps.ir.ranksys.metrics.rank.NoDiscountModel;
import es.uam.eps.ir.ranksys.metrics.rank.RankingDiscountModel;
import es.uam.eps.ir.ranksys.metrics.rel.BinaryRelevanceModel;
import es.uam.eps.ir.ranksys.metrics.rel.NoRelevanceModel;
import es.uam.eps.ir.ranksys.metrics.rel.RelevanceModel;
import es.uam.eps.ir.ranksys.nn.item.sim.ItemSimilarity;
import es.uam.eps.ir.ranksys.novdiv.distance.CosineFeatureItemDistanceModel;
import es.uam.eps.ir.ranksys.novdiv.distance.ItemDistanceModel;
import es.uam.eps.ir.ranksys.novelty.temporal.TimestampCalculator;
import es.uam.eps.ir.ranksys.rec.Recommender;
import es.uam.eps.ir.seqawareev.utils.CBBaselinesUtils.ITEMTRANSFORMATION;
import es.uam.eps.ir.seqawareev.utils.CBBaselinesUtils.USERTRANSFORMATION;
import es.uam.eps.ir.seqawareev.wrappers.SimpleFastTemporalFeaturePreferenceDataWrapper;
import es.uam.eps.ir.sr.data.POIProcessData;
import es.uam.eps.ir.sr.data.ProcessData;
import es.uam.eps.ir.sr.data.imputation.ImputedSimpleFastPreferenceDataAvgPreferencesUsers;
import es.uam.eps.ir.sr.data.imputation.ImputedSimpleFastPreferenceDataRecommender;
import es.uam.eps.ir.sr.data.imputation.ImputedSimpleFastPreferenceDataTestUsers;
import es.uam.eps.ir.sr.data.imputation.ImputedSimpleFastPreferenceDataTrainUsers;
import es.uam.eps.ir.sr.metrics.system.GiniRelevantUsers;
import es.uam.eps.ir.sr.metrics.system.RealAggregateDiversityRelevantUsers;
import es.uam.eps.ir.sr.metrics.system.RelevantUsersRecommended;
import es.uam.eps.ir.sr.metrics.system.UserCoverageRelevant;
import es.uam.eps.ir.sr.utils.PredicatesStrategies;
import es.uam.eps.ir.sr.utils.SequentialRecommendersUtils;



/***
 * Main class to launch the different recommenders / evaluation / data processing
 * @author Pablo Sanchez (psperez@icai.comillas.edu)
 *
 */
public class Experiment {
	/********************************/
	/**All options of the arguments**/
	/********************************/
	//
	/**Option of the case**/
	private static final String OPT_CASE = "option";

	/**Train and test files**/
	private static final String OPT_TRAIN_FILE = "train file";
	private static final String OPT_FULL_FILE = "full file";
	private static final String OPT_TRAIN_FILE_2 = "train file 2";
	private static final String OPT_TEST_FILE = "test file";
	private static final String OPT_TEST_FILE_2 = "testFile 2";


	/**For LCS**/
	private static final String OPT_LCS_NORMALIZATION = "LCS normalization";
	private static final String OPT_LCS_TRAIN_MODEL = "LCS trainModel";
	// Factor to multiply rating or ID for LCS
	private static final String OPT_IDENTIFICATION_FACTOR = "identification factor";
	private static final String OPT_RATING_FACTOR = "rating factor";
	
	/**For recommendation (general recommendation)**/
	
	/*Generic*/
	//Complete indexes of training and test sets (e.g. considering items and users in the test set)
	private static final String OPT_COMPLETE_INDEXES = "complete indexes";
	private static final String OPT_ITEMS_RECOMMENDED = "items recommended";
	private static final String OPT_RECOMMENDATION_STRATEGY = "recommendation strategy";
	private static final String OPT_RANKSYS_REC = "ranksysRecommender";	
	private static final String OPT_RANKSYS_REC2 = "ranksysRecommender2";
	
	//Some parameters
	private static final String OPT_LAST_N_SESSIONS = "last n sessions for the user";
	private static final String OPT_NEIGH = "neighbours";
	private static final String OPT_INVERSE = "inverse";
	private static final String OPT_NEGANTI = "negative anti";


	private static final String OPT_RANKSYS_SIM = "ranksys similarity";
	private static final String OPT_RANKSYS_SIM_2 = "ranksys similarity 2";
	private static final String OPT_SIMFILE = "similarity file";
	private static final String OPT_SIMFILE_2 = "similarity file 2";
	private static final String OPT_AGGR_TRAJ_SIM = "aggregation trajectory sim";
	private static final String OPT_SOFTEN_SIM_TRA = "soften similarity";
	
	private static final String OPT_COMPOSED_SIM = "composed similarity";
	
	private static final String OPT_MARKOV_CHAIN_ORDER = "Markov chain order";
	private static final String OPT_CONFIDENCE = "confidencesimilarity";
	private static final String OPT_CONFIDENCE_ESTIMATION = "confidence stimation";
	private static final String OPT_FILTER = "apply a filter";
	private static final String OPT_PREFERENCE_RATING = "preference rating";
	private static final String OPT_SOC_INFLUENCE_FILE = "social influence file";
	private static final String OPT_TEMPORAL_WEIGHT = "temporal weight used in temporal recommenders";
	private static final String OPT_LAMBDA = "lambda variable";
	private static final String OPT_LAMBDA_2 = "lambda variable (2)";
	private static final String OPT_NORMALIZE = "normalize";
	private static final String OPT_TEMPORAL_LAMBDA = "temporal lambda";
	
	//Hybrid
	private static final String OPT_HYBRID_WEIGHTS = "HybridWeights";

	
	// Indexes of BF recommender
	private static final String OPT_INDEX_BACKWARDS = "index backwards";
	private static final String OPT_INDEX_FORWARDS = "index forwards";
	private static final String OPT_IN_INDEX_FILE = "input indexes file";
	
	// Parameters for Ranksys Recommender MatrixFactorization
	private static final String OPT_K_FACTORIZER = "k factorizer value";
	private static final String OPT_ALPHA_FACTORIZER = "alpha factorizer value";
	private static final String OPT_ALPHA_FACTORIZER2 = "alpha factorizer (2) value";
	private static final String OPT_THETA_FACTORIZER = "theta factorizer value";
	private static final String OPT_USE_SIGMOID = "useSigmoid";

	


	private static final String OPT_LAMBDA_FACTORIZER = "lambda factorizer value";
	private static final String OPT_NUM_INTERACTIONS = "num interactions value";
	private static final String OPT_ETA = "eta"; 		


	// Parameters of Aggregation library
	private static final String OPT_NORM_AGGREGATE_LIBRARY = "norm of aggregate library";
	private static final String OPT_COMB_AGGREGATE_LIBRARY = "comb of aggregate library";
	private static final String OPT_WEIGHT_AGGREGATE_LIBRARY = "weight of aggregate library";
	
	// Regularization for SLIM recommender
	private static final String OPT_REGULARIZATION_L1 = "Regularization 1";
	private static final String OPT_REGULARIZATION_L2 = "Regularization 2";
	
	// Parameters for SVD/MFs
	private static final String OPT_SVD_NUMBER_BINS = "Bins for SVD recommender";
	private static final String OPT_SVD_BETA = "Beta value for SVD";
	private static final String OPT_SVD_TIME_STRATEGY = "Time strategy for SVD";
	private static final String OPT_SVD_REG_USER = "regularization for users for SVD";
	private static final String OPT_SVD_REG_ITEM = "regularization for items for SVD";
	private static final String OPT_SVD_REG_BIAS = "regularization for biases for SVD";
	private static final String OPT_SVD_LEARN_RATE = "learning rate for SVD";
	private static final String OPT_SVD_MAX_LEARN_RATE = "maximum learning rate for SVD";
	private static final String OPT_SVD_DECAY = "decay for SVD";
	private static final String OPT_SVD_IS_BOLD_DRIVER = "bold river for SVD";
	private static final String OPT_SVD_REG_IMP_ITEM = "regularization for implicit items for SVD";
	
	// Parameters for RankGeoFM
	private static final String OPT_EPSILON = "epsilon";
	private static final String OPT_C = "c";
	
	//Rerankers
	private static final String OPT_RERANKERS = "re-rankers";
	private static final String OPT_PROB_SMOOTH = "probability smoothing";
	private static final String OPT_PATTERN_LENGHT = "length of the pattern to search";
	private static final String OPT_FILL_STRATEGY = "opt filling strategy";
	private static final String OPT_RERANK_MAINTAIN_FIRST = "maintainFirstItemRecommended";

	
	//For caching Items (tour recommenders)
	private static final String OPT_CACHED_ITEMS = "cachedItems";
	private static final String OPT_NUM_ITEMS_COMPUTE = "number of items compute";
	
	//POI and tour recommenders
	private static final String OPT_POI_CITY_FILE = "poi city mapping file";
	private static final String OPT_CITY_SIM_FILE = "city similarity file";
	private static final String OPT_MATCHING_CATEGORY = "matching category for recommendation";
	private static final String OPT_POI_STIMATED_IME = "estimated poi time";
	private static final String OPT_CITY = "city selection";
	private static final String OPT_MAX_DISTANCE = "maximum distance";
	
	//TimePop parameters
	private static final String OPT_TIMEPOP_BETA = "Time pop beta";
	private static final String OPT_TIMEPOP_MINIMUM_CANDIDATES = "Minimum number of candidates to not use popularity";
	private static final String OPT_TIMEPOP_REFERING_TIMESTAMP = "Referring timestamp for TimePop";

	//Other
	private static final String OPT_MAX_POP_ITEMS = "max pop items";
	private static final String OPT_NUMBERRUNS = "number runs";
	private static final String OPT_REVERSE = "reverse";
	private static final String OPT_IMPUTE_ALL_USERS = "imputeAllUsers";
	private static final String OPT_IMPUTE_USERS_TESTS = "imputeUsersTest";
	private static final String OPT_BINARY_IMPUTATION = "binaryImputation";
	private static final String OPT_RANDOM_USERS_IMPUTATION = "random users imputation";
	private static final String OPT_GROUPING_IMPUTATION = "use grouping for imputation";
	private static final String OPT_SIZEGROUPING_IMPUTATION = "use size for imputation";
	
	
	
	/**For generating result files**/
	private static final String OPT_OUT_RESULT_FILE = "output result file";
	private static final String OPT_OUT_RESULT_FILE_2 = "output result file 2";
	private static final String OPT_OVERWRITE = "output result file overwrite";
	private static final String OPT_OUT_SIM_FILE = "output similarity file";
	private static final String OPT_EVALUATION_FILE = "evaluation file";
	private static final String OPT_PERCENTAGE = "percentage";
	private static final String OPT_USE_PERCENTAGE = "use percentage";
	private static final String OPT_NUMBER_BINS = "number bins";
	private static final String OPT_POP_EV_METRICS = "popularity evaluation";

	
	/**Feature variables (content based support and so on)**/
	private static final String OPT_FEATURES_ITEM_COL = "features item column";
	private static final String OPT_FEATURES_FEAT_COL = "features feature column";
	private static final String OPT_USER_FEATURE_FILE = "user feature file";
	private static final String OPT_FEATURE_DISTRIBUTION = "feature distribution";
	private static final String OPT_USER_FEATURE_SELECTED = "user feature selected";
	private static final String OPT_ITEM_FEATURE_SELECTED = "item feature selected";

	private static final String OPT_LST_USER_FEATURES = "list of user features for recommender";
	private static final String OPT_DIFFERENT_TRAININGS_USER_FEAUTURES = "see if we create a training set for each user feature";
	private static final String OPT_FEATURES_SEP = "features separator";
	private static final String OPT_FEATURES_HEADER = "features header";
	private static final String OPT_SAVING_FEATURES_SCHEME = "savingFeaturesScheme";
	private static final String OPT_ACTORS_SAVING_FILE = "actors saving file";
	private static final String OPT_DIRECTORSSAVINGFILE = "directorsSavingFile";
	private static final String OPT_GENRESSAVINGFILE = "genresSavingFile";
	private static final String OPT_FEATURE_READING_FILE = "featureFile";
	private static final String OPT_FEATURE_READING_FILE2 = "featureFile2";
	private static final String OPT_GENRESFILE = "genresFile";
	private static final String OPT_DIRECTORSFILE = "directorsFile";
	private static final String OPT_ACTORSFILE = "actorsFile";
	private static final String OPT_TAGSFILE = "tagsFile";
	private static final String OPT_FEATURE_CB = "feature cb similarity";
	private static final String OPT_CB_USER_TRANSFORMATION = "cb user transformation";
	private static final String OPT_CB_ITEM_TRANSFORMATION = "cb item transformation";
	private static final String OPT_CONSECUTIVE_SAME_CHECKINS = "remove consecutive checkins if they items are the same";
	private static final String OPT_MAX_SPEED_CHECKINS = "max speed between checkisn for not removing";
	
	/**Variables only used in evaluation**/
	private static final String OPT_RECOMMENDED_FILE = "recommendedFile";
	private static final String OPT_RECOMMENDED_FILE2 = "recommendedFile2";

	private static final String OPT_RECOMMENDED_FILES = "recommendedFiles";
	private static final String OPT_CUTOFF = "ranksysCutoff";
	private static final String OPT_MIN_ITEMS = "minRatingItems";
	private static final String OPT_MIN_USERS = "minRatingUsers";
	
	private static final String OPT_MIN_UNIQUE_RATINGS_ITEMS = "minUniqueRatingItems";
	private static final String OPT_MIN_UNIQUE_RATINGS_USERS = "minUniqueRatingUsers";
	
	private static final String OPT_THRESHOLD = "relevance or matching threshold";
	private static final String OPT_USER_MODEL_WGT = "user weight model weight";
	private static final String OPT_DEFAULT_SCORE = "default score";
	private static final String OPT_RANDOM_SEED = "random seed";
	private static final String OPT_STARTING_INDEX_RANKSIM = "index for ranksim";

	
	// Parameters for RankSys non accuracy evaluation
	private static final String OPT_RANKSYS_RELEVANCE_MODEL = "ranksys relevance model";
	private static final String OPT_RANKSYS_DISCOUNT_MODEL = "ranksys discount model";
	private static final String OPT_RANKSYS_BACKGROUND = "ranksys background for relevance model";
	private static final String OPT_RANKSYS_BASE = "ranksys base for discount model";
	
	// Files of items and times for our freshness/novelty metrics
	private static final String OPT_TIMES_AVERAGE_ITEMS = "file of average time items";
	private static final String OPT_TIMES_FIRST_ITEMS = "file of first time items";
	private static final String OPT_TIMES_LAST_ITEMS = "file of last time items";
	private static final String OPT_TIMES_MEDIAN_ITEMS = "file of median time items";
	private static final String OPT_RELEASE_DATE_ITEMS = "file of release dates of items";
	
	//For evaluation of novelty/diversity/freshness/POI metrics...
	private static final String OPT_COMPUTE_ANTI_METRICS = "Compute anti-metrics";
	private static final String OPT_COMPUTE_RATIOS = "compute ratios";
	private static final String OPT_ANTI_RELEVANCE_THRESHOLD = "Anti relevance threshold";
	private static final String OPT_COMPUTE_ONLY_ACC = "compute also dividing by the number of users in the test set";
	private static final String OPT_NOT_EXACT_METRIC = "execute the non exact relevance metric";

	private static final String OPT_SCORE_FREQ = "use simple score or the frequency for the AvgDis";
	private static final String OPT_LCS_EVALUATION = "LCS evaluation for taking into ccount the temporal order of the items";
	private static final String OPT_COMP_DISTANCES = "compute distances ev";
	
	private static final String OPT_COMPUTE_USER_FILTER = "compute user filter";
	private static final String OPT_PERUSER = "per User";


	/**For selecting best recommenders**/
	private static final String OPT_RECOMMENDERS_SEL = "recommendersSelected";
	private static final String OPT_METRIC_COLUMN = "metric column";
	private static final String OPT_METRIC_SELECTED = "metric selected";
	private static final String OPT_RESULT_COLUMN = "result column";
	private static final String OPT_RECOMMENDER_COLUMN = "recommender column";
	private static final String OPT_MATCHING_STRING = "matching string";

	/**Generate and filter new datasets**/
	//General
	private static final String OPT_IGNORE_FIRST_LINE = "ignore first line";
	private static final String OPT_POI_CATEGORY_FILE = "poi category file";
	private static final String OPT_MODE = "specific mode of an algorithm";
	// Mappings
	private static final String OPT_MAPPING_ITEMS = "file of mapping items";
	private static final String OPT_MAPPING_USERS = "file of mapping users";
	private static final String OPT_MAPPING_CATEGORIES = "file of mapping categories";
	// Generate datasets
	private static final String OPT_NEW_DATASET = "file for the new dataset";
	private static final String OPT_COORD_FILE = "poi coordinate file";
	private static final String OPT_WRAPPER_STRATEGY = "wrapper strategy for ranksys (removing timestamps)";
	private static final String OPT_WRAPPER_STRATEGY_TIMESTAMPS = "wrapper strategy for ranksys (strategy for the timestamps)";
	private static final String OPT_USING_WRAPPER = "using wrapper in ranksys recommenders";
	//Working with sessions
	private static final String OPT_MIN_SESSIONS = "min number of sessions";
	private static final String OPT_LIMIT_BETWEEN_ITEMS = "Limit condition between the items to be considered them as a session";
	private static final String OPT_SESSION_TYPE = "Type of session (distance, time etc)";
	private static final String OPT_MIN_DIFF_TIME = "min diff time";
	private static final String OPT_MAX_DIFF_TIME = "Max difference between timestamps";
	private static final String OPT_MIN_DIFF_TIME_2 = "min diff time2";
	private static final String OPT_MIN_CLOSE_PREF_BOT = "min close pref bot";
	private static final String OPT_PRINT_UTC = "print UTC";
	private static final String OPT_MIN_TIME = "min time";
	private static final String OPT_MAX_TIME = "max time";
	
	//Other
	private static final String OPT_LIST_NUMBERS = "list of numbers";
	private static final String OPT_PERCENTAGE_INCREMENT = "percentage increment";
	
	/***********************/
	/**Some default values**/
	/***********************/
	
	/**For Recommendation**/
	//General
	private static final String DEFAULT_SVDBETA = "1";
	private static final String DEFAULT_REG1 = "0.01"; // from mymedialite
	private static final String DEFAULT_REG2 = "0.001"; // from mymedialite
	private static final String DEFAULT_ETA = "0.05"; //from the USG recommender 
	private static final String DEFAULT_COMPDISTANCES = "false";
	private static final String DEFAULT_NORMALIZE = "true";
	
	//Some default values from MFs
	private static final String DEFAULT_FACTORS = "20";
	private static final String DEFAULT_ITERACTIONS = "20";
	
	//Some default values for HKV
	private static final String DEFAULT_ALPHA_HKV = "1"; 
	private static final String DEFAULT_LAMBDA_HKV = "0.1";
	
	//SVD parameters from librec: https://github.com/guoguibing/librec/blob/3ab54c4432ad573678e7b082f05785ed1ff31c5c/core/src/main/java/net/librec/recommender/MatrixFactorizationRecommender.java
	private static final String DEFAULT_SVD_LEARNING_RATE = "0.01f";
	private static final String DEFAULT_SVD_MAX_LEARNING_RATE = "1000f";
	private static final String DEFAULT_SVD_REG_USER = "0.01f";
	private static final String DEFAULT_SVD_REG_ITEM = "0.01f";
	private static final String DEFAULT_SVD_REG_BIAS = "0.01f";
	private static final String DEFAULT_IS_BOLD_DRIVER = "false";
	private static final String DEFAULT_DECAY = "1f";
	private static final String DEFAULT_REG_IMP_ITEM = "0.025f";
	
	// Some values for rank geoFM from Librec
	private static final String DEFAULT_EPSILON = "0.3";
	private static final String DEFAULT_C = "1.0";
	
	// Some default values for recommenders
	private static final String DEFAULT_COMPLETE_INDEXES = "false";
	private static final String DEFAULT_SCORE = "SIMPLE";
	private static final String DEFAULT_CB_USER_TRANSFORMATION = "WITHRATING";
	private static final String DEFAULT_CB_ITEM_TRANSFORMATION = "BINARY";
	private static final String DEFAULT_MAX_DISTANCE = "0.1";
	private static final String DEFAULT_MAX_POP_ITEMS= "300";

	
	
	// Generate evaluation
	private static final String DEFAULT_COMPUTE_ANTI_METRICS = "false";
	private static final String DEFAULT_ITEMS_RECOMMENDED = "100";
	private static final String DEFAULT_RECOMMENDATION_STRATEGY = "TRAIN_ITEMS";

	private static final String DEFAULT_OVERWRITE = "false";
	private static final String DEFAULT_OPT_COMPUTE_ONLY_ACC = "true";
	private static final String DEFAULT_OPT_NOT_EXACT_METRIC = "false";
	private static final String DEFAULT_COMPUTE_USER_FILTER = "false";

	// Default for metrics
	private static final String DEFAULT_OPT_COMPUTE_RATIOS = "false";
	private static final String DEFAULT_LCS_EVALUTATION = "false";
	private static final String DEFAULT_ANTI_RELEVANCE_THRESHOLD = "0";
	private static final String DEFAULT_MATCHING_THRESHOLD = "1";
	private static final String DEFAULT_USER_MODEL_WGT = "NORMAL";



	//Imputation
	private static final String DEFAULT_OPT_PERUSER = "perUser";
	private static final String DEFAULT_OPT_REVERSE = "false";
	
	private static final String DEFAULT_GROUPING_IMPUTATION = "false";
	private static final String DEFAULT_SIZEGROUPING_IMPUTATION = "false";
	
	
	
	public static void main(String[] args) throws Exception {

		String step = "";
		CommandLine cl = getCommandLine(args);
		if (cl == null) {
			System.out.println("Error in arguments");
			return;
		}

		// Obtain the arguments these 2 are obligatory
		step = cl.getOptionValue(OPT_CASE);
		String trainFile = cl.getOptionValue(OPT_TRAIN_FILE);

		System.out.println(step);
		switch (step) {


		/***
		 * Pure ranksys recommenders sections (using Preference data with NO timestamps) 
		 */
		// RankSys only (not reading similarities). Not working with time nor repetitions datasets
		// (working with FastPreferenceData only)
		case "ranksysOnlyComplete": {
			System.out.println(Arrays.toString(args));
			String trainFile2 = cl.getOptionValue(OPT_TRAIN_FILE_2);

			String outputFile = cl.getOptionValue(OPT_OUT_RESULT_FILE);
			String ranksysSimilarity = cl.getOptionValue(OPT_RANKSYS_SIM);
			String ranksysSimilarity2 = cl.getOptionValue(OPT_RANKSYS_SIM_2);
			Integer maxPopItems = Integer.parseInt(cl.getOptionValue(OPT_MAX_POP_ITEMS, DEFAULT_MAX_POP_ITEMS));

			
			String userFeatureFile = cl.getOptionValue(OPT_USER_FEATURE_FILE);
			String userFeaturesForRecommender = cl.getOptionValue(OPT_LST_USER_FEATURES);
			
			String combiner = cl.getOptionValue(OPT_COMB_AGGREGATE_LIBRARY);
			String normalizer = cl.getOptionValue(OPT_NORM_AGGREGATE_LIBRARY);

			String ranksysRecommender = cl.getOptionValue(OPT_RANKSYS_REC);
			String additionalRankysRecommenders = cl.getOptionValue(OPT_RANKSYS_REC2);
			String testFile = cl.getOptionValue(OPT_TEST_FILE);
			String testFile2 = cl.getOptionValue(OPT_TEST_FILE_2); // this is the test file for foursqr with the cross
																	// domain approach
			String socialInfluenceFile = cl.getOptionValue(OPT_SOC_INFLUENCE_FILE);
			String hybridWeights = cl.getOptionValue(OPT_HYBRID_WEIGHTS);
			
			//MFs Parameters
			Integer numIterations = Integer.parseInt(cl.getOptionValue(OPT_NUM_INTERACTIONS, DEFAULT_ITERACTIONS));
			Integer numFactors = Integer.parseInt(cl.getOptionValue(OPT_K_FACTORIZER, DEFAULT_FACTORS));
			Double thresholdAntiSim = Double.parseDouble(cl.getOptionValue(OPT_THRESHOLD, "3"));
			String mode = cl.getOptionValue(OPT_MODE);
			
			//SVD Parameters
			Double regUser = Double.parseDouble(cl.getOptionValue(OPT_SVD_REG_USER, DEFAULT_SVD_REG_USER));
			Double regItem = Double.parseDouble(cl.getOptionValue(OPT_SVD_REG_ITEM, DEFAULT_SVD_REG_ITEM));
			Double regBias = Double.parseDouble(cl.getOptionValue(OPT_SVD_REG_BIAS, DEFAULT_SVD_REG_BIAS));
			Double learnRate = Double.parseDouble(cl.getOptionValue(OPT_SVD_LEARN_RATE, DEFAULT_SVD_LEARNING_RATE));
			Double maxRate = Double.parseDouble(cl.getOptionValue(OPT_SVD_MAX_LEARN_RATE, DEFAULT_SVD_MAX_LEARNING_RATE));
			Boolean isboldDriver = Boolean.parseBoolean(cl.getOptionValue(OPT_SVD_IS_BOLD_DRIVER, DEFAULT_IS_BOLD_DRIVER));
			Double regImpItem = Double.parseDouble(cl.getOptionValue(OPT_SVD_REG_IMP_ITEM, DEFAULT_REG_IMP_ITEM));
			Double decay = Double.parseDouble(cl.getOptionValue(OPT_SVD_DECAY, DEFAULT_DECAY));
			Double beta = Double.parseDouble(cl.getOptionValue(OPT_SVD_BETA, DEFAULT_SVDBETA));
					
			
			Double maxDistance = Double.parseDouble(cl.getOptionValue(OPT_MAX_DISTANCE, DEFAULT_MAX_DISTANCE));
			
			// Rank Geo
			Double epsilon = Double.parseDouble(cl.getOptionValue(OPT_EPSILON, DEFAULT_EPSILON));
			Double c = Double.parseDouble(cl.getOptionValue(OPT_C, DEFAULT_C));

			
			//eta
			Double eta = Double.parseDouble(cl.getOptionValue(OPT_ETA, DEFAULT_ETA));
			
			//HKV Factorizer
			Double alphaFactorizer = Double.parseDouble(cl.getOptionValue(OPT_ALPHA_FACTORIZER, DEFAULT_ALPHA_HKV));
			Double alphaFactorizer2 = Double.parseDouble(cl.getOptionValue(OPT_ALPHA_FACTORIZER2, DEFAULT_ALPHA_HKV));
			Double theta = Double.parseDouble(cl.getOptionValue(OPT_THETA_FACTORIZER, DEFAULT_ALPHA_HKV));


			Double lambdaFactorizer = Double.parseDouble(cl.getOptionValue(OPT_LAMBDA_FACTORIZER, DEFAULT_LAMBDA_HKV));		

			boolean useSigmoid = Boolean.parseBoolean(cl.getOptionValue(cl.getOptionValue(OPT_USE_SIGMOID, "false")));
			
			// For candidates items matching a category
			String featureFile = cl.getOptionValue(OPT_FEATURE_READING_FILE);
			String realFeatureFile = cl.getOptionValue(OPT_FEATURE_READING_FILE2);

			String itemFeatureColumn = cl.getOptionValue(OPT_FEATURES_ITEM_COL);
			String featureColumn = cl.getOptionValue(OPT_FEATURES_FEAT_COL);
			String categoryMatch = cl.getOptionValue(OPT_MATCHING_CATEGORY);

			String poiCoordsFile = cl.getOptionValue(OPT_COORD_FILE);
			String poiCityFile = cl.getOptionValue(OPT_POI_CITY_FILE);
			String citySimFile = cl.getOptionValue(OPT_CITY_SIM_FILE);
			Boolean inverse = Boolean.parseBoolean(cl.getOptionValue(OPT_INVERSE, "true"));
			Boolean negativeAnti = Boolean.parseBoolean(cl.getOptionValue(OPT_NEGANTI, "true"));

			
			
			String regularization1S = cl.getOptionValue(OPT_REGULARIZATION_L1, DEFAULT_REG1);
			String regularization2S = cl.getOptionValue(OPT_REGULARIZATION_L2, DEFAULT_REG2);
			
			String scoreFreq = cl.getOptionValue(OPT_SCORE_FREQ, DEFAULT_SCORE);
			Boolean normalize = Boolean.parseBoolean(cl.getOptionValue(OPT_NORMALIZE, DEFAULT_NORMALIZE));

			Double regularization1 = Double.parseDouble(regularization1S);
			Double regularization2 = Double.parseDouble(regularization2S);
			
			final int numberItemsRecommend = Integer.parseInt(cl.getOptionValue(OPT_ITEMS_RECOMMENDED));
			int neighbours = Integer.parseInt(cl.getOptionValue(OPT_NEIGH));

			String recommendationStrategy = cl.getOptionValue(OPT_RECOMMENDATION_STRATEGY, DEFAULT_RECOMMENDATION_STRATEGY);
			String userCBTransformationS = cl.getOptionValue(OPT_CB_USER_TRANSFORMATION, DEFAULT_CB_USER_TRANSFORMATION);
			String itemCBTransformationS = cl.getOptionValue(OPT_CB_ITEM_TRANSFORMATION, DEFAULT_CB_ITEM_TRANSFORMATION);
			
			USERTRANSFORMATION ut = ExperimentUtils.obtCBUserTransformation(userCBTransformationS);
			ITEMTRANSFORMATION it = ExperimentUtils.obtCBItemTransformation(itemCBTransformationS);
			UserFeatureData<Long, String, Double> ufD = null;
			if (userFeatureFile != null) {
				ufD = SimpleUserFeatureData.load(SimpleUserFeaturesReader.get().read(userFeatureFile, lp, sp));
			}
			
			FeatureData<Long, String, Double> featureData = null;
			if (realFeatureFile != null) {
				featureData = SimpleFeatureData.load(SimpleFeaturesReader.get().read(realFeatureFile, lp, sp));
			}
			
			//In case we want to filter the training set by users/item features
			String userSelFeature = cl.getOptionValue(OPT_USER_FEATURE_SELECTED);
			String itemSelFeature = cl.getOptionValue(OPT_ITEM_FEATURE_SELECTED);


			
			SCORES_FREQUENCY scoref = ExperimentUtils.obtScoreFreq(scoreFreq);


			String[] outputsFiles = outputFile.split(",");
			String overwrite = cl.getOptionValue(OPT_OVERWRITE, DEFAULT_OVERWRITE);

			boolean allExist = true;
			for (String outputSingleFile : outputsFiles) {
				File f = new File(outputSingleFile);
				if (f.exists() && !f.isDirectory() && Boolean.parseBoolean(overwrite) == false) {
					System.out.println("Ignoring " + f + " because it already exists");
				} else {
					System.out.println(f + " does not exist. Computing.");
					allExist = false;
				}
			}
			if (allExist && Boolean.parseBoolean(overwrite) == false){
				 // If all the outputFiles exist and the overwrite is in false, we end
				return;
			}
			//This array must be 1 except for the HybridRecommender
			String completeOrNotS[] = cl.getOptionValue(OPT_COMPLETE_INDEXES, DEFAULT_COMPLETE_INDEXES).split(",");
			boolean completeOrNot = true;
			if (completeOrNotS.length == 1) {
				completeOrNot = Boolean.parseBoolean(completeOrNotS[0]);
			}
					

			String[] differentTests = testFile.split(",");

			FastPreferenceData<Long, Long> trainPrefData = null;
			FastPreferenceData<Long, Long> testPrefData = null;
			
			// Testfile2 should be a concatenation of all test files that we are using
			if (differentTests.length == 1) {
				trainPrefData = ExperimentUtils.loadTrainFastPreferenceData(trainFile, testFile, completeOrNot, true);
			}
			else {
				trainPrefData = ExperimentUtils.loadTrainFastPreferenceData(trainFile, testFile2, completeOrNot, true);
			}
			
			
			
			trainPrefData = SequentialRecommendersUtils.filterPreferenceDataByUserItemsFeatures(trainPrefData, 
					ufD, featureData, userSelFeature, itemSelFeature);
			
			System.out.println("Some stats: ");
			System.out.println("Users: " + trainPrefData.numUsers());
			System.out.println("Items: " + trainPrefData.numItems());
			System.out.println("Users with preferences: " + trainPrefData.numUsersWithPreferences());
			System.out.println("Items with preferences: " + trainPrefData.numItemsWithPreferences());
			System.out.println("Num preferences: " + trainPrefData.numPreferences());
			


			
			System.out.println("Not using wrapper (working with no timestamps)");

			Recommender<Long, Long> rankSysrec = SequentialRecommendersUtils.obtRankSysRecommeder(ranksysRecommender, additionalRankysRecommenders,
					ranksysSimilarity, ranksysSimilarity2, trainPrefData, neighbours, numFactors, alphaFactorizer, lambdaFactorizer, 
					numIterations, poiCoordsFile, poiCityFile, citySimFile, regularization1, regularization2, scoref, regUser, regItem, learnRate, maxRate, regBias, decay, 
					isboldDriver, regImpItem, socialInfluenceFile, eta, beta, normalize, ut, it, epsilon, c, featureData, ufD, combiner, normalizer, hybridWeights, completeOrNotS, 
					trainFile, testFile, maxDistance, mode, theta, alphaFactorizer2, useSigmoid, userFeaturesForRecommender, maxPopItems, thresholdAntiSim, inverse, negativeAnti);
			

			FastPreferenceData<Long, Long> ranksysTrainDataOriginal2 = trainPrefData;
			if (trainFile2 != null) {
				ranksysTrainDataOriginal2 = ExperimentUtils.loadTrainFastPreferenceData(trainFile2, testFile, completeOrNot, true);
			}
			System.out.println("TRAIN2: " + ranksysTrainDataOriginal2.numPreferences());
			System.out.println("TRAIN_NORMAL: " + trainPrefData.numPreferences());

			

			String[] candidatesFile = null;
			String[] categoriesMatch = null;
			
			if (featureFile != null) {
				candidatesFile = featureFile.split(",");
			}
			
			if (categoryMatch != null) {
				categoriesMatch = categoryMatch.split(",");
			}

			for (int i = 0; i < differentTests.length; i++) {
				System.out.println("Analyzing " + differentTests[i]);
				System.out.println("Recommender file " + outputsFiles[i]);
				File f = new File(outputsFiles[i]);
				if (f.exists() && !f.isDirectory() && Boolean.parseBoolean(overwrite) == false) {
					System.out.println(outputsFiles[i] + " exist, ignoring that recommender file");
					continue;
				}

				testPrefData = ExperimentUtils.loadTrainFastPreferenceData(trainFile, differentTests[i], true, false);
				
				testPrefData = SequentialRecommendersUtils.filterPreferenceDataByUserItemsFeatures(testPrefData, 
						ufD, featureData, userSelFeature, itemSelFeature);
				
				// Normal recommendation. Only recommend for test users items that have not been
				// seen in train
				if (featureFile == null || itemFeatureColumn == null || featureColumn == null || categoryMatch == null) {
					System.out.println("Writing recommended file. Not items candidates file provided. All candidates are the items not seen by that user in train.");
					PredicatesStrategies.ranksysWriteRanking(trainPrefData, ranksysTrainDataOriginal2, testPrefData, rankSysrec,
							outputsFiles[i], numberItemsRecommend,
							ExperimentUtils.obtRecommendationStrategy(recommendationStrategy));
				} else {
					System.out.println("Matching category " + categoriesMatch[i]);
					System.out.println("candidates " + candidatesFile[i]);
					System.out.println("Writing recommended file. Item candidates file provided. Candidates are items from that file that have not been rated by that user in train.");
					Stream<Long> candidates = SequentialRecommendersUtils.obtainCandidatesFile(candidatesFile[i],
							categoriesMatch[i], lp, Integer.parseInt(featureColumn),
							Integer.parseInt(itemFeatureColumn));
					PredicatesStrategies.ranksysWriteRankingFromCandidatesList(trainPrefData, testPrefData,
							rankSysrec, outputsFiles[i], numberItemsRecommend, candidates.collect(Collectors.toSet()), ExperimentUtils.obtRecommendationStrategy(recommendationStrategy) );
				}
			}

		}
			break;
			
		
			
			

		
			
		/***
		 * Skylines recommenders that read the test file and perform the test recommendations
		 */
			
		case "skylineRecommenders":
		case "testRecommenders": {
			// They need to receive the test set. We do not work with the train file
			System.out.println(Arrays.toString(args));
			String trainFile2 = cl.getOptionValue(OPT_TRAIN_FILE_2);
			Double antiRelTh = Double
					.parseDouble(cl.getOptionValue(OPT_ANTI_RELEVANCE_THRESHOLD, DEFAULT_ANTI_RELEVANCE_THRESHOLD));
			String outputFile = cl.getOptionValue(OPT_OUT_RESULT_FILE);
			String testFile = cl.getOptionValue(OPT_TEST_FILE);
			Double threshold = Double.parseDouble(cl.getOptionValue(OPT_THRESHOLD, DEFAULT_MATCHING_THRESHOLD));
			boolean completeOrNot = Boolean
					.parseBoolean(cl.getOptionValue(OPT_COMPLETE_INDEXES, DEFAULT_COMPLETE_INDEXES));
			String recommendationStrategy = cl.getOptionValue(OPT_RECOMMENDATION_STRATEGY,
					DEFAULT_RECOMMENDATION_STRATEGY);
			final int numberItemsRecommend = Integer
					.parseInt(cl.getOptionValue(OPT_ITEMS_RECOMMENDED, DEFAULT_ITEMS_RECOMMENDED));

			String rankSysRecommender = cl.getOptionValue(OPT_RANKSYS_REC);
			String overwrite = cl.getOptionValue(OPT_OVERWRITE, DEFAULT_OVERWRITE);

			File f = new File(outputFile);
			if (f.exists() && !f.isDirectory() && Boolean.parseBoolean(overwrite) == false) {
				System.out.println("Ignoring " + f + " because it already exists");
				return;
			}

			FastPreferenceData<Long, Long> ranksysTrainDataOriginal = ExperimentUtils
					.loadTrainFastPreferenceData(trainFile, testFile, completeOrNot, true);

			FastPreferenceData<Long, Long> ranksysTestDataOriginal = ExperimentUtils
					.loadTrainFastPreferenceData(trainFile, testFile, true, false);

			FastTemporalPreferenceDataIF<Long, Long> ranksysTestTemporal = ExperimentUtils
					.loadTrainFastTemporalFeaturePreferenceData(trainFile, testFile, completeOrNot, false);

			Recommender<Long, Long> rec = SequentialRecommendersUtils.obtRankSysSkylineRecommender(rankSysRecommender,
					ranksysTrainDataOriginal, ranksysTestDataOriginal, ranksysTestTemporal, threshold, antiRelTh);
			
			FastPreferenceData<Long, Long> ranksysTrainDataOriginal2 = ranksysTrainDataOriginal;
			if (trainFile2 != null) {
				ranksysTrainDataOriginal2 = ExperimentUtils.loadTrainFastPreferenceData(trainFile2, testFile, completeOrNot, true);
			}
			

			PredicatesStrategies.ranksysWriteRanking(ranksysTrainDataOriginal, ranksysTrainDataOriginal2, ranksysTestDataOriginal, rec, outputFile,
					numberItemsRecommend, ExperimentUtils.obtRecommendationStrategy(recommendationStrategy));

		}
		break;

		
				
		/***
		 * Data preprocessing section. K-core, new ids for the dataset, different splits, statistics...	
		 */
		
		//Apply a k-core to a dataset (careful, changed to use parameters)
		case "Kcore":
		case "DatasetReduction": {
			System.out.println(Arrays.toString(args));
			String outputFile = cl.getOptionValue(OPT_OUT_RESULT_FILE);
			Integer minUserRatings = Integer.parseInt(cl.getOptionValue(OPT_MIN_USERS));
			Integer minItemRatings = Integer.parseInt(cl.getOptionValue(OPT_MIN_ITEMS));
			
			Integer minUserUniqueRatings = Integer.parseInt(cl.getOptionValue(OPT_MIN_UNIQUE_RATINGS_USERS, "0"));
			Integer minItemUniqueRatings = Integer.parseInt(cl.getOptionValue(OPT_MIN_UNIQUE_RATINGS_ITEMS, "0"));
			
			ProcessData.DatasetReductionRatings(trainFile, outputFile, minUserRatings, minItemRatings, minUserUniqueRatings, minItemUniqueRatings);
		}
			break;

			

			
		//Global split (temporal)
		case "TemporalGlobalSplit": {
			System.out.println("-o TemporalGlobalSplit -trf completeFile dstPathTrain dstPathTest trainPercent");
			System.out.println(Arrays.toString(args));
			ProcessData.DatasetTemporalGlobalSplit(trainFile, args[4], args[5], Double.parseDouble(args[6]));
		}
			break;
		
			
		
			
		//Statistics of the datasets
		case "Statistics": {
			System.out.println("-o Statistics -trf completeFile statisticDestination");
			System.out.println(Arrays.toString(args));
			ProcessData.Stats(trainFile, args[4]);
		}
		break;
		
		case "SatatisticsPercentageCheckinsByMidPoint": {
			System.out.println(Arrays.toString(args));
			String originalPOICoordsFile = cl.getOptionValue(OPT_COORD_FILE);
			
			String scoreFreq = cl.getOptionValue(OPT_SCORE_FREQ, DEFAULT_SCORE);
			String listNumbers = cl.getOptionValue(OPT_LIST_NUMBERS);
			String resultFile = cl.getOptionValue(OPT_OUT_RESULT_FILE);
			
			SCORES_FREQUENCY scoref = ExperimentUtils.obtScoreFreq(scoreFreq);

			
            Map<Long, Tuple2<Double, Double>> mapCoordinates = SequentialRecommendersUtils.POICoordinatesMap(originalPOICoordsFile, lp);

		
			FastPreferenceData<Long, Long> data = ExperimentUtils.loadTrainFastPreferenceData(trainFile, trainFile, false, true);

			UsersMidPoints<Long,Long> usermid = new UsersMidPoints<>(data, mapCoordinates, scoref);
			

			List<String> listNumbersL = Arrays.asList(listNumbers.split(","));
			List<Double> listNumbersF = listNumbersL.stream().mapToDouble(s -> Double.parseDouble(s)).boxed().collect(Collectors.toList());
			
			Collections.sort(listNumbersF);
			
			POIProcessData.generatePercentageCheckinsPerUserByDistance(resultFile, data, usermid, mapCoordinates, listNumbersF);
            
		}
		break;


		
		
		case "GenerateNewDataset":
		case "ImputationDataset": {
			System.out.println(Arrays.toString(args));

			String outputFile = cl.getOptionValue(OPT_OUT_RESULT_FILE);
			String ranksysSimilarity = cl.getOptionValue(OPT_RANKSYS_SIM);
			String ranksysSimilarity2 = cl.getOptionValue(OPT_RANKSYS_SIM_2);
			Integer maxPopItems = Integer.parseInt(cl.getOptionValue(OPT_MAX_POP_ITEMS, DEFAULT_MAX_POP_ITEMS));
			Boolean imputeAllUsers = Boolean.parseBoolean(cl.getOptionValue(OPT_IMPUTE_ALL_USERS, "false"));
			Boolean imputeTestUsers = Boolean.parseBoolean(cl.getOptionValue(OPT_IMPUTE_USERS_TESTS, "false"));
			Boolean useGrouping = Boolean.parseBoolean(cl.getOptionValue(OPT_GROUPING_IMPUTATION, DEFAULT_GROUPING_IMPUTATION));
			Boolean useSizeGrouping = Boolean.parseBoolean(cl.getOptionValue(OPT_SIZEGROUPING_IMPUTATION, DEFAULT_SIZEGROUPING_IMPUTATION));

			Boolean binaryImputation = Boolean.parseBoolean(cl.getOptionValue(OPT_BINARY_IMPUTATION, "true"));


			
			String userFeatureFile = cl.getOptionValue(OPT_USER_FEATURE_FILE);
			String userFeaturesForRecommender = cl.getOptionValue(OPT_LST_USER_FEATURES);
			
			String combiner = cl.getOptionValue(OPT_COMB_AGGREGATE_LIBRARY);
			String normalizer = cl.getOptionValue(OPT_NORM_AGGREGATE_LIBRARY);

			String ranksysRecommender = cl.getOptionValue(OPT_RANKSYS_REC);
			String additionalRankysRecommenders = cl.getOptionValue(OPT_RANKSYS_REC2);
			String testFile = cl.getOptionValue(OPT_TEST_FILE);
			String testFile2 = cl.getOptionValue(OPT_TEST_FILE_2); // this is the test file for foursqr with the cross
																	// domain approach
			String socialInfluenceFile = cl.getOptionValue(OPT_SOC_INFLUENCE_FILE);
			String hybridWeights = cl.getOptionValue(OPT_HYBRID_WEIGHTS);
			
			//MFs Parameters
			Integer numIterations = Integer.parseInt(cl.getOptionValue(OPT_NUM_INTERACTIONS, DEFAULT_ITERACTIONS));
			Integer numFactors = Integer.parseInt(cl.getOptionValue(OPT_K_FACTORIZER, DEFAULT_FACTORS));
			Double thresholdAntiSim = Double.parseDouble(cl.getOptionValue(OPT_THRESHOLD, "3"));
			String mode = cl.getOptionValue(OPT_MODE);
			
			//SVD Parameters
			Double regUser = Double.parseDouble(cl.getOptionValue(OPT_SVD_REG_USER, DEFAULT_SVD_REG_USER));
			Double regItem = Double.parseDouble(cl.getOptionValue(OPT_SVD_REG_ITEM, DEFAULT_SVD_REG_ITEM));
			Double regBias = Double.parseDouble(cl.getOptionValue(OPT_SVD_REG_BIAS, DEFAULT_SVD_REG_BIAS));
			Double learnRate = Double.parseDouble(cl.getOptionValue(OPT_SVD_LEARN_RATE, DEFAULT_SVD_LEARNING_RATE));
			Double maxRate = Double.parseDouble(cl.getOptionValue(OPT_SVD_MAX_LEARN_RATE, DEFAULT_SVD_MAX_LEARNING_RATE));
			Boolean isboldDriver = Boolean.parseBoolean(cl.getOptionValue(OPT_SVD_IS_BOLD_DRIVER, DEFAULT_IS_BOLD_DRIVER));
			Double regImpItem = Double.parseDouble(cl.getOptionValue(OPT_SVD_REG_IMP_ITEM, DEFAULT_REG_IMP_ITEM));
			Double decay = Double.parseDouble(cl.getOptionValue(OPT_SVD_DECAY, DEFAULT_DECAY));
			Double beta = Double.parseDouble(cl.getOptionValue(OPT_SVD_BETA, DEFAULT_SVDBETA));
			Boolean reverse = Boolean.parseBoolean(cl.getOptionValue(OPT_REVERSE, DEFAULT_OPT_REVERSE));
			Boolean randomUsersImputation = Boolean.parseBoolean(cl.getOptionValue(OPT_RANDOM_USERS_IMPUTATION, "false"));		
			
			Double maxDistance = Double.parseDouble(cl.getOptionValue(OPT_MAX_DISTANCE, DEFAULT_MAX_DISTANCE));
			
			// Rank Geo
			Double epsilon = Double.parseDouble(cl.getOptionValue(OPT_EPSILON, DEFAULT_EPSILON));
			Double c = Double.parseDouble(cl.getOptionValue(OPT_C, DEFAULT_C));

			
			//eta
			Double eta = Double.parseDouble(cl.getOptionValue(OPT_ETA, DEFAULT_ETA));
			
			//HKV Factorizer
			Double alphaFactorizer = Double.parseDouble(cl.getOptionValue(OPT_ALPHA_FACTORIZER, DEFAULT_ALPHA_HKV));
			Double alphaFactorizer2 = Double.parseDouble(cl.getOptionValue(OPT_ALPHA_FACTORIZER2, DEFAULT_ALPHA_HKV));
			Double theta = Double.parseDouble(cl.getOptionValue(OPT_THETA_FACTORIZER, DEFAULT_ALPHA_HKV));


			Double lambdaFactorizer = Double.parseDouble(cl.getOptionValue(OPT_LAMBDA_FACTORIZER, DEFAULT_LAMBDA_HKV));		

			boolean useSigmoid = Boolean.parseBoolean(cl.getOptionValue(cl.getOptionValue(OPT_USE_SIGMOID, "false")));
			
			// For candidates items matching a category
			String realFeatureFile = cl.getOptionValue(OPT_FEATURE_READING_FILE2);



			String poiCoordsFile = cl.getOptionValue(OPT_COORD_FILE);
			String poiCityFile = cl.getOptionValue(OPT_POI_CITY_FILE);
			String citySimFile = cl.getOptionValue(OPT_CITY_SIM_FILE);
			Boolean inverse = Boolean.parseBoolean(cl.getOptionValue(OPT_INVERSE, "true"));
			Boolean negativeAnti = Boolean.parseBoolean(cl.getOptionValue(OPT_NEGANTI, "true"));

			
			
			String regularization1S = cl.getOptionValue(OPT_REGULARIZATION_L1, DEFAULT_REG1);
			String regularization2S = cl.getOptionValue(OPT_REGULARIZATION_L2, DEFAULT_REG2);
			
			String scoreFreq = cl.getOptionValue(OPT_SCORE_FREQ, DEFAULT_SCORE);
			Boolean normalize = Boolean.parseBoolean(cl.getOptionValue(OPT_NORMALIZE, DEFAULT_NORMALIZE));

			Double regularization1 = Double.parseDouble(regularization1S);
			Double regularization2 = Double.parseDouble(regularization2S);
			
			int neighbours = Integer.parseInt(cl.getOptionValue(OPT_NEIGH));

			String userCBTransformationS = cl.getOptionValue(OPT_CB_USER_TRANSFORMATION, DEFAULT_CB_USER_TRANSFORMATION);
			String itemCBTransformationS = cl.getOptionValue(OPT_CB_ITEM_TRANSFORMATION, DEFAULT_CB_ITEM_TRANSFORMATION);
			
			USERTRANSFORMATION ut = ExperimentUtils.obtCBUserTransformation(userCBTransformationS);
			ITEMTRANSFORMATION it = ExperimentUtils.obtCBItemTransformation(itemCBTransformationS);
			UserFeatureData<Long, String, Double> ufD = null;
			if (userFeatureFile != null) {
				ufD = SimpleUserFeatureData.load(SimpleUserFeaturesReader.get().read(userFeatureFile, lp, sp));
			}
			
			FeatureData<Long, String, Double> featureData = null;
			if (realFeatureFile != null) {
				featureData = SimpleFeatureData.load(SimpleFeaturesReader.get().read(realFeatureFile, lp, sp));
			}
			
			//In case we want to filter the training set by users/item features
			String userSelFeature = cl.getOptionValue(OPT_USER_FEATURE_SELECTED);
			String itemSelFeature = cl.getOptionValue(OPT_ITEM_FEATURE_SELECTED);


			
			SCORES_FREQUENCY scoref = ExperimentUtils.obtScoreFreq(scoreFreq);


			String[] outputsFiles = outputFile.split(",");
			String overwrite = cl.getOptionValue(OPT_OVERWRITE, DEFAULT_OVERWRITE);

			boolean allExist = true;
			for (String outputSingleFile : outputsFiles) {
				File f = new File(outputSingleFile);
				if (f.exists() && !f.isDirectory() && Boolean.parseBoolean(overwrite) == false) {
					System.out.println("Ignoring " + f + " because it already exists");
				} else {
					System.out.println(f + " does not exist. Computing.");
					allExist = false;
				}
			}
			if (allExist && Boolean.parseBoolean(overwrite) == false){
				 // If all the outputFiles exist and the overwrite is in false, we end
				return;
			}
			//This array must be 1 except for the HybridRecommender
			String completeOrNotS[] = cl.getOptionValue(OPT_COMPLETE_INDEXES, DEFAULT_COMPLETE_INDEXES).split(",");
			boolean completeOrNot = true;
			if (completeOrNotS.length == 1) {
				completeOrNot = Boolean.parseBoolean(completeOrNotS[0]);
			}
					

			String[] differentTests = testFile.split(",");

			FastPreferenceData<Long, Long> trainPrefData = null;
			
			// Testfile2 should be a concatenation of all test files that we are using
			if (differentTests.length == 1) {
				trainPrefData = ExperimentUtils.loadTrainFastPreferenceData(trainFile, testFile, completeOrNot, true);
			}
			else {
				trainPrefData = ExperimentUtils.loadTrainFastPreferenceData(trainFile, testFile2, completeOrNot, true);
			}
			
			
			
			trainPrefData = SequentialRecommendersUtils.filterPreferenceDataByUserItemsFeatures(trainPrefData, 
					ufD, featureData, userSelFeature, itemSelFeature);
			
			System.out.println("Some stats: ");
			System.out.println("Users: " + trainPrefData.numUsers());
			System.out.println("Items: " + trainPrefData.numItems());
			System.out.println("Users with preferences: " + trainPrefData.numUsersWithPreferences());
			System.out.println("Items with preferences: " + trainPrefData.numItemsWithPreferences());
			System.out.println("Num preferences: " + trainPrefData.numPreferences());
			


			
			System.out.println("Not using wrapper (working with no timestamps)");

			Recommender<Long, Long> rankSysrec = SequentialRecommendersUtils.obtRankSysRecommeder(ranksysRecommender, additionalRankysRecommenders,
					ranksysSimilarity, ranksysSimilarity2, trainPrefData, neighbours, numFactors, alphaFactorizer, lambdaFactorizer, 
					numIterations, poiCoordsFile, poiCityFile, citySimFile, regularization1, regularization2, scoref, regUser, regItem, learnRate, maxRate, regBias, decay, 
					isboldDriver, regImpItem, socialInfluenceFile, eta, beta, normalize, ut, it, epsilon, c, featureData, ufD, combiner, normalizer, hybridWeights, completeOrNotS, 
					trainFile, testFile, maxDistance, mode, theta, alphaFactorizer2, useSigmoid, userFeaturesForRecommender, maxPopItems, thresholdAntiSim, inverse, negativeAnti);
			

			Boolean perUser = Boolean.parseBoolean(cl.getOptionValue(OPT_PERUSER, DEFAULT_OPT_PERUSER));
			Double percentageIncrement = Double.parseDouble(cl.getOptionValue(OPT_PERCENTAGE_INCREMENT));
			
			
			Map<Long, Tuple2<Double, Double>> mapCoordinates = SequentialRecommendersUtils.POICoordinatesMap(cl.getOptionValue(OPT_COORD_FILE), lp);
            
            UsersMidPoints<Long,Long> usermid = new UsersMidPoints<>(trainPrefData, mapCoordinates, SCORES_FREQUENCY.FREQUENCY);
			double limitKm = Double.parseDouble(cl.getOptionValue(OPT_THRESHOLD));
			
			//Here is the imputation technique (HERE)
			
			System.out.println("limit KM is: " + limitKm);
			ImputedSimpleFastPreferenceDataRecommender<Long, Long> imputedPreferenceData = null;
			if (imputeAllUsers) {
				imputedPreferenceData = new ImputedSimpleFastPreferenceDataTrainUsers<>(trainPrefData, rankSysrec, perUser, percentageIncrement, PredicatesStrategies.recommendationStrategy.MIDPOINT_USER_TRAIN_ITEMS, mapCoordinates, usermid, limitKm, binaryImputation, useGrouping, useSizeGrouping);											
			}
			else if (imputeTestUsers) {
				FastPreferenceData<Long, Long> ranksysTestData = ExperimentUtils.loadTrainFastPreferenceData(trainFile, testFile, completeOrNot, false);
				imputedPreferenceData = new ImputedSimpleFastPreferenceDataTestUsers<>(trainPrefData, ranksysTestData, rankSysrec, perUser, percentageIncrement, PredicatesStrategies.recommendationStrategy.MIDPOINT_USER_TRAIN_ITEMS, mapCoordinates, usermid, limitKm, binaryImputation, useGrouping, useSizeGrouping);											
			}
			else
			{
				imputedPreferenceData = new ImputedSimpleFastPreferenceDataAvgPreferencesUsers<>(trainPrefData, rankSysrec, perUser, percentageIncrement, PredicatesStrategies.recommendationStrategy.MIDPOINT_USER_TRAIN_ITEMS, mapCoordinates, usermid, limitKm, reverse, binaryImputation, randomUsersImputation, useGrouping, useSizeGrouping);
			}
			imputedPreferenceData.writeimputeData(outputFile);	

		}
		break;
		
		//Method to remove from a rec file interactions that have already appeared in the training set
		case "filterRecFileByTrain": {
			System.out.println(Arrays.toString(args));

			String outputFile = cl.getOptionValue(OPT_OUT_RESULT_FILE);
		 	String recFile = cl.getOptionValue(OPT_RECOMMENDED_FILE);
		 	Integer limit = Integer.parseInt(cl.getOptionValue(OPT_THRESHOLD));
		 	
			FastTemporalPreferenceDataIF<Long, Long> ranksysTrainTemporal = ExperimentUtils.loadTrainFastTemporalFeaturePreferenceData(trainFile, trainFile, false, true);

			ProcessData.filterRecByTraining(ranksysTrainTemporal, recFile, limit, outputFile);

		}
		break;
		

			
		case "SpecificPOICoords": {
			System.out.println(
					"-o SpecificPOICoords -trf originalFullFile -coordFile coordFile -orf outputFileCoordsOfFile");
			System.out.println(Arrays.toString(args));

			String poisCoordsFile = cl.getOptionValue(OPT_COORD_FILE);
			String outputFile = cl.getOptionValue(OPT_OUT_RESULT_FILE);
			
			Map<Long, Tuple2<Double, Double>> coordinates = SequentialRecommendersUtils.POICoordinatesMap(poisCoordsFile, lp);

			POIProcessData.specificPOICoords(trainFile, coordinates, outputFile);
		}
			break;	


			
		case "generateNewCheckingFileWithTimeStamps": {
			System.out.println(Arrays.toString(args));

			String userMapFile = cl.getOptionValue(OPT_MAPPING_USERS);
			String itemsMapFile = cl.getOptionValue(OPT_MAPPING_ITEMS);
			String newDataset = cl.getOptionValue(OPT_NEW_DATASET);
			Boolean printUTC = Boolean.parseBoolean(cl.getOptionValue(OPT_PRINT_UTC, "true"));

			POIProcessData.generateNewCheckinFileWithTimeStamps(trainFile, userMapFile, itemsMapFile, newDataset, printUTC);
		}
			break;
		

			
		case "printClosestPOIs": {
			System.out.println(Arrays.toString(args));

			// String trainingFile, String poisCoordsFile, String resultFileNN, String
			// resultFileCoordsOfTrain, int nn
			POIProcessData.printClosestPOIs(trainFile, args[4], args[5], args[6], Integer.parseInt(args[7]));
		}
			break;


		
		/***
		 * Wrapper aggregation sections. Some utils for aggregating implicit datasets (when the user consumes the same item more than once)
		 */
		case "AggregateWithWrapper": {
			System.out.println(Arrays.toString(args));
			String newDataset = cl.getOptionValue(OPT_NEW_DATASET);
			String wrapperStrategy = cl.getOptionValue(OPT_WRAPPER_STRATEGY);

			FastTemporalPreferenceDataIF<Long, Long> ranksysTrainTemporal = ExperimentUtils.loadTrainFastTemporalFeaturePreferenceData(trainFile, trainFile, false, true);
			// Create the wrapper to be able to work with the ranksys datamodel
			SimpleFastTemporalFeaturePreferenceDataWrapper<Long, Long> wrpTrain = new SimpleFastTemporalFeaturePreferenceDataWrapper<>(
					ranksysTrainTemporal, ExperimentUtils.obtRepetitionsStrategy(wrapperStrategy));

			wrpTrain.writePreferences(newDataset);
		}
			break;

		case "AggregateWithWrapperTimeStamps": {
			System.out.println(Arrays.toString(args));
			String newDataset = cl.getOptionValue(OPT_NEW_DATASET);
			String wrapperStrategyTime = cl.getOptionValue(OPT_WRAPPER_STRATEGY_TIMESTAMPS);
			String wrapperStrategyPref = cl.getOptionValue(OPT_WRAPPER_STRATEGY);

			boolean completeOrNot = Boolean.parseBoolean(cl.getOptionValue(OPT_COMPLETE_INDEXES, DEFAULT_COMPLETE_INDEXES));
			
			FastTemporalPreferenceDataIF<Long, Long> ranksysTrainTemporal = ExperimentUtils.loadTrainFastTemporalFeaturePreferenceData(trainFile, trainFile, completeOrNot, true);
			// Create the wrapper to be able to work with the ranksys datamodel
			SimpleFastTemporalFeaturePreferenceDataWrapper<Long, Long> wrpTrain = new SimpleFastTemporalFeaturePreferenceDataWrapper<Long, Long>(
					ranksysTrainTemporal, ExperimentUtils.obtRepetitionsStrategy(wrapperStrategyPref), ExperimentUtils.obtRepetitionsStrategyTime(wrapperStrategyTime));

			wrpTrain.writePreferencesTimeStamps(newDataset);
		}
			break;
			
		
		

		
		case "ParseMyMediaLite": {
			System.out.println("-o ParseMyMediaLite -trf myMediaLiteRecommendation testFile newRecommendation");
			System.out.println(Arrays.toString(args));
			Tuple4<List<Long>, List<Long>, List<Long>, List<Long>> indexes = ExperimentUtils
					.retrieveTrainTestIndexes(args[4], args[4], false, lp, lp);

			List<Long> usersTest = indexes.v1;
			FastUserIndex<Long> userIndexTest = SimpleFastUserIndex.load(usersTest.stream());

			ProcessData.parseMyMediaLite(trainFile, userIndexTest, args[5]);
		}
			break;	
			

		

		case "CheckRecommendationsFile": {
			System.out.println(Arrays.toString(args));
			String featureFile = cl.getOptionValue(OPT_FEATURE_READING_FILE);
			String recFile = cl.getOptionValue(OPT_RECOMMENDED_FILE);
			String matchCat = cl.getOptionValue(OPT_MATCHING_CATEGORY);
			String testFile = cl.getOptionValue(OPT_TEST_FILE);
			POIProcessData.checkRecommendationCorrection(trainFile, featureFile, recFile, matchCat, testFile);
		}
			break;
			


		/**
		 * RankSys evaluation section or reranking
		 */		
		// Ranksys with non accuracy metrics
		case "ranksysNonAccuracyWithoutFeatureMetricsEvaluation":
		case "ranksysNonAccuracyMetricsEvaluation":
		case "ranksysNonAccuracyMetricsEvaluationPerUser":
		case "ranksysNonAccuracyWithoutFeatureMetricsEvaluationPerUser": {
			/*
			 * -Train file -Test file -Recommended file -Item feature file -Ranksys Metric
			 * -Output file -Threshold -Cutoff
			 */
			System.out.println(Arrays.toString(args));

			String recommendedFile = cl.getOptionValue(OPT_RECOMMENDED_FILE);
			String itemFeatureFile = cl.getOptionValue(OPT_FEATURE_READING_FILE);
			String userFeaturFile = cl.getOptionValue(OPT_USER_FEATURE_FILE);
			String userSelFeature = cl.getOptionValue(OPT_USER_FEATURE_SELECTED);
			Double confidenceSim = Double.parseDouble(cl.getOptionValue(OPT_CONFIDENCE, "0"));

			int threshold = Integer.parseInt(cl.getOptionValue(OPT_THRESHOLD));
			int antiRelevanceThreshold = Integer.parseInt(cl.getOptionValue(OPT_ANTI_RELEVANCE_THRESHOLD, "0"));
			int numberBins = Integer.parseInt(cl.getOptionValue(OPT_NUMBER_BINS, "0"));

			
			String cutoffs = cl.getOptionValue(OPT_CUTOFF);
			String overwrite = cl.getOptionValue(OPT_OVERWRITE, DEFAULT_OVERWRITE);
			String userWeightModel = cl.getOptionValue(OPT_USER_MODEL_WGT, DEFAULT_USER_MODEL_WGT);

			Boolean computeAntiMetrics = Boolean
					.parseBoolean(cl.getOptionValue(OPT_COMPUTE_ANTI_METRICS, DEFAULT_COMPUTE_ANTI_METRICS));
			Boolean computeOnlyAcc = Boolean
					.parseBoolean(cl.getOptionValue(OPT_COMPUTE_ONLY_ACC, DEFAULT_OPT_COMPUTE_ONLY_ACC));

			Boolean computeNonExactMetric = Boolean
					.parseBoolean(cl.getOptionValue(OPT_NOT_EXACT_METRIC, DEFAULT_OPT_NOT_EXACT_METRIC));
			Boolean computeUserFilter = Boolean
					.parseBoolean(cl.getOptionValue(OPT_COMPUTE_USER_FILTER, DEFAULT_COMPUTE_USER_FILTER));
			
			Boolean computePopEvMetrics = Boolean.parseBoolean(cl.getOptionValue(OPT_POP_EV_METRICS, "false"));

			String ranksysRelevanceModel = cl.getOptionValue(OPT_RANKSYS_RELEVANCE_MODEL);
			String ranksysDiscountModel = cl.getOptionValue(OPT_RANKSYS_DISCOUNT_MODEL);
			String ranksysBackground = cl.getOptionValue(OPT_RANKSYS_BACKGROUND);
			String ranksysBase = cl.getOptionValue(OPT_RANKSYS_BASE);

			String itemsTimeAvgFile = cl.getOptionValue(OPT_TIMES_AVERAGE_ITEMS);
			String itemsTimeFirstFile = cl.getOptionValue(OPT_TIMES_FIRST_ITEMS);
			String itemsTimeMedianFile = cl.getOptionValue(OPT_TIMES_MEDIAN_ITEMS);
			String itemsTimeLastFile = cl.getOptionValue(OPT_TIMES_LAST_ITEMS);
			String itemsReleaseDate = cl.getOptionValue(OPT_RELEASE_DATE_ITEMS);
			String mapCoordinates = cl.getOptionValue(OPT_COORD_FILE);
			String itemSim = cl.getOptionValue(OPT_RANKSYS_SIM);

			Boolean lcsEvaluation = Boolean
					.parseBoolean(cl.getOptionValue(OPT_LCS_EVALUATION, DEFAULT_LCS_EVALUTATION));
			Boolean computeRatios = Boolean
					.parseBoolean(cl.getOptionValue(OPT_COMPUTE_RATIOS, DEFAULT_OPT_COMPUTE_RATIOS));
			Boolean computeDistances = Boolean
					.parseBoolean(cl.getOptionValue(OPT_COMP_DISTANCES, DEFAULT_COMPDISTANCES));

			// we can filter out from test:
			// users with less than n ratings in training
			int minUserRatings = Integer.parseInt(cl.getOptionValue(OPT_MIN_USERS, "0"));
			// items with less than m ratings in training
			int minItemRatings = Integer.parseInt(cl.getOptionValue(OPT_MIN_ITEMS, "0"));

			Boolean isPerUser = step.contains("PerUser");

			// This ones can be more than one file
			String outputFileS = cl.getOptionValue(OPT_OUT_RESULT_FILE);
			String testFileS = cl.getOptionValue(OPT_TEST_FILE);

			//output and test files need to be separated by -
			String outputFileArr[] = outputFileS.split(",");
			String testFileArr[] = testFileS.split(",");
			
			if (outputFileArr.length == 1) {
				File f = new File(outputFileArr[0]);
				
				// If file of ranksys evaluation already exist then nothing
				if (f.exists() && !f.isDirectory() && Boolean.parseBoolean(overwrite) == false) {
					System.out.println("Ignoring " + f + " because it already exists");
					break;
				}
			}



			if (isPerUser) {
				System.out.println("Per user evaluation. We will compute the metrics for every user indvidually.");
			} else {
				System.out.println("Normal evaluation, computing the average of the results in the metrics.");
			}

			// This section is to use a temporal preference data in the test set, for a time
			// aware metric
			FastTemporalPreferenceDataIF<Long, Long> rankSysTemporalTestData = null;

			final PreferenceData<Long, Long> trainData = SimplePreferenceData
					.load(SimpleRatingPreferencesReader.get().read(trainFile, lp, lp));

			final FastPreferenceData<Long, Long> trainDataFast = ExperimentUtils
					.loadTrainFastPreferenceData(trainFile, trainFile, false, true);
			
			final PreferenceData<Long, Long> originalRecommendedData = SimplePreferenceData
					.load(SimpleRatingPreferencesReader.get().read(recommendedFile, lp, lp));
			
			
			
			FeatureData<Long, String, Double> featureData = null;
			ItemDistanceModel<Long> dist = null;
			if (!computeOnlyAcc && itemFeatureFile != null) {
				featureData = SimpleFeatureData.load(SimpleFeaturesReader.get().read(itemFeatureFile, lp, sp));
				// COSINE DISTANCE
				dist = new CosineFeatureItemDistanceModel<>(featureData);
			}
			
			double ranksysBackgroundD = 0.0;
			double ranksysBaseD = 0.0;

			if (ranksysBackground != null) {
				ranksysBackgroundD = Double.parseDouble(ranksysBackground);
			}

			if (ranksysBase != null) {
				ranksysBaseD = Double.parseDouble(ranksysBase);
			}
			
			RankingDiscountModel discModel = null;
			if (ranksysDiscountModel == null) {
				discModel = new NoDiscountModel();
			} else {
				discModel = SequentialRecommendersUtils.obtRankingDiscountModel(ranksysDiscountModel, ranksysBaseD);
			}
			
			//Obtain the similarities for non exact
			
			List<ItemSimilarity<Long>> itemSims = new ArrayList<>();
			List<String> nameSims = new ArrayList<>();
			
			if (computeNonExactMetric) {
				//Read full file
				String fullFileTrainTest = cl.getOptionValue(OPT_FULL_FILE);
				FastPreferenceData<Long, Long> fullDataset = ExperimentUtils.loadTrainFastPreferenceData(fullFileTrainTest, fullFileTrainTest, false, true);
				
				String [] pathsSims = itemSim.split(","); 
				for (String filePath: pathsSims) {
					File f = new File(filePath);
					String nameSim = f.getName();
					ItemIndexReadingSimilarity<Long> isim =  new ItemIndexReadingSimilarity<>(filePath, confidenceSim, lp, fullDataset);
					ItemSimilarity<Long> simItemSim = new ItemSim<>(fullDataset, isim);
					itemSims.add(simItemSim);
					nameSims.add(nameSim);
					System.out.println(nameSim+ " " + " having " + simItemSim.numItems());
				}
			}
			
			for (int i = 0; i < outputFileArr.length; i++) {
				String outputFile = outputFileArr[i];
				String testFile = testFileArr[i];
				
				File f = new File(outputFile);
				
				// If file of ranksys evaluation already exist then nothing
				if (f.exists() && !f.isDirectory() && Boolean.parseBoolean(overwrite) == false) {
					System.out.println("Ignoring " + f + " because it already exists");
					continue;
				}

				if (lcsEvaluation) {
					System.out.println("Reading test set with timestamps for LCS evaluation");

					rankSysTemporalTestData = ExperimentUtils.loadTrainFastTemporalFeaturePreferenceData(trainFile,
							testFile, true, false);
				}

				// End of test Temporal data

				
				
				final PreferenceData<Long, Long> testDataNoFilter = SimplePreferenceData
						.load(SimpleRatingPreferencesReader.get().read(testFile, lp, lp));
				
				final PreferenceData<Long, Long> testData = SequentialRecommendersUtils
						.filterTestBasedOnTraining(trainData, testDataNoFilter, minUserRatings, minItemRatings);

				final PreferenceData<Long, Long> totalData = new ConcatPreferenceData<>(trainData, testData);

				final Set<Long> testUsers = testData.getUsersWithPreferences().collect(Collectors.toSet());
				
				// recommended data has to be filtered to avoid evaluating users not in test
				final PreferenceData<Long, Long> recommendedData = SequentialRecommendersUtils
						.filterPreferenceData(originalRecommendedData, testUsers, null);

				////////////////////////
				// INDIVIDUAL METRICS //
				////////////////////////
				TimestampCalculator<Long, Long> its = ExperimentUtils.createTimeStampCalculator(trainFile,
						itemsTimeAvgFile, itemsTimeFirstFile, itemsTimeMedianFile, itemsTimeLastFile);


				IntentModel<Long, Long, String> intentModel = null;
				if (!computeOnlyAcc && itemFeatureFile != null) {
					// INTENT MODEL
					intentModel = new FeatureIntentModel<>(totalData, featureData);
				}

				// Binary relevance and anti relevance model for ranking metrics
				BinaryRelevanceModel<Long, Long> binRel = new BinaryRelevanceModel<>(false, testData, threshold);

				// Anti relevance model
				BinaryAntiRelevanceModel<Long, Long> selectedAntiRelevance = new BinaryAntiRelevanceModel<>(false,
						testData, antiRelevanceThreshold);
				
				// Relevance model for novelty/diversity (can be with or without relevance)
				RelevanceModel<Long, Long> selectedRelevance = null;

				if (ranksysRelevanceModel == null) {
					selectedRelevance = new NoRelevanceModel<>();
				} else {
					selectedRelevance = SequentialRecommendersUtils.obtRelevanceModelRanksys(ranksysRelevanceModel,
							testData, threshold, ranksysBackgroundD);
				}

				

				int numUsersTest = testData.numUsersWithPreferences();
				int numUsersRecommended = recommendedData.numUsersWithPreferences();
				int numItems = totalData.numItemsWithPreferences(); // Num items with preferences in the data
				int numItemsTr = trainData.numItemsWithPreferences(); // Num items with preferences in the trainingData

				
				System.out.println("\n\nNum users in the test set " + numUsersTest);
				System.out.println("\n\nNum users to whom we have made recommendations " + numUsersRecommended);
				System.out.println("Modified ratios normalization");

				Map<String, SystemMetric<Long, Long>> sysMetrics = new HashMap<>();
				// Ranking metrics (avg for recommendation)
				Map<String, RecommendationMetric<Long, Long>> recMetricsAvgRelUsers = new HashMap<>();
				Map<String, RecommendationMetric<Long, Long>> recMetricsAvgAntiRelUsers = new HashMap<>();
				Map<String, RecommendationMetric<Long, Long>> recMetricsAllRecUsers = new HashMap<>();
				Map<String, RecommendationMetric<Long, Long>> recMetricsAllRecUsersForStd = new HashMap<>();

				String[] differentCutoffs = cutoffs.split(",");
				
				
				Map<Long, Double> prefsItem = null;
				Map<Long, Double> usersItem = null;
				Map<Long, Double> popByRatingItems = null;
				Map<String, Double> featuresCounterItems = null;
				Map<String, Double> featuresCounterPrefs = null;
				Map<String, Double> featuresCounterPopByRatingsItems = null;

				
				if (computePopEvMetrics) {
					prefsItem = ProcessData.numberPreferencesItems(trainData);
					usersItem = ProcessData.numberUsersRatedItems(trainData);
					popByRatingItems = ProcessData.popByRatingItems(trainData);
					
					if (featureData != null) {
						featuresCounterItems = ProcessData.featureCountsByItems(featureData, trainData);
						featuresCounterPrefs = ProcessData.featureCountsByPreferences(featureData, trainData);
						featuresCounterPopByRatingsItems = ProcessData.featureCountsPopByRatings(featureData, trainData);
					}
				}

				
				UserFeatureData<Long, String, Double> ufD = null;
				if (!isPerUser) {
					if (userFeaturFile != null) {
						ufD = SimpleUserFeatureData.load(SimpleUserFeaturesReader.get().read(userFeaturFile, lp, sp));
					}
				}
				UserFeatureData<Long, String, Double> ufD2 = ufD;
				
				
				for (String cutoffS : differentCutoffs) {
					int cutoff = Integer.parseInt(cutoffS);
					ExperimentUtils.addMetrics(recMetricsAvgRelUsers, recMetricsAvgAntiRelUsers, recMetricsAllRecUsers,
							recMetricsAllRecUsersForStd, threshold, antiRelevanceThreshold, cutoff, trainData,
							trainDataFast, testData, its, itemsTimeAvgFile, itemsTimeFirstFile, itemsTimeMedianFile,
							itemsTimeLastFile, selectedRelevance, binRel, selectedAntiRelevance, discModel, featureData,
							dist, intentModel, itemsReleaseDate, step, computeAntiMetrics, computeOnlyAcc,
							rankSysTemporalTestData, lcsEvaluation, computeRatios, computeDistances, mapCoordinates,
							computeNonExactMetric, itemSims, nameSims, computePopEvMetrics, 
							numberBins, prefsItem, usersItem, popByRatingItems, featuresCounterItems, featuresCounterPrefs, featuresCounterPopByRatingsItems);

					// SYSTEM METRICS. Only for normal evaluation, not per user

					if (!isPerUser) {
						
		

						sysMetrics.put("aggrdiv@" + cutoff, new AggregateDiversityMetric<>(cutoff, selectedRelevance));
						sysMetrics.put("gini@" + cutoff, new GiniIndex<>(cutoff, numItems));
						sysMetrics.put("gini_rel@" + cutoff, new GiniRelevantUsers<>(cutoff, numItems, binRel));
						sysMetrics.put("giniITr@" + cutoff, new GiniIndex<>(cutoff, numItemsTr));
						sysMetrics.put("giniITr_rel@" + cutoff, new GiniRelevantUsers<>(cutoff, numItemsTr, binRel));
						sysMetrics.put("RealAD@" + cutoff, new RealAggregateDiversity<Long, Long>(cutoff));
						sysMetrics.put("RealAD_rel@" + cutoff, new RealAggregateDiversityRelevantUsers<Long, Long>(cutoff, binRel));
						sysMetrics.put("RelUsers@" + cutoff, new RelevantUsersRecommended<Long, Long>(cutoff, binRel));
					}
				}
				
				//User cov goes without cutoff
				sysMetrics.put("usercov", new UserCoverage<Long, Long>());
				sysMetrics.put("usercov_rel", new UserCoverageRelevant<Long, Long>(binRel));

				// Average of all only for normal evaluation, not per user
				if (!isPerUser) {

					UserMetricWeight weight = ExperimentUtils.obtUserMetricWeight(userWeightModel);
					System.out.println("Working with user weight " + weight.toString());

					WeightedModelUser<Long, Long, String, Double> wmu = new WeightedModelUser<>(trainData, weight,
							userSelFeature, ufD);

					recMetricsAvgRelUsers.forEach((name, metric) -> sysMetrics.put(name + "_rec",
							new WeightedAverageRecommendationMetricIgnoreNoRelevantUsersAndNaNs<>(metric, binRel,
									wmu)));
					recMetricsAvgAntiRelUsers.forEach((name, metric) -> sysMetrics.put(name + "_rec",
							new WeightedAverageRecommendationMetricIgnoreNoRelevantUsersAndNaNs<>(metric,
									selectedAntiRelevance, wmu)));

					// Previous to 12 of April of 2019 we used this. Non accuracy divided by all
					// users recommended
					// Now using the same as the rest of the metrics
					/*
					 * recMetricsAllRecUsers.forEach((name, metric) -> sysMetrics.put(name + "_rec",
					 * new AverageRecommendationMetric<>(metric, numUsersRecommended)));
					 */

					recMetricsAllRecUsers.forEach((name, metric) -> sysMetrics.put(name + "_rec",
							new WeightedAverageRecommendationMetricIgnoreNoRelevantUsersAndNaNs<>(metric, binRel,
									wmu)));


					RecommendationFormat<Long, Long> format = new SimpleRecommendationFormat<>(lp, lp);

					/*
					 * format.getReader(recommendedFile).readAll() .forEach(rec ->
					 * sysMetrics.values().forEach(metric -> metric.add(rec)));
					 */

					System.out.println("Computing user filter: " + computeUserFilter);

					format.getReader(recommendedFile).readAll().filter(
							rec -> ExperimentUtils.isValidUser(rec.getUser(), computeUserFilter, ufD2, userSelFeature))
							.forEach(rec -> sysMetrics.values().forEach(metric -> metric.add(rec)));

					PrintStream out = new PrintStream(new File(outputFile));
					sysMetrics.forEach((name, metric) -> out.println(name + "\t" + metric.evaluate()));
					out.close();
				}

				// If it is per user we will just compute the per user evaluations with no _rec
				// nor _test and no user coverages
				else {
					Map<String, Map<Long, Double>> perUserEvaluations = new HashMap<String, Map<Long, Double>>();
					Map<String, RecommendationMetric<Long, Long>> recMetricstoSend = new HashMap<>();
					
					for (String metric : recMetricsAvgRelUsers.keySet()) {
						recMetricstoSend.put(metric, recMetricsAvgRelUsers.get(metric));
					}
					for (String metric : recMetricsAvgAntiRelUsers.keySet()) {
						recMetricstoSend.put(metric, recMetricsAvgAntiRelUsers.get(metric));
					}
					for (String metric : recMetricsAllRecUsers.keySet()) {
						recMetricstoSend.put(metric, recMetricsAllRecUsers.get(metric));
					}
					ExperimentUtils.evaluateRecommenderFile(recommendedFile, recMetricstoSend, perUserEvaluations);
					PrintStream out = new PrintStream(new File(outputFile));
					for (String s : perUserEvaluations.keySet()) {
						double avg = 0;
						for (Long user : perUserEvaluations.get(s).keySet()) {
							//If the user feature is null or it is not null but it matches the selected user feature
							if (ufD == null || (ufD != null && ufD.getUserFeatures(user).map(t -> t.v1).anyMatch(uf -> uf.equals(userSelFeature)))) {
								avg += perUserEvaluations.get(s).get(user);
								out.println(s + "\t" + user + "\t" + perUserEvaluations.get(s).get(user));
							}
						}
						out.println(s + "\t" + "all" + "\t" + avg / perUserEvaluations.get(s).keySet().size());
					}
				}
			}

		}
			break;

		default:
			System.out.println("Option " + step + " not recognized");
			break;
		}
	}

	/**
	 * Method that will obtain a command line using the available options of the
	 * arguments
	 *
	 * @param args
	 *            the arguments that the program will receive
	 * @return a command line if the arguments are correct or null if an error
	 *         occurred
	 */
	private static CommandLine getCommandLine(String[] args) {
		Options options = new Options();

		// Number of the case
		Option caseIdentifier = new Option("o", OPT_CASE, true, "option of the case");
		caseIdentifier.setRequired(true);
		options.addOption(caseIdentifier);

		// Train file
		Option trainFile = new Option("trf", OPT_TRAIN_FILE, true, "input file train path");
		trainFile.setRequired(true);
		options.addOption(trainFile);

		// Train file (2)
		Option trainFile2 = new Option("trf2", OPT_TRAIN_FILE_2, true, "input file train path (2)");
		trainFile2.setRequired(false);
		options.addOption(trainFile2);
		
		// Here not required
		// TestFile file
		Option testFile = new Option("tsf", OPT_TEST_FILE, true, "input file test path");
		testFile.setRequired(false);
		options.addOption(testFile);

		// TestFile file (test file full, only for Foursqr when using cross domain)
		Option testFile2 = new Option("tsfFull", OPT_TEST_FILE_2, true, "input file test path containing all test Files");
		testFile2.setRequired(false);
		options.addOption(testFile2);

		// Similarity file
		Option similarityFile = new Option("sf", OPT_SIMFILE, true, "similarity file");
		similarityFile.setRequired(false);
		options.addOption(similarityFile);

		// Similarity file 2
		Option similarityFile2 = new Option("sf2", OPT_SIMFILE_2, true, "similarity file");
		similarityFile2.setRequired(false);
		options.addOption(similarityFile2);

		// Actors file
		Option actorsFile = new Option("af", OPT_ACTORSFILE, true, "actors file");
		actorsFile.setRequired(false);
		options.addOption(actorsFile);

		// Directors file
		Option directorsFile = new Option("df", OPT_DIRECTORSFILE, true, "directors file");
		directorsFile.setRequired(false);
		options.addOption(directorsFile);

		// Genres file
		Option genresFile = new Option("gf", OPT_GENRESFILE, true, "genres file");
		genresFile.setRequired(false);
		options.addOption(genresFile);

		// Tags file
		Option tagsfile = new Option("tf", OPT_TAGSFILE, true, "tags file");
		tagsfile.setRequired(false);
		options.addOption(tagsfile);

		// LCS normalization: true or false
		Option lcsNormalization = new Option("lcsn", OPT_LCS_NORMALIZATION, true, "lcs normalization");
		lcsNormalization.setRequired(false);
		options.addOption(lcsNormalization);

		// Neighbours
		Option neighbours = new Option("n", OPT_NEIGH, true, "neighbours");
		neighbours.setRequired(false);
		options.addOption(neighbours);

		// TrainModel
		Option lcsTrainModel = new Option("lcstm", OPT_LCS_TRAIN_MODEL, true, "lcs training model");
		lcsTrainModel.setRequired(false);
		options.addOption(lcsTrainModel);

		// threshold
		Option threshold = new Option("thr", OPT_THRESHOLD, true, "relevance or matching threshold");
		threshold.setRequired(false);
		options.addOption(threshold);
		
		Option matchingString = new Option("matchString", OPT_MATCHING_STRING, true, "Matching string");
		matchingString.setRequired(false);
		options.addOption(matchingString);

		// NumberItemsRecommended
		Option numberItemsRecommended = new Option("nI", OPT_ITEMS_RECOMMENDED, true, "Number of items recommended");
		numberItemsRecommended.setRequired(false);
		options.addOption(numberItemsRecommended);

		// numberMinRatItems
		Option numberMinRatItems = new Option("mri", OPT_MIN_ITEMS, true, "Minimum number of item ratings");
		numberMinRatItems.setRequired(false);
		options.addOption(numberMinRatItems);

		// numberMinRatUsers
		Option numberMinRatUsers = new Option("mru", OPT_MIN_USERS, true, "Minimum number of user ratings");
		numberMinRatUsers.setRequired(false);
		options.addOption(numberMinRatUsers);
		
		// numberMinRatItems
		Option numberMinUniqueRatItems = new Option("mUniqri", OPT_MIN_UNIQUE_RATINGS_ITEMS, true, "Minimum number of item unique ratings");
		numberMinUniqueRatItems.setRequired(false);
		options.addOption(numberMinUniqueRatItems);

		// numberMinRatUsers
		Option numberMinUniqueRatUsers = new Option("mUniqru", OPT_MIN_UNIQUE_RATINGS_USERS, true, "Minimum number of user unique ratings");
		numberMinUniqueRatUsers.setRequired(false);
		options.addOption(numberMinUniqueRatUsers);		

		// OutResultfile
		Option outfile = new Option("orf", OPT_OUT_RESULT_FILE, true, "output result file");
		outfile.setRequired(false);
		options.addOption(outfile);
		
		// OutResultfile 2
		Option outfile2 = new Option("orf2", OPT_OUT_RESULT_FILE_2, true, "output result file (2)");
		outfile2.setRequired(false);
		options.addOption(outfile2);

		// OutSimilarityfile
		Option outSimfile = new Option("osf", OPT_OUT_SIM_FILE, true, "output similarity file");
		outSimfile.setRequired(false);
		options.addOption(outSimfile);

		// Ranksys similarity
		Option rankSysSim = new Option("rs", OPT_RANKSYS_SIM, true, "ranksys similarity");
		rankSysSim.setRequired(false);
		options.addOption(rankSysSim);

		// Ranksys similarity 2
		Option rankSysSim2 = new Option("rs2", OPT_RANKSYS_SIM_2, true, "ranksys similarity 2");
		rankSysSim2.setRequired(false);
		options.addOption(rankSysSim2);

		// Ranksys recommender
		Option rankSysRecommender = new Option("rr", OPT_RANKSYS_REC, true, "ranksys recommeder");
		rankSysRecommender.setRequired(false);
		options.addOption(rankSysRecommender);
		
		// Ranksys recommender 2
		Option rankSysRecommender2 = new Option("rr2", OPT_RANKSYS_REC2, true, "additional ranksys recommeder");
		rankSysRecommender2.setRequired(false);
		options.addOption(rankSysRecommender2);

		// Overwrite result
		Option outputOverwrite = new Option("ovw", OPT_OVERWRITE, true, "overwrite");
		outputOverwrite.setRequired(false);
		options.addOption(outputOverwrite);

		// ItemColFeature
		Option itemColFeature = new Option("fic", OPT_FEATURES_ITEM_COL, true, "item column of feature file");
		itemColFeature.setRequired(false);
		options.addOption(itemColFeature);

		// FeatColFeature
		Option featColFeature = new Option("ffc", OPT_FEATURES_FEAT_COL, true, "feature column of feature file");
		featColFeature.setRequired(false);
		options.addOption(featColFeature);
		
		// metricColumn
		Option metricColumn = new Option("mcol", OPT_METRIC_COLUMN, true, "metric column");
		metricColumn.setRequired(false);
		options.addOption(metricColumn);
		
		//metricSelected
		Option metricSelected = new Option("metricSel", OPT_METRIC_SELECTED, true, "metric selected");
		metricSelected.setRequired(false);
		options.addOption(metricSelected);
		
		//resultColumn
		Option resultColumn = new Option("rescol", OPT_RESULT_COLUMN, true, "result column");
		resultColumn.setRequired(false);
		options.addOption(resultColumn);
		
		//resultColumn
		Option recommenderColumn = new Option("reccol", OPT_RECOMMENDER_COLUMN, true, "recommender column");
		recommenderColumn.setRequired(false);
		options.addOption(recommenderColumn);

		// sepFeature
		Option sepFeature = new Option("fs", OPT_FEATURES_SEP, true, "separator of feature file");
		sepFeature.setRequired(false);
		options.addOption(sepFeature);

		// headerFeature
		Option headerFeature = new Option("fh", OPT_FEATURES_HEADER, true, "does the feature file has a header?");
		headerFeature.setRequired(false);
		options.addOption(headerFeature);

		// SavingFeature
		Option savingFeature = new Option("sfs", OPT_SAVING_FEATURES_SCHEME, true, "saving feature scheme");
		savingFeature.setRequired(false);
		options.addOption(savingFeature);

		// Actors saving file
		Option actorsSavingPath = new Option("asf", OPT_ACTORS_SAVING_FILE, true, "actors saving file");
		actorsSavingPath.setRequired(false);
		options.addOption(actorsSavingPath);

		// Directors saving file
		Option directorsSavingPath = new Option("dsf", OPT_DIRECTORSSAVINGFILE, true, "directors saving file");
		directorsSavingPath.setRequired(false);
		options.addOption(directorsSavingPath);

		// Genres saving file
		Option genresSavingPath = new Option("gsf", OPT_GENRESSAVINGFILE, true, "genres saving file");
		genresSavingPath.setRequired(false);
		options.addOption(genresSavingPath);

		// Feature reading file
		Option featureReadingFile = new Option("ff", OPT_FEATURE_READING_FILE, true, "features file");
		featureReadingFile.setRequired(false);
		options.addOption(featureReadingFile);

		// Feature reading file2
		Option featureReadingFile2 = new Option("ff2", OPT_FEATURE_READING_FILE2, true, "features file (2)");
		featureReadingFile2.setRequired(false);
		options.addOption(featureReadingFile2);

		// First recommended file
		Option recommendedFile = new Option("rf", OPT_RECOMMENDED_FILE, true, "recommended file");
		recommendedFile.setRequired(false);
		options.addOption(recommendedFile);
		
		// Second recommended file
		Option recommendedFile2 = new Option("rf2", OPT_RECOMMENDED_FILE2, true, "recommended file 2");
		recommendedFile2.setRequired(false);
		options.addOption(recommendedFile2);

		// RankSysMetric
		Option ranksysCutoff = new Option("rc", OPT_CUTOFF, true, "ranksyscutoff");
		ranksysCutoff.setRequired(false);
		options.addOption(ranksysCutoff);

		// Confidence similarity (if similarity is lower than this value, the similarity
		// will consider )
		Option confidenceSimilarity = new Option("cs", OPT_CONFIDENCE, true, "confidence similarity");
		confidenceSimilarity.setRequired(false);
		options.addOption(confidenceSimilarity);

		// Preference rating (ratings with less that this value will not be consider)
		Option preferenceRating = new Option("pr", OPT_PREFERENCE_RATING, true, "preference rating");
		preferenceRating.setRequired(false);
		options.addOption(preferenceRating);

		// Recommenders files
		Option recommendedFiles = new Option("rfs", OPT_RECOMMENDED_FILES, true, "recommendenders files");
		recommendedFiles.setRequired(false);
		options.addOption(recommendedFiles);
		
		// Recommenders selected
		Option recommendersSelected = new Option("recSel", OPT_RECOMMENDERS_SEL, true, "recommendenders selected");
		recommendersSelected.setRequired(false);
		options.addOption(recommendersSelected);
		
		// Confidence estimation
		Option confidenceEstimation = new Option("ce", OPT_CONFIDENCE_ESTIMATION, true, "confidence estimation");
		confidenceEstimation.setRequired(false);
		options.addOption(confidenceEstimation);

		// Identification factor
		Option identificatorFactor = new Option("idLCSf", OPT_IDENTIFICATION_FACTOR, true, "identification factor");
		identificatorFactor.setRequired(false);
		options.addOption(identificatorFactor);

		// Identification factor
		Option ratingFactor = new Option("rLCSf", OPT_RATING_FACTOR, true, "rating factor");
		ratingFactor.setRequired(false);
		options.addOption(ratingFactor);

		// IndexBackWards
		Option indexBackward = new Option("indexb", OPT_INDEX_BACKWARDS, true, "index backwards");
		indexBackward.setRequired(false);
		options.addOption(indexBackward);

		// IndexForwards
		Option indexForward = new Option("indexf", OPT_INDEX_FORWARDS, true, "index forwards");
		indexForward.setRequired(false);
		options.addOption(indexForward);

		// Feature CB SImilarity
		Option featureCBSim = new Option("featureCB", OPT_FEATURE_CB, true, "feature content based");
		featureCBSim.setRequired(false);
		options.addOption(featureCBSim);

		// CB User transformation
		Option cbUserTransformation = new Option("userTransformationCB", OPT_CB_USER_TRANSFORMATION, true,
				"cb user transformation");
		cbUserTransformation.setRequired(false);
		options.addOption(cbUserTransformation);

		// CB Item transformation
		Option cbItemTransformation = new Option("itemTransformationCB", OPT_CB_ITEM_TRANSFORMATION, true,
				"cb item transformation");
		cbItemTransformation.setRequired(false);
		options.addOption(cbItemTransformation);

		// RankSys factorizers (k)
		Option kFactorizer = new Option("kFactorizer", OPT_K_FACTORIZER, true, "k factorizer");
		kFactorizer.setRequired(false);
		options.addOption(kFactorizer);

		// RankSys factorizers (alpha)
		Option alphaFactorizer = new Option("aFactorizer", OPT_ALPHA_FACTORIZER, true, "alpha factorizer");
		alphaFactorizer.setRequired(false);
		options.addOption(alphaFactorizer);
		
		// Factorizers (theta)
		Option thetaFactorizer = new Option("thetaFactorizer", OPT_THETA_FACTORIZER, true, "theta factorizer");
		thetaFactorizer.setRequired(false);
		options.addOption(thetaFactorizer);
		
		// factorizers (alpha2)
		Option alphaFactorizer2 = new Option("aFactorizer2", OPT_ALPHA_FACTORIZER2, true, "alpha factorizer (2)");
		alphaFactorizer2.setRequired(false);
		options.addOption(alphaFactorizer2);
		

		// RankSys factorizers (lambda)
		Option lambdaFactorizer = new Option("lFactorizer", OPT_LAMBDA_FACTORIZER, true, "lambda factorizer");
		lambdaFactorizer.setRequired(false);
		options.addOption(lambdaFactorizer);

		// RankSys factorizers (numInteractions)
		Option numInteractionsFact = new Option("nIFactorizer", OPT_NUM_INTERACTIONS, true,
				"numInteractions factorizer");
		numInteractionsFact.setRequired(false);
		options.addOption(numInteractionsFact);

		// IndexesFile
		Option inindexesFile = new Option("iindexFile", OPT_IN_INDEX_FILE, true, "indexes file (input)");
		inindexesFile.setRequired(false);
		options.addOption(inindexesFile);

		// Lambda parameter for ranksys Time recommender
		Option temporalLambda = new Option("tempLambda", OPT_TEMPORAL_LAMBDA, true, "temporal lambda");
		temporalLambda.setRequired(false);
		options.addOption(temporalLambda);

		// RankLibrary norm value
		Option normAggregate = new Option("normAgLib", OPT_NORM_AGGREGATE_LIBRARY, true,
				"normalization aggregate library");
		normAggregate.setRequired(false);
		options.addOption(normAggregate);

		// RankLibrary comb value
		Option combAggregate = new Option("combAgLib", OPT_COMB_AGGREGATE_LIBRARY, true, "combination aggregate library");
		combAggregate.setRequired(false);
		options.addOption(combAggregate);

		// RankLibrary norm value
		Option weightAggregate = new Option("weightAgLib", OPT_WEIGHT_AGGREGATE_LIBRARY, true,
				"weight aggregate library");
		weightAggregate.setRequired(false);
		options.addOption(weightAggregate);

		// RansksyLibrary NonAccuracyEvaluationParameters
		Option ranksysRelevanceModel = new Option("ranksysRelModel", OPT_RANKSYS_RELEVANCE_MODEL, true,
				"ranksys relevance model");
		ranksysRelevanceModel.setRequired(false);
		options.addOption(ranksysRelevanceModel);

		Option ranksysDiscountModel = new Option("ranksysDiscModel", OPT_RANKSYS_DISCOUNT_MODEL, true,
				"ranksys discount model");
		ranksysDiscountModel.setRequired(false);
		options.addOption(ranksysDiscountModel);

		Option ranksysBackground = new Option("ranksysBackRel", OPT_RANKSYS_BACKGROUND, true,
				"ranksys background for relevance model");
		ranksysBackground.setRequired(false);
		options.addOption(ranksysBackground);

		Option ranksysBase = new Option("ranksysBaseDisc", OPT_RANKSYS_BASE, true, "ranksys base for discount model");
		ranksysBase.setRequired(false);
		options.addOption(ranksysBase);

		// Freshness/novelty metrics
		Option itemsTimeAverage = new Option("ItimeAvgFile", OPT_TIMES_AVERAGE_ITEMS, true, "average item time file");
		itemsTimeAverage.setRequired(false);
		options.addOption(itemsTimeAverage);

		Option itemsTimeFirst = new Option("ItimeFstFile", OPT_TIMES_FIRST_ITEMS, true, "first item time file");
		itemsTimeFirst.setRequired(false);
		options.addOption(itemsTimeFirst);

		Option itemsTimeLast = new Option("ItimeLastFile", OPT_TIMES_LAST_ITEMS, true, "last item time file");
		itemsTimeLast.setRequired(false);
		options.addOption(itemsTimeLast);

		Option itemsTimeMedian = new Option("ItimeMedFile", OPT_TIMES_MEDIAN_ITEMS, true,
				"median base for discount model");
		itemsTimeMedian.setRequired(false);
		options.addOption(itemsTimeMedian);

		Option itemsReleaseDate = new Option("IReleaseDate", OPT_RELEASE_DATE_ITEMS, true,
				"release dates file for items");
		itemsReleaseDate.setRequired(false);
		options.addOption(itemsReleaseDate);

		// File containing 2 columns, old id and new id, for datasets that we have
		// changed the ids (for items)
		Option itemsMapping = new Option("IMapping", OPT_MAPPING_ITEMS, true, "mapping file of items");
		itemsMapping.setRequired(false);
		options.addOption(itemsMapping);

		// File containing 2 columns, old id and new id, for datasets that we have
		// changed the ids
		Option usersMapping = new Option("UMapping", OPT_MAPPING_USERS, true, "mapping file of users");
		usersMapping.setRequired(false);
		options.addOption(usersMapping);

		// File containing 2 columns, old id and new id, for datasets that we have
		// changed the ids
		Option mapCategories = new Option("CMapping", OPT_MAPPING_CATEGORIES, true, "mapping file of categories");
		mapCategories.setRequired(false);
		options.addOption(mapCategories);

		// Option for creating a new dataset when transforming the users and items
		Option newDatasetFile = new Option("newDataset", OPT_NEW_DATASET, true, "new dataset destination");
		newDatasetFile.setRequired(false);
		options.addOption(newDatasetFile);

		// Option for creating a new dataset when transforming the users and items
		Option poiCoordFile = new Option("coordFile", OPT_COORD_FILE, true, "poi coordinates file");
		poiCoordFile.setRequired(false);
		options.addOption(poiCoordFile);
		
		Option poiEstimatedTime = new Option("poiEstimatedTimeFile", OPT_POI_STIMATED_IME, true, "poi estimated time file");
		poiEstimatedTime.setRequired(false);
		options.addOption(poiEstimatedTime);
		

		Option poiCityFile = new Option("poiCityFile", OPT_POI_CITY_FILE, true, "poi city mapping file");
		poiCityFile.setRequired(false);
		options.addOption(poiCityFile);

		Option citySimFile = new Option("citySimFile", OPT_CITY_SIM_FILE, true, "city similarity file");
		citySimFile.setRequired(false);
		options.addOption(citySimFile);

		// Option of the wrapper strategy to parse datasets or to use it in ranksys
		Option wrapperStrategyPreferences = new Option("wStrat", OPT_WRAPPER_STRATEGY, true,
				"strategy to use in the wrapper");
		wrapperStrategyPreferences.setRequired(false);
		options.addOption(wrapperStrategyPreferences);

		// Option of the wrapper strategy to parse datasets or to use it in ranksys
		Option wrapperStrategyTimeStamps = new Option("wStratTime", OPT_WRAPPER_STRATEGY_TIMESTAMPS, true,
				"strategy to use in the wrapper (for timestamps)");
		wrapperStrategyTimeStamps.setRequired(false);
		options.addOption(wrapperStrategyTimeStamps);

		Option matchingCategory = new Option("matchCat", OPT_MATCHING_CATEGORY, true,
				"matching category for candidate items");
		matchingCategory.setRequired(false);
		options.addOption(matchingCategory);

		Option usingWrapper = new Option("usingWrapper", OPT_USING_WRAPPER, true, "using wrapper in ranksys");
		usingWrapper.setRequired(false);
		options.addOption(usingWrapper);

		Option useCompleteIndex = new Option("cIndex", OPT_COMPLETE_INDEXES, true,
				"use the complete indexes (train + test)");
		useCompleteIndex.setRequired(false);
		options.addOption(useCompleteIndex);

		Option recommendationStrategy = new Option("recStrat", OPT_RECOMMENDATION_STRATEGY, true,
				"recommendation strategy");
		recommendationStrategy.setRequired(false);
		options.addOption(recommendationStrategy);

		Option regularizationl1 = new Option("reg1", OPT_REGULARIZATION_L1, true, "regularization L1");
		regularizationl1.setRequired(false);
		options.addOption(regularizationl1);

		Option regularizationl2 = new Option("reg2", OPT_REGULARIZATION_L2, true, "regularization L2");
		regularizationl2.setRequired(false);
		options.addOption(regularizationl2);

		Option svdBins = new Option("svdBins", OPT_SVD_NUMBER_BINS, true, "bins for svd++ recommender");
		svdBins.setRequired(false);
		options.addOption(svdBins);

		Option svdBeta = new Option("svdBeta", OPT_SVD_BETA, true, "beta for svd++ recommender");
		svdBeta.setRequired(false);
		options.addOption(svdBeta);

		Option svdTimeStrategy = new Option("svdTimeStrat", OPT_SVD_TIME_STRATEGY, true, "svd time strategy");
		svdTimeStrategy.setRequired(false);
		options.addOption(svdTimeStrategy);
		
		Option svdRegUser = new Option("svdRegUser", OPT_SVD_REG_USER, true, "svd regularization for users");
		svdRegUser.setRequired(false);
		options.addOption(svdRegUser);
		
		Option svdRegItem = new Option("svdRegItem", OPT_SVD_REG_ITEM, true, "svd regularization for items");
		svdRegItem.setRequired(false);
		options.addOption(svdRegItem);
		
		Option svdRegBias = new Option("svdRegBias", OPT_SVD_REG_BIAS, true, "svd regularization for biases");
		svdRegBias.setRequired(false);
		options.addOption(svdRegBias);

		Option svdLearnRate = new Option("svdLearnRate", OPT_SVD_LEARN_RATE, true, "svd learning rate");
		svdLearnRate.setRequired(false);
		options.addOption(svdLearnRate);
		
		Option svdMaxLearnRate = new Option("svdMaxLearnRate", OPT_SVD_MAX_LEARN_RATE, true, "svd maximum learning rate");
		svdMaxLearnRate.setRequired(false);
		options.addOption(svdMaxLearnRate);
		
		Option svdDecay = new Option("svdDecay", OPT_SVD_DECAY, true, "svd decay");
		svdDecay.setRequired(false);
		options.addOption(svdDecay);
		
		Option svdIsboldDriver = new Option("svdIsboldDriver", OPT_SVD_IS_BOLD_DRIVER, true, "svd bold driver");
		svdIsboldDriver.setRequired(false);
		options.addOption(svdIsboldDriver);
		
		Option svdRegImpItem = new Option("svdRegImpItem", OPT_SVD_REG_IMP_ITEM, true, "reg imp item");
		svdRegImpItem.setRequired(false);
		options.addOption(svdRegImpItem);
		
				
		Option markovChainOrder = new Option("markovChainOrder", OPT_MARKOV_CHAIN_ORDER, true, "order for MC");
		markovChainOrder.setRequired(false);
		options.addOption(markovChainOrder);
		
		Option patternLenght = new Option("patternLenght", OPT_PATTERN_LENGHT, true, "length of the pattern to search");
		patternLenght.setRequired(false);
		options.addOption(patternLenght);

		Option antiMetrics = new Option("computeAntiMetrics", OPT_COMPUTE_ANTI_METRICS, true,
				"boolean in order to compute anti metrics");
		antiMetrics.setRequired(false);
		options.addOption(antiMetrics);

		Option antiRelevanceThreshold = new Option("antiRelTh", OPT_ANTI_RELEVANCE_THRESHOLD, true,
				"integer in order to indicate the antiRelevanceTh");
		antiRelevanceThreshold.setRequired(false);
		options.addOption(antiRelevanceThreshold);
		
		Option computeOnlyAccuracy = new Option("onlyAcc", OPT_COMPUTE_ONLY_ACC, true, "if we are computing just accuracy metrics or also novelty-diversity");
		computeOnlyAccuracy.setRequired(false);
		options.addOption(computeOnlyAccuracy);
		
		Option computeRatios = new Option("compRatio", OPT_COMPUTE_RATIOS, true, "if we want to compute the ratios");
		computeRatios.setRequired(false);
		options.addOption(computeRatios);
		
		Option scoreFreq = new Option("scoreFreq", OPT_SCORE_FREQ, true, "simple or frequency");
		scoreFreq.setRequired(false);
		options.addOption(scoreFreq);
		
		Option temporalWeight = new Option("temporalWeight", OPT_TEMPORAL_WEIGHT, true, "temporal weight for UB temporal similarities");
		temporalWeight.setRequired(false);
		options.addOption(temporalWeight);
		
		Option lcsEv = new Option("lcsEv", OPT_LCS_EVALUATION, true, "LCS evaluation for taking into ccount the temporal order of the items");
		lcsEv.setRequired(false);
		options.addOption(lcsEv);
		
		Option limit = new Option("limitBetweenItems", OPT_LIMIT_BETWEEN_ITEMS, true, "limit of time or distance between items");
		limit.setRequired(false);
		options.addOption(limit);
		
		Option sessionType = new Option("sessionType", OPT_SESSION_TYPE, true, "session stype");
		sessionType.setRequired(false);
		options.addOption(sessionType);
		
		Option TimePopBeta = new Option("tpBeta", OPT_TIMEPOP_BETA, true, "time pop beta");
		TimePopBeta.setRequired(false);
		options.addOption(TimePopBeta);
		
		Option TimePopMinimumCandidates = new Option("tpMinimumCand", OPT_TIMEPOP_MINIMUM_CANDIDATES, true, "time pop minimumcandidates");
		TimePopMinimumCandidates.setRequired(false);
		options.addOption(TimePopMinimumCandidates);
		
		Option TimePopReferingTimeStamp = new Option("tpRefTimeStamp", OPT_TIMEPOP_REFERING_TIMESTAMP, true, "time pop referring timestamp");
		TimePopReferingTimeStamp.setRequired(false);
		options.addOption(TimePopReferingTimeStamp);
		
		Option maxDiffBetweenTimestamps = new Option("maxDiffTime", OPT_MAX_DIFF_TIME, true, "maxDifftime");
		maxDiffBetweenTimestamps.setRequired(false);
		options.addOption(maxDiffBetweenTimestamps);
		
		Option minDiffBetweenTimestamps = new Option("minDiffTime", OPT_MIN_DIFF_TIME, true, "minDifftime");
		minDiffBetweenTimestamps.setRequired(false);
		options.addOption(minDiffBetweenTimestamps);
		
		Option minDiffBetweenTimestamps2 = new Option("minDiffTime2", OPT_MIN_DIFF_TIME_2, true, "minDifftime2");
		minDiffBetweenTimestamps2.setRequired(false);
		options.addOption(minDiffBetweenTimestamps2);
		
		Option minClosePrefBot = new Option("minClosePrefBot", OPT_MIN_CLOSE_PREF_BOT, true, "minDifftime");
		minClosePrefBot.setRequired(false);
		options.addOption(minClosePrefBot);
		
		Option printUTC = new Option("printUTC", OPT_PRINT_UTC, true, "printUTC");
		printUTC.setRequired(false);
		options.addOption(printUTC);
		
		Option minTime = new Option("minTime", OPT_MIN_TIME, true, "mintime");
		minTime.setRequired(false);
		options.addOption(minTime);
		
		Option maxTime = new Option("maxTime", OPT_MAX_TIME, true, "maxtime");
		maxTime.setRequired(false);
		options.addOption(maxTime);
		
		Option poiCategoryFile = new Option("poiCat", OPT_POI_CATEGORY_FILE, true, "poiCat");
		poiCategoryFile.setRequired(false);
		options.addOption(poiCategoryFile);
			
		Option ignoreFirstLine = new Option("igFst", OPT_IGNORE_FIRST_LINE, true, "ignore");
		ignoreFirstLine.setRequired(false);
		options.addOption(ignoreFirstLine);
		
		Option rerankers = new Option("rerankers", OPT_RERANKERS, true, "renrankers");
		rerankers.setRequired(false);
		options.addOption(rerankers);
		
		Option socialInfluenceFile = new Option("socialInfluenceFile", OPT_SOC_INFLUENCE_FILE, true, "socialInfluenceFile");
		socialInfluenceFile.setRequired(false);
		options.addOption(socialInfluenceFile);
		
		Option eta = new Option("eta", OPT_ETA, true, "eta");
		eta.setRequired(false);
		options.addOption(eta);
		
		Option probabilitySmoothing = new Option("probSmooth", OPT_PROB_SMOOTH, true, "probability smoothing");
		probabilitySmoothing.setRequired(false);
		options.addOption(probabilitySmoothing);
		
		Option minSessions = new Option("minSessions", OPT_MIN_SESSIONS, true, "minimum number sessions");
		minSessions.setRequired(false);
		options.addOption(minSessions);
		
		Option computeDistances = new Option("compDistances", OPT_COMP_DISTANCES, true, "compute the distances of the paths");
		computeDistances.setRequired(false);
		options.addOption(computeDistances);
		
		Option normalize = new Option("normalize", OPT_NORMALIZE, true, "apply a normalization function");
		normalize.setRequired(false);
		options.addOption(normalize);
		
		Option lambda = new Option("lambda", OPT_LAMBDA, true, "lambda variable");
		lambda.setRequired(false);
		options.addOption(lambda);
		
		Option lambda2 = new Option("lambda2", OPT_LAMBDA_2, true, "lambda variable (2)");
		lambda2.setRequired(false);
		options.addOption(lambda2);
		
		Option cachedItems = new Option("cachedI", OPT_CACHED_ITEMS, true, "cachedItems");
		cachedItems.setRequired(false);
		options.addOption(cachedItems);
		
		Option numItemsCompute = new Option("nIComp", OPT_NUM_ITEMS_COMPUTE, true, "number of items to compute");
		numItemsCompute.setRequired(false);
		options.addOption(numItemsCompute);
		
		Option fillRerankingStrategy = new Option("fillStrat", OPT_FILL_STRATEGY, true, "fill reranking strategy");
		fillRerankingStrategy.setRequired(false);
		options.addOption(fillRerankingStrategy);
		
		Option citySelected = new Option("city", OPT_CITY, true, "city selection");
		citySelected.setRequired(false);
		options.addOption(citySelected);
		
		Option epsilon = new Option("epsilon", OPT_EPSILON, true, "epsilon");
		epsilon.setRequired(false);
		options.addOption(epsilon);
		
		Option c = new Option("c", OPT_C, true, "c");
		c.setRequired(false);
		options.addOption(c);
		
		Option evaluationFile = new Option("evFile", OPT_EVALUATION_FILE, true, "evFile");
		evaluationFile.setRequired(false);
		options.addOption(evaluationFile);
		
		Option userModelWeight = new Option("userModelW", OPT_USER_MODEL_WGT, true, "user model weight");
		userModelWeight.setRequired(false);
		options.addOption(userModelWeight);
		
		Option userFeatureFile = new Option("uff", OPT_USER_FEATURE_FILE, true, "user feature file");
		userFeatureFile.setRequired(false);
		options.addOption(userFeatureFile);
		
		Option userFeatureSelected = new Option("ufs", OPT_USER_FEATURE_SELECTED, true, "user feature selected");
		userFeatureSelected.setRequired(false);
		options.addOption(userFeatureSelected);
		
		Option itemFeatureSelected = new Option("ifs", OPT_ITEM_FEATURE_SELECTED, true, "item feature selected");
		itemFeatureSelected.setRequired(false);
		options.addOption(itemFeatureSelected);
		
		Option notExactMetric = new Option("nonEM", OPT_NOT_EXACT_METRIC, true, "not exact metric");
		notExactMetric.setRequired(false);
		options.addOption(notExactMetric);
		
		
		Option computeUserFilter = new Option("compUF", OPT_COMPUTE_USER_FILTER, true, "compute user filter");
		computeUserFilter.setRequired(false);
		options.addOption(computeUserFilter);
		
		Option maintainFirst = new Option("rerankMaintainFst", OPT_RERANK_MAINTAIN_FIRST, true, "maintain first item test");
		maintainFirst.setRequired(false);
		options.addOption(maintainFirst);
		
		Option hybridWeights = new Option("hybridWeights", OPT_HYBRID_WEIGHTS, true, "hybrid weights");
		hybridWeights.setRequired(false);
		options.addOption(hybridWeights);
		
		Option composedSimilarity = new Option("composedSim", OPT_COMPOSED_SIM, true, "composed similarity");
		composedSimilarity.setRequired(false);
		options.addOption(composedSimilarity);
		
		
		Option applyFilter = new Option("applyFilter", OPT_FILTER, true, "apply filter");
		applyFilter.setRequired(false);
		options.addOption(applyFilter);
		
		Option defaultScore = new Option("defScore", OPT_DEFAULT_SCORE, true, "default score");
		defaultScore.setRequired(false);
		options.addOption(defaultScore);
		
		Option randomSeed = new Option("randomSeed", OPT_RANDOM_SEED, true, "random seed");
		randomSeed.setRequired(false);
		options.addOption(randomSeed);
		
		Option percentage = new Option("percentage", OPT_PERCENTAGE, true, "percentage");
		percentage.setRequired(false);
		options.addOption(percentage);
		
		Option usePercentage = new Option("usePercentage", OPT_USE_PERCENTAGE, true, "use percentage");
		usePercentage.setRequired(false);
		options.addOption(usePercentage);
		
		Option numberBins = new Option("numberBins", OPT_NUMBER_BINS, true, "number bins");
		numberBins.setRequired(false);
		options.addOption(numberBins);
		
		Option popEv = new Option("popEv", OPT_POP_EV_METRICS, true, "popularity evaluation metrics");
		popEv.setRequired(false);
		options.addOption(popEv);
		
		Option maxDist = new Option("maxDist", OPT_MAX_DISTANCE, true, "maximum distance");
		maxDist.setRequired(false);
		options.addOption(maxDist);
		
		Option mode = new Option("mode", OPT_MODE, true, "mode of an algorithm");
		mode.setRequired(false);
		options.addOption(mode);
		
		Option useSigmoid = new Option("useSigmoid", OPT_USE_SIGMOID, true, "use sigmoid");
		useSigmoid.setRequired(false);
		options.addOption(useSigmoid);
		
		Option lstUserFeatures = new Option("lstUserFeat", OPT_LST_USER_FEATURES, true, "list user features");
		lstUserFeatures.setRequired(false);
		options.addOption(lstUserFeatures);
		
		Option differentTrainingsByUserFeatures = new Option("diffTrUF", OPT_DIFFERENT_TRAININGS_USER_FEAUTURES, true, "different trainings from user features");
		differentTrainingsByUserFeatures.setRequired(false);
		options.addOption(differentTrainingsByUserFeatures);
		
		Option startingIndexRankSim = new Option("sIndexRankSim", OPT_STARTING_INDEX_RANKSIM, true, "different trainings from user features");
		startingIndexRankSim.setRequired(false);
		options.addOption(startingIndexRankSim);
		
		Option aggregationTrajectorySim = new Option("aggrTrajecSim", OPT_AGGR_TRAJ_SIM, true, "aggregationTrajectorySims");
		aggregationTrajectorySim.setRequired(false);
		options.addOption(aggregationTrajectorySim);
		
		Option softenSimTrajectories = new Option("softenSim", OPT_SOFTEN_SIM_TRA, true, "softenSimTrayectories");
		softenSimTrajectories.setRequired(false);
		options.addOption(softenSimTrajectories);
		
		Option lastNSessionsPerUser = new Option("lastNSessions", OPT_LAST_N_SESSIONS, true, "softenSimTrayectories");
		lastNSessionsPerUser.setRequired(false);
		options.addOption(lastNSessionsPerUser);
		
		Option fullFile = new Option("fullFile", OPT_FULL_FILE, true, "fullFile");
		fullFile.setRequired(false);
		options.addOption(fullFile);
		
		Option maxPopItems = new Option("maxPopItems", OPT_MAX_POP_ITEMS, true, "maxPopItems");
		maxPopItems.setRequired(false);
		options.addOption(maxPopItems);
		
		Option consecutiveSameCheckins = new Option("consecutiveSameCheckins", OPT_CONSECUTIVE_SAME_CHECKINS, true, "consecutiveSameCheckins");
		consecutiveSameCheckins.setRequired(false);
		options.addOption(consecutiveSameCheckins);
		
		Option maxSpeed = new Option("maxSpeed", OPT_MAX_SPEED_CHECKINS, true, "maxSpeed");
		maxSpeed.setRequired(false);
		options.addOption(maxSpeed);
		
		Option inverse = new Option("inverse", OPT_INVERSE, true, "inverse");
		inverse.setRequired(false);
		options.addOption(inverse);
		
		Option negAntiN = new Option("negAntiN", OPT_NEGANTI, true, "negative anti-neigh");
		negAntiN.setRequired(false);
		options.addOption(negAntiN);
		
		Option numberRuns = new Option("numberRuns", OPT_NUMBERRUNS, true, "number of runs");
		numberRuns.setRequired(false);
		options.addOption(numberRuns);
		
		Option featureDistribution = new Option("featDistr", OPT_FEATURE_DISTRIBUTION, true, "feature distribution");
		featureDistribution.setRequired(false);
		options.addOption(featureDistribution);
		
		Option listNumbers = new Option("listNumbers", OPT_LIST_NUMBERS, true, "list of numbers");
		listNumbers.setRequired(false);
		options.addOption(listNumbers);
		
		Option perUser = new Option("perUser", OPT_PERUSER, true, "perUser");
		perUser.setRequired(false);
		options.addOption(perUser);
		
		Option percentageIncrementDouble = new Option("perIncr", OPT_PERCENTAGE_INCREMENT, true, "percentageIncrement");
		percentageIncrementDouble.setRequired(false);
		options.addOption(percentageIncrementDouble);
		
		Option reverse = new Option("reverse", OPT_REVERSE, true, "reverse");
		reverse.setRequired(false);
		options.addOption(reverse);
		
		Option imputeAllUsers = new Option("imputeAllUsers", OPT_IMPUTE_ALL_USERS, true, "imputeAllUsers");
		imputeAllUsers.setRequired(false);
		options.addOption(imputeAllUsers);
		
		Option imputeUsersTest = new Option("imputeUsersTest", OPT_IMPUTE_USERS_TESTS, true, "imputeUsersTest");
		imputeUsersTest.setRequired(false);
		options.addOption(imputeUsersTest);
		
		Option binaryImputation = new Option("binaryImputation", OPT_BINARY_IMPUTATION, true, "binaryImputation");
		binaryImputation.setRequired(false);
		options.addOption(binaryImputation);
		
		Option randomUsersImputation = new Option("randomUsersImputation", OPT_RANDOM_USERS_IMPUTATION, true, "randomImputation");
		randomUsersImputation.setRequired(false);
		options.addOption(randomUsersImputation);
		
		Option useGroupingImputation = new Option("useGroupingImp", OPT_GROUPING_IMPUTATION, true, "useGroupingImp");
		useGroupingImputation.setRequired(false);
		options.addOption(useGroupingImputation);
		
		Option sizeGroupImputation = new Option("sizeGroupImputation", OPT_SIZEGROUPING_IMPUTATION, true, "sizeGroupImputation");
		sizeGroupImputation.setRequired(false);
		options.addOption(sizeGroupImputation);
		
		
		CommandLineParser parser = new DefaultParser();
		HelpFormatter formatter = new HelpFormatter();
		CommandLine cmd;

		try {
			cmd = parser.parse(options, args);
		} catch (ParseException e) {
			System.out.println(e.getMessage());
			formatter.printHelp("utility-name", options);

			return null;
		}
		return cmd;

	}

}
