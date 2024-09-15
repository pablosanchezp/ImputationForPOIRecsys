package es.uam.eps.ir.sr.data.imputation;

import java.io.PrintStream;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.Collectors;

import org.jooq.lambda.tuple.Tuple2;
import org.ranksys.core.util.tuples.Tuple2od;

import es.uam.eps.ir.crossdomainPOI.utils.UsersMidPoints;
import es.uam.eps.ir.ranksys.core.Recommendation;
import es.uam.eps.ir.ranksys.fast.preference.FastPreferenceData;
import es.uam.eps.ir.ranksys.rec.Recommender;
import es.uam.eps.ir.sr.utils.FrequencyGrouping;
import es.uam.eps.ir.sr.utils.PredicatesStrategies;
import es.uam.eps.ir.sr.utils.PredicatesStrategies.recommendationStrategy;
import es.uam.eps.ir.sr.utils.comparators.PreferenceComparators;


/***
 * Imputation preference class that will increase the number of preferences for the users until a percentage
 * @author Pablo Sanchez (psperez@icai.comillas.edu)
 *
 * @param <U>
 * @param <I>
 */
public class ImputedSimpleFastPreferenceDataAvgPreferencesUsers<U, I>
		extends ImputedSimpleFastPreferenceDataRecommender<U, I> {

	private double percentageIncrement; //Value between 0 and positive value
	private boolean reverse;
	private boolean random;
	
	private boolean userGrouping;
	private boolean userGropingUseSize;

	public ImputedSimpleFastPreferenceDataAvgPreferencesUsers(FastPreferenceData<U, I> prefData,
			Recommender<U, I> recommender, boolean perUser, double percentajeIncrement, recommendationStrategy imputationStrategyFilter, 
			Map<I, Tuple2<Double, Double>> coordinatesPOI, UsersMidPoints<U, I> usermid, double limitKm, boolean reverse, boolean binaryImputation, boolean random, 
			boolean userGrouping, boolean userGropingUseSize) {
		super(prefData, recommender, perUser, imputationStrategyFilter, coordinatesPOI, usermid, limitKm, binaryImputation);
		this.percentageIncrement = percentajeIncrement;
		this.reverse = reverse;
		this.random = random;
		this.userGrouping = userGrouping;
		this.userGropingUseSize = userGropingUseSize;
		
		System.out.println("ImputedSimpleFastPreferenceDataAvgPreferencesUsers");
		System.out.println("Binary imputation: " + binaryImputation);
		System.out.println("Reversed: " + reverse);
		System.out.println("Random:" + random);
		System.out.println("Increment: " + percentajeIncrement);
		System.out.println("User grouping: " + userGrouping);
		System.out.println("User Groping size: " + userGropingUseSize);

	}

	@Override
	protected void imputeData(PrintStream newPreferenceData) {
		double totalPreferences = (double) prefData.getUsersWithPreferences()
				.mapToLong(u -> prefData.getUserPreferences(u).count()).sum();
		
		double originalDenominator = ((double) prefData.numUsers() * (double) prefData.numItems()); // Original sparsity
		double currentSparsity = totalPreferences / originalDenominator;
		
		double currentAvgPrefUser = totalPreferences / (double) prefData.numUsers(); // Original preferences users

		double expectedSparsity = ((currentSparsity *  percentageIncrement) / 100.0) + currentSparsity;
		System.out.println("Expected sparsity: " + expectedSparsity);
		// This would be the new number of preferences for the user (basically that all
		// users should have this number of preferences)
		double finalAvgPrefUser = (expectedSparsity * currentAvgPrefUser) / currentSparsity;

		//Obtain the tuples ordered
		List<Tuple2<Integer, Integer>> lst = prefData.getAllUsers().map(u -> new Tuple2<Integer, Integer>(prefData.user2uidx(u), (int) prefData.getUserPreferences(u).count())).collect(Collectors.toList());
		
		if (this.random) {
			Collections.shuffle(lst);
		}
		else {
			lst.sort(PreferenceComparators.Tuple2ComparatorSecondFirstElement());
			
			if (reverse) {
		        Collections.reverse(lst);
			}
		}
		System.out.println("First user: " + prefData.uidx2user(lst.get(0).v1) + " " + lst.get(0).v2);
		System.out.println("Last user: " + prefData.uidx2user(lst.get(lst.size() - 1).v1) + " " + lst.get(lst.size() - 1).v2);

		
		//We now order the users from lower to higher sparsity. Hence 
		int incrementedInteractions = 0;
		boolean cont = true;
		for (Tuple2<Integer, Integer> tupleUidxPreferences: lst) {
			
				//No continuation if we achieved the expected sparsity
				if (!cont) {
					break;
				}
				U u = prefData.uidx2user(tupleUidxPreferences.v1);


				AtomicInteger countPrefs = new AtomicInteger(tupleUidxPreferences.v2);
	
				// Once we have finished with the REAL preferences, we need now to compute the
				// new ones
	
				int necessaryImputedValues = (int) (finalAvgPrefUser - countPrefs.get());
	
				// Now, we impute the next preferences for the user (the ones necessary to
				// achieve the stipulated percentage)
				if (necessaryImputedValues > 0) {
	
					// We generate the necessary imputed values (top N recommendation for that user)
					Recommendation<U, I> rec = recommender.getRecommendation(u, PredicatesStrategies.trainItems(u, this.prefData));
					
					FrequencyGrouping group = new FrequencyGrouping(necessaryImputedValues);
					if (userGrouping) {
						List<Double> listUser = this.prefData.getUserPreferences(u).map(t -> t.v2).collect(Collectors.toList());
						if (userGropingUseSize) {
							group.groupFrequenciesBySize(listUser);
						} else {
							group.groupFrequenciesByValues(listUser);
						}
						
					}
					
					int counter = 0;
					for (Tuple2od<I> t : rec.getItems()) {
						double achievedSparsity = (totalPreferences + incrementedInteractions) / originalDenominator;

						//We managed to obtain the expected sparsity, we finish
						if (achievedSparsity >= expectedSparsity) {
							cont = false;
							System.out.println("Total (original) preferences: " + totalPreferences);
							System.out.println("Incremented Interactions: " + incrementedInteractions);
							System.out.println("Original |U| x |I| " + originalDenominator);
							System.out.println("Achieved desired Sparsity: " + achievedSparsity);
							break;
						}

						// Repeat until we end or if the achieved sparsity is higher or equal the expected sparsity
						if (necessaryImputedValues < 0) {
							break;
						}

						double val = this.userGrouping ? group.giveValue(counter) : t.v2; 

						
						// We binarize the imputation
						boolean printed = super.printImputedPreference(newPreferenceData, u, t.v1, val);
						
						//The preference might not pass the test (i.e. for example, by distance. In that case, we do not consider it)
						if (printed) {
							
							necessaryImputedValues--;
							incrementedInteractions++;
							counter++;
						}
					}
	
				}
		}
		return;

	}
}


