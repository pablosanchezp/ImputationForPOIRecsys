package es.uam.eps.ir.sr.data.imputation;

import java.io.PrintStream;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import org.jooq.lambda.tuple.Tuple2;
import org.ranksys.core.util.tuples.Tuple2od;

import es.uam.eps.ir.crossdomainPOI.utils.UsersMidPoints;
import es.uam.eps.ir.ranksys.core.Recommendation;
import es.uam.eps.ir.ranksys.fast.preference.FastPreferenceData;
import es.uam.eps.ir.ranksys.rec.Recommender;
import es.uam.eps.ir.sr.utils.FrequencyGrouping;
import es.uam.eps.ir.sr.utils.PredicatesStrategies;
import es.uam.eps.ir.sr.utils.PredicatesStrategies.recommendationStrategy;


/***
 * Imputation preference class that will increase the number of preferences for EVERY USER in test appearing in training until a percentage.
 * For example, a value of 10 of percentage will denote that every user in the trainFile will increment their number of preferences by a 10%..
 * 
 * @author Pablo Sanchez (psperez@icai.comillas.edu)
 *
 * @param <U>
 * @param <I>
 */
public class ImputedSimpleFastPreferenceDataTestUsers<U, I>
		extends ImputedSimpleFastPreferenceDataRecommender<U, I> {

	private double percentageIncrement; //Value between 0 and positive value
	private FastPreferenceData<U, I> prefTestDataTest;
	private boolean userGrouping;
	private boolean userGropingUseSize;

	public ImputedSimpleFastPreferenceDataTestUsers(FastPreferenceData<U, I> prefData, FastPreferenceData<U, I> prefTestDataTest,
			Recommender<U, I> recommender, boolean perUser, double percentajeIncrement, recommendationStrategy imputationStrategyFilter, 
			Map<I, Tuple2<Double, Double>> coordinatesPOI, UsersMidPoints<U, I> usermid, double limitKm, boolean binaryImputation, 
			boolean userGrouping, boolean userGropingUseSize) {
		super(prefData, recommender, perUser, imputationStrategyFilter, coordinatesPOI, usermid, limitKm, binaryImputation);
		this.percentageIncrement = percentajeIncrement;
		this.binaryImputation = binaryImputation;
		this.prefTestDataTest = prefTestDataTest;		
		this.userGrouping = userGrouping;
		this.userGropingUseSize = userGropingUseSize;
		
		System.out.println("ImputedSimpleFastPreferenceDataTestUsers");
		System.out.println("Binary imputation: " + binaryImputation);

		System.out.println("Increment: " + percentajeIncrement);
		System.out.println("User grouping: " + userGrouping);
		System.out.println("User Groping size: " + userGropingUseSize);

	}

	@Override
	protected void imputeData(PrintStream newPreferenceData) {

		//We now order the users from lower to higher sparsity. Hence 
		this.prefTestDataTest.getAllUsers().forEach(u -> {
			if (this.prefData.containsUser(u)) {
			
				long numberPreferences = this.prefData.getUserPreferences(u).count();
				
				int necessaryImputedValues = (int) Math.ceil(numberPreferences * this.percentageIncrement / 100.0);
				
				FrequencyGrouping group = new FrequencyGrouping(necessaryImputedValues);
				if (userGrouping) {
					List<Double> listUser = this.prefData.getUserPreferences(u).map(t -> t.v2).collect(Collectors.toList());
					if (userGropingUseSize) {
						group.groupFrequenciesBySize(listUser);
					} else {
						group.groupFrequenciesByValues(listUser);
					}
					
				}
			
				//We generate a prediction of that number of required prefences.
	
				// We generate the necessary imputed values (top N recommendation for that user)
				Recommendation<U, I> rec = recommender.getRecommendation(u, PredicatesStrategies.trainItems(u, this.prefData));
	
				int counter = 0;
				for (Tuple2od<I> t : rec.getItems()) {
	
					// Repeat until we end or if the achieved sparsity is higher or equal the expected sparsity
					if (necessaryImputedValues <= 0) {
						break;
					}
					
					double val = this.userGrouping ? group.giveValue(counter) : t.v2; 

	
					// We binarize the imputation
					boolean printed = super.printImputedPreference(newPreferenceData, u, t.v1, val);
					
					//The preference might not pass the test (i.e. for example, by distance. In that case, we do not consider it)
					if (printed) {
						necessaryImputedValues--;
						counter++;
					}
				}
			}

			
		});
			
				
		
		return;

	}
}


