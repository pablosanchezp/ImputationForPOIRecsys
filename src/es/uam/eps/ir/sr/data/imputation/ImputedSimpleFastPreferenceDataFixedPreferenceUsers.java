package es.uam.eps.ir.sr.data.imputation;

import java.io.PrintStream;
import java.util.Map;
import java.util.function.Predicate;

import org.jooq.lambda.tuple.Tuple2;
import org.ranksys.core.util.tuples.Tuple2od;

import es.uam.eps.ir.crossdomainPOI.utils.UsersMidPoints;
import es.uam.eps.ir.ranksys.core.Recommendation;
import es.uam.eps.ir.ranksys.fast.preference.SimpleFastPreferenceData;
import es.uam.eps.ir.ranksys.rec.Recommender;
import es.uam.eps.ir.sr.utils.PredicatesStrategies;
import es.uam.eps.ir.sr.utils.PredicatesStrategies.recommendationStrategy;

/***
 * Imputed preference class that input preferences by a specific number of imputed preferences (if the value of the predicted preference if higher than the confidence)
 * @author Pablo Sanchez (psperez@icai.comillas.edu)
 *
 * @param <U>
 * @param <I>
 */
public class ImputedSimpleFastPreferenceDataFixedPreferenceUsers<U, I>
		extends ImputedSimpleFastPreferenceDataRecommender<U, I> {

	protected int numberImputedPreferences;

	public ImputedSimpleFastPreferenceDataFixedPreferenceUsers(SimpleFastPreferenceData<U, I> prefData,
			Recommender<U, I> recommender, boolean perUser, int numberImputedPreferences, recommendationStrategy imputationStrategyFilter, Map<I, Tuple2<Double, Double>> coordinatesPOI, UsersMidPoints<U, I> usermid, double limitKm, boolean binaryImputation) {
		super(prefData, recommender, perUser, imputationStrategyFilter, coordinatesPOI, usermid, limitKm, binaryImputation);
		this.numberImputedPreferences = numberImputedPreferences;

	}

	@Override
	public void imputeData(PrintStream newPreferenceData) {

		prefData.getAllUsers().forEach(u -> {

			// Now, we impute the next preferences for the user (the ones necessary to
			// achieve the stipulated percentage)
			if (this.numberImputedPreferences > 0) {

				// We generate the necessary imputed values (top N recommendation for that user)
				Recommendation<U, I> rec = recommender.getRecommendation(u, numberImputedPreferences, PredicatesStrategies.trainItems(u, this.prefData));
				for (Tuple2od<I> t : rec.getItems()) {
					super.printImputedPreference(newPreferenceData, u, t.v1, t.v2);
				}

			}

		});

	}
}
