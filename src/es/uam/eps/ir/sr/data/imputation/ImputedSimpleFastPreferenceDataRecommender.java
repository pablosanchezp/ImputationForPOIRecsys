package es.uam.eps.ir.sr.data.imputation;

import java.io.PrintStream;
import java.util.Map;
import java.util.function.Predicate;

import org.jooq.lambda.tuple.Tuple2;

import es.uam.eps.ir.crossdomainPOI.utils.UsersMidPoints;
import es.uam.eps.ir.ranksys.fast.preference.FastPreferenceData;
import es.uam.eps.ir.ranksys.fast.preference.SimpleFastPreferenceData;
import es.uam.eps.ir.ranksys.rec.Recommender;
import es.uam.eps.ir.sr.utils.PredicatesStrategies;
import es.uam.eps.ir.sr.utils.PredicatesStrategies.recommendationStrategy;

/***
 * ImputedSimpleFastPreferenceDataRecommender class.
 * 
 * Extended definition of the ImputedSimpleFastPreferenceData
 * 
 * It adds a confidence value and a recommender
 * @author Pablo Sanchez (psperez@icai.comillas.edu)
 *
 * @param <U> the users
 * @param <I> the items
 */
public abstract class ImputedSimpleFastPreferenceDataRecommender<U, I> extends ImputedSimpleFastPreferenceData<U, I>{
	protected Recommender<U, I> recommender;
	
	
	protected recommendationStrategy imputationStrategyFilter;

	//Geographical imputation
	protected Map<I, Tuple2<Double, Double>> coordinatesPOI;
	protected UsersMidPoints<U, I> usermid;
	protected double limitKm;

	
	public ImputedSimpleFastPreferenceDataRecommender(FastPreferenceData<U, I> prefData, Recommender<U, I> recommender, boolean perUser, 
			recommendationStrategy imputationStrategyFilter, Map<I, Tuple2<Double, Double>> coordinatesPOI, UsersMidPoints<U, I> usermid, double limitKm, boolean binaryImputation) {
		super(prefData, perUser, binaryImputation);
		this.recommender = recommender;
		this.imputationStrategyFilter = imputationStrategyFilter;
		
		//Details of the geographical imputation
		this.coordinatesPOI =coordinatesPOI;
		this.usermid = usermid;
		this.limitKm = limitKm;
		this.binaryImputation = binaryImputation;
	}
	
	@Override
	protected boolean printImputedPreference(PrintStream newPreferenceData, U user, I item, double rating) {
		Predicate<I> pred = PredicatesStrategies.selectGeographicalPredicate(this.prefData, imputationStrategyFilter, user, coordinatesPOI, usermid, item, null, limitKm);
		if (pred.test(item)) {
			return super.printImputedPreference(newPreferenceData, user, item, rating);
		}	
		return false;
	}
	
	
}
