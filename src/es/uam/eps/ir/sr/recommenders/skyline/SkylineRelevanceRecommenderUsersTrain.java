package es.uam.eps.ir.sr.recommenders.skyline;

import es.uam.eps.ir.ranksys.fast.preference.FastPreferenceData;
import es.uam.eps.ir.ranksys.rec.fast.FastRankingRecommender;
import it.unimi.dsi.fastutil.ints.Int2DoubleMap;
import it.unimi.dsi.fastutil.ints.Int2DoubleOpenHashMap;

/***
 * Skyline recommender (working with the test set).
 * 
 * It only works with the test file returning the same preferences. As it it a
 * FastRanking recommender, the higher the score (if the score is higher than a
 * specific threshold), the better. In this case it only returns recommendations for the users in the training set.
 * 
 * @author Pablo Sanchez (psperez@icai.comillas.edu)
 *
 * @param <U> type of the users
 * @param <I> type of the items
 */
public class SkylineRelevanceRecommenderUsersTrain<U, I> extends FastRankingRecommender<U, I> {
	protected final FastPreferenceData<U, I> dataTest;
	protected final FastPreferenceData<U, I> dataTrain;

	private double relTh;

	public SkylineRelevanceRecommenderUsersTrain(FastPreferenceData<U, I> dataTrain, FastPreferenceData<U, I> dataTest, double relTh) {
		super(dataTrain, dataTrain);
		this.dataTest = dataTest;
		this.dataTrain = dataTrain;
		this.relTh = relTh;
	}

	@Override
	public Int2DoubleMap getScoresMap(int uidx) {
		Int2DoubleOpenHashMap scoresMap = new Int2DoubleOpenHashMap();
		scoresMap.defaultReturnValue(0.0);
		
		//The indexes of the train and test set may not be the same
		
		int uidxTest = this.dataTest.user2uidx(this.dataTrain.uidx2user(uidx));

		this.dataTest.getUidxPreferences(uidxTest).forEach(i -> {
			if (i.v2 >= relTh) {
				
				int iidxTrain = this.dataTrain.item2iidx(this.dataTest.iidx2item(i.v1));
				
				if (iidxTrain != -1) {		
					scoresMap.put(iidxTrain, i.v2);
				}
			}
		});

		return scoresMap;
	}
}