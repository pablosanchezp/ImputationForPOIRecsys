package es.uam.eps.ir.sr.recommenders.POI.FMFMGM;

import java.util.Map;

import org.jooq.lambda.tuple.Tuple2;

import es.uam.eps.ir.ranksys.fast.preference.FastPreferenceData;
import es.uam.eps.ir.ranksys.rec.fast.FastRankingRecommender;
import it.unimi.dsi.fastutil.ints.Int2DoubleMap;
import it.unimi.dsi.fastutil.ints.Int2DoubleOpenHashMap;

/***
 * Multicenter gaussian model + Poisson factor model.
 * 
 * Proposed in:
 *
 * [1] Fused Matrix Factorization with Geographical and Social Influence in Location-Based Social Networks (2012)
 *
 * @author Pablo Sanchez (psperez@icai.comillas.edu)
 *
 * @param <U>
 * @param <I>
 */
public class FMFMGM <U, I extends Comparable<I>> extends FastRankingRecommender<U, I>{

	private MultiCenterGaussianModel<U, I> MGM;
	private PoissonFactorModel<U, I> PFM;
	
	public FMFMGM(FastPreferenceData<U, I> fastPreferenceData, Map<I, Tuple2<Double, Double>> coordinatesItems, double alpha, double theta, double dmax, int iterations, int factors, double alpha2, double beta, 
			double learningRate, boolean useSigmoid) {
		super(fastPreferenceData, fastPreferenceData);
		
		this.MGM = new MultiCenterGaussianModel<>(fastPreferenceData, coordinatesItems, alpha, theta, dmax);
		this.PFM = new PoissonFactorModel<>(fastPreferenceData, iterations, factors, alpha2, beta, learningRate, useSigmoid);
	}

	@Override
	public Int2DoubleMap getScoresMap(int uidx) {
		Int2DoubleOpenHashMap scoresMap = new Int2DoubleOpenHashMap();
		scoresMap.defaultReturnValue(0.0);

		Int2DoubleMap uidxMGM = this.MGM.getScoresMap(uidx);
		Int2DoubleMap uidxPFM = this.PFM.getScoresMap(uidx);

		
		this.iIndex.getAllIidx().forEach(itemIndex -> {
			scoresMap.put(itemIndex, uidxMGM.get(itemIndex) * uidxPFM.get(itemIndex));
		});
		return scoresMap;
	}

}
