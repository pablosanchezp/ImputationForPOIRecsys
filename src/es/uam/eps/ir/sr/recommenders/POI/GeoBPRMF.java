package es.uam.eps.ir.sr.recommenders.POI;

import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.Collectors;

import org.jooq.lambda.tuple.Tuple2;


import cern.colt.matrix.DoubleMatrix1D;

import es.uam.eps.ir.ranksys.fast.preference.FastPreferenceData;
import es.uam.eps.ir.seqawareev.utils.CernMatrixUtils;
import es.uam.eps.ir.sr.data.POIProcessData;
import es.uam.eps.ir.sr.recommenders.cf.mf.BPRMF;

/***
 * BPR with geographical influence from paper:
 * 
 * -[1] Joint Geo-Spatial Preference and Pairwise Ranking for Point-of-Interest
 * Recommendation(2016) : https://ieeexplore.ieee.org/document/7814578
 * 
 * Based on the code of:
 * 
 * -[2] Geographic Diversification of Recommended POIs in Frequently Visited
 * Areas (2019) : https://dlnext.acm.org/doi/abs/10.1145/3362505
 * 
 * @author Pablo Sanchez (psperez@icai.comillas.edu)
 *
 * @param <U>
 * @param <I>
 */
public class GeoBPRMF<U, I> extends BPRMF<U, I> {

	// Arguments
	private final Map<I, Tuple2<Double, Double>> mapCoordinates;
	private final double maxDistNeighKm;

	// Learned
	private Map<Integer, Set<Integer>> poiNeighs;

	public GeoBPRMF(FastPreferenceData<U, I> fastPreferenceData, double regUser, double regItem, double regBias,
			int iterations, int factors, double learningRate, Map<I, Tuple2<Double, Double>> mapCoordinates, 
			double distNeigh, double initMean, double initStd, int maxTrials) {
		super(fastPreferenceData, regUser, regItem, regBias, iterations, factors, learningRate, initMean, initStd, true, maxTrials);

		this.mapCoordinates = mapCoordinates;
		this.maxDistNeighKm = distNeigh;
		this.poiNeighs = new HashMap<>();
		computeGeoNeighPOIs();
	}

	public GeoBPRMF(FastPreferenceData<U, I> fastPreferenceData, double regUser, double regItem, double regBias,
			int iterations, int factors, double learningRate, Map<I, Tuple2<Double, Double>> mapCoordinates,
			double distNeigh) {
		this(fastPreferenceData, regUser, regItem, regBias, iterations, factors, learningRate, mapCoordinates, distNeigh, 0.0, 0.1, 100);
	}

	/***
	 * Method to fit the model according to its parameters
	 */
	@Override
	public void fit() {

		double lastLoss = Double.MAX_VALUE;

		for (int i = 0; i < this.iterations; i++) {
			double actualLoss = 0;

			// For each user
			for (int uidx = 0; uidx < this.numUsers; uidx++) {
				List<Integer> itemsRatedUserRandom = this.data.getUidxPreferences(uidx).mapToInt(t -> t.v1).boxed().collect(Collectors.toList());
				Collections.shuffle(itemsRatedUserRandom);
								
				for (Integer iPosItemIdx: itemsRatedUserRandom) {
					Tuple2<Integer, Integer> tupleItemIndexSampled = getSamplesUser(uidx, iPosItemIdx);
	

					Set<Integer> itemsRatedUserSet = new HashSet<>(itemsRatedUserRandom);
	
					int gGeoItemIdx = tupleItemIndexSampled.v1;
					int jNegItemIdx = tupleItemIndexSampled.v2;
	
					if (gGeoItemIdx == -1 || jNegItemIdx == -1) {
						// Could not found any geoNeigh, according to [2], apply normal BPR
						
						if (gGeoItemIdx == -1) {
							//System.out.println("Could not find geographical neighbour, entering normal BPR");
						}
						
						actualLoss += optimizeBPRUserPosItemNegItem(uidx, iPosItemIdx, jNegItemIdx);
	
					} else {
	
						DoubleMatrix1D userFactors = this.userFactors.viewRow(uidx);
						DoubleMatrix1D iPosItemIdxFactors = this.itemFactors.viewRow(iPosItemIdx);
						DoubleMatrix1D jNegItemIdxFactors = this.itemFactors.viewRow(jNegItemIdx);
						DoubleMatrix1D gGeoItemIdxFactors = this.itemFactors.viewRow(gGeoItemIdx);
	
						double posItemBias = this.itemBias.getQuick(iPosItemIdx);
						double negItemBias = this.itemBias.getQuick(jNegItemIdx);
						double geoItemBias = this.itemBias.getQuick(gGeoItemIdx);
	
						double yui = userFactors.zDotProduct(iPosItemIdxFactors) + posItemBias;
						double yuj = userFactors.zDotProduct(jNegItemIdxFactors) + negItemBias;
						double yug = userFactors.zDotProduct(gGeoItemIdxFactors) + geoItemBias;
	
						double weightIg = getWeight(gGeoItemIdx, itemsRatedUserSet);
	
						double cig = CernMatrixUtils.logistic(-weightIg * (yui - yug)) * weightIg;
						double cgj = CernMatrixUtils.logistic(-(yug - yuj));
	
						actualLoss += -Math.log(CernMatrixUtils.logistic(weightIg * (yui - yug)));
						actualLoss += -Math.log(CernMatrixUtils.logistic((yug - yuj)));
	
						// update the factors
	
						for (int k = 0; k < this.factors; k++) {
							double userFactorVal = userFactors.get(k);
							double iPosItemFactorVal = iPosItemIdxFactors.get(k);
							double gGeoItemFactorVal = gGeoItemIdxFactors.get(k);
							double jNegItemFactorVal = jNegItemIdxFactors.get(k);
	
							double newuserFactorVal = userFactorVal + this.learningRate * (cig * (iPosItemFactorVal - gGeoItemFactorVal) + cgj * (gGeoItemFactorVal - jNegItemFactorVal) - this.regUser * userFactorVal);
							double newiPosItemFactorVal = iPosItemFactorVal + this.learningRate * (cig * userFactorVal - this.regItem * iPosItemFactorVal);
							double newgGeoItemFactorVal = gGeoItemFactorVal + this.learningRate * ((cgj - cig) * userFactorVal - this.regItem * gGeoItemFactorVal);
							double newjNegItemFactorVal = jNegItemFactorVal + this.learningRate * (-cgj * userFactorVal - this.regItem * jNegItemFactorVal);
	
							this.userFactors.setQuick(uidx, k, newuserFactorVal);
							this.itemFactors.setQuick(iPosItemIdx, k, newiPosItemFactorVal);
							this.itemFactors.setQuick(gGeoItemIdx, k, newgGeoItemFactorVal);
							this.itemFactors.setQuick(jNegItemIdx, k, newjNegItemFactorVal);
	
							// Loss of everything except bias
							actualLoss += this.regUser * userFactorVal * userFactorVal
									+ this.regItem * iPosItemFactorVal * iPosItemFactorVal
									+ this.regItem * gGeoItemFactorVal * gGeoItemFactorVal
									+ this.regItem * jNegItemFactorVal * jNegItemFactorVal;
	
						}
	
						this.itemBias.setQuick(iPosItemIdx, posItemBias + this.learningRate * (cig - this.regBias * posItemBias));
						this.itemBias.setQuick(gGeoItemIdx, geoItemBias + this.learningRate * (cgj - cig - this.regBias * geoItemBias));
						this.itemBias.setQuick(jNegItemIdx, negItemBias + this.learningRate * (-1.0 * cgj - this.regBias * negItemBias));
	
						// Loss of bias
						actualLoss += this.regBias * posItemBias * posItemBias;
						actualLoss += this.regBias * geoItemBias * geoItemBias;
						actualLoss += this.regBias * negItemBias * negItemBias;
	
					}
				
				}

			}
			System.out.println("Iter: " + (i + 1) + " ActualLoss: " + actualLoss + " LastLoss: " + lastLoss);
			lastLoss = actualLoss;

		}

	}

	/***
	 * Method to obtain the weight associated to a geographical neighbour
	 * 
	 * @param gGeoItemIdx
	 *            the geo neigh
	 * @param itemsRatedUserSet
	 *            the set of items rated of the user
	 * @return the weight
	 */
	private double getWeight(int gGeoItemIdx, Set<Integer> itemsRatedUserSet) {
		Set<Integer> neighs = this.poiNeighs.get(gGeoItemIdx);
		Set<Integer> union = new HashSet<>(itemsRatedUserSet);
		union.retainAll(neighs);

		return 1.0 / (1.0 + (double) union.size());
	}

	/***
	 * Method that will retrieve a sample of items for the user u
	 * 
	 * @param uidx
	 *            the user
	 * @return a tuple containing the samples of the geoItem and
	 *         the negative item
	 */
	private Tuple2<Integer, Integer> getSamplesUser(int uidx, int iidx) {

		List<Integer> itemsRatedUser = this.data.getUidxPreferences(uidx).mapToInt(t -> t.v1).boxed().collect(Collectors.toList());
		
		Set<Integer> itemsRatedUserSet = new HashSet<>(itemsRatedUser);
		itemsRatedUser = itemsRatedUserSet.stream().collect(Collectors.toList());
		
		int randomPosItemG = -1;
		int randomGIdx = -1;
		List<Integer> neighsItem = this.poiNeighs.get(iidx).stream().collect(Collectors.toList());

		// find geo neigh of I that have not been rated by U
		int count = 0;
		if (neighsItem != null && neighsItem.size() != 0) {
			do {
				
				if (count >= this.maxTrials) {
					randomGIdx = -1;
					break;
				}
				count++;
				randomPosItemG = CernMatrixUtils.uniformRandomNumber(0, neighsItem.size() - 1);
				randomGIdx = neighsItem.get(randomPosItemG);
			} while (itemsRatedUserSet.contains(randomGIdx));
		}

		int randomNegIdx = iidx;

		count = 0;
		do {
			if (count >= this.maxTrials) {
				randomNegIdx = -1;
				break;
			}
			count++;
			randomNegIdx = CernMatrixUtils.uniformRandomNumber(0, this.numItems - 1);
			
		} while (itemsRatedUserSet.contains(randomNegIdx));

		return new Tuple2<>(randomGIdx, randomNegIdx);

	}

	/**
	 * Method to obtain all the POIs that are closer to a specific distance to all
	 * POIs in the system. Will fill the poiNeighsMap
	 */
	private void computeGeoNeighPOIs() {
		AtomicInteger count = new AtomicInteger(0);
		
		data.getIidxWithPreferences().forEach(iidx -> {
			I item1 = this.data.iidx2item(iidx);
			
			count.addAndGet(1);
			
			if (count.get() % 100 == 0) {
				//System.out.print(count.get() + " out of " + this.data.getIidxWithPreferences().count() + " ");
			}
			
			
			if (!poiNeighs.containsKey(iidx)) {
				poiNeighs.put(iidx, new HashSet<>());
			}
			
			data.getIidxWithPreferences().filter(iidx2 -> iidx2 != iidx).forEach(iidx2 -> {
				I item2 = this.data.iidx2item(iidx2);
				
				if (!poiNeighs.containsKey(iidx2)) {
					poiNeighs.put(iidx2, new HashSet<>());
				}

				Tuple2<Double, Double> coordItem1 = this.mapCoordinates.get(item1);
				Tuple2<Double, Double> coordItem2 = this.mapCoordinates.get(item2);
				
				//If they are already in the map, do not compute the distance again
				if (!poiNeighs.get(iidx).contains(iidx2) && !poiNeighs.get(iidx2).contains(iidx)) {
					
					double dist = POIProcessData.haversine(coordItem1.v1, coordItem1.v2, coordItem2.v1, coordItem2.v2, false);
					if (dist <= this.maxDistNeighKm) {
						poiNeighs.get(iidx).add(iidx2);
						poiNeighs.get(iidx2).add(iidx);
					}
					
				}

			});
		});
		//System.out.println();
	}

}
