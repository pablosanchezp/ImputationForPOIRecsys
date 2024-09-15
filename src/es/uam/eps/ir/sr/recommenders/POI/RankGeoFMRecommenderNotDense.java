/*******************************************************************************
 * Copyright (C) 2018 Pablo Sánchez, Information Retrieval Group at Universidad Autónoma de Madrid, http://ir.ii.uam.es
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 ******************************************************************************/
package es.uam.eps.ir.sr.recommenders.POI;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import org.jooq.lambda.tuple.Tuple2;

import com.google.common.collect.HashBasedTable;
import com.google.common.collect.Table;
import com.google.common.util.concurrent.AtomicDouble;

import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import es.uam.eps.ir.ranksys.fast.preference.FastPreferenceData;
import es.uam.eps.ir.ranksys.fast.preference.IdxPref;
import es.uam.eps.ir.ranksys.rec.fast.FastRankingRecommender;
import es.uam.eps.ir.seqawareev.comparators.WeightComparatorTuple2;
import es.uam.eps.ir.seqawareev.utils.CernMatrixUtils;
import es.uam.eps.ir.sr.data.POIProcessData;
import it.unimi.dsi.fastutil.ints.Int2DoubleMap;
import it.unimi.dsi.fastutil.ints.Int2DoubleOpenHashMap;

/***
 * RankGeoFM recommender from paper: 
 * 	-Rank-GeoFM: A Ranking based Geographical FactorizationMethod for Point of Interest Recommendation
 * Implementation based on the one provided in LibRec:
 * https://github.com/guoguibing/librec/blob/38c3d3325a9dc5e3e71a3c0d60e4db9d4ea8c814/core/src/main/java/net/librec/recommender/poi/RankGeoFMRecommender.java
 * 
 * @author Pablo Sanchez (psperez@icai.comillas.edu)
 *
 * @param <U>
 * @param <I>
 */
public class RankGeoFMRecommenderNotDense<U, I> extends FastRankingRecommender<U, I> {

	// Final variables
	private final Map<Long, Tuple2<Double, Double>> mapCoordinates;
	private final FastPreferenceData<U, I> prefData;
	private final int iterations;
	private final boolean boldDriver;
	private final double decay;
	private final double maxLearingRate;
	private final int numberUsers;
	private final int numberPois;
	private final int numberFactors;


	// User Factors
	private DenseDoubleMatrix2D userFactors;

	// User Factors for geographical influence score
	private DenseDoubleMatrix2D geoUserFactors;

	// The poiFactors
	private DenseDoubleMatrix2D poiFactors;

	// knn influence matrix for geographical influence score
	private DenseDoubleMatrix2D geoInfluenceMatrix;

	// Simple variables
	// margin for ranking
	private double epsilon;

	// Regularization radius
	private double C;

	// weight of the radious
	private double alpha;

	// number of neighbours of the pois
	private int knnPois;

	// array for converting i into E[i] for loss
	double[] E;

	private double learnRate;

	private List<Set<Integer>> usersPoisSet;

	private Map<Integer, List<Integer>> poisNeighs;
	private Table<Integer, Integer, Double> poisNeighsWeights;

	public RankGeoFMRecommenderNotDense(FastPreferenceData<U, I> prefData, Map<Long, Tuple2<Double, Double>> mapCoordinates, int numberFactors,
			int knnPois, double alpha, double C, double epsilon, int numberIterations, double learnRate,
			double maxLearningRate, boolean boldDriver, double decay) {
		super(prefData, prefData);
		this.mapCoordinates = mapCoordinates;
		this.numberUsers = prefData.numUsers();
		this.numberPois = prefData.numItems();
		this.numberFactors = numberFactors;

		this.prefData = prefData;
		this.iterations = numberIterations;
		this.learnRate = learnRate;
		this.boldDriver = boldDriver;
		this.maxLearingRate = maxLearningRate;
		this.decay = decay;

		this.epsilon = epsilon;
		this.C = C;
		this.alpha = alpha;
		this.knnPois = knnPois;

		double initStd = 0.1;
		double initMean = 0.0; // from Matrix Factorization Recommender form librec

		this.geoInfluenceMatrix = new DenseDoubleMatrix2D(this.numberPois, numberFactors);
		this.userFactors = new DenseDoubleMatrix2D(this.numberUsers, numberFactors);
		this.geoUserFactors = new DenseDoubleMatrix2D(this.numberUsers, numberFactors);
		this.poiFactors = new DenseDoubleMatrix2D(this.numberPois, numberFactors);

		CernMatrixUtils.initilizeRandomGaussian(this.userFactors, initMean, initStd);
		CernMatrixUtils.initilizeRandomGaussian(this.geoUserFactors, initMean, initStd);
		CernMatrixUtils.initilizeRandomGaussian(this.poiFactors, initMean, initStd);

		this.usersPoisSet = getUserPoisSet();
		this.poisNeighs = new HashMap<>();
		this.poisNeighsWeights = getPoiKNNWeightMatrix();

		E = new double[numberPois + 1];
		for (int i = 1; i <= numberPois; i++) {
			E[i] = E[i - 1] + 1.0 / i;
		}
		
		System.out.println("Version Reset GeoInfluenceMatrix");
		System.out.println("Version random sampling the users and the POIs visited by each user per iteraction");
		this.trainModel();
	}

	@Override
	public Int2DoubleMap getScoresMap(int uidx) {
		Int2DoubleOpenHashMap scoresMap = new Int2DoubleOpenHashMap();
		scoresMap.defaultReturnValue(0.0);
		if (uidx == -1) {
			return scoresMap;
		}
		this.iIndex.getAllIidx().forEach(itemIndex -> {
			scoresMap.put(itemIndex, this.userFactors.viewRow(uidx).zDotProduct(this.poiFactors.viewRow(itemIndex))
					+ this.geoUserFactors.viewRow(uidx).zDotProduct(this.geoInfluenceMatrix.viewRow(itemIndex)));

		});
		return scoresMap;

	}

	private void trainModel() {
		AtomicDouble actualLoss = new AtomicDouble(0.0);
		AtomicDouble lastLoss = new AtomicDouble(Double.MAX_VALUE);

		for (int iter = 1; iter <= this.iterations; iter++) {
			updateGeoInfluenceMatrix();
			actualLoss.set(0.0);
			DoubleMatrix2D tempUserFactors = this.userFactors.copy();
			DoubleMatrix2D tempGeoUserFactors = this.geoUserFactors.copy();
			DoubleMatrix2D tempPoiFactors = this.poiFactors.copy();

			// for each user in the system (shuffling the users)
			List<Integer> lstUidx = this.prefData.getUidxWithPreferences().boxed().collect(Collectors.toList());
			Collections.shuffle(lstUidx);			
			
			lstUidx.stream().forEach(uidx -> {
				
				// for each preference of user in the system (shuffling the preferences)
				List<IdxPref> lstIidxPrefUser = this.prefData.getUidxPreferences(uidx).collect(Collectors.toList());
				Collections.shuffle(lstIidxPrefUser);
				
				lstIidxPrefUser.stream().forEach(pref -> {
					int posPoiIdx = pref.v1;
					double posRealRating = pref.v2;

					int sampleCount = 0;
					double posPredictRating = tempUserFactors.viewRow(uidx)
							.zDotProduct(tempPoiFactors.viewRow(posPoiIdx))
							+ tempGeoUserFactors.viewRow(uidx).zDotProduct(this.geoInfluenceMatrix.viewRow(posPoiIdx));

					int negPoiIdx;
					double negPredictRating;
					double negRealRating;
					int incompatibility;

					while (true) {
						negPoiIdx = CernMatrixUtils.uniformRandomNumber(0, this.numberPois - 1);
						negPredictRating = tempUserFactors.viewRow(uidx).zDotProduct(tempPoiFactors.viewRow(negPoiIdx))
								+ tempGeoUserFactors.viewRow(uidx)
										.zDotProduct(this.geoInfluenceMatrix.viewRow(negPoiIdx));

						Set<Integer> poisSet = this.usersPoisSet.get(uidx);
						int negPoiIdx2 = negPoiIdx;

						if (poisSet.contains(negPoiIdx)) {
							negRealRating = this.prefData.getUidxPreferences(uidx)
									.filter(pref2 -> pref2.v1 == negPoiIdx2).findFirst().get().v2;
						} else {
							negRealRating = 0.0;
						}

						sampleCount++;
						incompatibility = indicator(posRealRating, negRealRating)
								* indicator(negPredictRating + this.epsilon, posPredictRating);
						if (incompatibility == 1 || sampleCount > this.numberPois) {
							break;
						}
					}

					if (incompatibility == 1) {
						int lowerBound = this.numberPois / sampleCount;
						double s = CernMatrixUtils.logistic(negPredictRating + this.epsilon - posPredictRating);
						actualLoss.addAndGet(E[lowerBound] * s);
						double uij = s * (1 - s);
						double ita = E[lowerBound] * uij;

						// update userFactors and geoUserFactors
						DoubleMatrix1D updateUserVec = this.poiFactors.viewRow(negPoiIdx).copy();
						updateUserVec.assign(this.poiFactors.viewRow(posPoiIdx), (x, y) -> x - y);
						updateUserVec.assign(x -> x * this.learnRate * ita);

						this.userFactors.viewRow(uidx).assign(updateUserVec, (x, y) -> x - y);

						DoubleMatrix1D updateGeoUserVec = this.geoInfluenceMatrix.viewRow(negPoiIdx).copy();
						updateGeoUserVec.assign(this.geoInfluenceMatrix.viewRow(posPoiIdx), (x, y) -> x - y);
						updateGeoUserVec.assign(x -> x * this.learnRate * ita);

						this.geoUserFactors.viewRow(uidx).assign(updateGeoUserVec, (x, y) -> x - y);

						// update poiFactors
						DoubleMatrix1D updatePoiVec = this.userFactors.viewRow(uidx).copy();
						updatePoiVec.assign(x -> x * this.learnRate * ita);
						this.poiFactors.viewRow(posPoiIdx).assign(updatePoiVec, (x, y) -> x + y);
						this.poiFactors.viewRow(negPoiIdx).assign(updatePoiVec, (x, y) -> x - y);

						// regularize userFactors and geoUserFactors
						double userVectorNorm = CernMatrixUtils.normalizeVector(this.userFactors.viewRow(uidx), 2.0);
						if (userVectorNorm > this.C) {
							this.userFactors.viewRow(uidx).assign(x -> x * (this.C / userVectorNorm));
						}
						double geoUserVectorNorm = CernMatrixUtils.normalizeVector(this.geoUserFactors.viewRow(uidx),
								2.0);
						if (geoUserVectorNorm > this.alpha * this.C) {
							this.geoUserFactors.viewRow(uidx)
									.assign(x -> x * (this.alpha * this.C / geoUserVectorNorm));
						}

						// regularize poiFactors
						double posPoiVectorNorm = CernMatrixUtils.normalizeVector(this.poiFactors.viewRow(posPoiIdx),
								2.0);
						if (posPoiVectorNorm > C) {
							this.poiFactors.viewRow(posPoiIdx).assign(x -> x * (this.C / posPoiVectorNorm));
						}
						double negPoiVectorNorm = CernMatrixUtils.normalizeVector(this.poiFactors.viewRow(negPoiIdx),
								2.0);
						if (negPoiVectorNorm > C) {
							this.poiFactors.viewRow(negPoiIdx).assign(x -> x * (this.C / negPoiVectorNorm));
						}
					}

				});
			});

			if (isConverged(iter, actualLoss.get(), lastLoss.get())) {
				break;
			}

			updateLearnRate(iter, actualLoss.get(), lastLoss.get());

			lastLoss.set(actualLoss.get());
		}
	}

	private boolean isConverged(int iter, double actualLoss, Double lastloss) {
		System.out.println("Iter " + iter + " actualLoss: " + actualLoss + " lastLoss " + lastloss);
		if (Math.abs(lastloss - actualLoss) <= 1e-4 /* || actualLoss > lastloss */) {
			System.out.println("Converged");
			return true;
		}
		return false;
	}

	/**
	 * Method to update the learn rate if needed
	 * 
	 * @param iter
	 *            the actual iteration
	 * @param actualLoss
	 *            the actual loss
	 * @param lastLoss
	 *            the previous loss
	 */
	private void updateLearnRate(int iter, double actualLoss, Double lastLoss) {
		if (this.learnRate < 0) {
			return;
		}
		if (iter > 1 && boldDriver) {
			this.learnRate = Math.abs(lastLoss) > Math.abs(actualLoss) ? this.learnRate * 1.05f : this.learnRate * 0.5f;
		} else {
			this.learnRate *= this.decay;
		}

		if (this.maxLearingRate > 0 && this.learnRate > this.maxLearingRate) {
			this.learnRate = this.maxLearingRate;
		}
	}

	public void updateGeoInfluenceMatrix() {
		this.geoInfluenceMatrix = new DenseDoubleMatrix2D(this.numberPois, numberFactors);

		
		for (int poiIdx = 0; poiIdx < this.numberPois; poiIdx++) {
			List<Integer> neigsIdx = poisNeighs.get(poiIdx);
			DoubleMatrix1D geoInfluenceItem = this.geoInfluenceMatrix.viewRow(poiIdx);

			Integer poiIdx2 = poiIdx; // Created because variable need to be final
			for (Integer neighx : neigsIdx) {
				DoubleMatrix1D neighxFactors = this.poiFactors.viewRow(neighx);
				geoInfluenceItem.assign(neighxFactors, (x, y) -> x + y);
				geoInfluenceItem.assign(x -> x * this.poisNeighsWeights.get(poiIdx2, neighx));
			}
		}
	}

	private Table<Integer, Integer, Double> getPoiKNNWeightMatrix() {
		Table<Integer, Integer, Double> dataTable = HashBasedTable.create();

		for (int poiIdx = 0; poiIdx < this.numberPois; poiIdx++) {
			List<Tuple2<Integer, Double>> locationNeighbors = new ArrayList<>();
			
			Tuple2<Double, Double> coordinatesOriginal = this.mapCoordinates.get(this.prefData.iidx2item(poiIdx));
			for (int neighborItemIdx = 0; neighborItemIdx < this.numberPois; neighborItemIdx++) {
				if (poiIdx != neighborItemIdx) {
					Tuple2<Double, Double> coordinatesNeigh = this.mapCoordinates.get(this.prefData.iidx2item(neighborItemIdx));
					
					locationNeighbors
							.add(new Tuple2<>(neighborItemIdx, POIProcessData.haversine(coordinatesOriginal.v1, coordinatesOriginal.v2, coordinatesNeigh.v1, coordinatesNeigh.v2, false)));
				}
			}
			// Sorting from the lower to the highest distance
			Collections.sort(locationNeighbors, new WeightComparatorTuple2());
			locationNeighbors = locationNeighbors.subList(0, this.knnPois);
			this.poisNeighs.put(poiIdx, locationNeighbors.stream().map(t -> t.v1).collect(Collectors.toList()));

			for (Tuple2<Integer, Double> neighAndDistance : locationNeighbors) {
				int neighborItemIdx = neighAndDistance.v1;
				double weight;
				if (neighAndDistance.v2 < 0.5) {
					weight = 1.0 / 0.5;
				} else {
					weight = 1.0 / (neighAndDistance.v2);
				}

				dataTable.put(poiIdx, neighborItemIdx, weight);
			}

		}

		Table<Integer, Integer, Double> normalizedDataTable = HashBasedTable.create();
		for (int itemIdx = 0; itemIdx < this.numberPois; itemIdx++) {
			double rowSum = sumRowNeighWeights(itemIdx, dataTable);
			List<Integer> neighsIdx = this.poisNeighs.get(itemIdx);

			for (Integer neighIdx : neighsIdx) {
				normalizedDataTable.put(itemIdx, neighIdx, dataTable.get(itemIdx, neighIdx) / rowSum);
			}
		}
		return normalizedDataTable;

	}

	private int indicator(double i, double j) {
		return i > j ? 1 : 0;
	}

	private List<Set<Integer>> getUserPoisSet() {
		List<Set<Integer>> userPoisSet = new ArrayList<>();

		for (int userIdx = 0; userIdx < this.numberUsers; ++userIdx) {
			Set<Integer> itemList = this.prefData.getUidxPreferences(userIdx).map(pref -> pref.v1)
					.collect(Collectors.toSet());
			userPoisSet.add(new HashSet<>(itemList));
		}
		return userPoisSet;
	}

	private double sumRowNeighWeights(int idx, Table<Integer, Integer, Double> table) {
		List<Integer> neighsIdx = this.poisNeighs.get(idx);
		double acc = 0;
		for (Integer neigh : neighsIdx) {
			acc += table.get(idx, neigh);
		}
		return acc;
	}

}
