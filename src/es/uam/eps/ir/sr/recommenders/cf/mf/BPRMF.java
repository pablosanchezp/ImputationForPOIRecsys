package es.uam.eps.ir.sr.recommenders.cf.mf;

import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import org.jooq.lambda.tuple.Tuple3;

import com.google.common.util.concurrent.AtomicDouble;

import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.impl.DenseDoubleMatrix1D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import es.uam.eps.ir.ranksys.fast.preference.FastPreferenceData;
import es.uam.eps.ir.ranksys.rec.fast.FastRankingRecommender;
import es.uam.eps.ir.seqawareev.utils.CernMatrixUtils;
import it.unimi.dsi.fastutil.ints.Int2DoubleMap;
import it.unimi.dsi.fastutil.ints.Int2DoubleOpenHashMap;
/***
 * BPRMF recommender (once it has been created the object, it is necessary to
 * call fit() method)
 * 
 * - [1] BPR: Bayesian Personalized Ranking from Implicit Feedback (2009) :
 * https://dl.acm.org/citation.cfm?id=1795167
 * 
 * Based on libRec and MyMedialite:
 * 
 * Librec:
 * https://github.com/guoguibing/librec/blob/38c3d3325a9dc5e3e71a3c0d60e4db9d4ea8c814/core/src/main/java/net/librec/recommender/cf/ranking/BPRRecommender.java
 * 
 * Mymedialite:
 * https://github.com/zenogantner/MyMediaLite/blob/d5e09470ee45ad2c19c60b6c67d9807f83b55830/src/MyMediaLite/ItemRecommendation/BPRMF.cs
 * 
 * 
 * @author Pablo Sanchez (psperez@icai.comillas.edu)
 *
 * @param <U>
 * @param <I>
 */
public class BPRMF<U, I> extends FastRankingRecommender<U, I> {

	// Final parameters of the model (arguments that will not be updated)
	protected final double regUser;
	protected final double regItem;
	protected final int iterations;
	protected final int factors;
	protected final double learningRate;
	protected final double regBias;
	protected final boolean withoutReplacementUniformPair;
	protected final int maxTrials;

	

	// Other (variables obtained from arguments that wont be modified)
	protected final int numUsers;
	protected final int numItems;
	protected final int numberPreferences;
	protected final FastPreferenceData<U, I> data;

	// Learned (updated every iteration)
	protected DenseDoubleMatrix2D userFactors;
	protected DenseDoubleMatrix2D itemFactors;
	protected DenseDoubleMatrix1D itemBias;

	public BPRMF(FastPreferenceData<U, I> fastPreferenceData, double regUser, double regItem, double regBias,
			int iterations, int factors, double learningRate, double initMean, double initStd, boolean withoutReplacementUniformPair, int maxTrials) {
		super(fastPreferenceData, fastPreferenceData);

		this.numUsers = fastPreferenceData.numUsers();
		this.numItems = fastPreferenceData.numItems();

		this.regUser = regUser;
		this.regItem = regItem;

		this.iterations = iterations;
		this.factors = factors;
		this.learningRate = learningRate;
		this.data = fastPreferenceData;
		this.regBias = regBias;
		this.withoutReplacementUniformPair = withoutReplacementUniformPair;
		this.maxTrials = maxTrials;
		

		this.userFactors = new DenseDoubleMatrix2D(this.numUsers, this.factors);
		this.itemFactors = new DenseDoubleMatrix2D(this.numItems, this.factors);
		this.itemBias = new DenseDoubleMatrix1D(this.numItems);

		this.numberPreferences = fastPreferenceData.numPreferences();

		CernMatrixUtils.initilizeRandomGaussian(this.userFactors, initMean, initStd);
		CernMatrixUtils.initilizeRandomGaussian(this.itemFactors, initMean, initStd);
		this.itemBias.assign(CernMatrixUtils.gaussianRandomGenerator(initMean, initStd));
	}

	public BPRMF(FastPreferenceData<U, I> fastPreferenceData, double regUser, double regItem, double regBias,
			int iterations, int factors, double learningRate) {
		this(fastPreferenceData, regUser, regItem, regBias, iterations, factors, learningRate, 0.0, 0.1, true, 100);
	}

	/***
	 * Method to fit the model according to its parameters
	 */
	public void fit() {

		double lastLoss = Double.MAX_VALUE;

		for (int i = 0; i < this.iterations; i++) {
			double actualLoss = 0;
			
			if (this.withoutReplacementUniformPair)
				actualLoss = optimizationBPRIterateWithoutReplacementUniformPair();
			else {
				actualLoss = optimizationBPRByPreferences();
			}
			System.out.println("Iter: " + (i + 1) + " ActualLoss: " + actualLoss + " LastLoss: " + lastLoss);
			lastLoss = actualLoss;
		}

	}

	/**
	 * Optimize the BPR taking uniform sampling from the preferences
	 * 
	 * @return the loss
	 */
	protected double optimizationBPRByPreferences() {
		double actualLoss = 0.0;
		for (int j = 0; j < this.numberPreferences; j++) {
			Tuple3<Integer, Integer, Integer> tuple = getSamplesUserPosItemNegItem();
			int uidx = tuple.v1;
			int iPosItemIdx = tuple.v2;
			int jNegItemIdx = tuple.v3;

			actualLoss += optimizeBPRUserPosItemNegItem(uidx, iPosItemIdx, jNegItemIdx);

		}
		return actualLoss;
	}
	
	/**
	 * Optimize the BPR using the preferences (only sampling neg item)
	 * 
	 * @return the loss
	 */
	protected double optimizationBPRIterateWithoutReplacementUniformPair() {
		AtomicDouble actualLoss = new AtomicDouble(0.0);
		
		//Randomize the users and their preferences
		List<Integer> uidxs = this.data.getUidxWithPreferences().boxed().collect(Collectors.toList());
		Collections.shuffle(uidxs);
		
		uidxs.stream().forEach(uidx -> {
			List<Integer> itemsRatedUserList = this.data.getUidxPreferences(uidx).map(t -> t.v1).collect(Collectors.toList());
			Collections.shuffle(itemsRatedUserList);

			Set<Integer> itemsRatedUserSet = new HashSet<>(itemsRatedUserList);
			
			itemsRatedUserList.stream().forEach(posItemIdx -> {	
				int jNegItemIdx = randomSampleNegItem(itemsRatedUserSet);
				actualLoss.addAndGet(optimizeBPRUserPosItemNegItem(uidx, posItemIdx, jNegItemIdx));
			});
			
		});
			
		return actualLoss.get();
	}
	
	

	/***
	 * Method to optimize the BPR once having the user index, the positive item
	 * index and the negative one
	 * 
	 * @param uidx
	 *            the user index
	 * @param iPosItemIdx
	 *            the positive item index
	 * @param jNegItemIdx
	 *            the negative item index
	 * @return the loss
	 */
	protected double optimizeBPRUserPosItemNegItem(int uidx, int iPosItemIdx, int jNegItemIdx) {
		if (jNegItemIdx == -1) {
			//System.out.println("Could not find negative item for user " + this.data.uidx2user(uidx));
			return 0;
		}
		
		double actualLoss = 0.0;

		DoubleMatrix1D userFactors = this.userFactors.viewRow(uidx);
		DoubleMatrix1D iPosItemIdxFactors = this.itemFactors.viewRow(iPosItemIdx);
		DoubleMatrix1D jNegItemIdxFactors = this.itemFactors.viewRow(jNegItemIdx);

		double yui = userFactors.zDotProduct(iPosItemIdxFactors) + this.itemBias.getQuick(iPosItemIdx);
		double yuj = userFactors.zDotProduct(jNegItemIdxFactors) + this.itemBias.getQuick(jNegItemIdx);

		double diff = yui - yuj;

		actualLoss += -Math.log(CernMatrixUtils.logistic(diff));

		double deriValue = CernMatrixUtils.logistic(-diff);

		// Adjust factors
		for (int k = 0; k < this.factors; k++) {
			double userFactorVal = this.userFactors.getQuick(uidx, k);
			double posItemFactorVal = this.itemFactors.getQuick(iPosItemIdx, k);
			double negItemFactorVal = this.itemFactors.getQuick(jNegItemIdx, k);

			this.userFactors.setQuick(uidx, k, userFactorVal + this.learningRate * (deriValue * (posItemFactorVal - negItemFactorVal) - this.regUser * userFactorVal));
			this.itemFactors.setQuick(iPosItemIdx, k, posItemFactorVal + this.learningRate * (deriValue * userFactorVal - this.regItem * posItemFactorVal));
			this.itemFactors.setQuick(jNegItemIdx, k, negItemFactorVal + this.learningRate * (deriValue * (-userFactorVal) - this.regItem * negItemFactorVal));

			actualLoss += this.regUser * userFactorVal * userFactorVal
					+ this.regItem * posItemFactorVal * posItemFactorVal
					+ this.regItem * negItemFactorVal * negItemFactorVal;

		}

		double posItemBiasVal = this.itemBias.getQuick(iPosItemIdx);
		double negItemBiasVal = this.itemBias.getQuick(jNegItemIdx);

		// Adjust Bias
		this.itemBias.setQuick(iPosItemIdx, posItemBiasVal + this.learningRate * (deriValue - this.regBias * posItemBiasVal));
		this.itemBias.setQuick(jNegItemIdx, negItemBiasVal + this.learningRate * (-deriValue - this.regBias * negItemBiasVal));

		actualLoss += this.regBias * posItemBiasVal * posItemBiasVal;
		actualLoss += this.regBias * negItemBiasVal * negItemBiasVal;

		return actualLoss;

	}

	@Override
	public Int2DoubleMap getScoresMap(int uidx) {
		Int2DoubleOpenHashMap scoresMap = new Int2DoubleOpenHashMap();
		scoresMap.defaultReturnValue(0.0);

		this.iIndex.getAllIidx().forEach(itemIndex -> {
			scoresMap.put(itemIndex, this.predictScore(uidx, itemIndex));
		});
		return scoresMap;
	}

	/***
	 * Prediction score of a classic MF algorithm
	 * 
	 * @param uidx
	 *            user index
	 * @param iidx
	 *            item index
	 * @return the score associated to that user and item
	 */
	protected double predictScore(int uidx, int iidx) {
		DoubleMatrix1D pu = this.userFactors.viewRow(uidx);
		DoubleMatrix1D qi = this.itemFactors.viewRow(iidx);
		return qi.zDotProduct(pu) + this.itemBias.getQuick(iidx);
	}

	/***
	 * Method to obtain random samples of the users, items rated and negative rated
	 * items
	 * 
	 * @return a tuple of indexes of the suer, positive item and negative item
	 */
	private Tuple3<Integer, Integer, Integer> getSamplesUserPosItemNegItem() {
		// Select a random user and then a random item that the user HAS rated
		// (positive) and that HAS NOT rated (negative)
		int randomUser = CernMatrixUtils.uniformRandomNumber(0, this.numUsers - 1);

		List<Integer> itemsRatedUser = this.data.getUidxPreferences(randomUser).mapToInt(t -> t.v1).boxed()
				.collect(Collectors.toList());
		Set<Integer> itemsRatedUserSet = new HashSet<>(itemsRatedUser);

		itemsRatedUser = itemsRatedUserSet.stream().collect(Collectors.toList());

		int randomPosItem = CernMatrixUtils.uniformRandomNumber(0, itemsRatedUser.size() - 1);
		int randomItemIdx = itemsRatedUser.get(randomPosItem);

		int randomNegIdx = randomSampleNegItem(itemsRatedUserSet);

		return new Tuple3<>(randomUser, randomItemIdx, randomNegIdx);
	}
	
	
	/***
	 * Method to uniform sampling a random negative item (item not rated in the set of items) 
	 * @param itemsRatedUserSet
	 * @return
	 */
	private int randomSampleNegItem(Set<Integer> itemsRatedUserSet) {
		int randomNegIdx = -1;
		
		int count = 0;
		do {
			
			if (count >= this.maxTrials) {
				randomNegIdx = -1;
				break;
			}
			
			randomNegIdx = CernMatrixUtils.uniformRandomNumber(0, this.numItems - 1);
			count ++;
		} while (itemsRatedUserSet.contains(randomNegIdx));
		
		return randomNegIdx;
	}
	

	/**
	 * Method to update the learn rate if needed
	 * @param iter the actual iteration
	 * @param actualLoss the actual loss
	 * @param lastLoss the previous loss
	 */
	/*
	protected void updateLearnRate(int iter, double actualLoss, Double lastLoss) {
		if (this.learningRate < 0) {
			return;
		}
		if (iter > 1 && this.isboldDriver) {
			this.learningRate = Math.abs(lastLoss) > Math.abs(actualLoss) ? this.learningRate * 1.05f : this.learningRate * 0.5f;
		} else {
			this.learningRate *= this.decay;
		}
		
		if (this.maxLearnRate > 0 && this.learningRate > this.maxLearnRate) {
			this.learningRate = this.maxLearnRate;
		}
	}
	*/

}
