package es.uam.eps.ir.sr.recommenders.POI.FMFMGM;

import org.apache.commons.math3.distribution.AbstractRealDistribution;
import org.apache.commons.math3.distribution.GammaDistribution;

import com.google.common.util.concurrent.AtomicDouble;

import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import cern.colt.matrix.impl.SparseDoubleMatrix2D;
import cern.colt.matrix.linalg.Algebra;
import es.uam.eps.ir.ranksys.fast.preference.FastPreferenceData;
import es.uam.eps.ir.ranksys.rec.fast.FastRankingRecommender;
import it.unimi.dsi.fastutil.ints.Int2DoubleMap;
import it.unimi.dsi.fastutil.ints.Int2DoubleOpenHashMap;


/***
 * Poisson Factor Model (PMF)
 * 
 * Part of the model proposed in:
 * 
 * [1] Fused Matrix Factorization with Geographical and Social Influence in Location-Based Social Networks (2012)
 * 
 * Based on the code of:
 * 
 * [2] An Experimental Evaluation of Point-of-interest Recommendation in Location-based Social Networks (2017)
 * 
 * @author Pablo Sanchez (psperez@icai.comillas.edu)
 *
 * @param <U>
 * @param <I>
 */
public class PoissonFactorModel <U, I> extends FastRankingRecommender<U, I>{

	private final int iterations;
	private final int factors;
	private final double alpha;
	private final double beta;
	private final double learningRate;
	private final boolean useSigmoid;
	
	private final AbstractRealDistribution gammaDistribution;
	
	private final int numUsers;
	private final int numItems;
	private final FastPreferenceData<U, I> data;


	private DenseDoubleMatrix2D uMatrix;
	private DenseDoubleMatrix2D iMatrix;
	private SparseDoubleMatrix2D fMatrix;
	

	
	
	public PoissonFactorModel(FastPreferenceData<U, I> fastPreferenceData, int iterations, int factors, double alpha, double beta, 
			double learningRate, boolean useSigmoid) {
		super(fastPreferenceData, fastPreferenceData);
		
		this.gammaDistribution = new GammaDistribution(alpha, beta);
		this.data = fastPreferenceData;
		this.iterations = iterations;
		this.factors = factors;
		this.numUsers = fastPreferenceData.numUsers();
		this.numItems = fastPreferenceData.numItems();
		this.alpha = alpha;
		this.beta = beta;
		this.learningRate = learningRate;
		this.useSigmoid = useSigmoid;
		
		this.uMatrix = new DenseDoubleMatrix2D(this.numUsers, this.factors);
		this.iMatrix = new DenseDoubleMatrix2D(this.numItems, this.factors);
		this.fMatrix = new SparseDoubleMatrix2D(this.numUsers, this.numItems);
		
		
		this.initFMatrix();
		
		
		
		this.uMatrix.assign(0.5 * Math.sqrt(gammaRandomGeneration(this.alpha, this.beta)) / this.factors);
		this.iMatrix.assign(0.5 * Math.sqrt(gammaRandomGeneration(this.alpha, this.beta)) / this.factors);

		this.train();

		
	}
	
	/**
	 * Train the model (basically the same as MF approach)
	 */
	private void train() {
		double tau = 10;
		double lastLoss = Double.MAX_VALUE;
		
		for (int i = 0; i < this.iterations; i++) {
			DoubleMatrix2D fMatrixCpy = fMatrix.copy();
			
			//For every item index
			this.data.getUidxWithPreferences().forEach(uidx -> {
				this.data.getUidxPreferences(uidx).forEach(pref -> {
					double newVal = (1.0 * fMatrix.get(uidx, pref.v1)) / this.uMatrix.viewRow(uidx).zDotProduct(this.iMatrix.viewRow(pref.v1)) - 1.0;
					fMatrixCpy.set(uidx, pref.v1, newVal);
				});
			});	
			
			double learningRatek = learningRate * tau / (tau + i);
			//Update U matrix
			DoubleMatrix2D fMatrixCpyDotL = fMatrixCpy.zMult(this.iMatrix, null);	
			this.uMatrix.assign(fMatrixCpyDotL, (x, y) -> x + learningRatek * (y + ((this.alpha - 1) / x) - (1 / this.beta )));


			//Update I matrix
			DoubleMatrix2D fMatrixCpyTransposeDotU = Algebra.DEFAULT.transpose(fMatrixCpy);
			fMatrixCpyTransposeDotU = fMatrixCpyTransposeDotU.zMult(this.uMatrix, null);
			this.iMatrix.assign(fMatrixCpyTransposeDotU, (x, y) -> x + learningRatek * (y + ((this.alpha - 1) / x) - (1 / this.beta )));

			AtomicDouble loss = new AtomicDouble(0.0);
			
			//For every item index
			this.data.getUidxWithPreferences().forEach(uidx -> {
				this.data.getUidxPreferences(uidx).forEach(pref -> {
					loss.addAndGet(Math.pow(fMatrix.getQuick(uidx, pref.v1) - this.uMatrix.viewRow(uidx).zDotProduct(this.iMatrix.viewRow(pref.v1)), 2));
				});
			});	
			System.out.println("Iteration: " + (i + 1) + ". Loss: " + loss.get());
			if (loss.get() > lastLoss) {
				System.out.println("Converged");
				return;
			}
			
			lastLoss = loss.get();
			
		}
	}
	
	/***
	 * Init the frequency matrix
	 */
	private void initFMatrix() {
		this.data.getUidxWithPreferences().forEach(uidx -> {
			this.data.getUidxPreferences(uidx).forEach(pref -> {
				this.fMatrix.set(uidx, pref.v1, pref.v2);
			});
		});
	}
	
	/***
	 * Method to obtain a random number based on the gamma distribution
	 * @param shape the first parameter
	 * @param scale the second parameter
	 * @return
	 */
	private double gammaRandomGeneration(double shape, double scale) {
		return this.gammaDistribution.sample();
	}

	@Override
	public Int2DoubleMap getScoresMap(int uidx) {
		Int2DoubleOpenHashMap scoresMap = new Int2DoubleOpenHashMap();
		scoresMap.defaultReturnValue(0.0);

		this.iIndex.getAllIidx().forEach(itemIndex -> {
			double score = this.uMatrix.viewRow(uidx).zDotProduct(this.iMatrix.viewRow(itemIndex));
			
			if (this.useSigmoid) {
				score = 1.0 / (1 + Math.exp(-score));
			}
			
			scoresMap.put(itemIndex, score);
		});
		
		return scoresMap;
	}

}
