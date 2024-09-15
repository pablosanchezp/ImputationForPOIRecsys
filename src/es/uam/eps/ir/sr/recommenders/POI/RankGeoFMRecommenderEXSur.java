package es.uam.eps.ir.sr.recommenders.POI;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.SingularValueDecomposition;
import org.jooq.lambda.tuple.Tuple2;


import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import cern.colt.matrix.impl.SparseDoubleMatrix2D;
import es.uam.eps.ir.ranksys.fast.preference.FastPreferenceData;
import es.uam.eps.ir.ranksys.fast.preference.IdxPref;
import es.uam.eps.ir.ranksys.rec.fast.FastRankingRecommender;
import es.uam.eps.ir.seqawareev.comparators.WeightComparatorTuple2;
import es.uam.eps.ir.seqawareev.utils.CernMatrixUtils;
import es.uam.eps.ir.sr.data.POIProcessData;
import it.unimi.dsi.fastutil.ints.Int2DoubleMap;
import it.unimi.dsi.fastutil.ints.Int2DoubleOpenHashMap;

/***
 * RankGeoFM recommender: 
 * 
 * Implementation based on the RankGeoFM recommender from paper:
 *
 * "An Experimental Evaluation of Point-of-interest Recommendation in Location-based Social Networks"
 * 
 * @author Pablo Sanchez (psperez@icai.comillas.edu)
 *
 * @param <U>
 * @param <I>
 */
public class RankGeoFMRecommenderEXSur<U, I> extends FastRankingRecommender<U, I> {

	private final Map<Long, Tuple2<Double, Double>> mapCoordinates;
	private final FastPreferenceData<U, I> prefData;
	private final int iterations;
	private final int numberUsers;
	private final int numberPois;
	private final int numberFactors;
	
	private DenseDoubleMatrix2D A;
	private int[][] aIndex;

	// User Factors
	private DenseDoubleMatrix2D U_1;

	// User Factors for geographical influence score
	private DenseDoubleMatrix2D U_2;

	// The poiFactors
	private DenseDoubleMatrix2D L_1;
	
	private double FactorInf;

	// knn influence matrix for geographical influence score
	private DenseDoubleMatrix2D F_G;
	
	//Checkin matrix
	private SparseDoubleMatrix2D B;

	// Simple variables
	// margin for ranking (maginal in Matlab code)
	private final double epsilon;

	// Regularization radius
	private final double C;

	// weight of the radious
	private final double alpha;

	// number of neighbours of the pois
	private int knnPois;

	// array for converting i into E[i] for loss
	double[] LossWeight;


	public RankGeoFMRecommenderEXSur(FastPreferenceData<U, I> prefData,
			Map<Long, Tuple2<Double, Double>> mapCoordinates, int numberFactors, int knnPois, double factorInf, double C,
			double epsilon, double alpha, int numberIterations) {
		super(prefData, prefData);
		this.mapCoordinates = mapCoordinates;
		this.numberUsers = prefData.numUsers();
		this.numberPois = prefData.numItems();
		this.numberFactors = numberFactors;
		this.knnPois = knnPois;

		this.prefData = prefData;
		this.iterations = numberIterations;
		
		//FACTORINF (ver que es, creo que es el anterior alpha, porque lo ponen a 0.2)
		this.FactorInf = factorInf;
		this.C = C;
		this.epsilon = epsilon;
		this.alpha = alpha; //me da que este alfa es el learning rate
		

		
		generateAIndexAMatrix();
		
		this.U_1 = new DenseDoubleMatrix2D(this.numberUsers, this.numberFactors);
		this.U_2 = new DenseDoubleMatrix2D(this.numberUsers, this.numberFactors);
		this.L_1 = new DenseDoubleMatrix2D(this.numberPois, this.numberFactors);
		
		//Create the factor matrices
		CernMatrixUtils.initilizeRandomGaussian(this.U_1, 0.0, 0.01);
		CernMatrixUtils.initilizeRandomGaussian(this.U_2, 0.0, 0.01);
		CernMatrixUtils.initilizeRandomGaussian(this.L_1, 0.0, 0.01);
		
		
		//Generate the checkin-matrix (B in the code of the experimental survey)
		this.B = new SparseDoubleMatrix2D(this.numberUsers, this.numberPois);
		this.prefData.getUidxWithPreferences().forEach(uidx -> {
			this.prefData.getUidxPreferences(uidx).forEach(iipref -> {
				this.B.setQuick(uidx, iipref.v1, iipref.v2);
			});
		});
		
		
		this.LossWeight = new double[this.numberPois];
		double total = 0;
		for (int i = 0; i < numberPois; i++) {
			this.LossWeight[i] = total + 1.0 / (i + 1);
			total += 1.0 / (i+1);
		}
		
		System.out.println("Version random sampling the users and the POIs visited by each user per iteraction (experimental survey)");

		trainModel();
	}

	

	
	
	private void trainModel() {	
		//For iterations
		for (int it = 0; it < this.iterations; it++) {
			
			for (int i = 0; i < this.numberUsers; i++) {
				double norm = CernMatrixUtils.normalizeVector(this.U_2.viewRow(i), 2);
				if (norm > this.C * this.FactorInf) {
					this.U_2.viewRow(i).assign(x -> (this.FactorInf * this.C * x) / norm);
				}	
			}
			
			//Reset the geoInfluenceMatrix
			this.F_G = new DenseDoubleMatrix2D(this.numberPois, this.numberFactors);
			for (int i = 0; i < this.numberPois; i++) {
				for (int j = 0; j< this.knnPois; j++) {
					int iI = i;
					int iJ = j;
					int indexNeigh = this.aIndex[i][j];
					DoubleMatrix1D l_1_n = this.L_1.viewRow(indexNeigh).copy();
					this.F_G.viewRow(i).assign(l_1_n, (x, y) -> x + this.A.getQuick(iI, iJ) * y);

				}				
			}
			
			DoubleMatrix2D U_1_pre = this.U_1.copy();
			DoubleMatrix2D U_2_pre = this.U_2.copy();
			DoubleMatrix2D L_1_pre = this.L_1.copy();
			
			
			DoubleMatrix2D UL = this.U_1.zMult(this.L_1.viewDice(), null);
			DoubleMatrix2D UFG = this.U_2.zMult(this.F_G.viewDice(), null);
			
			//Randomize the users
			List<Integer> lstUsers = this.prefData.getUidxWithPreferences().boxed().collect(Collectors.toList());
			Collections.shuffle(lstUsers);
			
			lstUsers.stream().forEach(uidx -> {
				int myu = uidx;
				
				//Randomize the preferences of the user 
				List<IdxPref> lstIidxPrefUser = this.prefData.getUidxPreferences(uidx).collect(Collectors.toList());
				Collections.shuffle(lstIidxPrefUser);

				lstIidxPrefUser.stream().forEach(iidxpref -> {
					int mya = iidxpref.v1;
					
					Tuple2<Integer, Integer> NandMyb = getN(myu, mya, UL, UFG);
	
					Integer N = NandMyb.v1;
					Integer Myb = NandMyb.v2;
					int myb = Myb;
					
					if (N != -1 && myb != -1 && N <= this.numberPois - 1) {
						N = (int) Math.floor( (this.numberPois - 1) / (N + 1));
						double delta = LossWeight[N - 1];
						
						//U_1(myu,:)=U_1(myu,:)+alpha*(delta*(L_1(mya,:)-L_1(myb,:)));
						DoubleMatrix1D L1_mya_acc = this.L_1.viewRow(mya).copy();
						DoubleMatrix1D L1_myb = this.L_1.viewRow(myb);
						L1_mya_acc.assign(L1_myb, (x, y) -> this.alpha * delta * (x - y));
						this.U_1.viewRow(myu).assign(L1_mya_acc, (x, y) -> x + y);
						
						//U_2(myu,:)=U_2(myu,:)+alpha*(delta*(F_G(mya,:)-F_G(myb,:)));
						DoubleMatrix1D FG_mya_acc = this.F_G.viewRow(mya).copy();
						DoubleMatrix1D FG_myb = this.F_G.viewRow(myb);
						FG_mya_acc.assign(FG_myb, (x, y) -> this.alpha * delta * (x - y));
						this.U_2.viewRow(myu).assign(FG_mya_acc, (x, y) -> x + y);
						
						
			            //   L_1(mya,:)=L_1(mya,:)+alpha*(delta*(U_1(myu,:)));
						this.L_1.viewRow(mya).assign(this.U_1.viewRow(myu), (x, y) -> x + this.alpha * delta * y);
	
						
			            //   L_1(myb,:)=L_1(myb,:)+alpha*(-delta*(U_1(myu,:)));
						this.L_1.viewRow(myb).assign(this.U_1.viewRow(myu), (x, y) -> x + this.alpha * -delta * y);
						
						/*
						if(norm(U_1(myu,:))>C)
			                   U_1(myu,:)=C*U_1(myu,:)/norm(U_1(myu,:));
			               end
			             */
						double normU1_myu = CernMatrixUtils.normalizeVector(this.U_1.viewRow(myu), 2);
						if (normU1_myu > C) {
							this.U_2.viewRow(myu).assign(x -> C * x / normU1_myu);
						}
						/*
						if(norm(U_2(myu,:))>C*FactorInf)
			                   U_2(myu,:)=C*U_2(myu,:)/norm(U_2(myu,:));
			               end
			             */
						double normU2_myu = CernMatrixUtils.normalizeVector(this.U_2.viewRow(myu), 2);
						if (normU2_myu > C * this.FactorInf) {
							this.U_2.viewRow(myu).assign(x -> C * x / normU2_myu);
						}
						
						/*
			               if(norm(L_1(mya,:))>C)
			                   L_1(mya,:)=C*L_1(mya,:)/norm(L_1(mya,:));
			               end
			               */
						double normL1_mya = CernMatrixUtils.normalizeVector(this.L_1.viewRow(mya), 2);
						if (normL1_mya > C) {
							this.L_1.viewRow(mya).assign(x -> C * x / normL1_mya);
						}
						
						/*
			               if(norm(L_1(myb,:))>C)
			                   L_1(myb,:)=C*L_1(myb,:)/norm(L_1(myb,:));
			               end
						*/
						double normL1_myb = CernMatrixUtils.normalizeVector(this.L_1.viewRow(myb), 2);
						if (normL1_myb > C) {
							this.L_1.viewRow(myb).assign(x -> C * x / normL1_myb);
						}
					}
					
					
					
				});
			});
			
			//Now, check diff and those things
			/*
			iteration=iteration+1;
		       diff=norm(U_1_pre-U_1)+norm(U_2_pre-U_2)+norm(L_1_pre-L_1);
		       fprintf('\niteration=%d,diff=%f\n',iteration,diff);
		       if(iteration>1000)
		           break;
		       end 
		     */
			double normU1PrevMinusU1 = normalizeMatrix(U_1_pre.assign(this.U_1, (x, y) -> x - y));
			double normU2PrevMinusU2 = normalizeMatrix(U_2_pre.assign(this.U_2, (x, y) -> x - y));
			double normL1PrevMinusL1 = normalizeMatrix(L_1_pre.assign(this.L_1, (x, y) -> x - y));
			
			double finalNorm = normU1PrevMinusU1 + normU2PrevMinusU2 + normL1PrevMinusL1;
			
			System.out.println("Iteration " + (it+1) + " Diff: " + finalNorm);
			
		}

		
		
	}
	
	private double normalizeMatrix(DoubleMatrix2D m) {

		//Norm 2 we need to cal the singular value decomposition from apache, so we need to copy the matrix
		Array2DRowRealMatrix m2 = new Array2DRowRealMatrix(m.rows(), m.columns());
		for (int i = 0; i < m.rows(); i++) {
			for (int j = 0; j < m.columns(); j++) {
				m2.setEntry(i, j, m.getQuick(i, j));
			}
		}
		
		double [] singularVal = new SingularValueDecomposition(m2).getSingularValues();
		double max = Arrays.stream(singularVal).max().getAsDouble();
		return max;
		
	}
	

	private Tuple2<Integer, Integer> getN(int myu, int mya, DoubleMatrix2D UL, DoubleMatrix2D UFG) {
		double S_ua = UL.getQuick(myu, mya);
		double S_a_gInf = UFG.getQuick(myu, mya);
		
		
		for (int i = 0 ; i < this.numberPois; i++) {
			int myb = CernMatrixUtils.uniformRandomNumber(0, this.numberPois - 1);
			
			if (this.B.getQuick(myu, mya) > this.B.getQuick(myu, myb) ) {
				double S_ub = UL.getQuick(myu, myb);
				double S_b_gInf = UFG.getQuick(myu, myb);
				double temp = S_ua + S_a_gInf - S_ub - S_b_gInf;
				
				if (temp < this.epsilon) {
					Tuple2<Integer, Integer> t = new Tuple2<Integer, Integer>(i, myb);
					return t;
				}

			}			
		}

		return new Tuple2<Integer, Integer>(-1, -1);
	}



	@Override
	public Int2DoubleMap getScoresMap(int uidx) {
		Int2DoubleOpenHashMap scoresMap = new Int2DoubleOpenHashMap();
		scoresMap.defaultReturnValue(0.0);
		if (uidx == -1) {
			return scoresMap;
		}
		this.iIndex.getAllIidx().forEach(itemIndex -> {
			scoresMap.put(itemIndex, this.predict(uidx, itemIndex));

		});
		return scoresMap;
	}



	private double predict(int uidx, int itemIndex) {
		//(U_1(i, :) * L_1' + U_2(i, :) * F_G')
		return this.U_1.viewRow(uidx).zDotProduct(this.L_1.viewRow(itemIndex)) + this.U_2.viewRow(uidx).zDotProduct(this.F_G.viewRow(itemIndex));
	}
	
	/**
	 * Equiparable to getDist from the experimental survey
	 * It will compute the A matrix and obtain the matrix index of neighoburs
	 */
	private void generateAIndexAMatrix() {
		this.A = new DenseDoubleMatrix2D(this.numberPois, this.knnPois);
		this.aIndex = new int[this.numberPois] [this.knnPois];
		
		//Generate the A matrix
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
			A.viewRow(poiIdx).assign(locationNeighbors.stream().mapToDouble(t -> t.v2).toArray());
			aIndex[poiIdx] = locationNeighbors.stream().mapToInt(t -> t.v1).toArray();
		}
		
		//Now, recompute the A matrix normalizing it

		for (int poiIdx = 0; poiIdx < this.numberPois; poiIdx++) {
			for (int indexNeigh = 0; indexNeigh < this.knnPois; indexNeigh ++) {
				if (A.getQuick(poiIdx, indexNeigh) < 0.5) {
					A.setQuick(poiIdx, indexNeigh, 0.5);
				}
			}
			double totalRow = A.viewRow(poiIdx).copy().assign(x -> 1.0 / x).zSum();
			A.viewRow(poiIdx).assign(x -> (1/x) / totalRow);
		}		
		
	}

}
