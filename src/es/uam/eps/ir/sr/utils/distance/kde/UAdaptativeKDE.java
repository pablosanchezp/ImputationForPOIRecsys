package es.uam.eps.ir.sr.utils.distance.kde;

import java.util.List;
import java.util.Map;

import org.jooq.lambda.tuple.Tuple2;

import es.uam.eps.ir.ranksys.core.preference.IdPref;

/***
 * Adaptative KDE from GeoSoca model. 
 * 
 * Implementation based on:
 * 
 * - [2] Geographic Diversification of Recommended POIs in Frequently Visited
 * Areas (2019) : https://dlnext.acm.org/doi/abs/10.1145/3362505
 * 
 * - [3] An Experimental Evaluation of Point-of-interest Recommendation in Location-based Social Networks (2017) -> http://www.vldb.org/pvldb/vol10/p1010-liu.pdf
 * @author Pablo Sanchez (psperez@icai.comillas.edu)
 *
 */
public class UAdaptativeKDE<U, I> {
	
	//Learned parameters
	private double H1;
	private double H2;
	private double N;
	private double geometricMeanUser;
	
	//Final parameters
	private final double alpha;
	private final List<IdPref<I>> userPreferences;
	private final Map<I, Tuple2<Double, Double>> coordinatesItems;
	
	public UAdaptativeKDE(List<IdPref<I>> userPreferences, Map<I, Tuple2<Double, Double>> coordinatesItems, double alpha){
		this.userPreferences = userPreferences;
		this.coordinatesItems = coordinatesItems;
		this.alpha = alpha;
		this.N = 0;
		this.H1 = 0;
		this.H2 = 0;
		this.geometricMeanUser = 0;
		this.computeH1AndH2();

		this.geometricMeanUser = geometricMeanUser();

	}
	
	/***
	 * Method to get the pilot density of an item (eq 1 Geosoca)
	 * @param item
	 * @return
	 */
	private double pilotDensity(I item) {
		double res = 0;
		for (IdPref<I> pref: this.userPreferences) {
			res += normalKernelFunction(item, pref.v1) * pref.v2;
		}
		res /= this.N;
		return res;
	}
	
	/***
	 * Method to get the kernelDensity (eq 3 Geosoca)
	 * @param item1 the first item
	 * @param item2 the second item
	 * @return
	 */
	private double normalKernelFunction(I item1, I item2) {
		double ret = 0.0;
		Tuple2<Double, Double> latLong1 = this.coordinatesItems.get(item1);
		Tuple2<Double, Double> latLong2 = this.coordinatesItems.get(item2);

		
		double power = -1.0 * (Math.pow((latLong1.v1 - latLong2.v1), 2.0) / (2 * this.H1 * this.H1)) - 1.0 * (Math.pow((latLong1.v2 - latLong2.v2), 2.0) / (2 * this.H2 * this.H2));
		double expPow = Math.exp(power);
		ret = expPow / (2.0 * Math.PI * this.H1 * this.H2);
		
		return ret;
		
	}
	

	
	
	/**
	 * Computation of H1 and H2 (eqs 4 and 5 in GeoSoca. 
	 * Careful. original formulas of the paper are wrong. The real formula of the standard deviation 
	 * is https://www.itl.nist.gov/div898/software/dataplot/refman2/ch2/weightsd.pdf)
	 * THis implementation is based in the one provided in paper [3]
	 */
	private void computeH1AndH2() {
		this.N = 0;
		
		double weightedMeanLatitude = 0.0;
		double weightedMeanLongitude = 0.0;
		
		double accH1 = 0.0;
		double accH2 = 0.0;
		
		//Compute weighted means
		for (IdPref<I> pref: this.userPreferences) {
			Tuple2<Double, Double> latLong = this.coordinatesItems.get(pref.v1);
			weightedMeanLatitude += latLong.v1 * pref.v2;
			weightedMeanLongitude += latLong.v2 * pref.v2;
			this.N += pref.v2;
		}
		weightedMeanLatitude /= this.N;
		weightedMeanLongitude /= this.N;
		
		
		for (IdPref<I> pref: this.userPreferences) {
			Tuple2<Double, Double> latLong = this.coordinatesItems.get(pref.v1);
			accH1 += pref.v2 * Math.pow((latLong.v1 - weightedMeanLatitude), 2);
			accH2 += pref.v2 * Math.pow((latLong.v2 - weightedMeanLongitude), 2);
		}
		accH1 = Math.sqrt(accH1 / this.N);
		accH2 = Math.sqrt(accH2 / this.N);	
		
		this.H1 = 1.06 * Math.pow(this.userPreferences.size(), -0.2) * accH1;
		this.H2 = 1.06 * Math.pow(this.userPreferences.size(), -0.2) * accH2;
	}
	
	/***
	 * Method to compute the hi formula (eq 6 GeoSoca)
	 * @return
	 */
	private double hi(I item) {
		double res = (1.0 / this.geometricMeanUser) * pilotDensity(item);
		return Math.pow(res, - this.alpha);		
	} 
	
	/***
	 * GeometricMean for user. GeoSoca version (eq 7)
	 * Modified to compute using the log
	 * @return
	 */
	private double geometricMeanUser() {
		double acc = 0;
		for (IdPref<I> prefsUser: this.userPreferences) {
			acc += Math.log(pilotDensity(prefsUser.v1));
		}
		return Math.exp(acc / this.userPreferences.size());
	}
	
	/***
	 * Method to get the final density of an item (eq 8 Geosoca)
	 * @param item the item to compute the 
	 * @return
	 */
	public double finalDensity(I item) {
		double res = 0;
		for (IdPref<I> pref: this.userPreferences) {
			res += kernelDensityWithHi(item, pref.v1) * pref.v2;
		}
		res /= this.N;
		return res;
	}
	
	/***
	 * Method to get the kernelDensity with hi (eq 9 GeoSoca)
	 * @param item1 the first item
	 * @param item2 the second item
	 * @return
	 */
	private double kernelDensityWithHi(I item1, I item2) {
		double ret = 0.0;
		Tuple2<Double, Double> latLong1 = this.coordinatesItems.get(item1);
		Tuple2<Double, Double> latLong2 = this.coordinatesItems.get(item2);
		double hi = hi(item2);
		hi = hi * hi;
		
		double power = -1.0 * (Math.pow((latLong1.v1 - latLong2.v1), 2.0) / (2 * this.H1 * this.H1 * hi)) - 1.0 * (Math.pow((latLong1.v2 - latLong2.v2), 2.0) / (2 * this.H2 * this.H2 * hi));
		double expPow = Math.exp(power);
		ret = expPow / (2.0 * Math.PI * this.H1 * this.H2 * hi);
		
		return ret;
	}
	

}
