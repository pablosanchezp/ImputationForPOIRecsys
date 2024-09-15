package es.uam.eps.ir.sr.recommenders.POI.FMFMGM;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.math3.distribution.MultivariateNormalDistribution;
import org.apache.commons.math3.distribution.MultivariateRealDistribution;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.SingularMatrixException;
import org.apache.commons.math3.stat.correlation.Covariance;
import org.jooq.lambda.tuple.Tuple2;

import es.uam.eps.ir.ranksys.fast.preference.FastPreferenceData;
import es.uam.eps.ir.ranksys.rec.fast.FastRankingRecommender;
import es.uam.eps.ir.sr.data.POIProcessData;
import it.unimi.dsi.fastutil.ints.Int2DoubleMap;
import it.unimi.dsi.fastutil.ints.Int2DoubleOpenHashMap;

/***
 * 
 * Multi Center Gaussian Model for POI recommendation
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
public class MultiCenterGaussianModel <U, I extends Comparable<I>> extends FastRankingRecommender<U, I>{
	
	//
	private Map<U, List<Location>> locations;
	private Map<U, List<Center>> centerUsers;
	
	private final Map<I, Tuple2<Double, Double>> coordinatesItems;
	private final FastPreferenceData<U, I> fastPreferenceData;
	private final double alpha;
	private final double theta;
	private final double dmax;


	public MultiCenterGaussianModel(FastPreferenceData<U, I> fastPreferenceData, Map<I, Tuple2<Double, Double>> coordinatesItems, double alpha, double theta, double dmax) {
		super(fastPreferenceData, fastPreferenceData);
		
		this.locations = new HashMap<>();
		this.centerUsers = new HashMap<>();
		
		this.fastPreferenceData = fastPreferenceData;
		this.coordinatesItems = coordinatesItems;
		this.alpha = alpha;
		this.theta = theta;
		this.dmax = dmax;
		
		this.buildUserProfiles();
		this.multiCenterDiscovery();
	}
	
	/***
	 * Method to build the list of locations from the preference data
	 */
	public void buildUserProfiles() {
		this.fastPreferenceData.getUsersWithPreferences().forEach(u -> {
			List<Location> locationsUser = new ArrayList<>();
			this.fastPreferenceData.getUserPreferences(u).forEach(pref -> {
				Tuple2<Double, Double> coords = this.coordinatesItems.get(pref.v1);
				locationsUser.add(new Location(pref.v1, coords.v1, coords.v2, pref.v2, -1));
			});
			this.locations.put(u, locationsUser);
		});
		
	}
	
	/**
	 * Method to build the centers of all the users in the system
	 */
	public void multiCenterDiscovery() {
		for (U user: locations.keySet()) {
			List<Center> centersUser = discoverUserCenters(locations.get(user));
			
			List<Center> finalCenterUser = new ArrayList<>();
			for (Center c: centersUser) {
				try {
					c.buildGaussian();
					finalCenterUser.add(c);
				} catch (SingularMatrixException ex) {
					//System.out.println("Singular matrix for user " + user);
				}
			}
			this.centerUsers.put(user, finalCenterUser);
		}		
	}
	
	/***
	 * Method to compute the user centers from the list of locations
	 * @param locationsUser the list of locations of the user
	 * @return the list of centers of that user
	 */
	public List<Center> discoverUserCenters(List<Location> locationsUser){
		List<Center> centerList = new ArrayList<>();
		double centerMinimumFreq = Math.max(locationsUser.stream().mapToDouble(loc -> loc.freq).sum() * this.theta, 2);
		Collections.sort(locationsUser, Collections.reverseOrder());
		
		int centerCount = 0;
		for (int i = 0; i < locationsUser.size(); i++) {
			Location locI = locationsUser.get(i);
			if (locI.centerIdx == -1) {
				centerCount ++;
				Center center = new Center();
				center.addLocation(locI);
				locI.centerIdx = centerCount;
				
				for (int j = i + 1; j < locationsUser.size(); j++) {
					Location locJ = locationsUser.get(j);
					double distanceBetween = POIProcessData.distance(this.coordinatesItems.get(locI.id).v1, 
							this.coordinatesItems.get(locI.id).v2, this.coordinatesItems.get(locJ.id).v1, this.coordinatesItems.get(locJ.id).v2, false);
					
					if (locJ.centerIdx == -1 && distanceBetween <= this.dmax) {
						locJ.centerIdx = centerCount;
						center.addLocation(locJ);
					}
				}
				if (center.totalFreq >= centerMinimumFreq) {
					centerList.add(center);
				}
			}
		}		
		return centerList;
		
	}


	@Override
	public Int2DoubleMap getScoresMap(int uidx) {
		Int2DoubleOpenHashMap scoresMap = new Int2DoubleOpenHashMap();
		scoresMap.defaultReturnValue(0.0);
		
		U user = this.fastPreferenceData.uidx2user(uidx);

		this.iIndex.getAllIidx().forEach(itemIndex -> {
			double prob = 0.0;
			
			I item = this.fastPreferenceData.iidx2item(itemIndex);
			Tuple2<Double, Double> coordinatesItem = this.coordinatesItems.get(item);
			
			List<Center> centersUser = centerUsers.get(user);
			if (centersUser != null && coordinatesItem != null) {
				double allCentersFreqs = centersUser.stream().mapToDouble(center -> Math.pow(center.totalFreq, this.alpha)).sum();
				double allCentersPdfs = centersUser.stream().mapToDouble(c -> c.pdf(coordinatesItem.v1, coordinatesItem.v2)).sum();
				
				if (allCentersPdfs != 0) {
					for (Center c: centersUser) {
						prob += (1.0 / (POIProcessData.distance(coordinatesItem.v1, coordinatesItem.v2, c.latitude, c.longitude, false) + 1.0)) * 
								(Math.pow(c.totalFreq, this.alpha) / allCentersFreqs) * 
								(c.pdf(coordinatesItem.v1, coordinatesItem.v2) / allCentersPdfs);
						
					}
				}
				
			}
			

			scoresMap.put(itemIndex, prob);
		});
		return scoresMap;
		
	}
	
	/***
	 * Inner location class. Each location has the id, the frequency of the visitsm the center associated and the latitude and the longitude
	 * @author Pablo Sanchez (pablo.sanchezp@uam.es)
	 *
	 */
	class Location implements Comparable<Location>{
	   private I id;
	   private Double latitude;
	   private Double longitude;
	   private Double freq;
	   private int centerIdx;

	   public Location(I id, Double latitude, Double longitude, Double freq, int center){
		   this.id = id;
		   this.latitude = latitude;
		   this.longitude = longitude;
		   this.freq = freq;
		   this.centerIdx = center;
	   }

		@Override
		public int compareTo(MultiCenterGaussianModel<U, I>.Location o) {
			int res = Double.compare(this.freq, o.freq);
			if (res == 0) {
				return this.id.compareTo(o.id);
			}
			return res;
		}
	   
	}
	
	
	/***
	 * Inner Center class. Each center has the list of locations, the total Frequency, the latitude, the longitude, and the multivariate normal distribution
	 * @author Pablo Sanchez (pablo.sanchezp@uam.es)
	 *
	 */
	class Center {
		private List<Location> locations;
		private int totalFreq;
		private double [] mu;
		private double [][] covariances;
		private MultivariateRealDistribution distribution;
		private Double latitude;
		private Double longitude;
		
		public Center() {
			this.locations = new ArrayList<>();
			this.totalFreq = 0;
			this.distribution = null;
			this.mu = null;
			this.covariances = null;
			this.latitude = null;
			this.longitude = null;
		}
		
		public void addLocation(Location loc) {
			this.locations.add(loc);
			this.totalFreq += loc.freq;
		}
		
		public void buildGaussian() throws SingularMatrixException {
			List<Double []> coordinatesSequence = new ArrayList<>();
			
			for (Location loc : locations) {
				for (int i = 0; i < loc.freq; i++) {
					coordinatesSequence.add(new Double[] {loc.latitude, loc.longitude});
				}
			}
			this.mu = new double[2];
			this.mu[0] = coordinatesSequence.stream().mapToDouble(arr -> arr[0]).sum() / coordinatesSequence.size();
			this.mu[1] = coordinatesSequence.stream().mapToDouble(arr -> arr[1]).sum() / coordinatesSequence.size();
			
			Double [] arrayLatsObj = coordinatesSequence.stream().map(arr -> arr[0]).toArray(Double[]::new);
			Double [] arrayLongsObj = coordinatesSequence.stream().map(arr -> arr[1]).toArray(Double[]::new);

			double [] arrayLats = Arrays.stream(arrayLatsObj).mapToDouble(Double::doubleValue).toArray();
			double [] arrayLongs = Arrays.stream(arrayLongsObj).mapToDouble(Double::doubleValue).toArray();
			
			/***
			 * In the python approach is like this but here we do not need to make the transpose
				double [][] toCov = new double [2][arrayLats.length];
				toCov[0] = arrayLats;
				toCov[1] = arrayLongs;
			***/
			double [][] toCov = new double [arrayLats.length][2];
			for (int i = 0; i < arrayLats.length; i++) {
				toCov[i][0] = arrayLats[i];
				toCov[i][1] = arrayLongs[i];

			}

			
			RealMatrix mx = MatrixUtils.createRealMatrix(toCov);
			RealMatrix cov = new Covariance(mx).getCovarianceMatrix();
			this.covariances = cov.getData();
						
			this.distribution = new MultivariateNormalDistribution(this.mu, this.covariances);
			this.latitude = this.mu[0];
			this.longitude = this.mu[1];


		}
		
		public double pdf(Double latitude, Double longitude) {
			return this.distribution.density(new double [] {latitude, longitude});
		}
		
	}
	
}
