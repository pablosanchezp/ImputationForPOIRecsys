package es.uam.eps.ir.sr.utils.distance.kde;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.jooq.lambda.tuple.Tuple2;

import es.uam.eps.ir.crossdomainPOI.utils.KernelDensityEstimation;
import es.uam.eps.ir.sr.data.POIProcessData;

/***
 * Kernel Density Estimation by distance (the distribution will be by distance)
 * 
 * Geographical influence component from
 * 
 * - iGSLR: Personalized Geo-Social Location Recommendation - A Kernel Density Estimation Approach
 * 
 * @author Pablo Sanchez (psperez@icai.comillas.edu)
 *
 * @param <U>
 */
public class KernelDensityEstimationDistance <U> extends KernelDensityEstimation<U> {

	private Map<U, List<Double>> mapUserDistanceUsers;
	private Map<U, Double> hUsers;
	
	public KernelDensityEstimationDistance(Map<U, List<Tuple2<Double, Double>>> usersCoordinates) {
		super(usersCoordinates);
		this.mapUserDistanceUsers = new HashMap<>();
		this.hUsers = new HashMap<>();
		for (U u : this.usersCoordinates.keySet()) {
			List<Tuple2<Double, Double>> userCoordinates = this.usersCoordinates.get(u);
			List<Double> distancesUser = new ArrayList<>();
			for (int i = 0; i < userCoordinates.size(); i++) {
				Tuple2<Double, Double> coorPOIi = userCoordinates.get(i);
				for (int j = i + 1; j < userCoordinates.size(); j++) {
					Tuple2<Double, Double> coorPOIj = userCoordinates.get(j);
					distancesUser.add(POIProcessData.haversine(coorPOIi.v1, coorPOIi.v2, coorPOIj.v1, coorPOIj.v2, false));
				}
			}
			double mean = mean(distancesUser);
			double std = std(mean, distancesUser);
			hUsers.put(u, h(distancesUser, std));
			mapUserDistanceUsers.put(u, distancesUser);
			
		}
	}

	@Override
	public double probability(U user, Tuple2<Double, Double> targetCoordinatesItem) {
		double hUser = hUsers.get(user);
		List<Double> distancesUser = mapUserDistanceUsers.get(user);
		List<Tuple2<Double, Double>> userCoordinates = this.usersCoordinates.get(user);
		double res = 0;
		for (Tuple2<Double, Double> coordinate: userCoordinates) {
			 res += f(POIProcessData.haversine(targetCoordinatesItem.v1, targetCoordinatesItem.v2, coordinate.v1, coordinate.v2, false), hUser, distancesUser);
		}
		
		return res;
	}
	
	private double f(double distance, double h, List<Double> distancesUser) {
		double ret = 0;
		for (Double dist: distancesUser) {
			ret += k((distance - dist)/h);
		}
		
		
		return ret/(h * distancesUser.size());
	}
	
	private double k(double x) {
		return (1.0 / Math.sqrt(2 * Math.PI)) * Math.pow(Math.E, -(x*x)/2);
	}
	
	private double h(List<Double> distancesUser, double std) {

		
		return 1.06 * std * Math.pow(distancesUser.size(), -0.2);
	}
	
	private double std (double mean, List<Double> distancesUser) {
		double res = 0;
		for (Double dist : distancesUser) {
			res += (dist - mean) * (dist - mean);
		}
		return Math.sqrt(res/distancesUser.size());
	}
	
	private double mean(List<Double> distancesUser) {
		double res = 0;
		for (Double dist : distancesUser) {
			res += dist;
		}
		return res/distancesUser.size();
	}
	

}
