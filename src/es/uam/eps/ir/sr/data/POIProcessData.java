package es.uam.eps.ir.sr.data;

import es.uam.eps.ir.crossdomainPOI.datamodel.SimpleFastTemporalPreferenceData;
import es.uam.eps.ir.crossdomainPOI.datamodel.temporal.interfaces.FastTemporalPreferenceDataIF;
import es.uam.eps.ir.crossdomainPOI.datamodel.temporal.preferences.IdTimePref;
import es.uam.eps.ir.crossdomainPOI.utils.UsersMidPoints;
import es.uam.eps.ir.crossdomainPOI.datamodel.temporal.interfaces.TemporalPreferenceDataIF;
import es.uam.eps.ir.ranksys.core.feature.FeatureData;
import es.uam.eps.ir.ranksys.core.feature.SimpleFeatureData;
import es.uam.eps.ir.ranksys.fast.index.FastItemIndex;
import es.uam.eps.ir.ranksys.fast.index.FastUserIndex;
import es.uam.eps.ir.ranksys.fast.index.SimpleFastItemIndex;
import es.uam.eps.ir.ranksys.fast.index.SimpleFastUserIndex;
import es.uam.eps.ir.ranksys.fast.preference.FastPreferenceData;
import es.uam.eps.ir.ranksys.fast.preference.SimpleFastPreferenceData;
import es.uam.eps.ir.seqawareev.comparators.WeightComparatorTuple2;
import es.uam.eps.ir.seqawareev.comparators.WeightComparatorTuple2String;
import es.uam.eps.ir.seqawareev.comparators.WeightComparatorTuple3;
import es.uam.eps.ir.seqawareev.utils.CernMatrixUtils;
import es.uam.eps.ir.seqawareev.utils.UsersSessions;
import es.uam.eps.ir.seqawareev.utils.UsersSessions.UserSession;
import es.uam.eps.ir.sr.mains.ExperimentUtils;
import es.uam.eps.ir.sr.utils.SequentialRecommendersUtils;
import es.uam.eps.ir.sr.utils.TimeStampUtils;
import es.uam.eps.ir.sr.utils.comparators.PreferenceComparators;

import static org.ranksys.formats.parsing.Parsers.lp;
import static org.ranksys.formats.parsing.Parsers.sp;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Calendar;
import java.util.Collections;
import java.util.Comparator;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Optional;
import java.util.Set;
import java.util.TimeZone;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.AtomicLong;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

import org.apache.commons.lang.ArrayUtils;
import org.jooq.lambda.tuple.Tuple2;
import org.jooq.lambda.tuple.Tuple3;
import org.jooq.lambda.tuple.Tuple4;
import org.ranksys.formats.feature.SimpleFeaturesReader;
import org.ranksys.formats.parsing.Parser;
import org.ranksys.formats.preference.SimpleRatingPreferencesReader;

import com.google.common.util.concurrent.AtomicDouble;

/**
 * *
 * Static class to process the data of Foursqr. Contains different methods to work with the
 * Global checkin dataset of https://sites.google.com/site/yangdingqi/home/foursquare-dataset
 * The final density is in the range of the expect density - tolerance.
 * 
 * For augmenting the density remove the sessions that have the items with less popularity
 *
 * @author Pablo Sanchez (psperez@icai.comillas.edu)
 *
 */
public final class POIProcessData {
	
	
	public static <U, I> void generatePercentageCheckinsPerUserByDistance(String resultFile, FastPreferenceData<U, I> data, UsersMidPoints<U, I> usermid, Map<Long, Tuple2<Double, Double>> mapCoordinates, List<Double> listDistances) {
		try {
			PrintStream outFileToWrite = new PrintStream(resultFile);
			
			//Write the  header
			outFileToWrite.print("User\t");
			for (Double d: listDistances) {
				outFileToWrite.print(d + "\t");
			}
			outFileToWrite.println();
			
			//For every user, we compute the distance of every visited venue with respect to its midpoint
			data.getAllUsers().forEach(u -> {
				Map<Double, Integer> mapDistances = new TreeMap<>();
				for (Double d: listDistances) {
					mapDistances.put(d, 0);
				}
				
				Tuple2<Double, Double> t = usermid.getUserAvgCoordinates(u);
				
				AtomicInteger count = new AtomicInteger(0);
				
				//Compute haversine distance for every visited venue with respect to the user midpoint
				data.getUserPreferences(u).forEach(pref -> {
					Tuple2<Double, Double> t2 = mapCoordinates.get(pref.v1);
					
					Double dist = haversine(t.v1, t.v2, t2.v1, t2.v2, false);
					//Increment in 1 if th edistance is lower than the bound sent in the list
					for (Double d: mapDistances.keySet()) {
						if (dist < d) {
							//add it to the map
							mapDistances.put(d, mapDistances.get(d) + 1);
						}
						
					}
					count.incrementAndGet();
				});
				
				//We have this for the user, then we print it
				outFileToWrite.print(u + "\t");
				//Write the result (number of checkins with lower distance in each )
				for (Double d: listDistances) {
					outFileToWrite.print((float) mapDistances.get(d) + "\t");
				}
				outFileToWrite.println();
			});
			
			outFileToWrite.close();
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}
	
	
	/***
	 * In case the home city is in the completeLocalTravelersFile
	 * @param cities
	 * @param completeLocalTravelersFile
	 * @param pathOutput
	 */
	public static void generateTravelersLocals(String cities, String completeLocalTravelersFile, String pathOutput) {
		Map<String, String> userHomeCity = new HashMap<>();
		Map<String, String> localcluster = new HashMap<>();
		Map<String, String> travelercluster = new HashMap<>();
		Set<String> allUsers = new HashSet<>();

		
		try (Stream<String> stream = Files.lines(Paths.get(completeLocalTravelersFile))) {			
			//Format of the file user-id,home-city.x,travel_cluster,local_cluster
			//As we have a header we remove the first line
			stream.skip(1).forEach(line -> {
				//72718,US_Austin,non-traveler,local-2
				String [] data = line.split(",", Integer.MAX_VALUE);
				String userID = data [0];
				String homeCity = data [1];
				String travelCluster = data[2];
				allUsers.add(userID);
				
				String localCluster = data[3];
				
				if (homeCity.length() > 1) {
					userHomeCity.put(userID, homeCity);
				}
				if (travelCluster.length() > 1 && !travelCluster.equals("non-traveler") && !travelCluster.equals("unclear_home")) {
					travelercluster.put(userID, travelCluster);
				}
				if (localCluster.length() > 1) {
					localcluster.put(userID, localCluster);		
				}
			});
			stream.close();		
			
			//Now, write the local and traveler files for that city
			
			String [] citiesS = cities.split(",");
			for (String city: citiesS) {
				PrintStream outLocals = new PrintStream(pathOutput+city+"_locals.txt");
				PrintStream outTravelers = new PrintStream(pathOutput+city+"_travelers.txt");
				
				for (String user: allUsers) {
					String homeCity = userHomeCity.get(user);
					String localCluster = localcluster.get(user);
					String tCluster = travelercluster.get(user);
					
					if (homeCity != null && homeCity.equals(city) && localCluster != null) {
						outLocals.println(user + "\t" + localcluster.get(user));
					} else {
						if (tCluster != null) {
							outTravelers.println(user + "\t" + travelercluster.get(user));
						}
					}
					
				}
				outLocals.close();
				outTravelers.close();
				
			}
			
			
			
			
		} catch (Exception e) {
			System.out.println(e);
		}
		
	}
	
	/***
	 * In case the home city is in the completeLocalTravelersFile
	 * @param cities
	 * @param completeLocalTravelersFile
	 * @param pathOutput
	 */
	public static void generateTravelersLocalsHomeDifferentFile(String cities, String completeLocalTravelersFile, String homeCityFile, String pathOutput) {
		Map<String, String> userHomeCity = new HashMap<>();
		Map<String, String> localcluster = new HashMap<>();
		Map<String, String> travelercluster = new HashMap<>();
		Set<String> allUsers = new HashSet<>();
		
		

		
		try  {
			Stream<String> streamHome = Files.lines(Paths.get(homeCityFile));
			streamHome.skip(1).forEach(line -> {
				String [] data = line.split(",");
				String userID = data[0];
				String homeCityRead = data[20];
				
				if (homeCityRead.length() > 1) {
					userHomeCity.put(userID, homeCityRead);
				}
			});
			
			streamHome.close();
			
			Stream<String> stream = Files.lines(Paths.get(completeLocalTravelersFile));
			//Format of the file user-id,home-city.x,travel_cluster,local_cluster
			//As we have a header we remove the first line
			stream.skip(1).forEach(line -> {
				//72718,US_Austin,non-traveler,local-2
				String [] data = line.split(",", Integer.MAX_VALUE);
				String userID = data [0];
				String travelCluster = data[2];
				allUsers.add(userID);
				
				String localCluster = data[1];
				
	
				if (travelCluster.length() > 1 && !travelCluster.equals("non-traveler") && !travelCluster.equals("unclear_home")) {
					travelercluster.put(userID, travelCluster);
				}
				if (localCluster.length() > 1) {
					localcluster.put(userID, localCluster);		
				}
			});
			stream.close();		
			
			//Now, write the local and traveler files for that city
			
			String [] citiesS = cities.split(",");
			for (String city: citiesS) {
				PrintStream outLocals = new PrintStream(pathOutput+"\\"+city+"_locals9D.txt");
				PrintStream outTravelers = new PrintStream(pathOutput+"\\"+city+"_travelers9D.txt");
				
				for (String user: allUsers) {
					String homeCity = userHomeCity.get(user);
					String localCluster = localcluster.get(user);
					String tCluster = travelercluster.get(user);
					
					if (homeCity != null && homeCity.equals(city) && localCluster != null) {
						outLocals.println(user + "\t" + localcluster.get(user));
					} else {
						if (tCluster != null) {
							outTravelers.println(user + "\t" + travelercluster.get(user));
						}
					}
					
				}
				outLocals.close();
				outTravelers.close();
				
			}
			
			
			
			
		} catch (Exception e) {
			System.out.println(e);
		}
		
	}
	
	
	public static void filterCheckins(String completeFile, String coordFile, boolean removeConsecutiveIfSameVenue, double maxSpeed, long minDiffforBots, int minPrefsToBeConsideredBot, String outputFile) {
		FastTemporalPreferenceDataIF<Long, Long> rankSysData = ExperimentUtils.loadTrainFastTemporalFeaturePreferenceData(completeFile, completeFile, false, true);
		
		
		Map<Long, Tuple2<Double, Double>> coordinates = SequentialRecommendersUtils.POICoordinatesMap(coordFile, lp);

		try {
			
			
			PrintStream out = new PrintStream(outputFile);
			rankSysData.getUsersWithPreferences().forEach(u -> {
				//preferences user ordered by timestamp
				List<IdTimePref<Long>> preferencesUser = rankSysData.getUserPreferences(u).sorted(PreferenceComparators.timeComparatorIdTimePref).collect(Collectors.toList());
				
				

				
				//If is not bot. Now for each preference filter by max speed and if remove Consecutive 
				if (!ProcessData.isBot(preferencesUser.stream().mapToLong(pref -> pref.v3).boxed().collect(Collectors.toList()), minDiffforBots, minPrefsToBeConsideredBot)) {
					
					List<IdTimePref<Long>> preferencesUserFilterIfSame = filterListConsecutiveIfSameVenue(preferencesUser);
					List<IdTimePref<Long>> preferencesUserFilterByMaxSpeed = filterListByMaxSpeed(preferencesUserFilterIfSame, coordinates, maxSpeed);

					for (IdTimePref<Long> temporalPref: preferencesUserFilterByMaxSpeed) {
						out.println(u + "\t" + temporalPref.v1 + "\t" + temporalPref.v2 + "\t" + temporalPref.v3);
					}
					
				}
				
			});
			
			
			out.close();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} finally {
		}

	}
	
	private static List<IdTimePref<Long>> filterListByMaxSpeed(List<IdTimePref<Long>> preferencesUser, Map<Long, Tuple2<Double, Double>> coordinates, double maxSpeed){
		List<IdTimePref<Long>> result = new ArrayList<>();
		if (preferencesUser.size() == 1) {
			return preferencesUser;
		}
		
		boolean ignoreNext = false;
		
		
		for (int i = 0; i < preferencesUser.size() - 1; i++) {
			IdTimePref<Long> temporalActualPref = preferencesUser.get(i);
			IdTimePref<Long> temporalNextPref = preferencesUser.get(i + 1);
			
			//If not ignored, we continue
			if (!ignoreNext) {
				
				
				Tuple2<Double, Double> coords1 = coordinates.get(temporalActualPref.v1); 
				Tuple2<Double, Double> coords2 = coordinates.get(temporalNextPref.v1); 
				
				double diffTimeSecs = temporalNextPref.v3 - temporalActualPref.v3;
				double distanceMeters = haversine(coords1.v1, coords1.v2, coords2.v1, coords2.v2, true);
				
				//Filter by max Speed
				if (!(distanceMeters/diffTimeSecs > maxSpeed)) {
					result.add(temporalActualPref);
				} else {
					ignoreNext = true; //next is also ignored
					System.out.println("Not considering " + temporalActualPref.v1 + "\t" + temporalActualPref.v2 + "\t" + temporalActualPref.v3);
					System.out.println("Next one (also ignored) " + temporalNextPref.v1 + "\t" + temporalNextPref.v2 + "\t" + temporalNextPref.v3);
					System.out.println("--------------------------MAXSPEED-------------------------------------------");
	
				}
			} else {
				ignoreNext = false;
				System.out.println("Here we ignore " + temporalActualPref.v1 + "\t" + temporalActualPref.v2 + "\t" + temporalActualPref.v3);

			}


		}
		
		return result;
	}

	
	/***
	 * Ignore consecutive checkins performed in the same venue and take just the last
	 * @param preferencesUser
	 * @return
	 */
	private static List<IdTimePref<Long>> filterListConsecutiveIfSameVenue (List<IdTimePref<Long>> preferencesUser){
		List<IdTimePref<Long>> result = new ArrayList<>();
		if (preferencesUser.size() == 1) {
			return preferencesUser;
		}
	
		for (int i = 0; i < preferencesUser.size() - 1; i++) {
			IdTimePref<Long> temporalActualPref = preferencesUser.get(i);
			IdTimePref<Long> temporalNextPref = preferencesUser.get(i + 1);
			if (!temporalActualPref.v1.equals(temporalNextPref.v1)) {
				result.add(temporalActualPref);
			} else {
				System.out.println("Not considering " + temporalActualPref.v1 + "\t" + temporalActualPref.v2 + "\t" + temporalActualPref.v3);
				System.out.println("Next one is "  + temporalNextPref.v1 + "\t" + temporalNextPref.v2 + "\t" + temporalNextPref.v3);
				System.out.println("---------------------------------------------------------------------");

			}

		}
		
		
		return result;
	}

	
	/***
	 * Filter checking file by venue city (only consider the checkins of a specific city)
	 * @param completeFile
	 * @param POICityFIle
	 * @param selected
	 */
	public static void filterCheckingFileByItemCity(String completeFile, String POICityFIle, String selectedCity, String outputFile) {
		Map<String, String> mapPOICityFile = new HashMap<>();
		
		try {
			PrintStream outputStream = new PrintStream(outputFile);
			
			Stream<String> stream = Files.lines(Paths.get(POICityFIle));
			stream.forEach(line -> {
				String [] data = line.split("\t");
				String poiId = data [0];
				String cityPoi = data[1];
				
				mapPOICityFile.put(poiId, cityPoi);
			});
			stream.close();

			//We now filter the checkins of this file if they match the city
			Stream<String> streamTrain = Files.lines(Paths.get(completeFile));
			
			streamTrain.forEach(line -> {
				String [] data = line.split("\t");
				String poiID = data [1];
				
				String city =  mapPOICityFile.get(poiID);
				
				if (city == null) {
					System.out.println(poiID + " has no city");
				}
				
				//Only print the checkins of the city
				//if (city!=null && city.equals(selectedCity)) {
				if (city!=null && city.contains(selectedCity)) {
					outputStream.println(line);
				}
				
			});

			outputStream.close();
			streamTrain.close();

		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}
	
	/***
	 * Method to obtain an output file of the coincidences between Gowalla and Foursquare
	 * @param sourceCoordinatesFs the source coordinates file from Foursquare
	 * @param sourceCoordinatesGowalla the source coordinates file from Gowalla 
	 * @param output the outputfile
	 * @param maxDistanceDiffMeters the maximum value of the distance discrepancies (e.g 5,10 meters)
	 */
	public static void findMixGowallaFoursquare(String sourceCoordinatesFs, String sourceCoordinatesGowalla, String output, double maxDistanceDiffMeters) {
		Map<String, Tuple2<Double, Double>> poisCoordsFs = new HashMap<String, Tuple2<Double,Double>>();
		Map<String, Tuple2<Double, Double>> poisCoodsGowalla = new HashMap<String, Tuple2<Double,Double>>();
        
		try {
			//Read Foursquare
			Stream<String> stream = Files.lines(Paths.get(sourceCoordinatesFs));
        	stream.forEach(line -> {
        		String data[] = line.split("\t");
        		String poiID = data[0];
        		Double latF = Double.parseDouble(data[1]);
        		Double longF  = Double.parseDouble(data[2]);
        		poisCoordsFs.put(poiID, new Tuple2<>(latF, longF));        		
        	});
        	stream.close();
			//Read Gowalla
        	
			Stream<String> stream2 = Files.lines(Paths.get(sourceCoordinatesGowalla));
			
			stream2.forEach(line -> {
        		String data[] = line.split("\t");
        		String poiID = data[4];
        		Double latF = Double.parseDouble(data[2]);
        		Double longF  = Double.parseDouble(data[3]);
        		poisCoodsGowalla.put(poiID, new Tuple2<>(latF, longF));        		
        	});
			stream2.close();
			
			/*
			stream2.skip(1).forEach(line -> {
        		String data[] = line.split(",");
        		String poiID = data[0];
        		Double latF = Double.parseDouble(data[2]);
        		Double longF  = Double.parseDouble(data[3]);
        		poisCoodsGowalla.put(poiID, new Tuple2<>(latF, longF));        		
        	});
			stream2.close();
			*/
			
    		PrintStream outCoincidenceFile = new PrintStream(output);
    		
			//Find one in other
			poisCoordsFs.keySet().forEach(poiFs-> {
				Tuple2<Double, Double> coordinatesFs = poisCoordsFs.get(poiFs);
				String selected = "";
				Double minDistance = Double.MAX_VALUE;
				for (String idGowalla: poisCoodsGowalla.keySet()){
					Tuple2<Double, Double> coordinatesGowalla = poisCoodsGowalla.get(idGowalla);
					
					//In this case we want the distance in meters
					double distance = haversine(coordinatesFs.v1, coordinatesFs.v2, coordinatesGowalla.v1, coordinatesGowalla.v2, true);
					if (distance == 0){
						selected = idGowalla;
						break;
					}
					if (minDistance > distance && distance < maxDistanceDiffMeters) {
						selected = idGowalla;
						minDistance = distance;
					}	
				}
				if (selected != "" && minDistance != Double.MAX_VALUE) {
					outCoincidenceFile.println(poiFs + "\t" + selected + "\t" + coordinatesFs.v1 + "\t" + coordinatesFs.v2 + "\t" + minDistance);
				}
			});
    		
			outCoincidenceFile.close();
		} catch (IOException e) {
	        e.printStackTrace();
	    }

	}

    private final static double EARTH_RADIUS = 6371.0; // Approx Earth radius in KM
    
    /***
     * Read a file according to this https://snap.stanford.edu/data/loc-gowalla.html
     * @param <U>
     * @param <I>
     * @param sourceFileGowalla
     * @param outputFile
     */
    public static void processGowallaDataset(String sourceFileGowalla, String outputFile, String utputCoordFile) {

        try (Stream<String> stream = Files.lines(Paths.get(sourceFileGowalla))) {
    		PrintStream outcheckinFile = new PrintStream(outputFile);
    		PrintStream outCoordFile = new PrintStream(outputFile);

    		Set<String> pois = new HashSet<>();
    		
        	stream.forEach(line -> {
        		String [] data = line.split("\n");
        		String userS = data[0];
        		String itemS = data[4];
        		String timeS = data[1]; //"Z means UTC"
        		Double latitude = Double.parseDouble(data[2]);
        		Double longitude = Double.parseDouble(data[3]);
        		
        		Long timestamp = parseGowallaDate(timeS);          		
        		if (!pois.contains(itemS)) {
        			outCoordFile.println(itemS + "\t" + latitude + "\t" + longitude);
        		}
        		
        		outcheckinFile.println(userS + "\t" + itemS + "\t" + "1.0" + "\t" + timestamp);

        	});
        	outcheckinFile.close();
        	outCoordFile.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
    
    
    public static Long parseGowallaDate(String utcDate) {
		SimpleDateFormat formatter = new SimpleDateFormat("yyyy-MM-dd'T'hh:mm:ss'Z'",  Locale.ENGLISH);	
		formatter.setTimeZone(TimeZone.getTimeZone("UTC"));
		
		Date date = null;
        try {
            date = formatter.parse(utcDate);
            Calendar cal = Calendar.getInstance();
            cal.setTime(date);
            date = cal.getTime();
        } catch (ParseException e) {
            System.out.println("Date "+ date + " not supported");
            throw new IllegalArgumentException("Incorrect date " + date);
        }

        return date.getTime() / 1000;
    }
    
    
    public static <U, I> void computeUsersWithSameCheckins(TemporalPreferenceDataIF<U, I> prefData, String outputFile, int maxSecondsDifference) {
    	try {
    		PrintStream out = new PrintStream(outputFile);
    		prefData.getItemsWithPreferences().forEach(i -> {
    			List<IdTimePref<U>> itemPrefs = prefData.getItemPreferences(i).collect(Collectors.toList());
    			Set<U> alreadyProcessed = new HashSet<>();
    			for (IdTimePref<U> itemPref1: itemPrefs) {
    				U user1 = itemPref1.v1;
        			for (IdTimePref<U> itemPref2: itemPrefs) {
        				U user2 = itemPref2.v1;
        				//For different users whose checkins time IN THE SAME POI is lower than the max, we print it in the outputFile
        				if (!user1.equals(user2) && !alreadyProcessed.contains(user2)) {
        					double diff = Math.abs(itemPref1.v3 - itemPref2.v3);
        					if (diff <= maxSecondsDifference) {
        						out.println(user1 + "\t" + user2 + "\t" + i + "\t" + itemPref1.v3 + "\t" + itemPref2.v3 + "\t" + diff);
        					}
        				}

        			}
        			alreadyProcessed.add(user1);
    			}
    		});
    		out.close();
    	} catch (IOException e) {
    		
    	}
    }
    
    /***
     * Method that generates a dummy coordinates file
     * @throws FileNotFoundException 
     */
    public static <U, I> void generateDummyFileCoordinates(String itemMappingFile, String outputFile) throws FileNotFoundException{
    	Map<String, Long> itemMap = new HashMap<>();    	
    	readMap(itemMappingFile, itemMap);
    	
    	PrintStream out = new PrintStream(outputFile);
    	for (String oldId: itemMap.keySet()) {
    		out.println(itemMap.get(oldId) + "\t" + 0.0 + "\t" + 0.0);
    	} 
    	out.close();

    	
    }

    
    private static <U, I> void printSessions(Map<U, List<List<IdTimePref<I>>>> sessions, String outputFile) throws FileNotFoundException {
    	PrintStream out = new PrintStream(outputFile);
    	int sessionID = 1;
    	
		for (U user: sessions.keySet()) {
	    	List<List<IdTimePref<I>>> listUsession = sessions.get(user);

			for (List<IdTimePref<I>> lst : listUsession) {
				
	    		for (IdTimePref<I> pref: lst) {
	    			out.println(user + "\t" + pref.v1 + "\t" + pref.v2 + "\t" + pref.v3 + "\t" + sessionID);
	    		}
	    		sessionID++;
			}
		}
		
		out.close();
    }
    
    public static <U, I extends Comparable<I>> void subsetDensityByInversePopularitySessionUsers(double expectedDensity, double maxDiffBetweenItems, UsersSessions<U, I> userSessions, FastTemporalPreferenceDataIF<U, I> dataWithUsersThatCannotBeRemoved, String outputFile) throws FileNotFoundException {
    	double tolerance = 0.005;
    	Double density = Double.MIN_VALUE;

    	FastTemporalPreferenceDataIF<U, I> data = userSessions.getData();
    	Map<U, List<List<IdTimePref<I>>>> sessions = new HashMap<>();
    	data.getUsersWithPreferences().forEach(u -> {
    	    sessions.put(u, userSessions.getUserSession(u, maxDiffBetweenItems, Double.MAX_VALUE, UserSession.TIME));
    	});
    	List<U> usersByLengthSessions = data.getUsersWithPreferences().map(u -> new Tuple3<>(data.user2uidx(u), (double) sessions.get(u).size(), (double) data.getUserPreferences(u).count())).sorted(new WeightComparatorTuple3()).map(t3 -> data.uidx2user(t3.v1)).collect(Collectors.toList());    	
    	
    	do {
    		density = computeDensitySessions(sessions);

			U candidate = usersByLengthSessions.remove(0);

    		if (dataWithUsersThatCannotBeRemoved == null || !dataWithUsersThatCannotBeRemoved.containsUser(candidate)) {
    			sessions.remove(candidate);
    		}


    		
    	}while (density + tolerance < expectedDensity && usersByLengthSessions.size() > 0);
    	
    	//Final stats
    	density = computeDensitySessions(sessions); 
    	
    	printSessions(sessions, outputFile);
    }
    
    
    
    /****
     * Method to obtain a subset with higher density from other dataset 
     * 
     * @param expectedDensity the expected density that the new dataset should have
     * @param userSessions the users sessions
     * @param outputFile the outputFile
     * @throws FileNotFoundException
     */
    public static <U, I extends Comparable<I>> void subsetDensityByInversePopularityItems(double expectedDensity, double maxDiffBetweenItems, UsersSessions<U, I> userSessions, FastTemporalPreferenceDataIF<U, I> dataWithUsersThatCannotBeRemoved, String outputFile) throws FileNotFoundException {
    	double tolerance = 0.005;
    	Double density = Double.MIN_VALUE;

    	FastTemporalPreferenceDataIF<U, I> data = userSessions.getData();
    	Map<U, List<List<IdTimePref<I>>>> sessions = new HashMap<>();
    	List<I> itemsByPopularity = data.getItemsWithPreferences().map(i -> new Tuple2<>(data.item2iidx(i), (double) data.getItemPreferences(i).count())).sorted(new WeightComparatorTuple2()).map(t -> data.iidx2item(t.v1)).collect(Collectors.toList());
    	
    	data.getUsersWithPreferences().forEach(u -> {
    		sessions.put(u, userSessions.getUserSession(u, maxDiffBetweenItems, Double.MAX_VALUE, UserSession.TIME));
    	});
    	
    	do {
    		density = computeDensitySessions(sessions);    	
    		I candidate = itemsByPopularity.remove(0);
    		for (U user: sessions.keySet()) {
    			if (dataWithUsersThatCannotBeRemoved == null || !dataWithUsersThatCannotBeRemoved.containsUser(user)) {
	    	    	List<List<IdTimePref<I>>> listUsession = sessions.get(user);
	    	    	
	    	    	for (Iterator<List<IdTimePref<I>>> iter = listUsession.listIterator(); iter.hasNext(); ) {    			
	    				Optional<IdTimePref<I>> op = iter.next().stream().filter(pref -> pref.v1.equals(candidate)).findFirst();
	    				if (op.isPresent()) {
	    					iter.remove();
	    				}
	    			}
	    	    	
    			}
    			
    		}

    		data.getUsersWithPreferences().forEach(u -> {
    		    if (sessions.containsKey(u) && sessions.get(u).isEmpty()) {
    		    	sessions.remove(u);
    		    }	
    		});

    		
    	}while (density + tolerance < expectedDensity && itemsByPopularity.size() > 0);
    	
    	//Final stats
    	density = computeDensitySessions(sessions);    	

    	
    	printSessions(sessions, outputFile);

    }

    
    /****
     * Method to obtain a subset from other dataset with a similar density to the one indicated in the expectedDensity parameter.
     * The final density is in the range of the expect density - tolerance.
     * 
     * For augmenting the density it will remove sessions selected randomly
     *  
     * @param expectedDensity
     * @param userSessions
     * @param outputFile
     * @throws FileNotFoundException
     */
    public static <U, I extends Comparable<I>> void subsetDensityRandom(double expectedDensity, double maxTimeBetweenItems, UsersSessions<U, I> userSessions, FastTemporalPreferenceDataIF<U, I> dataWithUsersThatCannotBeRemoved, String outputFile) throws FileNotFoundException {
    	double tolerance = 0.005;
    	Double density = Double.MIN_VALUE;
    	FastTemporalPreferenceDataIF<U, I> data = userSessions.getData();
    	Map<U, List<List<IdTimePref<I>>>> sessions = new HashMap<>();
    	
    	data.getUsersWithPreferences().forEach(u -> {
    		sessions.put(u, userSessions.getUserSession(u, maxTimeBetweenItems, Double.MAX_VALUE, UserSession.TIME));
    	});
    	
    	do {
    		density = computeDensitySessions(sessions);    		
    		int randomSession = -1;
    		do {
    			int randomUser = CernMatrixUtils.uniformRandomNumber(0, data.numUsers() - 1);
    			U user = data.uidx2user(randomUser);
    			if (dataWithUsersThatCannotBeRemoved == null || !dataWithUsersThatCannotBeRemoved.containsUser(user)) {
	    			if (sessions.containsKey(user)) {
	    				randomSession = CernMatrixUtils.uniformRandomNumber(0, sessions.get(user).size() - 1);
	    	    		sessions.get(user).remove(randomSession);
	    	    		if (sessions.get(user).isEmpty()) {
	    	    			sessions.remove(user);
	    	    		}
	    			}
    			}
    			
    		} while (randomSession == -1);
    	}while (density + tolerance < expectedDensity);
    	
    	printSessions(sessions, outputFile);

    }
    
    private static <U, I> Double computeDensitySessions(Map<U, List<List<IdTimePref<I>>>> sessions) {
    	double numUsers = sessions.keySet().size();
    	double numItems = 0;
    	double numPreferences = 0;
    	Set<I> setItems = new HashSet<>();
    	for (U user: sessions.keySet()) {
    		List<List<IdTimePref<I>>> sessionsUser = sessions.get(user);
    		for (List<IdTimePref<I>> lstSess: sessionsUser) {
    			for (IdTimePref<I> pref: lstSess) {
    				setItems.add(pref.v1);
    				numPreferences++;
    			}
    		}
    		
    	}
    	numItems = setItems.size();
    	double density = numPreferences / (numUsers * numItems);
    	System.out.println("Users: " + numUsers);
    	System.out.println("Items: " + numItems);
    	System.out.println("Preferences: " + numPreferences);
    	System.out.println("Density: " + density);
		return density;
	}



	/***
     * Method to obtain a map of categories. This categories are hand-made
     * @param fileCategories the file of categories
     * @return a map of grouped categories
     */
    public static Map<String, List<String>> groupingFoursquareCategories(String fileCategories) {
        Map<String, List<String>> result = new TreeMap<String, List<String>>();

        try {
            Stream<String> stream = Files.lines(Paths.get(fileCategories));
            stream.forEach(line -> {
                String key = getKeyAssociated(line);
                if (result.get(key) == null) {
                    result.put(key, new ArrayList<String>());
                }
                result.get(key).add(line);
            });
            stream.close();
        } catch (IOException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
        return result;
    }

    private static String getKeyAssociated(String line) {
        String line2 = line.toLowerCase();
        if (line2.contains("Restaurant".toLowerCase()) || line2.contains("Bar".toLowerCase())) {
            return "Restaurants and Bars";
        } else if (line2.contains("School".toLowerCase()) || line2.contains("College".toLowerCase()) || line2.contains("University".toLowerCase())) {
            return "Education";
        } else if (line2.contains("Airport".toLowerCase())) {
            return "Airport";
        } else if (line2.contains("Museum".toLowerCase()) || line2.contains("Art".toLowerCase()) || line2.contains("Movie".toLowerCase())) {
            return "Art & Museums";
        } else if (line2.contains("Soccer".toLowerCase()) || line2.contains("Football".toLowerCase())
                || line2.contains("Basketball".toLowerCase()) || line2.contains("Baseball".toLowerCase())
                || line2.contains("Tennis".toLowerCase()) || line2.contains("Track".toLowerCase()) || line2.contains("Volleyball".toLowerCase())) {
            return "Sports";
        } else if (line2.contains("Temple".toLowerCase()) || line2.contains("Church".toLowerCase())
                || line2.contains("Synagogue".toLowerCase()) || line2.contains("Mosque".toLowerCase())) {
            return "Temples";
        } else if (line2.contains("Store".toLowerCase()) || line2.contains("Shop".toLowerCase()) || line2.contains("Market".toLowerCase())) {
            return "Stores";
        }

        return "Other";
    }
    
    /***
     * Method to compute the statistics of the train files of the foursqr dataset
     * @param checkingTrainFileRep
     * @param checkingTrainFileNoRep
     * @param poiCitiesFile
     * @param resultStats
     */
    public static void statisticsFoursquare(String checkingTrainFileRep, String checkingTrainFileNoRep, String testFile, String poiCitiesFile, String resultStats) {
		String characterSplit = "\t";
    	try {
    		Map<String, String> mapPoisCities= new HashMap<>();
    		Stream<String> streamCitiesFile = Files.lines(Paths.get(poiCitiesFile));
    		streamCitiesFile.forEach(line -> {
    			String [] poiCity = line.split(characterSplit);
    			mapPoisCities.put(poiCity[0], poiCity[1]);
    		});
    		streamCitiesFile.close();
    		
            PrintStream resultFile = new PrintStream(resultStats);

            if (checkingTrainFileRep != null) {
            	processDataForStatistics(checkingTrainFileRep,mapPoisCities,resultFile);
            }
            if (checkingTrainFileRep != null)  {
            	processDataForStatistics(checkingTrainFileNoRep,mapPoisCities,resultFile);
            }
            if (testFile != null) {
            	processDataForStatistics(testFile,mapPoisCities,resultFile);
            }
            
			resultFile.close();
			
			
			
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
    	
    }
    
    
    private static void processDataForStatistics(String file, Map<String, String> mapPoisCities,PrintStream resultFile) {
    	int itemIndexRating = 1;
		int userIndexRating = 0;
		int timesStampIndex = 3;
		String characterSplit = "\t";
		try {
	    	Stream<String> lines = Files.lines(Paths.get(file));
			Map<String, Integer> users = new HashMap<>();
			Map<String, Integer> items = new HashMap<>();
			HashSet<String> mapRepetitions= new HashSet<>();
			Set<String> cities= new HashSet<>();
			AtomicInteger ratings = new AtomicInteger(0);
			AtomicInteger repetitions = new AtomicInteger(0);
			
			AtomicLong maxTimestamp= new AtomicLong(0L);
			AtomicLong minTimestamp= new AtomicLong(Long.MAX_VALUE);
			lines.forEach(line -> {
				String [] data = line.split(characterSplit);
				if (users.get(data[userIndexRating]) == null) {
					users.put(data[userIndexRating], 0);
				}
				
				if (items.get(data[itemIndexRating]) == null) {
					items.put(data[itemIndexRating], 0);
				}
				users.put(data[userIndexRating], users.get(data[userIndexRating]) + 1);
				items.put(data[itemIndexRating], items.get(data[itemIndexRating]) + 1);
				cities.add(mapPoisCities.get(data[itemIndexRating]));
				ratings.getAndIncrement();
				String key = data[userIndexRating] + "_"+data[itemIndexRating];
				if (mapRepetitions.contains(key)) {
					repetitions.incrementAndGet();
				}
				else {
					mapRepetitions.add(key);
				}
				if (data.length > timesStampIndex) {
					Long timeParsed = Long.parseLong(data[timesStampIndex]);
					if (maxTimestamp.get() < timeParsed) {
						maxTimestamp.set(timeParsed);
					}
					if (minTimestamp.get() > timeParsed) {
						minTimestamp.set(timeParsed);
					}
				}
				
				
			});
			lines.close();
			File f = new File(file);
			String filename= f.getName();
			SimpleDateFormat formatter = new SimpleDateFormat("EEE MMM dd HH:mm:ss Z yyyy",  Locale.ENGLISH);		

			resultFile.println("File: "+filename);
			resultFile.println("Number of ratings: " + ratings.get());
			resultFile.println("Number of users: " + users.keySet().size());
			resultFile.println("Number of items: " + items.keySet().size());
			resultFile.println("Average ratings per user " + ( (double) ratings.get()) / ( (double) users.keySet().size()));
			resultFile.println("Average ratings per items " + ( (double) ratings.get()) / ( (double) items.keySet().size()));
			resultFile.println("Total repetitions " + repetitions.get());
			resultFile.println("Different cities: ");
			if (maxTimestamp.get() != 0) {
				resultFile.println("Max timestamp: " + maxTimestamp.get() + "Date: " + formatter.format(new Date(maxTimestamp.get()*1000)));
				resultFile.println("Min timestamp: " + minTimestamp.get() + "Date: " + formatter.format(new Date(minTimestamp.get()*1000)));
			}
			for (String city: cities) {
				resultFile.print(city+" ");
			}
			resultFile.println("--------------------------");
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
    }
    
    
    
    /***
     * Method to generate a file of the center of each city by reading all its Pois associated
     * @param poiCoordsFile the coordinates file of the pois
     * @param poiCityFile the city file of the pois
     * @param resultfilepath the result file
     */
    public static void generateCityCenter(String poiCoordsFile, String poiCityFile, String resultfilepath) {
		try {
			String characterSplit="\t";
			Map<String,String> poiCityMap= new HashMap<>();
			Map<String, List<Tuple2<Double, Double>>> mapCity = new HashMap<String, List<Tuple2<Double, Double>>>();
			Stream<String> streamPoiCityFile = Files.lines(Paths.get(poiCityFile));
			streamPoiCityFile.forEach(line -> {
				String [] poiCity= line.split(characterSplit);
				poiCityMap.put(poiCity[0], poiCity[1]); // map POI -> city
			});
			streamPoiCityFile.close();
			
			//We have read all the POIs 
			Stream<String> streamPoiCoords = Files.lines(Paths.get(poiCoordsFile));
			streamPoiCoords.forEach(line -> {
				String [] poiCoor = line.split(characterSplit);
				String idPoi = poiCoor[0];
				String cityOfPoi = poiCityMap.get(idPoi);
				if (mapCity.get(cityOfPoi) == null)
					mapCity.put(cityOfPoi, new ArrayList<Tuple2<Double,Double>>());
				mapCity.get(cityOfPoi).add(new Tuple2<>(Double.parseDouble(poiCoor[1]), Double.parseDouble(poiCoor[2])));
			});
			streamPoiCoords.close();
			
            PrintStream resultFile = new PrintStream(resultfilepath);
			//Now we have for each city, a list of coordinates
			for (String city: mapCity.keySet()) {
				Tuple2<Double, Double> midPointCity = SequentialRecommendersUtils.midPointCoordinates(mapCity.get(city));
				resultFile.println(city+ characterSplit+midPointCity.v1+characterSplit+midPointCity.v2);				
			}
			resultFile.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
    }
    /**
     * Same method of get nnPois but with cities
     * @param coorCities
     * @param nn
     * @param result
     */
    public static void getNNCities(String coorCities, int nn, String fileResult) {
    	List<Tuple3<String,Double,Double>> lstCities = new ArrayList<>();
		String characterSplit="\t";

		try {
			Stream<String> streamPoiCoords = Files.lines(Paths.get(coorCities));
			streamPoiCoords.forEach(line -> {
				String [] tokens = line.split(characterSplit);
				lstCities.add(new Tuple3<>(tokens[0],Double.parseDouble(tokens[1]),Double.parseDouble(tokens[2])));
			});
			streamPoiCoords.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
    	
        getClosestPOIs(fileResult, lstCities, nn);
    }
    
    /**
     * Method to filter the checkins file by a country code
     * @param checkingFile the checking file
     * @param countryCodeFilter the country code 
     * @param filePoiCountryCode the file of the poi and the city
     * @param resultfilepath the result file
     */
    public static void filterCheckingFileByCountryCode(String checkingFile, String countryCodeFilter, String filePoiCountryCode, String resultfilepath) {
		try {
			String characterSplit="\t";
			String characterSplit2="_";			
			Set<String> validPois = new HashSet<String>(); //only pois matching that country code
			Stream<String> stream = Files.lines(Paths.get(filePoiCountryCode));
			stream.forEach(string ->{
				String [] data = string.split(characterSplit);
				String [] cCode = data[1].split(characterSplit2);
				if (cCode[0].equals(countryCodeFilter))
					validPois.add(data[0]);
			});
			stream.close();
			//Now filter the original checkin file
			
            PrintStream resultFile = new PrintStream(resultfilepath);
            Stream<String> stream2 = Files.lines(Paths.get(checkingFile));
			stream2.forEach(fullData -> {
				String idPoi = fullData.split(characterSplit)[1];
				if (validPois.contains(idPoi))
					resultFile.println(fullData);
			});
			stream2.close();
			resultFile.close();
			
			
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
    	
    }
    
    public static void filterByNNCities(String checkingFile, String city, String nnFile,int nnCitiesToConsider, String filePoiCountryCode, String resultfilepath) {
    	String characterSplit="\t";
		Set<String> validPois = new HashSet<String>(); //only pois matching the nn cities
		Set<String> validCities = new HashSet<String>(); //set to store the cities
    	try {
			Stream<String> stream = Files.lines(Paths.get(nnFile));
			Optional<String> matched = stream.filter(line -> line.split(characterSplit)[0].equals(city)).findAny();
			if (!matched.isPresent()) {
				System.out.println("No city in the file matchin "+ city);
				stream.close();
				return;
			}
			String [] tokens= matched.get().split(characterSplit);
			//The limit will be the lower in the nn cities to consider and the length of columns 
			int limit = (nnCitiesToConsider > tokens.length - 1) ? tokens.length - 1 : nnCitiesToConsider;
			System.out.println("NN cities of "+ city + " are:");
			for (int i = 0; i < limit+1; i++) {
				validCities.add(tokens[i]);
				System.out.println(tokens[i]);
			}
			stream.close();
			
			//Now we have all the valid cities. We must read the filePoiCountryCode and store only the pois of the cities we want
			
			Stream<String> stream2 = Files.lines(Paths.get(filePoiCountryCode));
			stream2.forEach(string ->{
				String [] data = string.split(characterSplit);
				if (validCities.contains(data[1]))
					validPois.add(data[0]);
			});
			stream2.close();
			
			PrintStream resultFile = new PrintStream(resultfilepath);
            Stream<String> stream3 = Files.lines(Paths.get(checkingFile));
			stream3.forEach(fullData -> {
				String idPoi = fullData.split(characterSplit)[1];
				if (validPois.contains(idPoi))
					resultFile.println(fullData);
			});
			stream3.close();
			resultFile.close();
			
			
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
    	
    	
		
	}
    
    
    
    /**
     * Method to analyze the distances between the most popular items of a dataset
     * @param trainFile the training file
     * @param poiCoordsFile the 3 column file cotaining the identifiers of the items and the latitudes and longitudes
     * @param nMostPopular the n items most popular we want to compute
     */
    public static void analyzeDistanceMostPopular(String trainFile, String poiCoordsFile, int nMostPopular) {
    	try {
    		String characterSplit="\t";
    		int indexPoi=1;
    		Map<String, Tuple2<Double,Double>> mapCoordinates = new HashMap<>();
			Stream<String> stream = Files.lines(Paths.get(poiCoordsFile));
			stream.forEach(line ->{
				String [] split = line.split(characterSplit);
				mapCoordinates.put(split[0], new Tuple2<>(Double.parseDouble(split[1]),Double.parseDouble(split[2])));
			});
			stream.close();
			//We have now stored all the coordinates, now we compute the mos popular
			Map<String, Integer> mapPreferences = new HashMap<>();
			Stream<String> stream2 = Files.lines(Paths.get(trainFile));
			stream2.forEach(line -> {
				String [] split = line.split(characterSplit);
				String poi = split[indexPoi];
				if (mapPreferences.get(poi) == null)
					mapPreferences.put(poi, 0);
				mapPreferences.put(poi, mapPreferences.get(split[indexPoi]) + 1);
			});			
			stream2.close();
			List<Tuple2<String, Integer>> lstOfPreferences = new ArrayList<>();
			for (String poi: mapPreferences.keySet())
				lstOfPreferences.add(new Tuple2<>(poi, mapPreferences.get(poi)));
			Comparator<Tuple2<String,Integer>> comparatorTuples= (o1,o2) ->  Double.compare(o1.v2, o2.v2);
			lstOfPreferences.sort(comparatorTuples.reversed());
			//Now we have the list of most popular computed, we take the nn sublist
			List<Tuple2<String, Integer>> lstOfPreferences2 = lstOfPreferences.subList(0, nMostPopular);
			AtomicDouble total = new AtomicDouble(0);
			//This could be optimized to only compute half of the matrix
			lstOfPreferences2.stream().forEach(pref -> {
				lstOfPreferences2.stream().forEach(pref2 -> {
					if (!pref.v1.equals(pref2.v1)) {
						double distance = haversine(mapCoordinates.get(pref.v1).v1, mapCoordinates.get(pref.v1).v2, mapCoordinates.get(pref2.v1).v1, mapCoordinates.get(pref2.v1).v2, false);
						System.out.println(distance);
						total.set(total.get() + distance);
					}
						
				});
			});
			System.out.println("Average distance of top most popular "+ total.get()/(lstOfPreferences2.size()* lstOfPreferences2.size() - lstOfPreferences2.size()));
	
			AtomicDouble total2 = new AtomicDouble(0);
			AtomicInteger count = new AtomicInteger(0);
			lstOfPreferences.stream().forEach(pref -> {
				count.getAndIncrement();
				if ((double) count.get() / lstOfPreferences.size() > 0.1) {
					System.out.println("10%");
					count.set(0);
				}
				lstOfPreferences.stream().forEach(pref2 -> {
					if (!pref.v1.equals(pref2.v1)) {
						double distance = haversine(mapCoordinates.get(pref.v1).v1, mapCoordinates.get(pref.v1).v2, mapCoordinates.get(pref2.v1).v1, mapCoordinates.get(pref2.v1).v2, false);
						total2.set(total2.get() + distance);
					}
						
				});
			});
			System.out.println("Average distance of full dataset "+ total2.get()/(lstOfPreferences.size()* lstOfPreferences.size() - lstOfPreferences.size()));

			
			
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
    	
    	
    }
    


    

    public static Map<String, String[]> readFoursquareCategories(String inputFile) throws IOException {
        Map<String, String[]> categoriesWithLevels = new HashMap<>();
        Stream<String> stream = Files.lines(Paths.get(inputFile), StandardCharsets.ISO_8859_1);
        // this file has a header, so we can skip it
        stream.skip(1).forEach(line -> {
            // level1_name,level1_id,level2_name,level2_id,level3_name,level3_id,level4_name,level4_id
            String[] data = line.split(",", -1);
            //System.out.println(data.length + "<->" + line);
            String level1 = data[0];
            String level2 = data[2];
            String level3 = data[4];
            String level4 = data[6];
            List<String> values = new ArrayList<>();
            String key = level4;
            if (key.isEmpty()) {
                key = level3;
                if (key.isEmpty()) {
                    key = level2;
                    if (key.isEmpty()) {
                        key = level1;
                        // this level should never be empty!
                        values.add(level1);
                    } else {
                        values.add(level1);
                        values.add(level2);
                    }
                } else {
                    values.add(level1);
                    values.add(level2);
                    values.add(level3);
                }
            } else {
                values.add(level1);
                values.add(level2);
                values.add(level3);
                values.add(level4);
            }
            String[] valuesAsArray = new String[values.size()];
            values.toArray(valuesAsArray);
            categoriesWithLevels.put(key, valuesAsArray);
        });
        stream.close();
        return categoriesWithLevels;
    }

    private static String parseFoursquareCheckinCategory(String category) {
        switch (category) {
            case "Light Rail":
                return "Light Rail Stations";
            case "College & University":
                return "Colleges & Universities";
            case "Yogurt":
                return "Frozen Yogurt";
            case "Car Dealership":
                return "Auto Dealerships";
            case "Deli / Bodega":
                return "Delis / Bodegas";
            case "Athletic & Sport":
                return "Athletics & Sports";
            case "Subway":
                return "Metro Stations";
            case "Mall":
                return "Shopping Malls";
            case "Spa / Massage":
                return "Spas";
            case "Home (private)":
                return "Homes (private)";
            case "Gym / Fitness Center":
                return "Gyms or Fitness Centers";
            case "Gas Station / Garage":
                return "Gas Stations";
            case "Shop & Service":
                return "Shops & Services";
            case "Tennis":
                return "Tennis Stadiums";
            case "Residential Building (Apartment / Condo)":
                return "Residential Buildings (Apartments / Condos)";
            case "Hiking Trail":
                return "Trails";
            case "Boat or Ferry":
                return "Boats or Ferries";
            case "Embassy / Consulate":
                return "Embassies / Consulates";
            case "Bike Rental / Bike Share":
                return "Bike Rentals / Bike Shares";
            case "General College & University":
                return "General Colleges & Universities";
            case "Drugstore / Pharmacy":
                return "Pharmacies";
            case "Ferry":
                return "Boats or Ferries";
            case "Salon / Barbershop":
                return "Salons / Barbershops";
            case "Malaysian Restaurant":
                return "Malay Restaurants";
            case "Harbor / Marina":
                return "Harbors / Marinas";
            case "Theme Park Ride / Attraction":
                return "Theme Park Rides/Attractions";
            case "Ramen /  Noodle House":
                // this could also be "Noodle Houses"
                return "Ramen Restaurants";
            case "Monument / Landmark":
                return "Monuments / Landmarks";
        }
        if (category.startsWith("Caf") && category.length() < 6) {
            return "Cafés";
        }
        return category;
    }

    public static String matchingFoursquareCheckinCategory(String category, int level, Map<String, String[]> categoriesWithLevels) {
        // some plurals are not in the file
        String[] categoriesToMatch = new String[]{category, category + "s", category + "es", category.substring(0, category.length() - 1) + "ies", parseFoursquareCheckinCategory(category)};
        for (String c : categoriesToMatch) {
            String[] categories = categoriesWithLevels.get(c);
            if (categories != null) {
                // there is a matching: let's find the closest level
                while (level > categories.length) {
                    level--;
                }
                String matchedCategory = categories[level - 1];
                return matchedCategory;
            }
        }
        return null;
    }

    /**
     * *
     * Method to print the nn closest POIS to each POI of the dataset
     * The result file will contain the nn closest pois of each poi in the training dataset
     * @param trainingFile the file of total POIs (used to filter)
     * @param poisCoordsFile the poisCoordsFile
     * @param resultFileNN the result file
     * @param nn the number of neighbours to print
     */
    public static void printClosestPOIs(String trainingFile, String poisCoordsFile, String resultFileNN, String resultFileCoordsOfTrain, int nn) {

        String charactersSplit = "\t";
        int indexPoiTrainFile = 1;

        int indexPoiCoordsFile = 0;
        int indexLat = 1;
        int indexLong = 2;
        Set<String> trainingPOIs = new HashSet<String>();
        try {
            PrintStream resultFileCoordsFiltered = new PrintStream(resultFileCoordsOfTrain);
            //Get all training POIs. We will not compute distances of more POIs
            Stream<String> stream = Files.lines(Paths.get(trainingFile));
            stream.forEach(line -> {
                String[] data = line.split(charactersSplit);
                trainingPOIs.add(data[indexPoiTrainFile]);
            });
            stream.close();

            //Now we read the POIs file
            List<Tuple3<String, Double, Double>> poisTuple = new ArrayList<>();
            Stream<String> stream2 = Files.lines(Paths.get(poisCoordsFile));
            stream2.forEach(line -> {
                String[] data = line.split(charactersSplit);
                String poiIdentifier = data[indexPoiCoordsFile];
                if (trainingPOIs.contains(poiIdentifier)) {
                    poisTuple.add(new Tuple3<>(poiIdentifier, Double.parseDouble(data[indexLat]), Double.parseDouble(data[indexLong])));
                    resultFileCoordsFiltered.println(line);
                }
            });
            stream2.close();
            resultFileCoordsFiltered.close();
            //We have the list, we retrieve the map
            getClosestPOIs(resultFileNN, poisTuple, nn);

        } catch (IOException e) {

        }
    }

    /**
     * Receiving a list of POIs, will return a map of the closest nn pois
     *
     * @param lst the list of POIS represented with a tuple (ID, lat and long)
     * @param nn the number of POIS to be returned
     * @return
     */
    public static void getClosestPOIs(String resultFileNN, List<Tuple3<String, Double, Double>> lst, int nn) {
        try {
            PrintStream resultFile = new PrintStream(resultFileNN);
            lst.stream().forEach(tuple3Poi -> {
                List<Tuple2<String, Double>> lstReturn = new ArrayList<>();

                lst.stream().filter(tuple3Poi2 -> tuple3Poi2.v1 != tuple3Poi.v1).forEach(tuple3Poi2 -> {
                    lstReturn.add(new Tuple2<>(tuple3Poi2.v1, haversine(tuple3Poi.v2, tuple3Poi.v3, tuple3Poi2.v2, tuple3Poi2.v3, false)));
                });
                lstReturn.sort(closerPois); //Each element of the list will contain id poi and the distance. The lower the distance the closest the poi

                List<Tuple2<String, Double>> rest = lstReturn.subList(0, nn);
                resultFile.print(tuple3Poi.v1);
                for (Tuple2<String, Double> neigh : rest) {
                    resultFile.print("\t" + neigh.v1);
                }
                resultFile.println();
            });
            resultFile.close();
        } catch (Exception e) {

        }
    } 
    
    public static void computeAvgDistanceCities(String midPoincities, String sevenNcities, String resultFile) {
        Map<String, Tuple2<Double, Double>> citiesMidPoints = new HashMap<>();
    	
        try {
            PrintStream resultStr = new PrintStream(resultFile);

        	Stream<String> stream = Files.lines(Paths.get(midPoincities));
        	 stream.forEach(line -> {
                 String[] data = line.split("\t");
                 citiesMidPoints.put(data[0], new Tuple2<>(Double.parseDouble(data[1]), Double.parseDouble(data[2])));
             });
             stream.close();
             
             //Now, read the n closest cities
             Stream<String> stream2 = Files.lines(Paths.get(sevenNcities));
        	 stream2.forEach(line -> {
        		 double distance = 0;
                 String[] data = line.split("\t");
                 String dataOrigin = data [0];
                 Tuple2<Double, Double> city1 = citiesMidPoints.get(dataOrigin);
                 for (int i = 1 ; i < data.length; i++) {
                     Tuple2<Double, Double> city2 = citiesMidPoints.get(data[i]);
                	 distance += haversine(city1.v1, city1.v2, city2.v1, city2.v2, false);
                 }
                 resultStr.println(dataOrigin + "\t" + distance/data.length);

        	 });
             stream2.close();
             resultStr.close();

        } catch (IOException e) {

        }
    }

    

    private static Comparator<Tuple2<String, Double>> closerPois = (o1, o2) -> Double.compare(o1.v2, o2.v2);



    /**
     * *
     * Method to generate a feature file (IdPoi, idFeature, Rating)
     *
     * @param originalPoisFile the original pois file from foursqr
     * @param categoriesMapFile the map of categories
     * @param poisMapFile the map of pois
     * @param featureDestFile the destination path of the file
     */
    public static void generateFeatureFile(String originalPoisFile, String poisMapFile, String categoriesMapFile,
            String featureDestFile) {

        Map<String, Long> categoriesMap = new HashMap<>();
        readMap(categoriesMapFile, categoriesMap);

        Map<String, Long> poisMap = new HashMap<>();
        readMap(poisMapFile, poisMap);

        try {
            PrintStream resultFile = new PrintStream(featureDestFile);

            String characterSplit = "\t";

            int columnVenueId = 0;
            int columnvenueCategory = 1;

            Stream<String> stream = Files.lines(Paths.get(originalPoisFile), StandardCharsets.ISO_8859_1);
            stream.forEach(line -> {
                String[] originalData = line.split(characterSplit);
                String venueCategory = originalData[columnvenueCategory];
                String venueid = originalData[columnVenueId];
                resultFile.println(poisMap.get(venueid) + "\t" + categoriesMap.get(venueCategory) + "\t" + 1.0);

            });
            stream.close();
            resultFile.close();

        } catch (FileNotFoundException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        } catch (IOException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }

    }

    /**
     * *
     * Method to generate a mapping files of categories of POIS
     *
     * @param originalPoisFile the original pois file
     * @param destinationFile the destination file
     */
    public static void generateMapOfVenuesCategories(String originalPoisFile, String destinationFile) {
        Map<String, Long> mapCategoriesVenues = new HashMap<>();
        try {
            PrintStream resultFile = new PrintStream(destinationFile);

            String characterSplit = "\t";

            int columnvenueCategory = 1;

            AtomicLong counterTypesOfCategories = new AtomicLong(1L);
            Stream<String> stream = Files.lines(Paths.get(originalPoisFile), StandardCharsets.ISO_8859_1);
            stream.forEach(line -> {
                String[] originalData = line.split(characterSplit);
                String venueCategory = originalData[columnvenueCategory];
                if (mapCategoriesVenues.get(venueCategory) == null) {
                    mapCategoriesVenues.put(venueCategory, counterTypesOfCategories.getAndIncrement());
                }

            });
            stream.close();
            for (String key : mapCategoriesVenues.keySet()) {
                resultFile.println(key + "\t" + mapCategoriesVenues.get(key));
            }

            resultFile.close();

        } catch (FileNotFoundException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        } catch (IOException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }

    }

    /**
     * *
     * Method to obtain a file of POIs coordinates
     *
     * @param originalPoisFile the file of original pois (information associated
     * of the POIS/VENUES)
     * @param fileMapItems file containing the map of olderId -> New id
     * @param destinationFile the destination file of the latitudes and
     * longitudes of the pois
     */
    public static void generatePoisCoords(String originalPoisFile, String fileMapItems, String destinationFile) {
        final Map<String, Long> mapItems = new HashMap<>();

        readMap(fileMapItems, mapItems);

        BufferedReader br;

        try {
            PrintStream resultFile = new PrintStream(destinationFile);

            br = new BufferedReader(new FileReader(originalPoisFile));

            String characterSplit = "\t";

            int columnVenueId = 0;
            int columnLatitude = 1;
            int columnLongitude = 2;

            String line = br.readLine();
            while (line != null) {
                if (line != null) {
                    String[] originalData = line.split(characterSplit);
                    // New id, latitude and longitude
                    resultFile.println(mapItems.get(originalData[columnVenueId]) + "\t" + originalData[columnLatitude]
                            + "\t" + originalData[columnLongitude]);
                }
                line = br.readLine();

            }
            br.close();
            resultFile.close();

        } catch (FileNotFoundException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        } catch (IOException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
    }
    
    /****
     * Method to print the coordinates of the POIs of a specific training file using the coordinates of a higher file of POIs
     * @param trainFile the trainFile
     * @param coordinates the coordinates
     * @param outputFile the outputFile
     */
    public static void specificPOICoords(String trainFile, Map<Long, Tuple2<Double, Double>> coordinates, String outputFile) {
    	 try {
             PrintStream resultFile = new PrintStream(outputFile);
             String characterSplit = "\t";
             int columnVenue = 1;
             Set<Long> procesedVenues = new HashSet<Long>();

            
             Stream<String> stream = Files.lines(Paths.get(trainFile));
             stream.forEach(line -> {
                 String[] originalData = line.split(characterSplit);
                 Long idVenue = Long.parseLong(originalData[columnVenue]);
                 
                 if (!procesedVenues.contains(idVenue)) {
                	 procesedVenues.add(idVenue);
                	 resultFile.println(idVenue + "\t" + coordinates.get(idVenue).v1 + "\t" + coordinates.get(idVenue).v2);
                 }
             });
             stream.close();
             resultFile.close();

         } catch (FileNotFoundException e) {
             // TODO Auto-generated catch block
             e.printStackTrace();
         } catch (IOException e) {
             // TODO Auto-generated catch block
             e.printStackTrace();
         }
	}
    
    
    
    /****
     * Method to filter a datamodel by the POIs of a city 
     * @param original the original datamodel
     * @param city the city to filter
     * @param cityMapping the mapping of items-cities
     * @return a new datamodel containing only the pois of a city
     */
    public static <U, I> SimpleFastTemporalPreferenceData<U, I> cityFilter(SimpleFastTemporalPreferenceData<U, I> original, String city, Map<I, String> cityMapping) {

        Stream<Tuple4<U, I, Double, Long>> prev = original.getAsTemporalTuples();
        Stream<Tuple4<U, I, Double, Long>> prevFiltered = prev.filter(t -> cityMapping.get(t.v2).equals(city));
        List<Tuple4<U, I, Double, Long>> lst = prevFiltered.collect(Collectors.toList());
        Set<U> usersFiltered = new HashSet<>();
        Set<I> itemsFiltered = new HashSet<>();
        lst.stream().forEach(t -> {
            usersFiltered.add(t.v1);
            itemsFiltered.add(t.v2);
        });
        FastUserIndex<U> userIndex = SimpleFastUserIndex.load(usersFiltered.stream());
        FastItemIndex<I> itemIndex = SimpleFastItemIndex.load(itemsFiltered.stream());

        SimpleFastTemporalPreferenceData<U, I> result = SimpleFastTemporalPreferenceData.loadTemporal(lst.stream(), userIndex, itemIndex);

        return result;
    }
    /***
     *  Method to filter a datamodel by the POIs blonging of different cities
     * @param original the original datamodel
     * @param cities the cities
     * @param cityMapping the mapping of item-cities
     * @return a new datamodel with POIs belongin to that cities
     */
    public static <U, I> SimpleFastTemporalPreferenceData<U, I> citiesFilter(SimpleFastTemporalPreferenceData<U, I> original, Set<String> cities, Map<I, String> cityMapping) {

        Stream<Tuple4<U, I, Double, Long>> prev = original.getAsTemporalTuples();
        Stream<Tuple4<U, I, Double, Long>> prevFiltered = prev.filter(t -> cities.contains(cityMapping.get(t.v2)));
        List<Tuple4<U, I, Double, Long>> lst = prevFiltered.collect(Collectors.toList());
        Set<U> usersFiltered = new HashSet<>();
        Set<I> itemsFiltered = new HashSet<>();
        lst.stream().forEach(t -> {
            usersFiltered.add(t.v1);
            itemsFiltered.add(t.v2);
        });
        FastUserIndex<U> userIndex = SimpleFastUserIndex.load(usersFiltered.stream());
        FastItemIndex<I> itemIndex = SimpleFastItemIndex.load(itemsFiltered.stream());

        SimpleFastTemporalPreferenceData<U, I> result = SimpleFastTemporalPreferenceData.loadTemporal(lst.stream(), userIndex, itemIndex);

        return result;
    }

    /**
     *
     * Creates different splittings depending on the constraints received as
     * parameters. Repetitions are not allowed except in the first splitting. It
     * will also control for repetitions between the splits.
     *
     * @param constraints array of initial and end timestamps for each split
     * that should be created
     * @param fullData the original check-in file as PreferenceData
     * @param filePoiIdCity file with mapping between POIs and their
     * corresponding cities. It can be null if splitting per city is not needed
     * @param kCore constraint to be applied for the kCore
     * @param acceptableCities a list of the cities to be considered when
     * creating the datamodels, or null (or empty) if all the check-ins should
     * be used
     *
     * @return the different splittings created as temporal datamodels
     */
    public static Map<String, SimpleFastTemporalPreferenceData<Long, Long>[]> createFoursquareSplitting(final Long[][] constraints, final SimpleFastTemporalPreferenceData<Long, Long> fullData, final int kCore, final String filePoiIdCity, final Set<String> acceptableCities, boolean filterTestByTrain) {
        // read poi (with new id) <-> city mapping
        Map<Long, String> mapping = new HashMap<>();
        try {
            if (filePoiIdCity != null) {
                String characterSplit = "\t";
                Stream<String> streamMapping = Files.lines(Paths.get(filePoiIdCity));
                streamMapping.forEach(line -> {
                    String[] poiCity = line.split(characterSplit);
                    if (poiCity.length == 3) {
                        mapping.put(Long.parseLong(poiCity[0]), poiCity[1] + "_" + poiCity[2].replaceAll("\\s+", ""));
                    } else {
                        mapping.put(Long.parseLong(poiCity[0]), poiCity[1].replaceAll("\\s+", ""));
                    }
                });
                streamMapping.close();
            }
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }
        Map<String, SimpleFastTemporalPreferenceData<Long, Long>[]> cityModels = new HashMap<>();
        // perform temporal filtering
        SimpleFastTemporalPreferenceData<Long, Long> timeFilteredData = SequentialRecommendersUtils.temporalFilter(fullData, constraints[0][0], constraints[constraints.length - 1][constraints[constraints.length - 1].length - 1]);
        // run k-core
        SimpleFastTemporalPreferenceData<Long, Long> kCoreData = SequentialRecommendersUtils.KCore(timeFilteredData, kCore, kCore);
        // split into training and test
        Boolean[] allowedRepetitions = new Boolean[constraints.length];
        Arrays.fill(allowedRepetitions, false);
        allowedRepetitions[0] = true;
        final SimpleFastTemporalPreferenceData<Long, Long>[] splits = SequentialRecommendersUtils.splitData(kCoreData, constraints, allowedRepetitions, filterTestByTrain);
        cityModels.put("complete", splits);
        // for each city, filter training and test
        final Set<String> cities = new HashSet<>();
        if (filePoiIdCity != null) {
            if ((acceptableCities == null) || (acceptableCities.isEmpty())) {
                cities.addAll(mapping.values());
            } else {
                cities.addAll(acceptableCities);
            }
        }
        final SimpleFastTemporalPreferenceData<Long, Long>[] cityCompleteSplits = new SimpleFastTemporalPreferenceData[constraints.length];
        IntStream.range(0, constraints.length).forEach(i -> cityCompleteSplits[i] = citiesFilter(splits[i], cities, mapping));
        cityModels.put("complete_cities", cityCompleteSplits);
        for (String city : cities) {
            final SimpleFastTemporalPreferenceData<Long, Long>[] citySplits = new SimpleFastTemporalPreferenceData[constraints.length];
            IntStream.range(0, constraints.length).forEach(i -> citySplits[i] = cityFilter(splits[i], city, mapping));
            cityModels.put(city, citySplits);
        }
        return cityModels;
    }

    public static Map<String, SimpleFastTemporalPreferenceData<Long, Long>[]> createCityPercentageSplitting(final SimpleFastTemporalPreferenceData<Long, Long> fullData, final int kCore, final String filePoiIdCity, final Set<String> acceptableCities) {
        // read poi (with new id) <-> city mapping
        Map<Long, String> mapping = new HashMap<>();
        try {
            if (filePoiIdCity != null) {
                String characterSplit = "\t";
                Stream<String> streamMapping = Files.lines(Paths.get(filePoiIdCity));
                streamMapping.forEach(line -> {
                    String[] poiCity = line.split(characterSplit);
                    if (poiCity.length == 3) {
                        mapping.put(Long.parseLong(poiCity[0]), poiCity[1] + "_" + poiCity[2].replaceAll("\\s+", ""));
                    } else {
                        mapping.put(Long.parseLong(poiCity[0]), poiCity[1].replaceAll("\\s+", ""));
                    }
                });
                streamMapping.close();
            }
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }
        Map<String, SimpleFastTemporalPreferenceData<Long, Long>[]> cityModels = new HashMap<>();
        // perform city filtering
        final Set<String> cities = new HashSet<>();
        if (filePoiIdCity != null) {
            if ((acceptableCities == null) || (acceptableCities.isEmpty())) {
                cities.addAll(mapping.values());
            } else {
                cities.addAll(acceptableCities);
            }
        }
        SimpleFastTemporalPreferenceData<Long, Long> citiesFilteredData = citiesFilter(fullData, cities, mapping);
        // run k-core
        SimpleFastTemporalPreferenceData<Long, Long> kCoreData = SequentialRecommendersUtils.KCore(citiesFilteredData, kCore, kCore);
        // split into training and test
        Boolean[] allowedRepetitions = new Boolean[2];
        Arrays.fill(allowedRepetitions, false);
        allowedRepetitions[0] = true;
        final SimpleFastTemporalPreferenceData<Long, Long>[] splits = SequentialRecommendersUtils.splitPercentageData(kCoreData, 0.8, allowedRepetitions);
        cityModels.put("complete_cities", splits);
        // for each city, filter training and test
        for (String city : cities) {
            final SimpleFastTemporalPreferenceData<Long, Long>[] citySplits = new SimpleFastTemporalPreferenceData[splits.length];
            IntStream.range(0, splits.length).forEach(i -> citySplits[i] = cityFilter(splits[i], city, mapping));
            cityModels.put(city, citySplits);
        }
        return cityModels;
    }
    
    /***
     * Method that will create a new dataset from the original, following these steps:
     * -Filtering out the ratings that does not belong to the acceptable cities
     * -Removing duplicates.
     * -KCore
     * -Temporal filter
     * @param fullData the full dataset
     * @param kCore minimum number of ratings to filter both users and items
     * @param filePoiIdCity the file of pois IDs and cities
     * @param acceptableCities the cities to filter
     * @return
     */
    public static Map<String, SimpleFastTemporalPreferenceData<Long, Long>[]> createCityPercentageSplittingRemoveDuplicatesAndSplit(final SimpleFastTemporalPreferenceData<Long, Long> fullData, final int kCore, final String filePoiIdCity, final Set<String> acceptableCities) {
        // read poi (with new id) <-> city mapping
        Map<Long, String> mapping = new HashMap<>();
        try {
            if (filePoiIdCity != null) {
                String characterSplit = "\t";
                Stream<String> streamMapping = Files.lines(Paths.get(filePoiIdCity));
                streamMapping.forEach(line -> {
                    String[] poiCity = line.split(characterSplit);
                    if (poiCity.length == 3) {
                        mapping.put(Long.parseLong(poiCity[0]), poiCity[1] + "_" + poiCity[2].replaceAll("\\s+", ""));
                    } else {
                        mapping.put(Long.parseLong(poiCity[0]), poiCity[1].replaceAll("\\s+", ""));
                    }
                });
                streamMapping.close();
            }
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }
        Map<String, SimpleFastTemporalPreferenceData<Long, Long>[]> cityModels = new HashMap<>();
        // perform city filtering
        final Set<String> cities = new HashSet<>();
        if (filePoiIdCity != null) {
            if ((acceptableCities == null) || (acceptableCities.isEmpty())) {
                cities.addAll(mapping.values());
            } else {
                cities.addAll(acceptableCities);
            }
        }
        SimpleFastTemporalPreferenceData<Long, Long> citiesFilteredData = citiesFilter(fullData, cities, mapping);
        Boolean[] allowedRepetitions = new Boolean[2];
        Arrays.fill(allowedRepetitions, false);
        
        final SimpleFastTemporalPreferenceData<Long, Long>[] splits = SequentialRecommendersUtils.splitPercentageRemoveDuplicatesAndSplit(citiesFilteredData, 0.8, allowedRepetitions,kCore, kCore);
        
        // split into training and test
        cityModels.put("complete_cities", splits);
        System.out.println(splits[0].numPreferences());
        System.out.println(splits[1].numPreferences());
        // for each city, filter training and test
        for (String city : cities) {
            final SimpleFastTemporalPreferenceData<Long, Long>[] citySplits = new SimpleFastTemporalPreferenceData[splits.length];
            IntStream.range(0, splits.length).forEach(i -> citySplits[i] = cityFilter(splits[i], city, mapping));
            cityModels.put(city, citySplits);
        }
        return cityModels;
    }

    
    public static void simulateSplitting(final Long[][] constraints, final SimpleFastTemporalPreferenceData<Long, Long> data, final int kCore, final String filePoiCity, final Set<String> acceptableCities, boolean filterTestByTrain) {
        Map<String, SimpleFastTemporalPreferenceData<Long, Long>[]> splits = createFoursquareSplitting(constraints, data, kCore, filePoiCity, acceptableCities, filterTestByTrain);
        for (int i = 0; i < constraints.length; i++) {
            for (String city : splits.keySet()) {
                SimpleFastTemporalPreferenceData<Long, Long> s = splits.get(city)[i];
                System.out.println(city + "\tSplit between " + constraints[i][0] + " and " + constraints[i][1] + ": " + s.numPreferences());
            }
        }
    }

    public static void generateSplitting(final String folder, final String[] suffix, final Long[][] constraints, final SimpleFastTemporalPreferenceData<Long, Long> data, final int kCore, final String filePoiCity, final Set<String> acceptableCities, boolean filterTestByTrain) {
        Map<String, SimpleFastTemporalPreferenceData<Long, Long>[]> splits = createFoursquareSplitting(constraints, data, kCore, filePoiCity, acceptableCities, filterTestByTrain);
        for (int i = 0; i < constraints.length; i++) {
            for (String city : splits.keySet()) {
                SimpleFastTemporalPreferenceData<Long, Long> s = splits.get(city)[i];
                String outfile = folder + "split_" + city + "__" + suffix[i];
                try {
                	File f = new File(outfile);
                	if (!f.exists()) {
                		System.out.println(outfile + "does not exist. Computing");
	                    PrintStream out = new PrintStream(outfile);
	                    s.getAsTemporalTuples().forEach(t -> {
	                        out.println(t.v1 + "\t" + t.v2 + "\t" + t.v3 + "\t" + t.v4);
	                    });
	                    out.close();
                	}
                } catch (FileNotFoundException e) {
                    e.printStackTrace();
                }
            }
        }
    }
    
    public static void generateSplittingTemporalFirstRemovingDuplicates(final String folder, final String[] suffix, final SimpleFastTemporalPreferenceData<Long, Long> data, final int kCore, final String filePoiCity, final Set<String> acceptableCities) {
        Map<String, SimpleFastTemporalPreferenceData<Long, Long>[]> splits = createCityPercentageSplittingRemoveDuplicatesAndSplit(data, kCore, filePoiCity, acceptableCities);
        for (int i = 0; i < suffix.length; i++) {
            for (String city : splits.keySet()) {
                SimpleFastTemporalPreferenceData<Long, Long> s = splits.get(city)[i];
                String outfile = folder + "split_" + city + "__" + suffix[i];
                try {
                    PrintStream out = new PrintStream(outfile);
                    s.getAsTemporalTuples().forEach(t -> {
                        out.println(t.v1 + "\t" + t.v2 + "\t" + t.v3 + "\t" + t.v4);
                    });
                    out.close();
                } catch (FileNotFoundException e) {
                    e.printStackTrace();
                }
            }
        }
    }
    
    
/****
 * Method to check if a recommendation file is correct depending on some constraints:
 * -For every user, there must not be any any item recommended that has been previously seen in train.
 * -Every item recommended must match the matching category.
 * 
 * @param trainFile the train file 
 * @param poiCityFile the poiCityFile
 * @param recommendedFile the recommended file
 * @param matchCategory the matching category
 */
    public static void checkRecommendationCorrection(String trainFile, String poiCityFile, String recommendedFile, String matchCategory, String testFile) {
        AtomicBoolean result = new AtomicBoolean(true);
        AtomicBoolean hasCityFile = new AtomicBoolean(false);
        
        //
        AtomicBoolean userRatedItemInTrain = new AtomicBoolean(false);
        AtomicBoolean itemNotRatedInTrain = new AtomicBoolean(false);
        AtomicBoolean itemNotInCity = new AtomicBoolean(false);
        AtomicBoolean userNotInTest = new AtomicBoolean(false);
        if (poiCityFile != null) {
        	hasCityFile.set(true);
        }
        try {
        	Tuple2<List<Long>, List<Long>> setIndexUsersItems = ExperimentUtils.getUserItemsFromFile(trainFile,lp,lp);
        	
        	FastUserIndex<Long> userIndex = SimpleFastUserIndex.load(setIndexUsersItems.v1.stream());
            FastItemIndex<Long> itemIndex = SimpleFastItemIndex.load(setIndexUsersItems.v2.stream());

            String characterSplit = "\t";
            FastPreferenceData<Long, Long> ranksysTrainDataOriginal = SimpleFastPreferenceData
                    .load(SimpleRatingPreferencesReader.get().read(trainFile, lp, lp), userIndex, itemIndex);
            
            Tuple2<List<Long>, List<Long>> setIndexUsersItemsTest = ExperimentUtils.getUserItemsFromFile(testFile,lp,lp);
        	
        	FastUserIndex<Long> userIndexTest = SimpleFastUserIndex.load(setIndexUsersItemsTest.v1.stream());
            FastItemIndex<Long> itemIndexTest = SimpleFastItemIndex.load(setIndexUsersItemsTest.v2.stream());

            FastPreferenceData<Long, Long> ranksysTestDataOriginal = SimpleFastPreferenceData
                    .load(SimpleRatingPreferencesReader.get().read(testFile, lp, lp), userIndexTest, itemIndexTest);

            
            Set<Long> poisCities = new HashSet<>();
            if (hasCityFile.get()) {
	            Stream<String> stream = Files.lines(Paths.get(poiCityFile));
	            stream.forEach(line -> {
	                String[] split = line.split(characterSplit);
	                if (split[1].contains(matchCategory)) {
	                    poisCities.add(Long.parseLong(split[0]));
	                }
	            });
	            stream.close();
            }
            
            
            //We got now the pois associated with the city. We now check that the recommended items are not in train and are from that city

            Stream<String> stream2 = Files.lines(Paths.get(recommendedFile));
            stream2.forEach(line -> {
                String[] recommended = line.split(characterSplit);
                Long user = Long.parseLong(recommended[0]);
                Long item = Long.parseLong(recommended[1]);
                
                if (!ranksysTestDataOriginal.containsUser(user)) {
                	//System.out.println("User "+ user +" is not in the test subset");
                	result.set(false);
                	userNotInTest.set(true);
                }

                if (ranksysTrainDataOriginal.containsUser(user) && ranksysTrainDataOriginal.getUserPreferences(user).filter(pref -> pref.v1.equals(item)).findAny().isPresent()) {
                    System.out.println("Item " + item + " has been rated in train by user " + user);
                    result.set(false);
                    userRatedItemInTrain.set(true);
                }
                
                if (hasCityFile.get() && !poisCities.contains(item)) {
                    //System.out.println("Item " + item + " is not in city");
                    result.set(false);
                    itemNotInCity.set(true);
                }
                
                if (ranksysTrainDataOriginal.containsItem(item)) {
	                if (!ranksysTrainDataOriginal.getItemPreferences(item).findFirst().isPresent()) {
	                	//System.out.println("Item " + item + "is not rated in train");
	                	itemNotRatedInTrain.set(true);
	                	result.set(false);
	                }
                }else {
                	//System.out.println("Item "+ item +" is not in the train subset");
                	itemNotRatedInTrain.set(true);
                	result.set(false);
                }

            });
            stream2.close();
            if (result.get()) {
                System.out.println("Recommendations correct");
                System.out.println("No item rated in train by same user. No items not in city " + matchCategory + ". All items are rated in train. No user in recommendation file that it is not in test.");
            } else {
            	if (itemNotRatedInTrain.get()) {
            		System.out.println("Some items are not rated in the train subset");
            	}
            	
            	if (itemNotInCity.get()) {
            		System.out.println("Some items are not in the destination city");
            	}
            	
            	if (userRatedItemInTrain.get()) {
            		System.out.println("Some users rated the item in the training set");
            	}
            	
            	if (userNotInTest.get()) {
            		System.out.println("Some users do not appear in the test set");
            	}
            	
            }

        } catch (IOException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }

    }

    /**
     * Creates a datamodel from the original check-in file.
     *
     * @param originalDataset the original check-in file
     * @param processData flag to indicate if the date is already processed (it
     * is a timestamp) or not, the same for users or items
     * @param userMappingFile file containing the user mapping
     * @param itemMappingFile file containing the item mapping
     *
     * @return the check-in file in memory as a temporal datamodel
     */
    public static SimpleFastTemporalPreferenceData<Long, Long> readFullDataset(final String originalDataset, final boolean processData, final String userMappingFile, final String itemMappingFile) {
        String characterSplit = "\t";
        try {
            Set<Long> users = new HashSet<>();
            AtomicLong userCounter = new AtomicLong(1L);
            AtomicLong itemCounter = new AtomicLong(1L);
            Map<String, Long> userMapping = new HashMap<>();
            if (processData && new File(userMappingFile).exists()) {
                readMap(userMappingFile, userMapping);
            }
            Map<String, Long> itemMapping = new HashMap<>();
            if (processData && new File(itemMappingFile).exists()) {
                readMap(itemMappingFile, itemMapping);
            }
            Set<Long> items = new HashSet<>();
            // check-in
            Map<Long, Map<Long, List<Tuple2<Double, Long>>>> userItemInteractions = new HashMap<>();
            Stream<String> stream = Files.lines(Paths.get(originalDataset));
            stream.forEach(line -> {
                String[] originalData = line.split(characterSplit);
                String user = originalData[0];
                Long u = userMapping.get(user);
                if (!processData)
                	u = Long.parseLong(user);
                else if (u == null) {
	                    u = userCounter.getAndIncrement();
	                    userMapping.put(user, u);
	                }
                
                String poi = originalData[1];
                Long item = itemMapping.get(poi);
                if (!processData)
                	item = Long.parseLong(poi);
                else if (item == null) {
                    item = itemCounter.getAndIncrement();
                    itemMapping.put(poi, item);
                }
                users.add(u);
                items.add(item);
                String date = (processData ? originalData[2] : originalData[3]);
                Double score = (processData ? 1.0 : Double.parseDouble(originalData[2]));
                try {
                    Long processedDate = (processData ? parseFoursqrDateLocalTime(date, originalData[3]) : Long.parseLong(date));
                    if (!userItemInteractions.containsKey(u)) {
                        userItemInteractions.put(u, new HashMap<>());
                    }
                    if (!userItemInteractions.get(u).containsKey(item)) {
                        userItemInteractions.get(u).put(item, new ArrayList<>());
                    }
                    userItemInteractions.get(u).get(item).add(new Tuple2<>(score, processedDate));
                } catch (IllegalArgumentException e) {
                    System.out.println("Ignoring " + line);
                }
            });
            stream.close();
            // create the datamodels
            FastUserIndex<Long> uIndex = SimpleFastUserIndex.load(users.stream());
            FastItemIndex<Long> iIndex = SimpleFastItemIndex.load(items.stream());
            List<Tuple4<Long, Long, Double, Long>> tuples = new ArrayList<>();
            for (Long u : userItemInteractions.keySet()) {
                for (Long item : userItemInteractions.get(u).keySet()) {
                    for (Tuple2<Double, Long> t : userItemInteractions.get(u).get(item)) {
                        tuples.add(new Tuple4<>(u, item, t.v1, t.v2));
                    }
                }
            }
            SimpleFastTemporalPreferenceData<Long, Long> model = SimpleFastTemporalPreferenceData.loadTemporal(tuples.stream(), uIndex, iIndex);
            return model;
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }
        return null;
    }

    /**
     * *
     * Method to obtain the minimum and maximum timestamps for each city based
     * on the check-ins
     *
     * @param originalDataset the original check-in file
     * @param filePoiCity file with mapping between POIs and their corresponding
     * cities
     * @param processDate flag to indicate if the date is already processed (it
     * is a timestamp) or not
     */
    public static void computeTimeIntervalsPerCity(final String originalDataset, final String filePoiCity, final boolean processDate) {
        String characterSplit = "\t";
        try {
            Map<String, Long> minDates = new HashMap<>();
            Map<String, Long> maxDates = new HashMap<>();
            // read poi <-> city mapping
            Map<String, String> mapping = new HashMap<>();
            Stream<String> streamMapping = Files.lines(Paths.get(filePoiCity));
            streamMapping.forEach(line -> {
                String[] poiCity = line.split(characterSplit);
                mapping.put(poiCity[0], poiCity[1] + "_" + poiCity[2]);
            });
            streamMapping.close();
            // check-in
            Stream<String> stream = Files.lines(Paths.get(originalDataset));
            stream.forEach(line -> {
                String[] originalData = line.split(characterSplit);
                String city = mapping.get(originalData[1]);
                String date = originalData[2];
                try {
                    Long processedDate = (processDate ? parseFoursqrDateLocalTime(date,originalData[3]) : Long.parseLong(date));
                    minDates.put(city, Long.min(processedDate, minDates.getOrDefault(city, Long.MAX_VALUE)));
                    maxDates.put(city, Long.max(processedDate, maxDates.getOrDefault(city, Long.MIN_VALUE)));
                } catch (IllegalArgumentException e) {
                    System.out.println("Ignoring " + line);
                }
            });
            stream.close();

            minDates.keySet().forEach(k -> System.out.println(k + "\t" + minDates.get(k) + "\t" + maxDates.get(k)));
        } catch (FileNotFoundException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        } catch (IOException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
    }

    /****
     * Method to parse a file containing the old ids of POIS and cities to a file containing the new ids of pois and cities
     * @param originalFile the original pois-cities file
     * @param fileMapItems the map of items
     * @param destFile the destination file
     */
    public static void parsePOISCities(String originalFile, String fileMapItems, String destFile) {
        final Map<String, Long> mapItems = new HashMap<>();
        readMap(fileMapItems, mapItems);
        String characterSplit = "\t";

        try {
            PrintStream resultFile = new PrintStream(destFile);
            Stream<String> stream = Files.lines(Paths.get(originalFile));
            stream.forEach(line -> {
                String[] originalData = line.split(characterSplit);
                try {
                	if (mapItems.get(originalData[0]) == null) {
                		System.out.println(originalData[0] + " WRONG");
                	}
                    resultFile.println(mapItems.get(originalData[0]) + "\t" + originalData[1]);
                } catch (IllegalArgumentException e) {
                    System.out.println("Ignoring " + line);
                }
            });
            stream.close();
            resultFile.close();

        } catch (FileNotFoundException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        } catch (IOException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
    }

    /**
     * *
     * Method to generate a new file of check-ins from the original dataset but
     * with new ids of users/items and the original time switched to timestamp
     *
     * @param originalDataset the original check-in file
     * @param fileMapUsers map of users (old id of user -> new id of user)
     * @param fileMapItems map of items (old id of item -> new id of item)
     * @param newDataset the dataset with new ids of users/items and timestamps
     */
    public static void generateNewCheckinFileWithTimeStamps(String originalDataset, String fileMapUsers,
            String fileMapItems, String newDataset, boolean printUTC) {
        // Old id -> new id
        final Map<String, Long> mapUsers = new HashMap<>();
        readMap(fileMapUsers, mapUsers);

        // Old id -> new Id
        final Map<String, Long> mapItems = new HashMap<>();
        readMap(fileMapItems, mapItems);

        String characterSplit = "\t";

        try {
            PrintStream resultFile = new PrintStream(newDataset);
            Stream<String> stream = Files.lines(Paths.get(originalDataset));
            stream.forEach(line -> {
                String[] originalData = line.split(characterSplit);
                try {
                    // to work with ratings we will need to concatenate an 1.0 as a rating
                    String utcDate = originalData[2];
                    String offsetStr = originalData[3];
                    if (printUTC) {
                    	resultFile.println(mapUsers.get(originalData[0]) + "\t" + mapItems.get(originalData[1]) + "\t" + 1.0
                            + "\t" + parseFoursqrDateLocalTime(utcDate, offsetStr)  + "\t" + parseFoursqrDateUTC(utcDate));
                    } else {
                    	resultFile.println(mapUsers.get(originalData[0]) + "\t" + mapItems.get(originalData[1]) + "\t" + 1.0
                                + "\t" + parseFoursqrDateLocalTime(utcDate, offsetStr));
                    }
                    
                } catch (IllegalArgumentException e) {
                    System.out.println("Ignoring " + line);
                }
            });
            stream.close();
            resultFile.close();

        } catch (FileNotFoundException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        } catch (IOException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }

    }

    /**
     * Method to read map of older ids and new ids
     *
     * @param fileMap
     * @return the map read
     */
    public static void readMap(String fileMap, Map<String, Long> result) {
        String characterSplit = "\t";
        try {
            Stream<String> stream = Files.lines(Paths.get(fileMap));
            stream.forEach(line -> {
                String[] oldNewId = line.split(characterSplit);
                result.put(oldNewId[0], Long.parseLong(oldNewId[1]));

            });
            stream.close();
        } catch (FileNotFoundException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        } catch (IOException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
    }
    
    /***
     * Method to build an inverse map from a map. It is necessary that the values and the keys are
     * unique
     * @param original the original map
     * @param inverse the inverse map
     */
    public static<K, V> void inverseMap(Map<K, V> original, Map<V, K> inverse) {
    	
    	for (Map.Entry<K, V> entryOriginal : original.entrySet()) {
    		inverse.put(entryOriginal.getValue(), entryOriginal.getKey());
    	}
    	
    }

    /**
     * Method to obtain files of maps of items, and users with the old and the
     * new ids
     *
     * @param originalDataset the original dataset (check-ins file)
     * @param destMapUsers the destination file to store the map of users
     * @param destMapItems the destination file to store the map of items
     */
    public static void generateMapOfUsersVenues(String originalDataset, String destMapUsers, String destMapItems) {
        // Old id -> new id
        Map<String, Long> mapUsers = new HashMap<>();

        // Old id -> new Id
        Map<String, Long> mapItems = new HashMap<>();

        AtomicLong counterUsers = new AtomicLong(1L);
        AtomicLong counterItems = new AtomicLong(1L);

        try {
            PrintStream usersStream = new PrintStream(destMapUsers);
            PrintStream itemsStream = new PrintStream(destMapItems);

            String characterSplit = ",";

            // Parameters to configure
            int columnVenueId = 1;
            int columnUserId = 0;

            // First pass of the file
            System.out.println("Reading CheckingsFile file");
            Stream<String> stream = Files.lines(Paths.get(originalDataset));
            stream.forEach(line -> {
                String[] data = line.split(characterSplit);
                String venueID = data[columnVenueId];
                String userId = data[columnUserId];

                if (mapUsers.get(userId) == null) {
                    mapUsers.put(userId, counterUsers.getAndIncrement());
                }

                if (mapItems.get(venueID) == null) {
                    mapItems.put(venueID, counterItems.getAndIncrement());
                }
            });
            stream.close();

            for (String oldUserId : mapUsers.keySet()) {
                usersStream.println(oldUserId + "\t" + mapUsers.get(oldUserId));
            }

            for (String oldItemId : mapItems.keySet()) {
                itemsStream.println(oldItemId + "\t" + mapItems.get(oldItemId));
            }

            itemsStream.close();
            usersStream.close();

        } catch (FileNotFoundException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        } catch (IOException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }

    }
    
    /**
     * Method to obtain files of maps of items, and users with the old and the
     * new ids
     *
     * @param originalDataset the original dataset (check-ins file)
     * @param destMapUsers the destination file to store the map of users
     * @param destMapItems the destination file to store the map of items
     */
    public static void generateNewDatasetMapUsersItems(String originalDataset, String destMapUsers, String destMapItems, String outputFile) {
        // Old id -> new id
        Map<String, Long> mapUsers = new HashMap<>();

        // Old id -> new Id
        Map<String, Long> mapItems = new HashMap<>();

        AtomicLong counterUsers = new AtomicLong(1L);
        AtomicLong counterItems = new AtomicLong(1L);
        

        try {
            PrintStream usersStream = new PrintStream(destMapUsers);
            PrintStream itemsStream = new PrintStream(destMapItems);
            PrintStream out = new PrintStream(outputFile);

            
            String characterSplit = "\t";

            // Parameters to configure
            int columnUserId = 0;
            int columnItemId = 1;
            int columnRating = 2;
            int columnTime = 3;

            // First pass of the file
            System.out.println("Reading CheckingsFile file");
            Stream<String> stream = Files.lines(Paths.get(originalDataset));
            stream.forEach(line -> {
                String[] data = line.split(characterSplit);
                String itemID = data[columnItemId];
                String userId = data[columnUserId];

                if (mapUsers.get(userId) == null) {
                    mapUsers.put(userId, counterUsers.getAndIncrement());
                }

                if (mapItems.get(itemID) == null) {
                    mapItems.put(itemID, counterItems.getAndIncrement());
                }
                
                out.println(mapUsers.get(userId) + "\t" + mapItems.get(itemID) + "\t" + data[columnRating] + "\t" + data[columnTime]);
            });
            stream.close();

            for (String oldUserId : mapUsers.keySet()) {
                usersStream.println(oldUserId + "\t" + mapUsers.get(oldUserId));
            }

            for (String oldItemId : mapItems.keySet()) {
                itemsStream.println(oldItemId + "\t" + mapItems.get(oldItemId));
            }

            itemsStream.close();
            usersStream.close();
            out.close();

        } catch (FileNotFoundException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        } catch (IOException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }

    }

    /**
     * *
     * Method to generate check-ins files per cites
     *
     * @param venuesFile the file of Pois/venues
     * @param citiesFile the cities file
     * @param checkinsData the original check-ins file (Original)
     * @param destPath the destinationPaths
     */
    public static void generateFilesPerCity(String venuesFile, String citiesFile, String checkinsData, String destPath,
            String pathOfVenuesAndCity) {
        // For venues
        Map<String, String> venuesCountryCode;
        Map<String, FoursqrVenue> venuesIDInformation;
        Tuple2<Map<String, String>, Map<String, FoursqrVenue>> tuplesVenues;

        // For cities
        Map<String, String> countryCodeAndCityToCountryCode = new HashMap<>();
        Map<String, FoursqrCity> countryCodeAndCitiesComplete = new HashMap<>();
        Map<String, Set<String>> countryCodeToCountryCodeAndCities = new HashMap<>();
        Tuple3<Map<String, String>, Map<String, Set<String>>, Map<String, FoursqrCity>> tupleCities;

        // Reading the files
        tuplesVenues = readFoursqrVenues(venuesFile);

        venuesCountryCode = tuplesVenues.v1;
        venuesIDInformation = tuplesVenues.v2;

        tupleCities = readFoursqrCities(citiesFile);
        countryCodeAndCityToCountryCode = tupleCities.v1;
        countryCodeToCountryCodeAndCities = tupleCities.v2;
        countryCodeAndCitiesComplete = tupleCities.v3;

        // Venues to cities
        Map<String, String> venuesToCity = new HashMap<String, String>();

        // All the files associated with the countries
        Map<String, PrintStream> mapPrint = new HashMap<>();
        try {
            if (destPath != null) {
                for (String city : countryCodeAndCityToCountryCode.keySet()) {
                    mapPrint.put(city, new PrintStream(destPath + city + ".txt"));
                }
            }

            // All print streams opened. Now read the check-in file and process it
            String characterSplit = "\t";

            // Parameters to configure
            int columnVenueId = 1;

            // First pass of the file
            System.out.println("Reading CheckingsFile file");

            Map<String, Set<String>> countryCodeToCountryCodeAndCities2 = countryCodeToCountryCodeAndCities;
            Map<String, FoursqrCity> countryCodeAndCitiesComplete2 = countryCodeAndCitiesComplete;
            Stream<String> stream = Files.lines(Paths.get(checkinsData));
            stream.forEach(line -> {
                String[] data = line.split(characterSplit);
                String venueID = data[columnVenueId];
                // Get the venue ID , get the country associated, write in that print stream the
                // checkings
                String countryCode = venuesCountryCode.get(venueID);
                Set<String> citiesOfCountry = countryCodeToCountryCodeAndCities2.get(countryCode);
                if (venuesToCity.get(venueID) == null) {
                    // We have not found the city of this venue
                    String contryCodeAndcity = getCityFromVenue(venuesIDInformation.get(venueID), citiesOfCountry,
                            countryCodeAndCitiesComplete2);
                    venuesToCity.put(venueID, contryCodeAndcity);
                }

                // We have now the city
                if (mapPrint.containsKey(venuesToCity.get(venueID))) {
                    mapPrint.get(venuesToCity.get(venueID)).println(line);
                }

            });
            stream.close();

            // Write venues to cities file
            PrintStream venuetoCity = new PrintStream(pathOfVenuesAndCity);
            for (String venueID : venuesToCity.keySet()) {
                venuetoCity.println(venueID + "\t" + venuesToCity.get(venueID));
            }
            venuetoCity.close();

        } catch (FileNotFoundException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        } catch (IOException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        } finally {
            for (String city : mapPrint.keySet()) {
                mapPrint.get(city).close();
            }
        }
    }

    /**
     * *
     * Method that receives a list of cities and a venue and returns the city
     * associated to that venue
     *
     * @param v the venue
     * @param cities the list of cities
     * @param citiesComplete the map of complete information of the city
     * @return the city associated to the venue
     */
    private static String getCityFromVenue(FoursqrVenue v, Set<String> cities,
            Map<String, FoursqrCity> citiesComplete) {
        double minDistance = Double.MAX_VALUE;
        String finalCity = null;
        for (String city : cities) {
            double heaverD = heaversine(v, citiesComplete.get(city));
            if (minDistance > heaverD) {
                minDistance = heaverD;
                finalCity = city;
            }
        }

        return finalCity;
    }


    
    /***
     * Method to parse a Foursqr date without offsets
     * @param utcDate
     * @param offsetStr
     * @return
     */
    public static Long parseFoursqrDateUTC(String utcDate) {
		SimpleDateFormat formatter = new SimpleDateFormat("EEE MMM dd HH:mm:ss Z yyyy",  Locale.ENGLISH);		
		Date date = null;
        try {
            date = formatter.parse(utcDate);
        } catch (ParseException e) {
            System.out.println("Date "+ date + " not supported");
            throw new IllegalArgumentException("Incorrect date " + date);
        }

        return date.getTime() / 1000;
    }
    
    /**
     * Method to parse the date considering the local time. We need to add the minutes to the date. 
     * @param utcDate
     * @param offseStr
     * @return
     */
    public static Long parseFoursqrDateLocalTime(String utcDate, String offseStr) {
		SimpleDateFormat formatter = new SimpleDateFormat("EEE MMM dd HH:mm:ss Z yyyy",  Locale.ENGLISH);		
		Date date = null;
        try {
            date = formatter.parse(utcDate);
            Calendar cal = Calendar.getInstance();
            cal.setTime(date);
            cal.add(Calendar.MINUTE, Integer.parseInt(offseStr));
            date = cal.getTime();
        } catch (ParseException e) {
            System.out.println("Date "+ date + " not supported");
            throw new IllegalArgumentException("Incorrect date " + date);
        }

        return date.getTime() / 1000;
    }
    

    public static void obtainCitiesMaximizingCategory(String POICity, String POICategory, String matchingCategory) {
    	Map<String, String> poisMatchingCat = new HashMap<>();
    	Map<String, Integer> cityPoisMatchingCat = new HashMap<>();
    	
    	try {
    		Stream<String> stream = Files.lines(Paths.get(POICategory));
    		stream.forEach(line -> {
    			String [] data = line.split("\t");
    			
    			if (data[3].toLowerCase().contains(matchingCategory.toLowerCase())) {
    				poisMatchingCat.put(data[0], data[3]);
    			}    			
    		});
    		stream.close();
    		
    		Stream<String> stream2 = Files.lines(Paths.get(POICity));
    		stream2.forEach(line -> {
    			String [] data = line.split("\t");
    			String POI = data[0];
    			String city = data[1];
    			
    			//The poi is stored previously
    			if (poisMatchingCat.get(POI) != null) {
    				if (cityPoisMatchingCat.get(city) == null) {
    					cityPoisMatchingCat.put(city, 0);
    				}
					cityPoisMatchingCat.put(city, cityPoisMatchingCat.get(city) + 1);
    			}
    		});
    		stream2.close();
    		
    		Map<String, Integer> result = SequentialRecommendersUtils.sortByValue(cityPoisMatchingCat, true);
    		
    		for (String k: result.keySet()) {
    			System.out.println(k + "\t" + result.get(k));
    		}
    		
    		
    		
    	} catch (Exception e) {
    		// TODO Auto-generated catch block
            e.printStackTrace();	
    	}
    	
    }
    
    

    /**
     * *
     * Method to generate files of checkings per country
     *
     * @param venuesFile
     * @param citiesFile
     * @param checkinsData
     * @param destPath
     */
    public static void generateFilesPerCountry(String venuesFile, String citiesFile, String checkinsData,
            String destPath) {
        // For venues
        Map<String, String> venuesCountryCode;
        Map<String, FoursqrVenue> venuesIDInformation;
        Tuple2<Map<String, String>, Map<String, FoursqrVenue>> tuplesVenues;

        // For cities
        Map<String, String> countryCodeAndCityToCountryCode = new HashMap<>();
        Map<String, FoursqrCity> countryCodeCitiesComplete = new HashMap<>();
        Map<String, Set<String>> countryCodeToCountryCodeCities = new HashMap<>();
        Tuple3<Map<String, String>, Map<String, Set<String>>, Map<String, FoursqrCity>> tupleCities;

        // Reading the files
        tuplesVenues = readFoursqrVenues(venuesFile);

        venuesCountryCode = tuplesVenues.v1;
        venuesIDInformation = tuplesVenues.v2;

        tupleCities = readFoursqrCities(citiesFile);
        countryCodeAndCityToCountryCode = tupleCities.v1;
        countryCodeToCountryCodeCities = tupleCities.v2;
        countryCodeCitiesComplete = tupleCities.v3;

        // All the files associated with the countries
        Map<String, PrintStream> mapPrint = new HashMap<>();
        try {
            for (String countryCode : countryCodeToCountryCodeCities.keySet()) {
                List<String> cities = new ArrayList<>(countryCodeToCountryCodeCities.get(countryCode)); // we are only
                // interested in
                // the full name
                // of the
                // country
                mapPrint.put(countryCode,
                        new PrintStream(destPath + countryCodeCitiesComplete.get(cities.get(0)).country + ".txt"));
            }

            // All print streams opened. Now read the checkin file and process it
            String characterSplit = "\t";

            // Parameters to configure
            int columnVenueId = 1;

            // First pass of the file
            System.out.println("Reading CheckingsFile file");
            Stream<String> stream = Files.lines(Paths.get(checkinsData));
            stream.forEach(line -> {
                String[] data = line.split(characterSplit);
                String venueID = data[columnVenueId];
                // Get the venue ID , get the country associated, write in that printstream the
                // checkins
                String countryCode = venuesCountryCode.get(venueID);
                mapPrint.get(countryCode).println(line);

            });
            stream.close();

        } catch (FileNotFoundException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        } catch (IOException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        } finally {
            for (String countryCode : countryCodeToCountryCodeCities.keySet()) {
                mapPrint.get(countryCode).close();
            }
        }

    }

    /**
     * Read foursqr venues
     *
     * @param venuesFile
     * @return
     */
    private static Tuple2<Map<String, String>, Map<String, FoursqrVenue>> readFoursqrVenues(String venuesFile) {
        // Assume structure of the file is: UserdId(String) ItemID(String)
        // rating(double) timestamp(long)

        // Map of venuesID -> country
        Map<String, String> venuesCountryCode = new HashMap<>();

        // Map of venuesID, and information
        Map<String, FoursqrVenue> venuesIDInformation = new HashMap<>();

        String characterSplit = "\t";

        // Parameters to configure
        int columnVenueId = 0;
        int columnLatitude = 1;
        int columnLongitude = 2;
        int columnCountryCode = 4;

        try {

            // First pass of the file
            System.out.println("Reading POIS/Venues file");
            Stream<String> stream = Files.lines(Paths.get(venuesFile));
            stream.forEach(line -> {
                String[] data = line.split(characterSplit);
                String venueID = data[columnVenueId];
                String countryCode = data[columnCountryCode];
                Double longitude = Double.parseDouble(data[columnLongitude]);
                Double latitude = Double.parseDouble(data[columnLatitude]);

                if (venuesCountryCode.get(venueID) == null) {
                    venuesCountryCode.put(venueID, countryCode);
                }

                if (venuesIDInformation.get(venueID) == null) {
                    venuesIDInformation.put(venueID, new FoursqrVenue(countryCode, latitude, longitude));
                }

            });
            stream.close();

            return new Tuple2<>(venuesCountryCode, venuesIDInformation);

        } catch (FileNotFoundException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        } catch (IOException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
        return null;

    }

    /**
     * Read foursqr cities file
     *
     * @param citiesFile
     * @return
     */
    private static Tuple3<Map<String, String>, Map<String, Set<String>>, Map<String, FoursqrCity>> readFoursqrCities(
            String citiesFile) {
        // Assume structure of the file is: UserdId(String) ItemID(String)
        // rating(double) timestamp(long)
        // For each city there is only one country
        Map<String, String> countryCodeAndCityToCountryCode = new HashMap<>();

        // Each city with its latitude and longitude
        Map<String, FoursqrCity> countryCodeAndCitiesComplete = new HashMap<>();

        // Each country may have multiple cities
        Map<String, Set<String>> countryCodeToCountryCodeAndCities = new HashMap<>();

        String characterSplit = "\t";

        // Parameters to configure
        int columnCityName = 0;
        int columnLatitude = 1;
        int columnLongitude = 2;
        int columnCountryCode = 3;
        int columnCountry = 4;
        int columnCityType = 5;

        try {

            // First pass of the file
            System.out.println("Reading cities file");
            Stream<String> stream = Files.lines(Paths.get(citiesFile));
            stream.forEach(line -> {
                String[] data = line.split(characterSplit);
                String cityName = data[columnCityName];
                String countryCode = data[columnCountryCode];
                cityName = cityName.replaceAll("\\s+", ""); //Remove spaces

                countryCodeAndCityToCountryCode.put(countryCode + "_" + cityName, countryCode);

                countryCodeAndCitiesComplete.put(countryCode + "_" + cityName,
                        new FoursqrCity(data[columnCityType], data[columnCountry],
                                Double.parseDouble(data[columnLatitude]), Double.parseDouble(data[columnLongitude])));

                if (countryCodeToCountryCodeAndCities.get(countryCode) == null) {
                    countryCodeToCountryCodeAndCities.put(countryCode, new TreeSet<>());
                }

                countryCodeToCountryCodeAndCities.get(countryCode).add(countryCode + "_" + cityName);
            });
            stream.close();

            return new Tuple3<>(countryCodeAndCityToCountryCode, countryCodeToCountryCodeAndCities,
                    countryCodeAndCitiesComplete);

        } catch (FileNotFoundException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        } catch (IOException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
        return null;

    }

    /**
     * Specific class to work with foursqr venues or POIS
     *
     * @author Pablo Sanchez
     *
     */
    public static class FoursqrVenue {

        String countryCode;
        double latitude;
        double longitude;

        public FoursqrVenue(String countryCode, double latitude, double longitude) {
            this.countryCode = countryCode;
            this.latitude = latitude;
            this.longitude = longitude;
        }

        public String getCountryCode() {
            return countryCode;
        }

        public void setCountryCode(String cityType) {
            this.countryCode = cityType;
        }

        public double getLatitude() {
            return latitude;
        }

        public void setLatitude(double latitude) {
            this.latitude = latitude;
        }

        public double getLongitude() {
            return longitude;
        }

        public void setLongitude(double longitude) {
            this.longitude = longitude;
        }

    }

    /**
     * Specific class for working with Foursqr data
     *
     * @author Pablo Sanchez
     *
     */
    public static class FoursqrCity {

        String cityType;
        String country;
        double latitude;
        double longitude;

        public FoursqrCity(String cityType, String country, double latitude, double longitude) {
            this.cityType = cityType;
            this.latitude = latitude;
            this.longitude = longitude;
            this.country = country;
        }

        public String getCityType() {
            return cityType;
        }

        public void setCityType(String cityType) {
            this.cityType = cityType;
        }

        public double getLatitude() {
            return latitude;
        }

        public void setLatitude(double latitude) {
            this.latitude = latitude;
        }

        public double getLongitude() {
            return longitude;
        }

        public void setLongitude(double longitude) {
            this.longitude = longitude;
        }

    }

    /**
     * *
     * Method to compute the heaversine distance
     *
     * @param v the foursqr venue
     * @param fc the foursqr city
     * @return the distance (in km)
     */
    public static double heaversine(FoursqrVenue v, FoursqrCity fc) {
        double lat1 = v.latitude;
        double long1 = v.longitude;
        double lat2 = fc.latitude;
        double long2 = fc.longitude;

        return haversine(lat1, long1, lat2, long2, false);

    }

    /***
     * Method to compute the haversine distance between two coordinates
     * @param lat1 the first latitude
     * @param long1 the first longitude
     * @param lat2 the second latitude
     * @param long2 the second longitude
     * @return the distance between the two coordinates in km
     */
    public static double haversine(double lat1, double long1, double lat2, double long2, boolean meters) {
        double dlat = Math.toRadians(lat2 - lat1);
        double dLong = Math.toRadians(long2 - long1);

        double startLat = Math.toRadians(lat1);
        double endLat = Math.toRadians(lat2);

        double a = haversin(dlat) + Math.cos(startLat) * Math.cos(endLat) * haversin(dLong);
        double c = 2 * Math.atan2(Math.sqrt(a), Math.sqrt(1 - a));

        double distanceKm = EARTH_RADIUS * c;
        return meters ? distanceKm * 1000 : distanceKm;
    }

    private static double haversin(double val) {
        return Math.pow(Math.sin(val / 2), 2);
    }
    
    /**
     * Earth distace used in the experimental survey for the RankGeoFM recommender
     * @param lat1
     * @param lng1
     * @param lat2
     * @param lng2
     */
    public static double getEarthDistance(double lat1,double lng1, double lat2, double lng2) {
    	double radLat1 = Math.toRadians(lat1);
	    double radLat2 = Math.toRadians(lat2);
	    double a = radLat1 - radLat2;
	    double b = Math.toRadians(lng1) - Math.toRadians(lng2);
	    
	    double s = 2 * Math.asin(Math.sqrt(Math.pow(Math.sin(a/2), 2) + Math.cos(radLat1) * Math.cos(radLat2) * Math.pow(Math.sin(b/2), 2)));
	    
	    s*=6371.004;
	    return s;
    }
    

    
    
    
    /**
     * calculate the spherical distance between location(lat1, long1) and location (lat2, long2) (distance from librecRankGeoFM)
     * @param lat1
     * @param long1
     * @param lat2
     * @param long2
     * @return
     */
    public static double getDistance(double lat1, double long1, double lat2, double long2) {
        double a, b, R;
        //earth radius
        R = 6378137;
        lat1 = lat1 * Math.PI / 180.0;
        lat2 = lat2 * Math.PI / 180.0;
        a = lat1 - lat2;
        b = (long1 - long2) * Math.PI / 180.0;
        double sina, sinb;
        sina = Math.sin(a / 2.0);
        sinb = Math.sin(b / 2.0);
        double distance;
        distance = 2 * R
                * Math.asin(Math.sqrt(sina * sina + Math.cos(lat1)
                * Math.cos(lat2) * sinb * sinb));
        return distance / 1000;
    }
    
    /***
     * Another method to compute the distance between two coordinates. This is the method implemented in the USG code provided in 
     * http://spatialkeyword.sce.ntu.edu.sg/eval-vldb17/ It is equivalent to the other
     * 
     * @param lat1 the first latitude
     * @param long1 the first longitude
     * @param lat2 the second latitude
     * @param long2 the second longitude
     * @return the distance between the two coordinates in km
     */
    public static double distance(double lat1, double long1, double lat2, double long2, boolean meters) {
    	double degreesToRadians = Math.PI/180.0;
    	double phi1 = (90 - lat1) * degreesToRadians;
    	double phi2 = (90 - lat2) * degreesToRadians;
	    double theta1 = long1 * degreesToRadians;
	    double theta2 = long2 * degreesToRadians;
	    double cos = (Math.sin(phi1) * Math.sin(phi2) * Math.cos(theta1 - theta2) + Math.cos(phi1) * Math.cos(phi2));
	    double arc = Math.acos(cos);
	    
        double distanceKm = arc * EARTH_RADIUS;
        return meters ? distanceKm * 1000 : distanceKm;
    }
    
    
    /***
     * Method to see if a user is a tourist from his preferences. The method is based in the paper:
     * Automatic construction of travel itineraries using social bread-crumbs
     * 
     * There they consider a tourist every user that have not a difference of 21 days between their first and last timestamp
     * @param prefsUser
     * @param maxDiff
     * @return
     */
    public static <U, I> boolean isTouristByDiffMinMaxTimeStamp(Stream <? extends IdTimePref<I>> prefsUser, double maxDiff) {
    	Tuple2<Long, Long> minMaxTimeStampFourUser = TimeStampUtils.getMinMaxTimeStampOfUser(prefsUser);
    	
    	if ((minMaxTimeStampFourUser.v2 - minMaxTimeStampFourUser.v1) > maxDiff) {
    		return true; 
    	}
    	return false;
    }
    
    /***
     * Method to parse the trajectories of the dataset of IJCAI15 (https://sites.google.com/site/limkwanhui/datacode)
     * @param originalFile the original file
     * @param outputFile the output file
     */
    public static void processCityTrajectoriesIJCAI15(String originalFile, String outputMappingUsers, String outputMappingItems, String outputFile) {
    	try {
    		AtomicInteger countUsers = new AtomicInteger(1);
    		AtomicInteger countItems = new AtomicInteger(1);

    		Map<String, Integer> oldNewUsers = new HashMap<>();
    		Map<String, Integer> oldNewItems = new HashMap<>();

    		
    		//Our format is user, item, rating, timestamp, seqID
			Stream<String> stream = Files.lines(Paths.get(originalFile));
			PrintStream output = new PrintStream(outputFile);
			PrintStream outputMapUsers = new PrintStream(outputMappingUsers);
			PrintStream outputMapItems = new PrintStream(outputMappingItems);
			
			//First line is the header
			stream.skip(1).forEach(line -> {
				String [] data = line.split(";");
				String user = data [1].replace("\"", "");
				String item = data [3];
				String timeStamp = data [2];
				String seqID = data [6];
				
				if (oldNewUsers.get(user) == null) {
					oldNewUsers.put(user, countUsers.get());
					countUsers.incrementAndGet();
				}
				
				if (oldNewItems.get(item) == null) {
					oldNewItems.put(item, countItems.get());
					countItems.incrementAndGet();
				}
				output.println(oldNewUsers.get(user) + "\t" + oldNewItems.get(item) + "\t" + "1.0" + "\t" + timeStamp + "\t" + seqID);
			}); 
			stream.close();
			output.close();
			
			for (String oldUser: oldNewUsers.keySet()) {
				outputMapUsers.println(oldUser + "\t" + oldNewUsers.get(oldUser));
			}
			
			for (String oldItem: oldNewItems.keySet()) {
				outputMapItems.println(oldItem + "\t" + oldNewItems.get(oldItem));
			}
			
			
			outputMapUsers.close();
			outputMapItems.close();
    	} catch (Exception e) {
            e.printStackTrace();
		}
    }
    
    public static void obtainCoordinatesNewIdsIJCAI15(String originalInformationPOI, String mapItemsFile, String outputFile) {
    	Map<String, Integer> mapItems = new HashMap<>();
    	try {
        	PrintStream out = new PrintStream(outputFile);
			Stream<String> stream = Files.lines(Paths.get(mapItemsFile));
			stream.forEach(line -> {
				String [] oldNew = line.split("\t");
				mapItems.put(oldNew[0], Integer.parseInt(oldNew[1]));
			});
			stream.close();
			//Now, read the full POI information 
			Stream<String> stream2 = Files.lines(Paths.get(originalInformationPOI));
			stream2.skip(1).forEach(line -> {
				String [] originalPOIInfo = line.split(";"); 
				String oldId = originalPOIInfo[0];
				Double latitude = Double.parseDouble(originalPOIInfo[2]);
				Double longitude = Double.parseDouble(originalPOIInfo[3]);
				if (mapItems.get(oldId) != null) {
					out.println(mapItems.get(oldId) + "\t" + latitude + "\t" + longitude);
				} else {
					System.out.println("Item (old id): " + mapItems.get(oldId) + " is not in MAP-POI file");
				}
			});
			stream2.close();
			out.close();
    	} catch (Exception e) {
			e.printStackTrace();
		}
    }
    
    
    
    
    
    
    
    /***
	 * Method to parse the trajectories of the tripBuilder dataset (https://github.com/igobrilhante/TripBuilder)
	 * @param originFile
	 * @param outputFile
	 */
	public static void processCityTrajectoriesTripBuilder(String originFile, String outputFile) {
		try {			
			Stream<String> stream = Files.lines(Paths.get(originFile));
			PrintStream out = new PrintStream(outputFile);
			
			AtomicInteger session = new AtomicInteger(1);
			
			stream.forEach(line -> {
				
				//Trajectory of a user
				String [] data = line.split("\t");
				Integer user = Integer.parseInt(data[0]);

				
				for (int i = 1; i < data.length; i++) {
					String [] poiSession = data[i].split(";");
					String item = poiSession[0];
					String startingTime = poiSession[2];
					
					//User, item, time and session
					out.println(user + "\t" + item + "\t" + "1.0" + "\t" + startingTime + "\t" + session.get());
				}
				session.incrementAndGet();
			});
			stream.close();
			out.close();

		} catch (Exception e) {
			
		}
	} 
	
	/***
	 * Method to compute the matrix distance matrix of a list of pois
	 * @param itemStream the stream of items
	 * @param itemIndexes the indexes
	 * @param coordinates the coordinates
	 * @return a matrix with the same indexes as the FastItemIndex
	 */
	public static <I> float [][] distanceMatrix(List<I> itemStream, FastItemIndex<I> itemIndexes, Map<I, Tuple2<Double, Double>> coordinates, boolean meters){
		int size = itemStream.size();
		
		float [][] result = new float[size][size];
		
		for (I item1: itemStream) {
			int indexX = itemIndexes.item2iidx(item1);
			double latx = coordinates.get(item1).v1; 
			double longx = coordinates.get(item1).v2;		
					
			for (I item2: itemStream) {
				if (!item1.equals(item2)) {
					int indexY = itemIndexes.item2iidx(item2);
					if (result[indexX][indexY] == 0.0 || result[indexY][indexX] == 0.0) {
						double laty = coordinates.get(item2).v1; 
						double longy = coordinates.get(item2).v2;		
						double dist = POIProcessData.haversine(latx, longx, laty, longy, meters);
						result[indexX][indexY] = (float) dist;
						result[indexY][indexX] = (float) dist;
					}
				}
			}
		}
		
		return result;
	}

	
	
	/***
	 * Method to read a categories grouped file composed by two columns.
	 * The first colum is the groupped category and the second column is the unique category (key)
	 * @param categoriesFile
	 * @return a map of the unique categories and their grouped one
	 */
	public static Map<String, Set<String>> readGroupedCategoriesTripBuilder (String categoriesFile){
		Map<String, Set<String>> result = new HashMap<>();
		Stream<String> stream;
		int columnGrouped = 0;
		int columnUniq = 1;
		try {
			stream = Files.lines(Paths.get(categoriesFile));
			//It has a header
			stream.skip(1).forEach(line -> {
				String [] data  = line.split("\t");
				if (result.get(data[columnGrouped]) == null) {
					result.put(data[columnGrouped], new HashSet<>());
				}
				result.get(data[columnGrouped]).add(data[columnUniq]);
			});
			stream.close();
			return result;
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return null;

	}
	
	public static Map<String, String> readUniqueCategoriesTripBuilder (String categoriesFile){
		Map<String, String> result = new HashMap<>();
		Stream<String> stream;
		int columnGrouped = 0;
		int columnUniq = 1;
		try {
			stream = Files.lines(Paths.get(categoriesFile));
			//It has a header
			stream.skip(1).forEach(line -> {
				String [] data  = line.split("\t");
				result.put(data[columnUniq], data[columnGrouped]);
			});
			stream.close();
			return result;
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return null;

	}
	
	
	/***
	 * This method will receive the mapping between the gruopped ones and their set of unique categories and will retrieve list of tuples of the 
	 * frequency of each grouped category ordered by its popularity
	 * @param mapGroupedCategories
	 * @return
	 */
	public static List<String> getGroupedCategoriesFreqTripbuilder(Map<String, Set<String>> mapGroupedCategories){
		List<Tuple2<String, Double>> result = new ArrayList<>();

		for (String groupped: mapGroupedCategories.keySet()) {
			result.add(new Tuple2<>(groupped, (double) mapGroupedCategories.get(groupped).size()));
		}
		
		Collections.sort(result, new WeightComparatorTuple2String().reversed());
		return result.stream().map(tuple2 -> tuple2.v1).collect(Collectors.toList());
	}
	
	/***
	 * 
	 * @param originalFile the original poi-clusters file from Tripbuilder
	 * @param groupedCategoriesFile
	 * @param destinationFile
	 */
	public static void generatePOIGrupedCategoryTripBuilderByPopularity(String originalFile, String groupedCategoriesFile, String destinationFile) {
		Map<String, Set<String>> groupedCategories = readGroupedCategoriesTripBuilder(groupedCategoriesFile);
		List<String> categoriesGroupOrdered = getGroupedCategoriesFreqTripbuilder(groupedCategories);
		
		PrintStream out;
		Stream<String> stream;
		int columnCategories = 5;
		int columnCluster = 0;
		try {
			out = new PrintStream(destinationFile);
			stream = Files.lines(Paths.get(originalFile));
			stream.forEach(line -> {
				String data [] = line.split(",");
				String categories [] = data[columnCategories].split(";");
				String finalCategory = "";
				for (String popCat: categoriesGroupOrdered) {
					for (String uniqCate: groupedCategories.get(popCat)) {
						if (ArrayUtils.contains(categories, uniqCate)) {
							finalCategory = popCat;
							break;
						}
					}
						
				}
				out.println(data[columnCluster] + "\t" + finalCategory);
			});
			
			out.close();
			stream.close();
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}
	
	/***
	 * 
	 * @param originalFile
	 * @param groupedCategoriesFile
	 * @param destinationFile
	 */
	public static void generatePOIGrupedCategoryTripBuilderByMax(String originalFile, String groupedCategoriesFile, String destinationFile) {
		List<String> possible = new ArrayList<>();
		possible.add("Religious");
		possible.add("Museum");
		possible.add("Architectural");
		possible.add("Shops & Services");
		possible.add("Travel & Transport");
		possible.add("Outdoors & Recreation");
		
		
		Map<String, String> groupedCategories = readUniqueCategoriesTripBuilder(groupedCategoriesFile);
		
		PrintStream out;
		Stream<String> stream;
		int columnCategories = 5;
		int columnCluster = 0;
		try {
			out = new PrintStream(destinationFile);
			stream = Files.lines(Paths.get(originalFile));
			stream.skip(1).forEach(line -> {
				String data [] = line.split(",");
				String categories [] = data[columnCategories].split(";");
				String finalCategory = "";
				Map<String, Integer> grupedCatsCont = new HashMap<>();
				
				//For finding the maximum
				for (String uniqCate: categories) {
					String grouppedAssociated = groupedCategories.get(uniqCate);
					if (grupedCatsCont.get(grouppedAssociated) == null) {
						grupedCatsCont.put(grouppedAssociated, 0);
					}
					grupedCatsCont.put(grouppedAssociated, grupedCatsCont.get(grouppedAssociated) + 1);
				}
				
				finalCategory = getFinalGroupCategory(possible, grupedCatsCont);
				out.println(data[columnCluster] + "\t" + finalCategory);
			});
			
			out.close();
			stream.close();
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}
	
	private static String getFinalGroupCategory(List<String> possible, Map<String, Integer> grouppedCont) {
		
		List<Tuple2<String, Double>> lstOrdered = new ArrayList<>();
		for (String un: grouppedCont.keySet()) {
			lstOrdered.add(new Tuple2<>(un, (double) grouppedCont.get(un)));
		}
		lstOrdered = lstOrdered.stream().sorted(new WeightComparatorTuple2String().reversed()).collect(Collectors.toList());
		
		List<String> result = new ArrayList<>();
		double maxVal = lstOrdered.get(0).v2;
		for (int i = 0; i < lstOrdered.size(); i++) {
			if (lstOrdered.get(i).v2 < maxVal) {
				break;
			}
			result.add(lstOrdered.get(i).v1);
		}
		if (result.size() == 1) {
			return result.get(0);
		}
		
		for (String p: possible) {
			if (result.contains(p)) {
				return p;
			}
		}
		
		return null;
	}
	
	/***
	 * Method to generate a file with the coordinates of the POIS from the original poi-cluster file
	 * @param originalFile the original file
	 * @param outputFile the output file
	 */
	public static void generateClustersCoordsTripbuilder(String originalFile, String outputFile) {
		try {
			PrintStream out = new PrintStream(outputFile);
			Stream<String> stream = Files.lines(Paths.get(originalFile));
			stream.skip(1).forEach(line -> {
				String [] data = line.split(",");
				String clusterId = data[0];
				Double lat = Double.parseDouble(data[2]);
				Double longitude = Double.parseDouble(data[3]);
				out.println(clusterId + "\t" + lat + "\t" + longitude);
			});
			stream.close();
			out.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
		
	}

	
	/***
	 * Method to generate a file containing
	 * @param trainFile
	 * @param cityPoiFile
	 * @param POIcoordFile
	 * @param poiFeature
	 * @param output
	 */
	public static void generateUserVariablesForClusteringFoursquare(String trainFile, String cityPoiFile,
			String POIcoordFile, String cityCoordFile, String poiFeature, String output) {
		
		Set<String> variablesToConsider = new LinkedHashSet<>();
		variablesToConsider.add("ID_User"); //done
		variablesToConsider.add("Number_Of_Countries_Visited"); //done
		variablesToConsider.add("Number_Of_Cities_Visited"); //done
		
		variablesToConsider.add("DistanceKm_Between_Cities_Home_Visited"); //done
		variablesToConsider.add("Average_DistanceKm_Between_Cities_Home_Visited"); //done
		variablesToConsider.add("DistanceKm_Between_Checkins_In_Same_Cities"); //done
		variablesToConsider.add("Average_DistanceKm_Between_Checkins_In_Same_Cities"); //done

		variablesToConsider.add("AverageDistanceKm_Between_Checkins");
		
		variablesToConsider.add("Number_User_Checkins"); //done
		variablesToConsider.add("Number_Unique_User_Checkins"); //done
		
		variablesToConsider.add("Days_Registered_Foursquare"); // done
		
		variablesToConsider.add("Total_Number_Different_Trayectories"); // done
		variablesToConsider.add("DistanceKm_Trayectories"); //done
		variablesToConsider.add("Average_DistanceKm_Trayectories"); //done
		
		variablesToConsider.add("Total_Number_Different_Trips"); // done
		variablesToConsider.add("Total_Trip_TimeDays"); //done
		variablesToConsider.add("Average_Trip_TimeDays"); //done

		
		variablesToConsider.add("Total_Different_Features"); //done
		variablesToConsider.add("Average_Days_In_City"); //done
		variablesToConsider.add("Average_Checkins_PerDay_In_Foursquare"); //done

		
		
		//Media de checkins al dia

		
		try {
			FeatureData<Long, String, Double> featureData =  SimpleFeatureData.load(SimpleFeaturesReader.get().read(poiFeature, lp, sp));
			FastTemporalPreferenceDataIF<Long, Long> temporalData = ExperimentUtils.loadTrainFastTemporalFeaturePreferenceData(trainFile, trainFile, false, true);		
			Map<Long, Tuple2<Double, Double>> coordinatesPOIs = SequentialRecommendersUtils.POICoordinatesMap(POIcoordFile, lp);
			Map<String, Tuple2<Double, Double>> coordinatesCities = SequentialRecommendersUtils.POICoordinatesMap(cityCoordFile, sp);

			Map<Long, String> poiCity = poiCityFile(cityPoiFile, lp);
			
			UsersSessions<Long, Long> userSessions = new UsersSessions<Long, Long>(temporalData, null);
			
			PrintStream out = new PrintStream(output);
			
			for (String variable: variablesToConsider) {
				out.print(variable + "\t");
			}
			out.println();

			
			temporalData.getUsersWithPreferences().forEach(u -> {
				
				//Ordered the interactions in a temporal way
				List<IdTimePref<Long>> userData = temporalData.getUserPreferences(u).sorted(PreferenceComparators.timeComparatorIdTimePref).collect(Collectors.toList());
				////////////////////////////////
				String homeCityUser = homeCityUser(userData, poiCity);
				
				
				List<Double> values = new ArrayList<>();
				
				Tuple2<Set<String>, Set<String>> countriesCitiesVisited = countriesCitiesVisited(userData, poiCity);
				
				//Number of Countries
				values.add((double) countriesCitiesVisited.v1.size());
				//Number of Cities
				values.add((double) countriesCitiesVisited.v2.size());
				
				Tuple2<Double, Double> distanceAndAverageDistanceBetweenCities = getDistanceAndAverageDistanceBetweenCitiesVisitedAndHome(userData, coordinatesCities, poiCity, homeCityUser);
				
				//Distance cities
				values.add(distanceAndAverageDistanceBetweenCities.v1);

				//Average distance cities
				values.add(distanceAndAverageDistanceBetweenCities.v2);
				
				//distance and average distance between check-ins in the same city
				Tuple2<Double, Double> distanceAndAverageBetweenCheckinsSameCities = getDistanceAndAverageBetweenCheckinsSameCities(userData, countriesCitiesVisited.v2, coordinatesPOIs, poiCity);
				values.add(distanceAndAverageBetweenCheckinsSameCities.v1);
				values.add(distanceAndAverageBetweenCheckinsSameCities.v2);
				
				//AverageDistance between all check-ins
				double averageDistance = getAverageDistanceBetweenAllCheckins(userData, coordinatesPOIs);
				values.add(averageDistance);
				
				//number of check-ins of the user
				Double numberUserCheckins = (double) userData.size();
				values.add(numberUserCheckins);
				
				Set<Long> uniquePois = new HashSet<>();
				userData.stream().forEach(pref -> uniquePois.add(pref.v1));
				
				//number of unique check-ins of the user
				Double numberUniqueUserCheckins = (double) uniquePois.size();
				values.add(numberUniqueUserCheckins);
				
				//Number of days in Foursquare
				double daysInFoursquare = daysInFoursquare(userData);
				values.add(daysInFoursquare);
				
				////////////////////////////////////////////////////////////////
				// 8 * 3600 = intervals of 8 hours
				List<List<IdTimePref<Long>>> trajectories = userSessions.getUserSession(u, 8 * 3600,  Double.MAX_VALUE, UserSession.TIME);
				double totalTrajectories = trajectories.size();
				
				//total trajectories
				values.add(totalTrajectories);
				
				Tuple2<Double, Double> averageNumberTrajectoriesLength = getDistanceAndAverageDistanceOfTrajectories(trajectories, coordinatesPOIs);
				
				//Total trajectories Length
				values.add(averageNumberTrajectoriesLength.v1);
				
				//average number trajectories length
				values.add(averageNumberTrajectoriesLength.v2);
				
				//numberTrips and AverageTimeof the trips (outside home City)
				Tuple3<Double, Double, Double> numberTripsAverageTime = numberTripsAverageTimeOutsideHomeCity(userData, poiCity, homeCityUser);
				double numberTrips = numberTripsAverageTime.v1;
				double totalTime = numberTripsAverageTime.v2;
				double averageTime = numberTripsAverageTime.v3;
				values.add(numberTrips);
				values.add(totalTime);
				values.add(averageTime);

				
				
				double totalFeaturesUser = totalFeaturesUser(userData, featureData);
				//total features user
				values.add(totalFeaturesUser);
				
				double averageDaysInCity = averageDaysInCities(userData, poiCity);
				
				//Average days in city
				values.add(averageDaysInCity);
				
				
				double averageCheckinsPerDay = (double) userData.size() / daysInFoursquare;
				
				//average checkins per day
				values.add(averageCheckinsPerDay);
				
				out.print(u + "\t");
				
				for (Double val: values) {
					out.print(val + "\t");
				}
				out.println();
				
			});
			out.close();
			
			
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}	

	}
	
	/***
	 * Method to obtain the home city of the user
	 * @param userData
	 * @param poiCity the map between POIs and the cities
	 * @return the homeCity of the user
	 */
	private static String homeCityUser(List<IdTimePref<Long>> userData, Map<Long, String> poiCity) {
		Map<String, Integer> prefsByCity = new HashMap<>();

		String homeCity="";
		Integer maxVisits=Integer.MIN_VALUE;
		
		for (IdTimePref<Long> pref: userData) {
			String city = poiCity.get(pref.v1);
			if (city != null) {
				if (!prefsByCity.containsKey(city)) {
					prefsByCity.put(city,  0);
				}
				prefsByCity.put(city, prefsByCity.get(city) + 1);
			}
		}
		
		for (String city: prefsByCity.keySet()) {
			Integer visits = prefsByCity.get(city);
			if (visits > maxVisits) {
				homeCity = city;
				maxVisits = visits;
			}
		}
		return homeCity;
	}
	
	private static Tuple3<Double, Double, Double> numberTripsAverageTimeOutsideHomeCity(List<IdTimePref<Long>> userData, Map<Long, String> poiCity, String homeCity) {
		List<List<IdTimePref<Long>>> tripsOutsideHome = new ArrayList<>();
		
		
		//We have home city. Now, we find the trips outside the homeCity
		List<IdTimePref<Long>> possibleTrip = new ArrayList<>();;
		
		boolean computingTrip = false;
		for (int i = 1; i < userData.size(); i++) {
			IdTimePref<Long> prevPref = userData.get(i - 1);
			IdTimePref<Long> actualPref = userData.get(i);
			
			String cityPrev = poiCity.get(prevPref.v1);
			String actualCity = poiCity.get(actualPref.v1);

			
			if (homeCity.equals(cityPrev) && homeCity.equals(actualCity)) {
				//home city equals both cities
				continue;
			}
			
			if (homeCity.equals(cityPrev) && !homeCity.equals(actualCity)) {
				//Possible NEW trip
				possibleTrip = new ArrayList<>();
				computingTrip = true;
				continue;
			}
			
			//Finishing the trip
			if (!cityPrev.equals(actualCity) && actualCity.equals(homeCity)) {
				possibleTrip.add(prevPref);
				tripsOutsideHome.add(possibleTrip);
				possibleTrip = new ArrayList<>();
				computingTrip = false;
			} else {
				computingTrip = true;
				possibleTrip.add(prevPref);
			}

			
		}
		
		//Only for trips that the user does not return to the home city
		if (computingTrip) {
			tripsOutsideHome.add(possibleTrip);
		}
		
		
		double numberTrips = 0;
		double averageTime = 0;

		//Now, analyzing the trips
		for (List<IdTimePref<Long>> possibleTripAux: tripsOutsideHome) {
			if (possibleTripAux.size() < 2)
				continue;
			Long fstTime = possibleTripAux.get(0).v3;
			Long LastTime = possibleTripAux.get(possibleTripAux.size() - 1).v3;
			
			for (IdTimePref<Long> pref: possibleTripAux) {
				String city = poiCity.get(pref.v1);
				if (city.equals(homeCity)) {
					System.out.println("ERROR");
				}
			}
			
			//86400 are the seconds per day and 7 implies a week
			//If the if is satisfied, there is a trip
			if ((LastTime - fstTime) / (86400 * 7) > 1) {
				numberTrips+=1;
				averageTime+=(LastTime - fstTime) / (86400);
			}
			
			
		}
		
		double totalTime = averageTime;
		if (averageTime > 0 ) {
			averageTime /=numberTrips;
		}
		
		
		return new Tuple3<>(numberTrips, totalTime, averageTime);
		
	}
	
	private static double totalFeaturesUser(List<IdTimePref<Long>> userData, FeatureData<Long, String, Double> featureData) {
		Set<String> totalFeatures = new HashSet<String>();
		
		userData.stream().forEach(pref -> {
			featureData.getItemFeatures(pref.v1).forEach(feat ->{
				totalFeatures.add(feat.v1);
			});
		});
		return totalFeatures.size();
	}
	
	/**
	 * Method to compute the average Days in cities
	 * @param userData
	 * @param poiCity
	 * @return
	 */
	private static double averageDaysInCities(List<IdTimePref<Long>> userData, Map<Long, String> poiCity) {
		double averageDaysCity = 0.0;
		int totalCitiesChange = 1;
		
		if (userData.size() > 1) {
			
			String prevCity = null;
			Long prevTimestamp = null;
			for (IdTimePref<Long> pref: userData) {
				
				//First interaction
				if (prevCity == null) {
					prevCity = poiCity.get(pref.v1);
					prevTimestamp = pref.v3;
				} else {
					//We have a previous interaction, we check if the user is in the same city
					String actualCity = poiCity.get(pref.v1);
					Long actualTime = pref.v3;
					if (prevCity.equals(actualCity)) {
						averageDaysCity+=Math.abs(actualTime-prevTimestamp);
					} else {
						//change over city
						totalCitiesChange++;
						prevCity = actualCity;						
					}
					prevTimestamp = actualTime;

					
				}
			}
			
		}
		if (averageDaysCity == 0) {
			return 1;
		}
			
		return Math.max((averageDaysCity / (3600 * 24)) / totalCitiesChange, 1);
	}
	
	/***
	 * Method to compute the average length of the trajectories followed by the user
	 * @param trajectories the total trajectories for a specific user
	 * @param coordinatesPOIs the coordinates of the POIs
	 * @return
	 */
	private static Tuple2<Double, Double> getDistanceAndAverageDistanceOfTrajectories(List<List<IdTimePref<Long>>> trajectories, Map<Long, Tuple2<Double, Double>> coordinatesPOIs) {
		
		double averageTrajectoriesLength = 0.0;
		
		for (List<IdTimePref<Long>> trajectory: trajectories) {
			
			if (trajectory.size() > 1) {
				for (int i = 1; i< trajectory.size(); i++) {
					Tuple2<Double, Double> coordinatesPOIFirst = coordinatesPOIs.get(trajectory.get(i - 1).v1);
					Tuple2<Double, Double> coordinatesPOISecond = coordinatesPOIs.get(trajectory.get(i).v1);
					
					averageTrajectoriesLength+= haversine(coordinatesPOIFirst.v1, coordinatesPOIFirst.v2, coordinatesPOISecond.v1, coordinatesPOISecond.v2, false);
				}
			}
			
		}

		
		
		
		return new Tuple2<>(averageTrajectoriesLength, averageTrajectoriesLength/ trajectories.size());
	}

	/***
	 * Method to compute the number of days in the Foursquare dataset
	 * @param userCheckins  the total checkins of the user
	 * @return
	 */
	private static double daysInFoursquare(List<IdTimePref<Long>> userCheckins) {
		userCheckins.sort(PreferenceComparators.timeComparatorIdTimePref);
		
		if (userCheckins.size() == 1) {
			return 1;
		}
		
		Long timestampFst = userCheckins.get(0).v3;
		Long timestampLast = userCheckins.get(userCheckins.size() - 1).v3;
		
		return Math.max(((timestampLast - timestampFst) / (3600 * 24)), 1);
		
	}
	
	/**
	 * Method to get the average distance of all checkins of the user
	 * @param userCheckins
	 * @param coordinatesPOIs
	 * @return
	 */
	private static double getAverageDistanceBetweenAllCheckins(List<IdTimePref<Long>> userCheckins, Map<Long, Tuple2<Double, Double>> coordinatesPOIs) {
		double averageDistance = 0;
		
		for (int i = 1; i < userCheckins.size(); i++) {
			Tuple2<Double, Double> coordinatesPOIFirst = coordinatesPOIs.get(userCheckins.get(i - 1).v1);
			Tuple2<Double, Double> coordinatesPOISecond = coordinatesPOIs.get(userCheckins.get(i).v1);
			
			averageDistance+=haversine(coordinatesPOIFirst.v1, coordinatesPOIFirst.v2, coordinatesPOISecond.v1, coordinatesPOISecond.v2, false);

		}
		return averageDistance / userCheckins.size();
	}
	
	/***
	 * Method to compute the total distance between the checkins that belong to the same city (if the user visited more than 1 city, we consider all that checkins separately)
	 * @param userCheckins
	 * @param citiesVisited
	 * @param coordinatesPOIs
	 * @param poiCity
	 * @return
	 */
	private static Tuple2<Double, Double> getDistanceAndAverageBetweenCheckinsSameCities(List<IdTimePref<Long>> userCheckins, Set<String> citiesVisited, Map<Long, Tuple2<Double, Double>> coordinatesPOIs, Map<Long, String> poiCity) {
		double totalDistance = 0.0;
		
		for (String city: citiesVisited) {
			//Only the POis in that city
			List<Long> poisAndTimestampsCity = userCheckins.stream().map(pref -> pref.v1).filter(poi -> poiCity.get(poi).equals(city)).collect(Collectors.toList());
			
			if (poisAndTimestampsCity.size() > 1) {
				for (int i = 1; i < poisAndTimestampsCity.size(); i++) {
					Tuple2<Double, Double> coordinatesPOIFirst = coordinatesPOIs.get(poisAndTimestampsCity.get(i - 1));
					Tuple2<Double, Double> coordinatesPOISecond = coordinatesPOIs.get(poisAndTimestampsCity.get(i));
					
					
					totalDistance+=haversine(coordinatesPOIFirst.v1, coordinatesPOIFirst.v2, coordinatesPOISecond.v1, coordinatesPOISecond.v2, false);

					
				}
			}
			
		}
		return new Tuple2<>(totalDistance, totalDistance/citiesVisited.size());
	}
	
	/***
	 * Method to compute the total distance and the average between the cities visited
	 * @param userData 
	 * @param coordinatesCities
	 * @param citiesVisited
	 * @return
	 */
	private static Tuple2<Double, Double> getDistanceAndAverageDistanceBetweenCitiesVisitedAndHome(List<IdTimePref<Long>> userData, Map<String, Tuple2<Double, Double>> coordinatesCities, Map<Long, String> poiCity, String homeCity) {
		double totalDistance = 0;
	
		Set<String> citiesVisited = new HashSet<>();
		userData.stream().forEach(p -> citiesVisited.add(poiCity.get(p.v1)));
		
		
		if (citiesVisited.size() == 1) {
			return new Tuple2<Double, Double>(0.0, 0.0);
		}
		
		//Removing the homeCity
		citiesVisited.remove(homeCity);
		Tuple2<Double, Double> coordinatesCitiesHome = coordinatesCities.get(homeCity);

		for (String cityVisited: citiesVisited) {
			Tuple2<Double, Double> coordinatesCitiesSecond = coordinatesCities.get(cityVisited);
			
			totalDistance += haversine(coordinatesCitiesHome.v1, coordinatesCitiesHome.v2, coordinatesCitiesSecond.v1, coordinatesCitiesSecond.v2, false);
			
		}
		
		//Second value is the average
		return new Tuple2<Double, Double>(totalDistance, totalDistance/(double)citiesVisited.size());
	}

	/***
	 * Method to obtain the sets of the countries and the cities visited by the user
	 * @param userCheckins
	 * @param poiCity
	 * @return
	 */
	private static Tuple2<Set<String>, Set<String>> countriesCitiesVisited(List<IdTimePref<Long>> userCheckins, Map<Long, String> poiCity){
		Set<String> cities = new HashSet<>();
		Set<String> country = new HashSet<>();
		
		userCheckins.stream().forEach(pref -> {
			Long poi = pref.v1;
			
			String city = poiCity.get(poi);

			if (city != null) {
				String [] data = city.split("_");
				String countryOfPoi = data [0];
				cities.add(city);
				country.add(countryOfPoi);
			}
			
		});
		
		return new Tuple2<>(country, cities);
	}

	private static <I> Map<I, String> poiCityFile (String cityPoiFile, Parser<I> parserPoi){
		Map<I, String> mapPOICityFile = new HashMap<>();

		try {			
			Stream<String> stream = Files.lines(Paths.get(cityPoiFile));
			stream.forEach(line -> {
				String [] data = line.split("\t");
				String poiId = data [0];
				String cityPoi = data[1];
				
				mapPOICityFile.put(parserPoi.apply(poiId), cityPoi);
			});
			stream.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		return mapPOICityFile;
	}

	public static void analizeTravelersLocals(String trainFile, FastTemporalPreferenceDataIF<Long, Long> rankSysData,
			String selectedCity, String homeCities) {
		Map<Long, String> homeCitiesUsers = new HashMap<>();
		Map<String, Set<Long>> clustersAndUsers = new HashMap<>();
		
		try {			
			Stream<String> stream = Files.lines(Paths.get(homeCities));
			stream.skip(1).forEach(line -> {
				String [] data = line.split("\t");
				Long user = Long.parseLong(data[0]);
				String city = data[1];
				homeCitiesUsers.put(user, city);				
			});
			stream.close();
			Stream<String> stream2 = Files.lines(Paths.get(trainFile));
			stream2.skip(1).forEach(line -> {
				String [] data = line.split(",", Integer.MAX_VALUE);
				Long user = Long.parseLong(data[0]);
				String local = data[1];
				String traveler = data[2];
				
				
				if (rankSysData.containsUser(user)) {
					
					if (local.length() > 1) {
						if (!clustersAndUsers.containsKey(local)) {
							clustersAndUsers.put(local, new HashSet<>());
						}
						
						if (homeCitiesUsers.get(user).equals(selectedCity)) {
							clustersAndUsers.get(local).add(user);
						}
						
					}
					
					if (traveler.length() > 1 && !traveler.equals("non-traveler") && !traveler.equals("unclear_home")) {
						if (!clustersAndUsers.containsKey(traveler)) {
							clustersAndUsers.put(traveler, new HashSet<>());
						}
						
						//then, she is a traveler
						if (!homeCitiesUsers.get(user).equals(selectedCity)) {
							clustersAndUsers.get(traveler).add(user);
						}
					}
					
					
					
				
				}
				
			});
			stream2.close();
			System.out.println("Number of users in the test set " + rankSysData.numUsers());
			for (String c: clustersAndUsers.keySet()) {
				System.out.println(c + " in city of " + selectedCity + " are: " + clustersAndUsers.get(c).size());
			}

			
		} catch (Exception e) {
			e.printStackTrace();
		}
		
	}
	
	


	
	
    

	

}
