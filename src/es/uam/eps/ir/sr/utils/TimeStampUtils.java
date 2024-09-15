package es.uam.eps.ir.sr.utils;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Calendar;
import java.util.Collections;
import java.util.Date;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.atomic.AtomicLong;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import org.jooq.lambda.tuple.Tuple2;
import org.jooq.lambda.tuple.Tuple3;

import es.uam.eps.ir.crossdomainPOI.datamodel.temporal.interfaces.FastTemporalPreferenceDataIF;
import es.uam.eps.ir.crossdomainPOI.datamodel.temporal.preferences.IdTimePref;



/***
 * Final class that contains methods working with temporal data.
 * 
 * @author Pablo Sanchez (psperez@icai.comillas.edu)
 *
 */
public final class TimeStampUtils {
	
	public enum TimestampStrategy {MAX_TIMESTAMP_PERUSER, MAX_TIMESTAMP_TRAIN, MIN_TIMESTAMP_PERUSER, MIN_TIMESTAMP_TRAIN};

	
	/***
	 * Method to obtain the hour, the minute and the second from a timestamp
	 * @param timestamp the timestamp
	 * @param multiply boolean indicating if we are going to multiply the timestamp (java work on miliseconds, not seconds)
	 * @return a tuple with the hour, the minute and the second of the timestamp
	 */
	public static Tuple3<Integer, Integer, Integer> getHourMinSecondfromTimestamp(Long timestamp, boolean multiply){
		Long finalTime = timestamp;
		if (multiply) {
			finalTime *= 1000;
		}
		
		Date date = new Date(finalTime); 
		Calendar calendar = Calendar.getInstance(); // creates a new calendar instance
		calendar.setTime(date);   // assigns calendar to given date 
		return new Tuple3<>(calendar.get(Calendar.HOUR_OF_DAY), calendar.get(Calendar.MINUTE), calendar.get(Calendar.SECOND));
	}
	
	/**
     * *
     * Method to get the largest timestamp associated with a datamodel
     *
     * @param dataModel the datamodel
     * @return a 
     */
    public static <U, I> Tuple2<Long,Long> getMinMaxTimestampOfDataset(FastTemporalPreferenceDataIF<U, I> dataModel) {
        AtomicLong minTime = new AtomicLong(Long.MAX_VALUE);
        AtomicLong maxTime = new AtomicLong(-1);

        
        dataModel.getUsersWithPreferences().forEach(u -> {
            long maxOfUser = dataModel.getUserPreferences(u).mapToLong(p -> p.v3).max().getAsLong();
            if (maxTime.get() < maxOfUser) {
                maxTime.set(maxOfUser);
            }
            long minOfUser = dataModel.getUserPreferences(u).mapToLong(p -> p.v3).min().getAsLong();
            if (minTime.get() > minOfUser) {
                minTime.set(minOfUser);
            }
        });

        return new Tuple2<>(minTime.get(), maxTime.get());
    }
        
    /***
     * Get the number of days between two timestamps. 2 second based timestamps
     * @param t1 the first timestamp
     * @param t2 the second timestamp
     * @return the number of days between those timestamps (this method is symmetric, it returns the number of days in absolute)
     */
    public static int daysBetweenTimestamps(long t1, long t2) {
		 Long diff = Math.abs(t1 - t2);
		 return Math.round(diff / (24 * 3600));
	}
    

    /**
     * Method to get the average timestamps of an item
     *
     * @param dataModel the datamodel
     * @param item the item
     * @return the average of timestamps associated
     */
    public static <U, I> double getAverageTimeStampOfItem(Stream<? extends IdTimePref<U>> prefsItem) {
        AtomicLong acc = new AtomicLong(0);
        AtomicLong count = new AtomicLong(0);

        prefsItem.forEach(timePref -> {
            acc.set(acc.get() + timePref.v3);
            count.getAndIncrement();
        });

        return (double) acc.get() / (double) count.get();
    }

    /**
     * *
     * Method to get the minimum timestamp associated with an item
     *
     * @param dataModel the datamodel
     * @param item the item to compute the minimum timestamp
     * @return the minimum timestamp associated
     */
    public static <U, I> Tuple2<Long, Long>  getMinMaxTimeStampOfItem(Stream<? extends IdTimePref<U>> prefsItem) {
        AtomicLong minTime = new AtomicLong(Long.MAX_VALUE);
        AtomicLong maxTime = new AtomicLong(-1);

        prefsItem.forEach(timePref -> {
            if (timePref.v3 < minTime.get()) {
            	minTime.set(timePref.v3);
            }
            
            if (timePref.v3 > maxTime.get()) {
            	maxTime.set(timePref.v3);
            }
        });
    	return new Tuple2<>(minTime.get(), maxTime.get());
    }

    
    /***
     * Method to get the minimum and the maximum timestamp of the user
     * @param prefsUser the user preferences
     * @return a tuple with the minimum and maximum timestamp
     */
    public static <U, I> Tuple2<Long, Long> getMinMaxTimeStampOfUser(Stream<? extends IdTimePref<I>> prefsUser){
    	AtomicLong minTime = new AtomicLong(Long.MAX_VALUE);
    	AtomicLong maxTime = new AtomicLong(-1);
    	prefsUser.forEach(pref -> {
    		if (pref.v3 > maxTime.get()) {
    			maxTime.set(pref.v3);
    		}
    		if (pref.v3 < minTime.get()) {
    			minTime.set(pref.v3);
    		}
    	});
    	return new Tuple2<>(minTime.get(), maxTime.get());
    }


    /**
     * Method to get the median of the timestamps associated with an item
     *
     * @param dataModel the datamodel
     * @param item the item
     * @return the median associated with that item's timestamp
     */
    public static <U, I> Number getMedianTimeStampOfItem(Stream<? extends IdTimePref<U>> prefsItem) {
        List<Number> storation = new ArrayList<Number>();
        prefsItem.forEach(timepref -> {
            storation.add(timepref.v3);
        });
        return getMedian(storation);
    }
    
    public static Number getMedian(List<Number> lstNumbers) {
    	 if (lstNumbers.size() == 0) {
    		 return Double.NaN;
    	 }
    	
    	 Number[] toArr = new Number[lstNumbers.size()];
         toArr = lstNumbers.toArray(toArr);
         Arrays.sort(toArr);
         if (toArr.length % 2 == 0) {
             return (toArr[toArr.length / 2].doubleValue() + toArr[toArr.length / 2 - 1].doubleValue()) / 2.0;
         } else {
             return toArr[toArr.length / 2];
         }

    }
    
    public static List<Number> getQuartiles(List<Number> lstNumbers){
    	List<Number> result = new ArrayList<>();
    	Number q2 = getMedian(lstNumbers);
    	List<Number> prevQ2 = new ArrayList<>();
    	List<Number> postQ2 = new ArrayList<>();
    	for (Number n : lstNumbers) {
    		if (n.doubleValue() == q2.doubleValue()) {
    			postQ2.add(n);
    			prevQ2.add(n);
    		}
    		else {
	    		if (n.doubleValue() > q2.doubleValue()) {
	    			postQ2.add(n);
	    		} else {
	    			prevQ2.add(n);
	    		}
    		}
    	}
    	Number q1 = getMedian(prevQ2);
    	Number q3 = getMedian(postQ2);
        result.add(q1);
        result.add(q2);
        result.add(q3);
    	return result;
    }
    
    
  
    /***
     * Method to obtain a single timeStamp for representing the specific user
     * @param u the user
     * @param strat the timestamp strategy
     * @param data the temporal preference data
     * @param minTimestamp the minimum timestamp
     * @param maxTimestamp the maximum timestamp
     * @return the timestamp
     */
    public static <U, I> long getSpecificTimeStampForUser(U u, TimestampStrategy strat, FastTemporalPreferenceDataIF<U, I> data, Long minTimestamp, Long maxTimestamp) {
		switch (strat) {
			case MAX_TIMESTAMP_PERUSER:
				return data.getUserPreferences(u).mapToLong(IdTimePref<I>::v3).max().getAsLong();
			case MAX_TIMESTAMP_TRAIN:
				return minTimestamp;
			case MIN_TIMESTAMP_PERUSER:
				return data.getUserPreferences(u).mapToLong(IdTimePref<I>::v3).min().getAsLong();
			case MIN_TIMESTAMP_TRAIN:
				return maxTimestamp;
			default:
				return -1L;		
		}
			 
	}
    
    /**
     * Method to obtain the hours of a full stream of timestamps ()
     * @param timestamps
     * @return hours between the first and last timestamps
     */
	public static double getHoursListTime(List<Long> timestamps) {
		Collections.sort(timestamps);
		return (timestamps.get(timestamps.size() - 1) - timestamps.get(0)) / (3600.0);
	}
	
	public static <U, I> Map<I, Long> stimateTimeTransitionPOI(FastTemporalPreferenceDataIF<U, I> prefData) {
		Map<I, Long> result = new HashMap<>();
		Map<I, List<Long>> itemTransitions = new HashMap<>();
		
		prefData.getUsersWithPreferences().forEach(u -> {
			List<IdTimePref<I>> userPrefs = prefData.getUserPreferences(u).collect(Collectors.toList());
			Long prev = userPrefs.get(0).v3;
			I actualId = userPrefs.get(0).v1;
			if (userPrefs.size() > 1) {
				for (int i = 1; i < userPrefs.size(); i++) {
					Long nextTime = userPrefs.get(i).v3;
					Long diff = Math.abs(nextTime - prev);
					I nextId = userPrefs.get(i).v1;
					
					if (itemTransitions.get(actualId) == null) {
						itemTransitions.put(actualId, new ArrayList<>());
					} else {
						itemTransitions.get(actualId).add(diff);
					}
					prev = nextTime;
					actualId = nextId;
				}
			}
		});
		// Now we need to compute the average
		for (I item: itemTransitions.keySet()) {
			Double avg = itemTransitions.get(item).stream().mapToDouble(i -> i).average().getAsDouble();
			result.put(item, avg.longValue());
		}
		
		return result;
	}




}
