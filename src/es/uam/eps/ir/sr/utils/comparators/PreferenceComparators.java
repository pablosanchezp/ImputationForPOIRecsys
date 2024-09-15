package es.uam.eps.ir.sr.utils.comparators;

import java.util.Comparator;

import org.jooq.lambda.tuple.Tuple2;
import org.jooq.lambda.tuple.Tuple3;
import org.ranksys.core.util.tuples.Tuple2id;
import org.ranksys.core.util.tuples.Tuple2od;

import es.uam.eps.ir.ranksys.core.preference.IdPref;
import es.uam.eps.ir.ranksys.fast.preference.IdxPref;
import es.uam.eps.ir.crossdomainPOI.datamodel.temporal.preferences.IdTimePref;
import es.uam.eps.ir.crossdomainPOI.datamodel.temporal.preferences.IdxTimePref;

/***
 * Class that contain the preference comparators
 * 
 * @author Pablo Sanchez (psperez@icai.comillas.edu)
 *
 */
public final class PreferenceComparators {
	

	
	//Tuple2 ComparatorSecondElement
	public static<T extends Comparable<T>, K extends Comparable<K>> Comparator<Tuple2<T, K>> comparatorTuple2Second() {
		return new Comparator<Tuple2<T, K>>() {		
			@Override
			public int compare(Tuple2<T, K> o1, Tuple2<T, K> o2) {
				int res = o1.v2.compareTo(o2.v2);
				if (res ==  0)
					return o1.v1.compareTo(o2.v1);
				return res;
			}
		};
	}
	
	//Tuple3 ComparatorSecondElement
	public static<T, K extends Comparable<K>, C extends Comparable<C>> Comparator<Tuple3<T, K, C>> comparatorTuple3Third() {
			return new Comparator<Tuple3<T, K, C>>() {		
				@Override
				public int compare(Tuple3<T, K, C> o1, Tuple3<T, K, C> o2) {
					int res = o1.v3.compareTo(o2.v3);
					if (res ==  0) {
						return o1.v2.compareTo(o2.v2);
					}
					return res;
			}
		};
	}	
	
	



	

	
	//Comparator of preferences (IdxPref)
	public static Comparator<IdxPref> preferenceComparatorIdxPref = (o1, o2) -> o1.v2 != o2.v2 ? Double.compare(o1.v2, o2.v2): Integer.compare(o1.v1, o2.v1);
	
	//Comparator of indexes (IdxPref)
	public static Comparator<IdxPref> indexComparatorIdxPref = (o1, o2) -> Integer.compare(o1.v1, o2.v1);
	
	///
	//Comparator of preferences (IdxTimePref)
	public static Comparator<IdxTimePref> preferenceComparatorIdxTimePref = (o1, o2) ->  o1.v2 != o2.v2 ? Double.compare(o1.v2, o2.v2): Integer.compare(o1.v1, o2.v1);

		
	//Comparator of indexes (IdxTimePref)
	public static Comparator<IdxTimePref> indexComparatorIdxTimePref = (o1, o2) -> Integer.compare(o1.v1, o2.v1);
	
	//Comparator of time (IdxTimePref)
	public static Comparator<IdxTimePref> timeComparatorIdxTimePref = (o1, o2) -> o1.v3 != o2.v3 ? Long.compare(o1.v3, o2.v3): Integer.compare(o1.v1, o2.v1);
	
	
	//Id and IdTime pref	
	
	//Comparator of preferences (IdPref)
	public static Comparator<IdPref<Long>> preferenceComparatorIdPref = (o1, o2) -> o1.v2 != o2.v2 ?  Double.compare(o1.v2, o2.v2): (o1.v1).compareTo(o2.v1);
	
	
	//Comparator of indexes (IdPref)
	public static Comparator<IdPref<Long>> idComparatorIdPref = (o1, o2) -> (o1.v1).compareTo(o2.v1);
	
	///
	//Comparator of preferences (IdTimePref)
	public static Comparator<IdTimePref<Long>> preferenceComparatorIdTimePref = (o1, o2) -> o1.v2 != o2.v2 ?  Double.compare(o1.v2, o2.v2): (o1.v1).compareTo(o2.v1);
	public static<T extends Comparable<T>> Comparator<IdTimePref<T>> preferenceComparatorIdTimePref () {
		return new Comparator<IdTimePref<T>>() {
			
			@Override
			public int compare(IdTimePref<T> o1, IdTimePref<T> o2) {
				return o1.v2 != o2.v2 ?  Double.compare(o1.v2, o2.v2): (o1.v1).compareTo(o2.v1);
			}
		};
	}
	
	public static<T extends Comparable<T>, U extends Comparable<U>> Comparator<Tuple2<T, U>> Tuple2ComparatorSecondFirstElement () {
		
		return new Comparator<Tuple2<T, U>>() {
		    @Override
		    public int compare(Tuple2<T, U> t1, Tuple2<T, U> t2) {
		    	int cmp = t1.v2.compareTo(t2.v2());	
			    if (cmp != 0) {
			        return cmp;
			    } else {
			        return t1.v1.compareTo(t2.v1);
			    }
		    }
		};
	}
	
	public static<T extends Comparable<T>> Comparator<IdPref<T>> preferenceComparatorIdPref () {
		return new Comparator<IdPref<T>>() {			
			@Override
			public int compare(IdPref<T> o1, IdPref<T> o2) {
				return o1.v2 != o2.v2 ?  Double.compare(o1.v2, o2.v2): (o1.v1).compareTo(o2.v1);
			}
		};
	}
	
	
	
	
	//Comparator of indexes (IdxTimePref)
	public static Comparator<IdTimePref<Long>> idComparatorIdTimePref = (o1, o2) -> (o1.v1).compareTo(o2.v1);
	public static<T extends Comparable<T>> Comparator<IdTimePref<T>> idComparatorIdTimePref () {
		return new Comparator<IdTimePref<T>>() {		
			@Override
			public int compare(IdTimePref<T> o1, IdTimePref<T> o2) {
				return o1.v1.compareTo(o2.v1);
			}
		};
	}

	//Comparator of time (IdxTimePref)
	public static Comparator<IdTimePref<Long>> timeComparatorIdTimePref = (o1, o2) -> o1.v3 != o2.v3 ? Long.compare(o1.v3, o2.v3): (o1.v1).compareTo(o2.v1);
	public static<T extends Comparable<T>> Comparator<IdTimePref<T>> timeComparatorIdTimePref () {
		return new Comparator<IdTimePref<T>>() {
			
			@Override
			public int compare(IdTimePref<T> o1, IdTimePref<T> o2) {
				return o1.v3 != o2.v3 ? Long.compare(o1.v3, o2.v3): (o1.v1).compareTo(o2.v1);
			}
		};
	}
	

	
	//////////////////////////////////
	//Other comparators
	/////////////////////////////////
	public static Comparator<Tuple2od<Comparable<Object>>> recommendationComparatorTuple2od = (o1,o2) -> o1.v2 != o2.v2 ? Double.compare(o1.v2, o2.v2): (o1.v1).compareTo(o2.v1);;
	public static<T extends Comparable<T>> Comparator<Tuple2od<T>> recommendationComparatorTuple2od () {
		return new Comparator<Tuple2od<T>>() {
			
			@Override
			public int compare(Tuple2od<T> o1, Tuple2od<T> o2) {
				return o1.v2 != o2.v2 ? Double.compare(o1.v2, o2.v2): ((Comparable<T>) o1.v1).compareTo(o2.v1);
			}
		};
	}
	

	
	public static Comparator<Tuple2id> recommendationComparatorTuple2id = (o1,o2) -> o1.v2 != o2.v2 ? Double.compare(o1.v2, o2.v2): Integer.compare(o1.v1,o2.v1);
	
	
	public static Comparator<Object> comparatorStringLong = (l1, l2) -> l1.toString().compareTo(l2.toString());

	
	public static<K, T extends Comparable<T>> Comparator<Tuple2<K, T>> genericTupleComparator () {
		return new Comparator<Tuple2<K, T>>() {		
			@Override
			public int compare(Tuple2<K, T> o1, Tuple2<K, T> o2) {
				return o1.v2.compareTo(o2.v2);
			}
		};
	}
	
	public static<K, T extends Comparable<T>> Comparator<Tuple3<K, T, T>> genericTuple3Comparator () {
		return new Comparator<Tuple3<K, T, T>>() {		
			@Override
			public int compare(Tuple3<K, T, T> o1, Tuple3<K, T, T> o2) {
				int val = o1.v2.compareTo(o2.v2);
				if (val == 0) {
					return o1.v3.compareTo(o2.v3);
				} else {
					return val;
				}
			}
		};
	}
	
	
}
