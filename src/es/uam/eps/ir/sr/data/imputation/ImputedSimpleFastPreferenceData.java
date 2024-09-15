package es.uam.eps.ir.sr.data.imputation;




import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.concurrent.atomic.AtomicInteger;

import es.uam.eps.ir.ranksys.fast.preference.FastPreferenceData;

/***
 * ImputedSimpleFastPreferenceData class.
 * 
 * Basic definition of the imputation class. 
 * 
 * It contains the original preference data and a boolean indicating if the data to impute is per users (e.g. adding the users a specific number
 * of ratings or according to the users number of ratings) or if the data to impute is per item (e.g adding the items a specific number of ratings)
 * 
 * @author Pablo Sanchez (psperez@icai.comillas.edu) 
 *
 * @param <U> the users
 * @param <I> the items
 */
public abstract class ImputedSimpleFastPreferenceData<U, I> {
	
	protected FastPreferenceData<U, I> prefData;
	protected boolean perUser;
	protected boolean binaryImputation;

	
	public ImputedSimpleFastPreferenceData(FastPreferenceData<U, I> prefData, boolean perUser, boolean binaryImputation) {
		this.prefData = prefData;
		this.perUser = perUser;
		this.binaryImputation = binaryImputation;
		System.out.println("Per user: "  + this.perUser);
		System.out.println("Binary imputation: "  + this.binaryImputation);

	}
	
	/***
	 * 
	 * @param prefData the original preferenceData
	 * @param rec the recommender
	 * @param percentajeIncrement
	 * @param confidence
	 * @param destPath
	 * @return 
	 */
	public void writeimputeData(String destPath) {
		PrintStream newPreferenceData;
		try {
			newPreferenceData = new PrintStream(destPath);
			prefData.getAllUsers().forEach(u -> {
				AtomicInteger countPrefs = new AtomicInteger(0);
				
				//Write the original preferences
				prefData.getUserPreferences(u).forEach(pref -> {
					if (perUser) {
						newPreferenceData.println(u + "\t" + pref.v1 + "\t" + pref.v2);
					} else {
						newPreferenceData.println(pref.v1 + "\t" + u + "\t" + pref.v2);
					}
					countPrefs.incrementAndGet();
				});
				
			});
			
			imputeData(newPreferenceData);
			newPreferenceData.close();

		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return;
	}
	
	protected abstract void imputeData(PrintStream newPreferenceData);
	
	protected boolean printImputedPreference(PrintStream newPreferenceData, U user, I item, double rating) {
		if (perUser) {
			if (binaryImputation) {
				newPreferenceData.println(user + "\t" + item + "\t" + 1.0);
			} else {
				newPreferenceData.println(user + "\t" + item + "\t" + rating);
			}
		} else {
			if (binaryImputation) {
				newPreferenceData.println(item +  "\t" + user + "\t" + 1.0);
			} else {
				newPreferenceData.println(item +  "\t" + user + "\t" + rating);
			}
		}
		return true;
	}
}
