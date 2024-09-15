package es.uam.eps.ir.sr.metrics.system;

import java.util.Set;
import java.util.TreeSet;

import es.uam.eps.ir.ranksys.core.Recommendation;
import es.uam.eps.ir.ranksys.metrics.SystemMetric;
import es.uam.eps.ir.ranksys.metrics.rel.BinaryRelevanceModel;

/***
 * Metric to analyze the user coverage of relevant users (users that have at leats 1 relevant item in the test set according to a relevance model)
 * 
 * @author Pablo Sanchez (psperez@icai.comillas.edu)
 *
 * @param <U>
 * @param <I>
 */
public class UserCoverageRelevant<U, I> implements SystemMetric<U, I>{

    private Set<U> users;
    private BinaryRelevanceModel<U, I> relevanceModel;
    
    public UserCoverageRelevant(BinaryRelevanceModel<U, I> relModel) {
    	users = new TreeSet<>();
    	relevanceModel= relModel;
	}
    
	
	@Override
	public void add(Recommendation<U, I> recommendation) {
		if (relevanceModel.getModel(recommendation.getUser()).getRelevantItems().size() > 0) {
            users.add(recommendation.getUser());
        }	
	}

	@Override
	public double evaluate() {
        return users.size();
	}

	@Override
	public void combine(SystemMetric<U, I> other) {		
	}

	@Override
	public void reset() {
        users.clear();		
	}
}