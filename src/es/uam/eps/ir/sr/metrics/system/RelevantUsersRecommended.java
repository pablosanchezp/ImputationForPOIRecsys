package es.uam.eps.ir.sr.metrics.system;

import java.util.Set;
import java.util.TreeSet;

import es.uam.eps.ir.ranksys.core.Recommendation;
import es.uam.eps.ir.ranksys.metrics.SystemMetric;
import es.uam.eps.ir.ranksys.metrics.rel.BinaryRelevanceModel;
import es.uam.eps.ir.ranksys.metrics.rel.IdealRelevanceModel.UserIdealRelevanceModel;

/***
 * Metric to obtain the number of users to whom we have retrieved at least 1 relevant item in the recommendations made with a cutoff
 * 
 * @author Pablo Sanchez (psperez@icai.comillas.edu)
 *
 * @param <U>
 * @param <I>
 */
public class RelevantUsersRecommended<U, I> implements SystemMetric<U, I> {

	private final int cutoff;
	private Set<U> usersRelecantRecommended;
    private BinaryRelevanceModel<U, I> relevanceModel;


	public RelevantUsersRecommended(int cutoff, BinaryRelevanceModel<U, I> relModel) {
	    this.cutoff = cutoff;
		usersRelecantRecommended = new TreeSet<>();
		this.relevanceModel = relModel;
	}
	
	
	@Override
	public void add(Recommendation<U, I> recommendation) {
		UserIdealRelevanceModel<U, I> urel = relevanceModel.getModel(recommendation.getUser());

		boolean isARelevantItem = recommendation.getItems().stream().limit(cutoff).map(p -> p.v1).filter(i -> urel.isRelevant(i)).findAny().isPresent();
		
		if (isARelevantItem) {
			usersRelecantRecommended.add(recommendation.getUser());
		}
	}

	@Override
	public double evaluate() {
		return usersRelecantRecommended.size();
	}

	@Override
	public void combine(SystemMetric<U, I> other) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void reset() {
		usersRelecantRecommended.clear();		
	}
	
}
