package es.uam.eps.ir.sr.metrics.system;


import es.uam.eps.ir.crossdomainPOI.metrics.system.RealAggregateDiversity;
import es.uam.eps.ir.ranksys.core.Recommendation;
import es.uam.eps.ir.ranksys.metrics.rel.BinaryRelevanceModel;

/***
 * Real aggregate diversity that only takes into account the users that have something relevant
 * @author Pablo Sanchez (pablo.sanchezp@uam.es)
 *
 * @param <U>
 * @param <I>
 */
public class RealAggregateDiversityRelevantUsers<U, I> extends RealAggregateDiversity<U, I> {

    private BinaryRelevanceModel<U, I> relevanceModel;

	public RealAggregateDiversityRelevantUsers(int cutOff, BinaryRelevanceModel<U, I> relevanceModel) {
		super(cutOff);
		this.relevanceModel = relevanceModel;
	}

	@Override
	public void add(Recommendation<U, I> recommendation) {
		if (relevanceModel.getModel(recommendation.getUser()).getRelevantItems().size() > 0) {
			 super.add(recommendation);
        }	
	}

}
