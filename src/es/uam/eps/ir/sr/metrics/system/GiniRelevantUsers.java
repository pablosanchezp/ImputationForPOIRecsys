package es.uam.eps.ir.sr.metrics.system;


import es.uam.eps.ir.ranksys.core.Recommendation;
import es.uam.eps.ir.ranksys.diversity.sales.metrics.GiniIndex;
import es.uam.eps.ir.ranksys.metrics.rel.BinaryRelevanceModel;

/***
 * Version of gini metric that ignores the recommendations provided for a user that do not have a relevance model.
 * @author Pablo Sanchez (psperez@icai.comillas.edu)
 *
 * @param <U>
 * @param <I>
 */
public class GiniRelevantUsers <U, I> extends GiniIndex<U, I> {

    private BinaryRelevanceModel<U, I> relevanceModel;


    /**
     * Constructor.
     *
     * @param cutoff maximum length of the recommendation lists that is evaluated
     * @param numItems total number of items in the catalog
     */
    public GiniRelevantUsers(int cutoff, int numItems, BinaryRelevanceModel<U, I> binaryRelModel) {
        super(cutoff, numItems);
        this.relevanceModel = binaryRelModel;
    }
    
    @Override
	public void add(Recommendation<U, I> recommendation) {
		if (relevanceModel.getModel(recommendation.getUser()).getRelevantItems().size() > 0) {
            super.add(recommendation);
        }	
	}
}

