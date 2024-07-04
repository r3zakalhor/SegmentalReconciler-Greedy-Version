#ifndef SEGMENTALRECONCILER_H
#define SEGMENTALRECONCILER_H

#include <iostream>

#include <map>
#include "util.h"
#include "newicklex.h"
#include "node.h"
#include "genespeciestreeutil.h"
#include "treeiterator.h"
#include "hashtable.h"

#include <vector>
#include <algorithm>
#include <numeric> // For std::accumulate
#include <random>
#include <ctime>
#include <utility>

using namespace std;




#define GNode Node
#define SNode Node
#define GSMap unordered_map<GNode*, SNode*>







enum MoveType { UpMove, DownMove, BulkUpMove, BulkDownMove };

/**
Represent a possible remapping move.  Depending on the type of move, different variables will be set.  
In the case of an up move, g and s will be set (meaning remap g to the species s). 
In the case of a bulk move, s_src and s_dest will be set, indicating that all 
top dups in s_src should be moved to s_dest.  In that case, delta_dup is the max increase in height of a gene.
Polymorphism might be more appropriate here, please formulate complaints into a text file and store it in /dev/null
**/
struct RemapMove{
	
	MoveType type;
	
	int delta_dup;
	int delta_loss;
	double delta_cost;
	double prob;
	
	GNode* g;
	SNode* s;
	
	SNode* s_src;
	SNode* s_dest;
	
	
	RemapMove(){
		g = nullptr;
		s = nullptr;
		s_src = nullptr;
		s_dest = nullptr;
		delta_dup = -1;
	}
	
	
	RemapMove(MoveType type, int delta_dup, int delta_loss, double delta_cost) : RemapMove(){
		this->type = type;
		this->delta_dup = delta_dup;
		this->delta_loss = delta_loss;
		this->delta_cost = delta_cost;
	}
	
	
	
	//constructor intended for a standard up/down move
	RemapMove(MoveType type, int delta_dup, int delta_loss, double delta_cost, GNode* g, SNode* s) : RemapMove(){
		this->type = type;
		this->delta_dup = delta_dup;
		this->delta_loss = delta_loss;
		this->delta_cost = delta_cost;
		this->g = g;
		this->s = s;
	}


	void debug_print() {
		if (type == UpMove) {
			std::cout << "Upmove: Remap g=" << g->GetLabel() << " to s=" << s->GetLabel() 
					  << "  delta_dup=" << this->delta_dup << "  delta_loss=" << delta_loss << "  delta_cost=" << delta_cost << endl;
		}
		else if (type == DownMove) {
			std::cout << "Downmove: Remap g=" << g->GetLabel() << " to s=" << s->GetLabel()
				<< "  delta_dup=" << this->delta_dup << "  delta_loss=" << delta_loss << "  delta_cost=" << delta_cost << endl;
		}
		else if (type == BulkUpMove) {
			std::cout << "Bulk up move: Remap s_src=" << s_src->GetLabel() << " to s_dest=" << s_dest->GetLabel()
					  << "  delta_dup=" << delta_dup << "  delta_loss=" << delta_loss << "  delta_cost=" << delta_cost << endl;
		}
		else if (type == BulkDownMove) {
			std::cout << "Bulk down move: Remap s_src=" << s_src->GetLabel() << " to s_dest=" << s_dest->GetLabel()
				<< "  delta_dup=" << delta_dup << "  delta_loss=" << delta_loss << "  delta_cost=" << delta_cost << endl;
		}
	}
	
};



/**
Stores a list of possible moves.  Maintains the move with minimum delta_cost
**/
class RemapMoveList{

public:

	size_t best_move_index = -1;
	size_t random_move_index = -1;

	double stochastic_temperature;	//initialized in constructor

	RemapMoveList() {
		best_move_index = -1;
		random_move_index = -1;
		stochastic_temperature = 1;
	}
	
	void Reset(){
		moves.clear();
		best_move_index = -1;
		random_move_index = -1;
	}
	
	
	void AddUpMove(int delta_dup, int delta_loss, double delta_cost, Node* g, Node* s)
	{
		RemapMove mv(UpMove, delta_dup, delta_loss, delta_cost, g, s);
		mv.prob = std::exp(- delta_cost / (stochastic_temperature)); // boltzman_distribution
		moves.push_back(mv);
		
		//update best move index
		if (best_move_index == -1 || delta_cost <= moves[best_move_index].delta_cost) {
			best_move_index = moves.size() - 1;
		}
	}
	
	
	
	void AddBulkUpMove(int delta_dup, int delta_loss, double delta_cost, Node* s_src, Node* s_dest)
	{
		RemapMove mv(BulkUpMove, delta_dup, delta_loss, delta_cost);
		mv.s_src = s_src;
		mv.s_dest = s_dest;
		mv.prob = std::exp(- delta_cost / (stochastic_temperature)); // boltzman_distribution
		moves.push_back(mv);
		
		
		//update best move index
		if (best_move_index == -1 || delta_cost <= moves[best_move_index].delta_cost){
			best_move_index = moves.size() - 1;
		}
	}
	
	
	
	
	void AddDownMove(int delta_dup, int delta_loss, double delta_cost, Node* g, Node* s)
	{
		RemapMove mv(DownMove, delta_dup, delta_loss, delta_cost, g, s);
		mv.prob = std::exp(-delta_cost / (stochastic_temperature)); // boltzman_distribution
		moves.push_back(mv);

		//update best move index
		if (best_move_index == -1 || delta_cost <= moves[best_move_index].delta_cost) {
			best_move_index = moves.size() - 1;
		}
	}



	void AddBulkDownMove(int delta_dup, int delta_loss, double delta_cost, Node* s_src, Node* s_dest)
	{
		RemapMove mv(BulkDownMove, delta_dup, delta_loss, delta_cost);
		mv.prob = std::exp(-delta_cost / (stochastic_temperature)); // boltzman_distribution
		mv.s_src = s_src;
		mv.s_dest = s_dest;
		moves.push_back(mv);


		//update best move index
		if (best_move_index == -1 || delta_cost <= moves[best_move_index].delta_cost) {
			best_move_index = moves.size() - 1;
		}
	}
	
	
	// Function to normalize the prob values to values between 0 and 1
    void NormalizeProbabilities() {
        if (moves.empty()) {
            return;
        }
		
		auto min_it = std::min_element(moves.begin(), moves.end(), 
                                   [](const RemapMove& a, const RemapMove& b) {
                                       return a.prob < b.prob;
                                   });
		auto max_it = std::max_element(moves.begin(), moves.end(), 
                                   [](const RemapMove& a, const RemapMove& b) {
                                       return a.prob < b.prob;
                                   });
								  
        double min_val = min_it->prob;
        double max_val = max_it->prob;
        double range = max_val - min_val;
		
		if (range == 0) {
			// If all probabilities are the same, normalize them to 0.5 (or any consistent value)
			for (auto& m : moves) {
				m.prob = 0.5;
			}
		} else {
			for (auto& m : moves) {
				m.prob = (m.prob - min_val) / range;
			}
		}
    }
	
	
	// Function to perform weighted random selection from the moves vector based on prob values
    size_t WeightedRandomSelection() {
        double totalWeight = std::accumulate(moves.begin(), moves.end(), 0.0, 
                           [](double sum, const RemapMove& move) {
                               return sum + move.prob;
                           });

        if (totalWeight == 0) {
            return -1;
        }

        static std::default_random_engine generator(static_cast<unsigned>(std::time(0)));
        std::uniform_real_distribution<double> distribution(0.0, totalWeight);
        double randomNumber = distribution(generator);

        double cumulativeWeight = 0.0;
        for (size_t i = 0; i < moves.size(); ++i) {
            cumulativeWeight += moves[i].prob;
            if (randomNumber <= cumulativeWeight) {
                return i;
            }
        }

        // This should never happen if the weights sum up to a non-zero value
        return -1; // Error
    }

	
	size_t GetNbMoves()
	{
		return moves.size();
	}
	
	RemapMove GetMove(size_t i)
	{
		return moves[i];
	}
	
	size_t GetBestMoveIndex()
	{
		return best_move_index;
	}
	
	size_t GetRandomMoveIndex()
	{
		NormalizeProbabilities();
		random_move_index = WeightedRandomSelection();
		return random_move_index;
	}
	
	
private:
	
	vector<RemapMove> moves;

	
};



/**
Utility class that stores info to maintain through the execution of the algorithm.  This includes the current gene-species map (curmap), 
and the hashtable data structure.  
This also includes the tables that are needed by the dynamic programming algorithm.  The intent is to pass a RemapperInfo object around 
for the various remapping calculations and updates.
Once an iteration is done and a remap move is applied, Reset() should be called.
**/
class RemapperInfo{
public:
	//these are maintained across all iterations
	GSMap curmap;
	vector<hashlist> hashtable;
	double dupcost; 
	double losscost;
	
	//these are for stochastic version to keep the best mapping
	GSMap bestmap;
	double bestcost;
	int bestdupHeightSum;
	int bestNblosses;
	
	
	//these are reset every iteration
	unordered_map<GNode*, unordered_map<SNode*, int> > chain_lengths;
	RemapMoveList moves;
	map< GNode*, map< SNode*, map< SNode*, int > > > delta_gst;	//delta_gst[g, s, t] = change in dup height of t if we remap to s (or vice-versa?)
	map< GNode*, map< SNode*, int > > delta_gs;	//delta_gs[g, s] = change in all dup heights if we remap to s
	map< GNode*, map< SNode*, int > > lambda_gs;	//lambda_gs[g, s] = change in all losses if we remap to s
	
	
	void Reset()
	{
		chain_lengths.clear();	//needed?  maybe not, could be useful to not clear and reuse memory
		moves.Reset();
		delta_gst.clear();	//needed?
		delta_gs.clear();	//needed?
		lambda_gs.clear();	//needed?
	}
	
	void UpdateBestMapping(double cost, int dupHeightSum, int nbLosses){
		bestmap = curmap;
		bestcost = cost;
		bestdupHeightSum = dupHeightSum;
		bestNblosses = nbLosses;
	}
	
	RemapMove GetBestMove()
	{
		return moves.GetMove( moves.GetBestMoveIndex() );
	}
	
	RemapMove GetRandomMove()
	{
		return moves.GetMove( moves.GetRandomMoveIndex() );
	}
	
	
};




//TODO: maybe merge with RemapperInfo?
class SegmentalReconcilerInfo
{
public:
    GSMap curmap;
    int nbLosses;
    int dupHeightSum;
    bool isBad;

    SegmentalReconcilerInfo()
    {
        isBad = false;
        nbLosses = 0;
        dupHeightSum = 0;
    }

    double GetCost(double dupcost, double losscost)
    {
        return dupcost * (double)dupHeightSum + losscost * (double)nbLosses;
    }
};



class SegmentalReconciler
{
public:

    SegmentalReconciler(vector<Node*>& geneTrees, Node* speciesTree, GSMap& leaf_species_map, double dupcost, double losscost, int nbspecies, int nbgenes, int numintnodes, bool stochastic, double stochastic_temperature, int nbstochasticLoops);

    /**
     * @brief Reconcile
     * Performs the reconciliation.  The return value contains the mapping, the sum of duplication heights and number of losses.
     * If isBad is true in the returned info, then it means there exists no solution.
     * @return
     */
    SegmentalReconcilerInfo Reconcile();
	SegmentalReconcilerInfo ReconcileWithLCAMap();
    SegmentalReconcilerInfo ReconcileWithSimphy();

	int GetExhaustiveDupHeight(GSMap& curmap);

    /**
     * @brief GetMappingCost
     * @param fullMapping A mapping of each node of each gene tree to the species tree.  We assume this mapping is valid without checking.
     * @return The total segmental dup + loss cost.
     */
    double GetMappingCost(unordered_map<Node*, Node*>& fullMapping);


    int GetNbLosses(GSMap& fullMapping);

    int GetDupHeightSum(vector<hashlist>& hashtable);
    
	
	/**
     * @brief IsDuplication Returns true iff g is a duplication under partialMapping
     * @param g Internal node from some gene tree
     * @param partialMapping The current partial mapping.  Can actually be complete.
     * @return true or false
     */
    bool IsDuplication(Node* g, unordered_map<Node*, Node*>& partialMapping);

private:

    vector<Node*> geneTrees;
    Node* speciesTree;
    unordered_map<Node*, Node*> leaf_species_map;
    double dupcost;
    double losscost;
    int nbspecies, numintnodes;
    int nbgenes;
	bool stochastic;
	double stochastic_temperature;
	int nbstochasticLoops;
	
	unordered_map< Node*, unordered_map<Node*, int> > speciesTreeDistances;
    


    //holds the current best solution, so that we can do some branch-and-bound early stop if we know we acnnot beat this in a recursion
    SegmentalReconcilerInfo currentBestInfo;

	//idk wtf this is
	int Getspecieslbl(Node* s, int numintnodes);
    
    GSMap GetLCAMapping();

    //true iff g is a key in partialMapping
    bool IsMapped(Node* g, GSMap& partialMapping);

    
    Node* GetSimphyMapping(Node* g, unordered_map<Node*, Node*>& partialMapping);

    //returns the duplication height at species of the subtree rooted at g
    int GetDuplicationHeightUnder(Node* g, Node* species, GSMap& curmap);

    //returns the list of nodes on which minimalNode can be mapped to
    vector<SNode*> GetPossibleUpRemaps(SNode* spnode);
	
	void ComputeUpMove(GNode* g, SNode* s, RemapperInfo &remapinfo);
	bool ComputeBulkUpMove(SNode* s_src, SNode* s_dest, RemapperInfo &remapinfo);
	
	void ApplyUpMove(GNode* g, SNode* s, RemapperInfo &remapinfo);
	void ApplyBulkUpMove(SNode* s_src, SNode* s_dest, RemapperInfo &remapinfo);
	

    //Returns the number of edges between x and y in the species tree.  x and y must be comparable.
    int GetSpeciesTreeDistance(Node* x, Node* y);

	RemapperInfo BuildRemapperInfo();

	void ComputeAllUpMoves(RemapperInfo& remapinfo);

	void ComputeAllBulkUpMoves(RemapperInfo& remapinfo);
	
	void ComputeAllDownMoves(RemapperInfo& remapinfo);

	vector<SNode*> GetPossibleDownRemaps(GNode* gnode, GSMap& gsmap);
	void ApplyDownMove(GNode* g, SNode* s, RemapperInfo& remapinfo);
	bool ComputeBulkDownMove(SNode* s_src, SNode* s_dest, vector<GNode*>& deepest_nodes, RemapperInfo& remapinfo);
	void ComputeAllBulkDownMoves(RemapperInfo& remapinfo);
	void ApplyBulkDownMove(SNode* s_src, SNode* s_dest, RemapperInfo& remapinfo);

	bool ComputeDownMove(GNode* g, SNode* s, RemapperInfo& remapinfo);

	int GetNbLossesOnParentBranch(GNode* g_child, GSMap& curmap);
	int GetDupHeight(SNode* s, RemapperInfo& remapinfo);

	//returns gene tree nodes that are dups mapped to s, at depth target_depth (no return, fills the vector instead)
	void GetDupsAtDepth(GNode* curnode, SNode* s, GSMap& curmap, int target_depth, int cur_depth, vector<GNode*>& vec_to_fill);
};





#endif // SEGMENTALRECONCILER_H
