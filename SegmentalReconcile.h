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


/**
 * @brief The MultiGeneReconcilerInfo class is a basic structure to hold
 * various variables related to a partial mapping.  It is mainly used to pass all these
 * variables around into a single structure.
 */
class SegmentalReconcileInfo
{
public:
    unordered_map<Node*, Node*> partialMapping;
    int nbLosses;
    int dupHeightSum;
    bool isBad;

    SegmentalReconcileInfo()
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

class Triple
{
public:

    int slbl;
    int dupHeight;
    Node* n;
};

class  SegmentalReconcile
{
public:

    /**
     * @brief MultiGeneReconciler
     * @param geneTrees The set of gene trees contained in the forest.
     * @param speciesTree The species tree.
     * @param geneSpeciesMapping A mapping from the leaves of the gene trees to the leaves of the species tree.
     * @param dupcost The cost for one level of duplication.
     * @param losscost The cost for each loss.
     * @param maxDupHeight The maximum allowable duplication height.
     */
    SegmentalReconcile(vector<Node*>& geneTrees, Node* speciesTree, unordered_map<Node*, Node*>& geneSpeciesMapping, double dupcost, double losscost, int maxDupHeight, int nbspecies, int nbgenes, int numintnodes, string algorithm);

    /**
     * @brief Reconcile
     * Performs the reconciliation.  The return value contains the mapping, the sum of duplication heights and number of losses.
     * If isBad is true in the returned info, then it means there exists no solution.
     * @return
     */
    SegmentalReconcileInfo Reconcile();
    SegmentalReconcileInfo lca_algorithm();
    SegmentalReconcileInfo simphy_algorithm();
    SegmentalReconcileInfo greedy_algorithm();
    SegmentalReconcileInfo fastgreedy_algorithm();
    SegmentalReconcileInfo ultragreedy_algorithm();
    SegmentalReconcileInfo stochastic_algorithm();
    /**
     * @brief GetMappingCost
     * @param fullMapping A mapping of each node of each gene tree to the species tree.  We assume this mapping is valid without checking.
     * @return The total segmental dup + loss cost.
     */
    double GetMappingCost(unordered_map<Node*, Node*>& fullMapping);

    int GetnbLossesexceptgenetree(unordered_map<Node*, Node*>& fullMapping, int genetreenum);

    int GetnbLosses(unordered_map<Node*, Node*>& fullMapping);

    int GetdupHeightSum(vector<hashlist>& hashtable);
    /**
     * @brief IsDuplication Returns true iff g is a duplication under partialMapping
     * @param g Internal node from some gene tree
     * @param partialMapping The current partial mapping.  Can actually be compelte.
     * @return true or false
     */
    bool IsDuplication(Node* g, unordered_map<Node*, Node*>& partialMapping);

private:

    vector<Node*> geneTrees;
    Node* speciesTree;
    unordered_map<Node*, Node*> geneSpeciesMapping;
    unordered_map<Node*, Node*> lcaMapping;
    double dupcost;
    double losscost;
    int maxDupHeight;
    int nbspecies, numintnodes;
    int nbgenes;
    int** DupChanges;
    int*** TDupChanges;
    int** LossChanges;
    int** Chain;
    bool* Visit;
    string algorithm;
    vector<hashlist> hashtable;

    //Main recursive function for the computation of a mapping.  Takes the partial mapping in info and tries to map additional
    //nodes.  The given mapping must be clean.  Returns a MultiGeneReconcilerInfo containing a complete mapping if one can be
    //found.  If not, the returned info object will have isBad = true.  The variable duplicationHeights has one entry for each species.
    SegmentalReconcileInfo ReconcileRecursive(SegmentalReconcileInfo& info, unordered_map<Node*, int>& duplicationHeights);

    //holds the current best solution, so that we can do some branch-and-bound early stop if we know we acnnot beat this in a recursion
    SegmentalReconcileInfo currentBestInfo;

    //key1 = species 1, key2 = species2, int = dist.  Not sure how well this performs, but let's try
    unordered_map< Node*, unordered_map<Node*, int> > speciesTreeDistances;

    //fills up the lcaMapping variables
    void ComputeLCAMapping();

    //true iff g is a key in partialMapping
    bool IsMapped(Node* g, unordered_map<Node*, Node*>& partialMapping);

    //returns the lsit of unmapped nodes whose two children are mapped
    vector<Node*> GetMinimalUnmappedNodes(unordered_map<Node*, Node*>& partialMapping);

    SegmentalReconcileInfo GreedyRemapping(unordered_map<Node*, Node*>& partialMapping, vector<hashlist>& hashtable, vector<Node*>& geneTrees, Node* speciesTree, double LCAcost, double dupcost, double losscost, bool* improve);

    SegmentalReconcileInfo UltraGreedyRemapping(unordered_map<Node*, Node*>& partialMapping, vector<hashlist>& hashtable, vector<Node*>& geneTrees, Node* speciesTree, double LCAcost, double dupcost, double losscost, bool *improve);
    SegmentalReconcileInfo UltraGreedyRemapping1(unordered_map<Node*, Node*>& partialMapping, vector<hashlist>& hashtable, vector<Node*>& geneTrees, Node* speciesTree, double LCAcost, double dupcost, double losscost, bool* improve);
    SegmentalReconcileInfo FastGreedyRemapping(unordered_map<Node*, Node*>& partialMapping, vector<hashlist>& hashtable, int** Chain, int*** TDupChanges,int** DupChanges, int** LossChanges, bool* Visit, vector<Node*>& geneTrees, Node* speciesTree, double LCAcost, double dupcost, double losscost, bool* improve);
    SegmentalReconcileInfo StochasticRemapping(unordered_map<Node*, Node*>& partialMapping, vector<hashlist>& hashtable, int** Chain, int*** TDupChanges, int** DupChanges, int** LossChanges, bool* Visit, vector<Node*>& geneTrees, Node* speciesTree, double LCAcost, double dupcost, double losscost, bool* improve);

    double CalculateCostChange_v2(vector<Node*> n, Node* s, unordered_map<Node*, Node*>& partialMapping, vector<hashlist>& hashtable, int** Chain, int*** TDupChanges, int** DupChanges, int** LossChanges, bool* Visit, vector<Node*>& geneTrees, Node* speciesTree, double dupcost, double losscost, bool &tuple);
    double CalculateCostChange(Node* n, Node* s, unordered_map<Node*, Node*>& partialMapping, vector<hashlist>& hashtable, int** Chain, int*** TDupChanges, int** DupChanges, int** LossChanges, vector<Node*>& geneTrees, Node* speciesTree, double dupcost, double losscost);
    int ApplyChange(Node* n, Node* s, unordered_map<Node*, Node*>& partialMapping, vector<hashlist>& hashtable);
    //returns the lowest node of the species tree on which g can be mapped, ie the lca of the mappings of the 2 children of g
    Node* GetLowestPossibleMapping(Node* g, unordered_map<Node*, Node*>& partialMapping);

    Node* GetSimphyMapping(Node* g, unordered_map<Node*, Node*>& partialMapping);

    //false iff mapping g anywhere valid makes it a duplication
    bool IsRequiredDuplication(Node* g, unordered_map<Node*, Node*>& partialMapping);

    //true iff g is unmapped but its children are
    bool IsMinimalUnmapped(Node* g, unordered_map<Node*, Node*>& partialMapping);

    //true iff mapping g to its lowest possible place does not incerase duplication heights
    bool IsEasyDuplication(Node* g, unordered_map<Node*, Node*>& partialMapping, unordered_map<Node*, int>& duplicationHeights);

    //returns the duplication height at species of the subtree rooted at g
    int GetDuplicationHeightUnder(Node* g, Node* species, unordered_map<Node*, Node*>& partialMapping);



    //returns a node in the minimal nodes whose species is the lowest possible (multiple choices are possible, returns the first)
    Node* GetLowestMinimalNode(vector<Node*>& minimalNodes, unordered_map<Node*, Node*>& partialMapping);

    //returns the list of nodes on which minimalNode can be mapped to
    vector<Node*> GetPossibleSpeciesMapping(Node* minimalNode, unordered_map<Node*, Node*>& partialMapping);


    //Returns the number of edges between x and y in the species tree.  x and y must be comparable.
    int GetSpeciesTreeDistance(Node* x, Node* y);

    //Applies the cleaning phase on a partialMapping by mapping easy nodes until none are left.
    //This mapping can undergo modifications.  Only the minimal nodes and their ancestors can be modified.
    int CleanupPartialMapping(unordered_map<Node*, Node*>& partialMapping, unordered_map<Node*, int>& duplicationheights, vector<Node*>& minimalNodes);

};



class StochasticVectors {
public:
    std::vector<int> gindex;
    std::vector<int> sindex;
    std::vector<bool> gindex_to_sindex;
    std::vector<double> prob;

public:
    // Constructor
    StochasticVectors() {}

    // Function to add elements to a specified vector
    void addElement(int gindex1, int sindex1, bool gindex_to_sindex1, double prob1) {
        gindex.push_back(gindex1);
        sindex.push_back(sindex1);
        gindex_to_sindex.push_back(gindex_to_sindex1);
        prob.push_back(prob1);
    }

    // Function to get the size of a specified vector
    int getSize(int vectorNumber) {
        return gindex.size();
    }

    // Function to display elements of a specified vector
    void displayVector(int vectorNumber) {
        switch (vectorNumber) {
        case 1:
            printVector(gindex);
            break;
        case 2:
            printVector(sindex);
            break;
        case 3:
            printVector(gindex_to_sindex);
            break;
        case 4:
            printVector(prob);
            break;
        default:
            std::cout << "Invalid vector number!" << std::endl;
            break;
        }
    }


    // Function to normalize the prob vector to values between 0 and 1
    void normalizeProbabilities() {
        if (prob.empty()) {
            return;
        }

        double min_val = *std::min_element(prob.begin(), prob.end());
        double max_val = *std::max_element(prob.begin(), prob.end());
        double range = max_val - min_val;

        if (range == 0) {
            range = 1;
        }

        for (auto& p : prob) {
            p = (p - min_val) / range;
        }
    }

    // Function to perform weighted random selection from the prob vector
    std::pair<int, double> weightedRandomSelection() {
        double totalWeight = std::accumulate(prob.begin(), prob.end(), 0.0);

        if (totalWeight == 0) {
            return { -2, -2.0 };
        }

        static std::default_random_engine generator(static_cast<unsigned>(std::time(0)));
        std::uniform_real_distribution<double> distribution(0.0, totalWeight);
        double randomNumber = distribution(generator);

        double cumulativeWeight = 0.0;
        for (size_t i = 0; i < prob.size(); ++i) {
            cumulativeWeight += prob[i];
            if (randomNumber <= cumulativeWeight) {
                return { static_cast<int>(i), prob[i] };
            }
        }

        // This should never happen if the weights sum up to a non-zero value
        return { -1, -1.0 }; // Error
    }

    void printValuesAtIndex(int i) {
        if (i < 0 || i >= gindex.size()) {
            std::cout << "Index out of bounds!" << std::endl;
            return;
        }
        std::cout << "[ " << gindex[i] << ", " << sindex[i] << ", " << gindex_to_sindex[i] << ", " << prob[i] << " ]" << std::endl;
    }

private:
    // Helper function to print a vector
    void printVector(const std::vector<int>& vec) {
        for (int element : vec) {
            std::cout << element << " ";
        }
        std::cout << std::endl;
    }
    // Helper function to print a vector
    void printVector(const std::vector<bool>& vec) {
        for (bool element : vec) {
            std::cout << element << " ";
        }
        std::cout << std::endl;
    }
    // Helper function to print a vector
    void printVector(const std::vector<double>& vec) {
        for (double element : vec) {
            std::cout << element << " ";
        }
        std::cout << std::endl;
    }
};

#endif // SEGMENTALRECONCILER_H
