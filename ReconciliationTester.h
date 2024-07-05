#pragma once

#include <vector>
#include <map>

#include "SegmentalReconciler.h"

using namespace std;



class ReconciliationTester
{
public:
	int GetDupHeightUnder(GNode* g, SNode* s, GSMap& map);
	int GetSpeciesDistance(SNode* s, SNode* t);
	map<SNode*, int> GetDupHeights(SNode* species_tree, vector<GNode*>& gene_trees, GSMap& curmap);
	bool IsDup(GNode* g, GSMap& map);
	bool IsValidMap(SNode* species_tree, vector<GNode*>& gene_trees, GSMap& map, GSMap& lcamap);
	GNode* FindGeneTreeNode(GNode* g_root, string name);
	GSMap GetLcaMap(vector<GNode*> gene_trees, SNode* species_tree);
	GSMap ReadSegrecFile(string filename);
	int GetNbLossesOnParentBranch(GNode* g, GSMap& map);
	int GetNbLosses(SNode* species_tree, vector<GNode*>& gene_trees, GSMap& curmap);
};

