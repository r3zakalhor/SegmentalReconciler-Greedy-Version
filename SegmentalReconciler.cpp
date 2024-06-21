#include "SegmentalReconciler.h"
#include <random>
#include <ctime>

SegmentalReconciler::SegmentalReconciler(vector<Node*>& geneTrees, Node* speciesTree, GSMap& leaf_species_map, double dupcost, double losscost, int nbspecies, int nbgenes, int numintnodes, bool stochastic)
{
    this->geneTrees = geneTrees;
    this->speciesTree = speciesTree;
    this->leaf_species_map = leaf_species_map;
    this->dupcost = dupcost;
    this->losscost = losscost;
    
	this->nbspecies = nbspecies;
    this->numintnodes = numintnodes;
    this->nbgenes = nbgenes;
    this->stochastic = stochastic;
    
}



vector<SNode*> SegmentalReconciler::GetPossibleUpRemaps(SNode* spnode)
{
	vector<SNode*> sp;
	
	while (!spnode->IsRoot())
	{
		spnode = spnode->GetParent();
		sp.push_back(spnode);
	}
	
	return sp;
}





/**
This function has two side effects: 
- it assumes that remapinfo stuff has been calculated at the parent, and will add entries to remapinfo for n
- will add a RemapMove to the "remapinfo.moves" object
**/
void SegmentalReconciler::ComputeUpMove(GNode* g, SNode* s, RemapperInfo &remapinfo) {

	double delta_dup = 0;
    int slbl = Util::ToInt(s->GetLabel());
	int slbl_parent = slbl;
	if (!g->IsRoot())
        slbl_parent = Util::ToInt(remapinfo.curmap[g->GetParent()]->GetLabel());

	
    SNode* mu = remapinfo.curmap[g];
	int mu_specieslbl = Getspecieslbl(mu, numintnodes);

    SNode* mpu = mu;
    if (!g->IsRoot()){
        mpu = remapinfo.curmap[g->GetParent()];
	}
    

	
    if (!g->IsRoot() && mu == mpu)
        remapinfo.chain_lengths[g][mu] = remapinfo.chain_lengths[g->GetParent()][mu] + 1;
    else
        remapinfo.chain_lengths[g][mu] = 1;
    

    if (g->IsRoot() || slbl < slbl_parent)
        remapinfo.chain_lengths[g][s] = 1;
    else 
        remapinfo.chain_lengths[g][s] = remapinfo.chain_lengths[g->GetParent()][s] + 1;
	
    // dup changes, implements recurrences from paper
    if (g->IsRoot() || slbl < slbl_parent) {
	    
        remapinfo.delta_gst[g][s][s] = std::max(remapinfo.chain_lengths[g][s], remapinfo.hashtable[slbl].size()) - remapinfo.hashtable[slbl].size();
        if (remapinfo.hashtable[mu_specieslbl].is_unique_max(g))
            remapinfo.delta_gst[g][s][mu] = -1;
        else
            remapinfo.delta_gst[g][s][mu] = 0;
        remapinfo.delta_gs[g][s] = remapinfo.delta_gst[g][s][s] + remapinfo.delta_gst[g][s][mu];
    }
    else {
        
        remapinfo.delta_gst[g][s][s] = std::max(remapinfo.chain_lengths[g][s], remapinfo.hashtable[slbl].size()) - remapinfo.hashtable[slbl].size();

		
		if (Util::ToInt(mu->GetLabel()) < Util::ToInt(mpu->GetLabel()))
            remapinfo.delta_gst[g->GetParent()][s][mu] = 0;			//WHY?  This sets something of the parent
        else	//WHY?  mu = mpu here, no?
            remapinfo.delta_gst[g->GetParent()][s][mu] = remapinfo.delta_gst[g->GetParent()][s][mpu];
		
        if (mu == mpu) {
            int hprim = remapinfo.hashtable[mu_specieslbl].size() + remapinfo.delta_gst[g->GetParent()][s][mu];
            int hmu = remapinfo.hashtable[mu_specieslbl].dupheight_at_subtree(g);
            int tmp = 0;
            if (hprim == hmu && remapinfo.hashtable[mu_specieslbl].is_unique_element(g, hprim))
                tmp = -1;
            remapinfo.delta_gst[g][s][mu] = remapinfo.delta_gst[g->GetParent()][s][mu] + tmp;
        }
        else {
            if (remapinfo.hashtable[mu_specieslbl].is_unique_max(g))
                remapinfo.delta_gst[g][s][mu] = -1;
            else
                remapinfo.delta_gst[g][s][mu] = 0;
        }
		
        if (s == mpu) {
            remapinfo.delta_gs[g->GetParent()][s] = 0;	//WHY?
            remapinfo.delta_gst[g->GetParent()][s][s] = 0;
            remapinfo.delta_gst[g->GetParent()][s][mu] = 0;
        }
        remapinfo.delta_gs[g][s] = remapinfo.delta_gs[g->GetParent()][s] - remapinfo.delta_gst[g->GetParent()][s][s] - remapinfo.delta_gst[g->GetParent()][s][mu] +
									remapinfo.delta_gst[g][s][s] + remapinfo.delta_gst[g][s][mu];

    }

	
    // loss changes
    int imu = 0;
    int impu = 0;
    int im_us_pu = 0;

    if (IsDuplication(g, remapinfo.curmap))
        imu = 0;
    else
        imu = 2;

    if (!g->IsRoot()) {
        if (IsDuplication(g->GetParent(), remapinfo.curmap))
            impu = 0;
        else
            impu = 2;

		//WHY?  tmporary remap
        remapinfo.curmap[g] = s;
        if (IsDuplication(g->GetParent(), remapinfo.curmap))
            im_us_pu = 0;
        else
            im_us_pu = 2;
        remapinfo.curmap[g] = mu;
    }

    int lambda = GetSpeciesTreeDistance(s, remapinfo.curmap[g->GetChild(0)]) + GetSpeciesTreeDistance(s, remapinfo.curmap[g->GetChild(1)]) -
				 GetSpeciesTreeDistance(mu, remapinfo.curmap[g->GetChild(0)]) - GetSpeciesTreeDistance(mu, remapinfo.curmap[g->GetChild(1)]) + imu;
    if (g->IsRoot()) {
        remapinfo.lambda_gs[g][s] = lambda;
    }
    else if (slbl <= slbl_parent) {
        remapinfo.lambda_gs[g][s] = lambda + GetSpeciesTreeDistance(s, mpu) - im_us_pu - GetSpeciesTreeDistance(mu, mpu) + impu;
    }
    else if (slbl > slbl_parent) {	//WHY not just else?
        remapinfo.lambda_gs[g][s] = remapinfo.lambda_gs[g->GetParent()][s] + lambda - GetSpeciesTreeDistance(s, mu);
    }


    double costchange = remapinfo.lambda_gs[g][s] * remapinfo.losscost + remapinfo.delta_gs[g][s] * remapinfo.dupcost;
	
	remapinfo.moves.AddUpMove(remapinfo.delta_gs[g][s], remapinfo.lambda_gs[g][s], costchange, g, s);

}



bool SegmentalReconciler::ComputeBulkMove(SNode* s_src, SNode* s_dest, RemapperInfo &remapinfo){
	
	unordered_set<GNode*> gnodes = remapinfo.hashtable[Getspecieslbl(s_src, numintnodes)].return_max_heights();
	
    if (gnodes.size() == 0) {
        cout << "No gene mapped to src" << endl;
        return false;
    }

    if (gnodes.size() == 1) {   //bulk move must include two genes
        return false;
    }

	int nb_losses = 0;
    int max_dup_height_increase = -999999;
	
	for (GNode* g : gnodes){
		nb_losses += remapinfo.lambda_gs[g][s_dest];
		
		/**
		That line below calculates the total change in dup height if we remap g to s_dest, 
		WITH the assumption that all top dups mapped to s_src get remapped at once.  
		This means that the height in s_src is reduced by exactly 1.
		Sincere mapinfo.delta_gs[g][s_dest] includes the term remapinfo.delta_gst[g][s_dest][s_src], we subtract it,
		and add "-1", which is the new change in height under our assumption.
		**/
		int dup_change = remapinfo.delta_gs[g][s_dest] - remapinfo.delta_gst[g][s_dest][s_src] - 1;
		
		max_dup_height_increase = max((int)max_dup_height_increase, (int)dup_change);
	}
	
	
	double costchange = nb_losses * remapinfo.losscost + max_dup_height_increase * remapinfo.dupcost;
	
	remapinfo.moves.AddBulkMove(max_dup_height_increase, nb_losses, costchange, s_src, s_dest);
	
    return true;
}





void SegmentalReconciler::ApplyBulkMove(SNode* s_src, SNode* s_dest, RemapperInfo &remapinfo) {
	
	unordered_set<GNode*> gnodes = remapinfo.hashtable[Getspecieslbl(s_src, numintnodes)].return_max_heights();
	
    if (gnodes.size() == 0) {
        cout << "Error: no gene is a dup in " << s_src->GetLabel() << endl;
    }

	for (Node* g : gnodes){
		ApplyUpMove(g, s_dest, remapinfo);
	}
}






void SegmentalReconciler::ApplyUpMove(GNode* g, SNode* s, RemapperInfo &remapinfo) {
	
	GNode* g_cur = g;
	vector<GNode*> nodes_to_remap;
	
	bool done = false;
	while (!done){
		
		nodes_to_remap.push_back(g_cur);
        
		if (g_cur->IsRoot())
			done = true;
		else{
			g_cur = g_cur->GetParent();
			
			//we stop when g_cur is mapped to s or above, since the chain reaction of remaps can stop
			if (Util::ToInt(remapinfo.curmap[g_cur]->GetLabel()) > Util::ToInt(s->GetLabel()))
				done = true;
		}
	}


    //important for hashtable data structure to do it in reverse
    for (int j = nodes_to_remap.size() - 1; j >= 0; j--) {
        GNode* g_cur = nodes_to_remap[j];
        
        if (IsDuplication(g_cur, remapinfo.curmap)) {
            SNode* s_cur = remapinfo.curmap[g_cur];
            int slbl = Getspecieslbl(s_cur, numintnodes);
            int dupheight = GetDuplicationHeightUnder(g_cur, s_cur, remapinfo.curmap);
            remapinfo.hashtable[slbl].remove(dupheight, g_cur);
        }

    }
	
	
	//now we do the actual remap, making sure that we update the hash table
	for (GNode* g : nodes_to_remap){
		remapinfo.curmap[g] = s;

		int dupheight = GetDuplicationHeightUnder(g, s, remapinfo.curmap);
		if (dupheight > 0) {
			int slbl = Getspecieslbl(s, numintnodes);
			remapinfo.hashtable[slbl].add_cell(dupheight, g);
		}
	}
}





RemapperInfo SegmentalReconciler::BuildRemapperInfo(){
	
	RemapperInfo remapinfo;
	remapinfo.curmap = GetLCAMapping();	//maybe slow
	remapinfo.dupcost = this->dupcost;
	remapinfo.losscost = this->losscost;
	
	//build initial hashtable
	remapinfo.hashtable.resize(this->nbspecies);
	for (int i = 0; i < geneTrees.size(); i++)
    {
        Node* g_root = geneTrees[i];

        TreeIterator* it = g_root->GetPostOrderIterator();
        while (GNode* g = it->next())
        {
			SNode* s = remapinfo.curmap[g];
            int dupheight = GetDuplicationHeightUnder(g, s, remapinfo.curmap);
            
            if (dupheight > 0) {
                
                int slbl = Getspecieslbl(s, numintnodes);
                
                remapinfo.hashtable[slbl].add_cell(dupheight, g);
            }
        }
        g_root->CloseIterator(it);
    }

    cout << "dh=" << GetDupHeightSum(remapinfo.hashtable) << endl;
	
	return remapinfo;	
}





void SegmentalReconciler::ComputeAllUpMoves(RemapperInfo& remapinfo) {
    //compute up moves for every gene tree node
    for (int i = 0; i < geneTrees.size(); i++) {
        GNode* g_root = geneTrees[i];

        TreeIterator* it = g_root->GetPreOrderIterator();	//needs to be preorder for DP to work
        while (GNode* g = it->next()) {
            if (!g->IsLeaf()) {
                SNode* mu = remapinfo.curmap[g];
                vector<SNode*> possibleremapping = GetPossibleUpRemaps(mu);

                for (SNode* s : possibleremapping) {
                    ComputeUpMove(g, s, remapinfo);

                    /*myfile << "gene index: " << g->GetLabel() << " from " << mu->GetLabel() << " -> " << s->GetLabel()
                        << " Dup change: " << Util::ToString(remapinfo.delta_gs[g][s]) << " loss change: "
                        << Util::ToString(remapinfo.lambda_gs[g][s])
                        << " cost change: " << remapinfo.moves.GetMove(remapinfo.moves.GetNbMoves()-1).delta_cost << endl;*/

                }
            }
        }
        g_root->CloseIterator(it);
    }
}




void SegmentalReconciler::ComputeAllBulkMoves(RemapperInfo& remapinfo) {
    //compute bulk moves for every species 
    TreeIterator* it_s = speciesTree->GetPostOrderIterator();
    while (SNode* s_src = it_s->next()) {
        if (remapinfo.hashtable[Getspecieslbl(s_src, numintnodes)].size() >= 1) {
            vector<SNode*> possibleremapping = GetPossibleUpRemaps(s_src);
            for (SNode* s_dest : possibleremapping) {
                bool added_move = ComputeBulkMove(s_src, s_dest, remapinfo);

                /*if (added_move) {
                    myfile << "set of gene size of : " << remapinfo.hashtable[Getspecieslbl(s_src, numintnodes)].return_max_heights().size() <<
                        " from " << s_src->GetLabel() << " -> " << s_dest->GetLabel()
                        << " cost change: " << remapinfo.moves.GetMove(remapinfo.moves.GetNbMoves() - 1).delta_cost << endl;
                }*/


            }
        }
    }
    speciesTree->CloseIterator(it_s);
}











SegmentalReconcilerInfo SegmentalReconciler::Reconcile()
{

	RemapperInfo remapinfo = BuildRemapperInfo();
	int dupHeightSum = GetDupHeightSum(remapinfo.hashtable);
	int nbLosses = GetNbLosses(remapinfo.curmap);
	double cost = dupHeightSum * remapinfo.dupcost + nbLosses * remapinfo.losscost;
	remapinfo.UpdateBestMapping(cost, dupHeightSum, nbLosses);
	    
	bool done = false;
    int stopStochastic = 0;
	
	ofstream costs;
    costs.open("stochastic_cost_changes.txt");
	
    while (!done) {
        
		RemapMove bestmove;
		
		if(stochastic){
			
			ComputeAllUpMoves(remapinfo);
			ComputeAllBulkMoves(remapinfo);
			
			//At this point, remapinfo.moves has all the moves.  If nothing can be done, we exit
			if (remapinfo.moves.GetNbMoves() == 0){
				done = true;
				break;	//evil break
			}
			
			bestmove = remapinfo.GetRandomMove();
			
			if (bestmove.delta_cost >= 0){
				stopStochastic++;
				if(stopStochastic == 100){ // If nothing improves after 100 moves, we stop
					done = true;
					break;	//evil break
				}
			}
			else{
				stopStochastic = 0;
			}
			costs << bestmove.delta_cost << endl;
		}
		

		else{
			
			ComputeAllUpMoves(remapinfo);
			
			//At this point, remapinfo.moves has all the moves.  If nothing can be done, we exit
			if (remapinfo.moves.GetNbMoves() == 0){
				done = true;
				break;	//evil break
			}
			
			bestmove = remapinfo.GetBestMove();
			
			if (bestmove.delta_cost > 0) {

				//if no move is good, try bulk moves
				ComputeAllBulkMoves(remapinfo);
				bestmove = remapinfo.GetBestMove();
				
				//if stil lnothing, give up
				if (remapinfo.moves.GetNbMoves() == 0 || bestmove.delta_cost > 0) {
					done = true;
					break;	//evil break
				}
			}
		}
		
		
        bestmove.debug_print();
        

		if (bestmove.type == UpMove)
			ApplyUpMove(bestmove.g, bestmove.s, remapinfo);
		else if (bestmove.type == BulkMove)
			ApplyBulkMove(bestmove.s_src, bestmove.s_dest, remapinfo);
		
		if(stochastic){
			dupHeightSum = GetDupHeightSum(remapinfo.hashtable);
			nbLosses = GetNbLosses(remapinfo.curmap);
			cost = dupHeightSum * remapinfo.dupcost + nbLosses * remapinfo.losscost;
			if(cost < remapinfo.bestcost)
				remapinfo.UpdateBestMapping(cost, dupHeightSum, nbLosses);
		}
		
		
		remapinfo.Reset();
    }
 
	SegmentalReconcilerInfo return_info;
 
	if(stochastic){
		
		return_info.dupHeightSum = remapinfo.bestdupHeightSum;
		return_info.nbLosses = remapinfo.bestNblosses;
		return_info.curmap = remapinfo.bestmap;
		
	}
	
	else{
		
		return_info.dupHeightSum = GetDupHeightSum(remapinfo.hashtable);
		int dupheightsum = GetExhaustiveDupHeight(remapinfo.curmap);

		if (dupheightsum != return_info.dupHeightSum) {
			cout << "ERROR: GetExhaustiveDupHeight != GetDupHeightSum, something is probably wrong with the H data structure" << endl;
		}

		//cout << "ri dh=" << return_info.dupHeightSum << "  exdh=" << dupheightsum << endl;
		return_info.nbLosses = GetNbLosses(remapinfo.curmap);
		return_info.curmap = remapinfo.curmap;
	}
			
	return return_info;
}





int SegmentalReconciler::GetDupHeightSum(vector<hashlist>& hashtable)
{
    int dupHeightSum = 0;
    for (int i = 0; i < hashtable.size(); i++) {
        dupHeightSum += hashtable[i].size();
    }
    return dupHeightSum;
}





int SegmentalReconciler::Getspecieslbl(Node* s, int numintnodes) {
    int str_lbl;
    if (s->IsLeaf()) {
        string str = s->GetLabel();
        str.erase(std::remove(str.begin(), str.end(), '\''), str.end());
        str_lbl = Util::ToInt(str);
        str_lbl = str_lbl + numintnodes;
        //cout << "str " << str << " strlbl " << str_lbl << endl;
    }
    else {
        str_lbl = Util::ToInt(s->GetLabel());
    }
    return str_lbl;
}







bool SegmentalReconciler::IsMapped(Node* g, GSMap& curmap){
	return (curmap.find(g) != curmap.end());
}





GSMap SegmentalReconciler::GetLCAMapping()
{
    GSMap lcamap;
    for (int i = 0; i < geneTrees.size(); i++)
    {
        Node* g = geneTrees[i];

        TreeIterator* it = g->GetPostOrderIterator();
        while (GNode* n = it->next())
        {
            if (n->IsLeaf()){
                lcamap[n] = this->leaf_species_map[n];
            }
            else{
                lcamap[n] = lcamap[n->GetChild(0)]->FindLCAWith(lcamap[n->GetChild(1)]);
            }
        }
        g->CloseIterator(it);
    }
	return lcamap;
}



int SegmentalReconciler::GetNbLosses(GSMap& mapping)
{
    int nblosses = 0;

    for (int i = 0; i < geneTrees.size(); i++){
        Node* genetree = geneTrees[i];

        TreeIterator* it = genetree->GetPostOrderIterator();
        while (GNode* g = it->next()){
            if (!g->IsLeaf()){
                bool isdup = this->IsDuplication(g, mapping);

                int d1 = GetSpeciesTreeDistance(mapping[g], mapping[g->GetChild(0)]);
                int d2 = GetSpeciesTreeDistance(mapping[g], mapping[g->GetChild(1)]);

                int losses_tmp = (double)(d1 + d2);
                if (!isdup){
                    losses_tmp -= 2;
                }

                nblosses += losses_tmp;
            }
        }

        genetree->CloseIterator(it);
    }

    return nblosses;
}




bool SegmentalReconciler::IsDuplication(Node* g, GSMap& mapping)
{
    if (g->IsLeaf())
        return false;

    Node* s = mapping[g];
    Node* s1 = mapping[g->GetChild(0)];
    Node* s2 = mapping[g->GetChild(1)];

    if (s1->HasAncestor(s2) || s2->HasAncestor(s1))
        return true;

    if (s != s1->FindLCAWith(s2))
        return true;

    return false;
}


int SegmentalReconciler::GetDuplicationHeightUnder(Node* g, Node* species, GSMap& curmap)
{
    if (!IsDuplication(g, curmap) || curmap[g] != species)
        return 0;

    int d1 = GetDuplicationHeightUnder(g->GetChild(0), species, curmap);
    int d2 = GetDuplicationHeightUnder(g->GetChild(1), species, curmap);

    return 1 + max(d1, d2);
}







int SegmentalReconciler::GetExhaustiveDupHeight(GSMap &curmap) {
    unordered_map<SNode*, int> duplicationHeights;
    TreeIterator* it_s = speciesTree->GetPostOrderIterator();
    while (Node* s = it_s->next()) {
        duplicationHeights[s] = 0;
    }
    speciesTree->CloseIterator(it_s);


    for (int i = 0; i < geneTrees.size(); i++) {
        Node* genetree = geneTrees[i];

        TreeIterator* it = genetree->GetPostOrderIterator();
        while (GNode* g = it->next()) {
            if (!g->IsLeaf()) {

                bool isdup = this->IsDuplication(g, curmap);

                if (isdup) {
                    SNode* s = curmap[g];
                    int dupheight = GetDuplicationHeightUnder(g, s, curmap);
                    duplicationHeights[s] = max(duplicationHeights[s], dupheight);
                }


            }
        }
        genetree->CloseIterator(it);
    }


    TreeIterator* it_s2 = speciesTree->GetPostOrderIterator();
    int dup_height_sum = 0;
    while (Node* s = it_s2->next()) {
        dup_height_sum += duplicationHeights[s];
    }
    speciesTree->CloseIterator(it_s2);

    return dup_height_sum;
}





//TODO: has duped code
SegmentalReconcilerInfo SegmentalReconciler::ReconcileWithLCAMap() {

    GSMap lcamap = GetLCAMapping();

    unordered_map<SNode*, int> duplicationHeights;
    TreeIterator* it_s = speciesTree->GetPostOrderIterator();
    while (Node* s = it_s->next()){
        duplicationHeights[s] = 0;
    }
    speciesTree->CloseIterator(it_s);
	
	
	int nblosses = 0;

	for (int i = 0; i < geneTrees.size(); i++){
        Node* genetree = geneTrees[i];

        TreeIterator* it = genetree->GetPostOrderIterator();
        while (GNode* g = it->next()){
            if (!g->IsLeaf()){
				
				bool isdup = this->IsDuplication(g, lcamap);

                int d1 = GetSpeciesTreeDistance(lcamap[g], lcamap[g->GetChild(0)]);
                int d2 = GetSpeciesTreeDistance(lcamap[g], lcamap[g->GetChild(1)]);

                int losses_tmp = (double)(d1 + d2);
                if (!isdup){
                    losses_tmp -= 2;
                }

                nblosses += losses_tmp;
				
				
				if (isdup){
					SNode* s = lcamap[g];
					int dupheight = GetDuplicationHeightUnder(g, s, lcamap);
					duplicationHeights[s] = max(duplicationHeights[s], dupheight);
				}
				
				
            }
        }
        genetree->CloseIterator(it);
	}
	
	
	TreeIterator* it_s2 = speciesTree->GetPostOrderIterator();
	int dup_height_sum = 0;
    while (Node* s = it_s2->next()){
        dup_height_sum += duplicationHeights[s];
    }
    speciesTree->CloseIterator(it_s2);
	
	
	SegmentalReconcilerInfo info;
    info.curmap = lcamap;
    info.dupHeightSum = dup_height_sum;
    info.nbLosses = nblosses;
    info.isBad = false;
    bool improve = true;
    int numofruns = 0;
	
    cout << "LCA reconciliation is finished!" << endl;
    return info;
}






SegmentalReconcilerInfo SegmentalReconciler::ReconcileWithSimphy() {


	GSMap curmap(this->leaf_species_map);

    unordered_map<SNode*, int> duplicationHeights;
    TreeIterator* it_s = speciesTree->GetPostOrderIterator();
    while (Node* s = it_s->next()){
        duplicationHeights[s] = 0;
    }
    speciesTree->CloseIterator(it_s);
	
	int nblosses = 0;
	
	for (int i = 0; i < geneTrees.size(); i++){
        Node* genetree = geneTrees[i];

        TreeIterator* it = genetree->GetPostOrderIterator();
        while (GNode* g = it->next()){
            if (!g->IsLeaf()){
				
				Node* s = GetSimphyMapping(g, curmap);
				curmap[g] = s;
				
				
				bool isdup = this->IsDuplication(g, curmap);

                int d1 = GetSpeciesTreeDistance(curmap[g], curmap[g->GetChild(0)]);
                int d2 = GetSpeciesTreeDistance(curmap[g], curmap[g->GetChild(1)]);

                int losses_tmp = (double)(d1 + d2);
                if (!isdup){
                    losses_tmp -= 2;
                }

                nblosses += losses_tmp;
				
				
				if (isdup){
					int dupheight = GetDuplicationHeightUnder(g, s, curmap);
					duplicationHeights[s] = max(duplicationHeights[s], dupheight);
				}
				
				
            }
        }
        genetree->CloseIterator(it);
	}
	
	
	TreeIterator* it_s2 = speciesTree->GetPostOrderIterator();
	int dup_height_sum = 0;
    while (Node* s = it_s2->next()){
        dup_height_sum += duplicationHeights[s];
    }
    speciesTree->CloseIterator(it_s2);
	
	
	SegmentalReconcilerInfo info;
    info.dupHeightSum = dup_height_sum;
    info.nbLosses = nblosses;
    info.isBad = false;
    bool improve = true;
    int numofruns = 0;
	
	info.curmap = curmap;
	
    cout << "Simphy reconciliation is finished!" << endl;
    return info;
	
}







SNode* SegmentalReconciler::GetSimphyMapping(GNode* g, GSMap& curmap)
{
    string err = "";
    if (g->IsLeaf())
    {
        err = "g is a leaf.";
    }
    if (curmap.count(g))
    {
        err = "g is already mapped.";
    }
    if (!curmap.count(g->GetChild(0)))
    {
        err = "g's child 0 is not mapped.";
    }
    if (!curmap.count(g->GetChild(1)))
    {
        err = "g's child 1 is not mapped.";
    }

    if (err != "")
    {
        cout << "Error in GetSimphyMapping: " << err << endl;
        throw "Error in GetSimphyMapping: " + err;
    }
    char separator = '-';
    int i = 0;
    string glbl = g->GetLabel();
    bool flag = false;
    bool flag2 = false;
    string s = glbl;
    

    char delimiter = '_';

    
    std::stringstream ss(s);

    std::vector<std::string> tokens;
    std::string token;

    // Tokenize the input string by the delimiter and store each token in the vector
    while (std::getline(ss, token, delimiter)) {
        tokens.push_back(token);
    }
    string species, event;
    species = "none";
    event = "none";
    std::string myString = s;

    std::istringstream iss(myString);
    std::vector<std::string> parts;

    // Extract all parts separated by underscores
    while (std::getline(iss, myString, '_')) {
        parts.push_back(myString);
    }

    // Check if there are at least three parts
    if (parts.size() >= 3) {
        // Output the second and third parts
        species = parts[1];
        event = parts[2];
    }
    else {
        std::cout << "Insufficient parts found." << std::endl;
    }


    if (event == "Dup") {
        g->SetDup(true);
    }
    else {
        g->SetDup(false);
    }
    
	
	//TODO: this is suboptimal, but ok that works I suppose
    SNode * look_species_ = curmap[g->GetChild(0)]->FindLCAWith(curmap[g->GetChild(1)]);

    while (!look_species_->IsRoot()) {

        string look_species = look_species_->GetLabel();

        if (species == look_species) {
            return look_species_;
        }
        if(!look_species_->IsRoot())
            look_species_ = look_species_->GetParent();
    }
    
    return look_species_;
}






int SegmentalReconciler::GetSpeciesTreeDistance(Node* x, Node* y)
{
    if (speciesTreeDistances.find(x) != speciesTreeDistances.end())
    {
        if (speciesTreeDistances[x].find(y) != speciesTreeDistances[x].end())
            return speciesTreeDistances[x][y];
    }


    int dist = 99999;
    Node* d = NULL;
    Node* a = NULL;   //d = descendant, y = ancestor
    if (x->HasAncestor(y))
    {
        d = x;
        a = y;
    }
    else if (y->HasAncestor(x))
    {
        d = y;
        a = x;
    }

    if (d && a)
    {
        dist = 0;

        while (d != a)
        {
            d = d->GetParent();
            dist++;
        }
    }

    speciesTreeDistances[x][y] = dist;
    speciesTreeDistances[y][x] = dist;

    return dist;
}









