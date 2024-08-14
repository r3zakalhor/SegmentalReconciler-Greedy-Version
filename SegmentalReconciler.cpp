#include "SegmentalReconciler.h"
#include <random>
#include <ctime>

SegmentalReconciler::SegmentalReconciler(vector<Node*>& geneTrees, Node* speciesTree, GSMap& leaf_species_map, double dupcost, double losscost, int nbspecies, int nbgenes, int numintnodes, bool stochastic, double stochastic_temperature, int nbstochasticLoops)
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
    this->stochastic_temperature = stochastic_temperature;
    this->nbstochasticLoops = nbstochasticLoops;

    this->max_remap_distance = 999999;

    this->debug_mode = false;
}



vector<SNode*> SegmentalReconciler::GetPossibleUpRemaps(SNode* spnode, SNode* lcamap_node)
{
    vector<SNode*> sp;

    bool still_within_dist = true;

    while (!spnode->IsRoot() && still_within_dist)
    {
        spnode = spnode->GetParent();

        if (lcamap_node && max_remap_distance < 99999) {
            if (GetSpeciesTreeDistance(spnode, lcamap_node) > max_remap_distance) {
                still_within_dist = false;
            }
        }

        if (still_within_dist)
            sp.push_back(spnode);
    }

    return sp;
}





/**
This function has two side effects:
- it assumes that remapinfo stuff has been calculated at the parent, and will add entries to remapinfo for n
- will add a RemapMove to the "remapinfo.moves" object
**/
void SegmentalReconciler::ComputeUpMove(GNode* g, SNode* s, RemapperInfo& remapinfo) {

    double delta_dup = 0;
    //int slbl = Util::ToInt(s->GetLabel());
    int sindex = s->GetIndex();
    int sindex_parent = sindex;
    if (!g->IsRoot())
        sindex_parent = remapinfo.curmap[g->GetParent()]->GetIndex();
        //slbl_parent = Util::ToInt(remapinfo.curmap[g->GetParent()]->GetLabel());


    SNode* mu = remapinfo.curmap[g];
    //int mu_specieslbl = Getspecieslbl(mu, numintnodes);
    int mu_speciesindex = mu->GetIndex();
    SNode* mpu = mu;
    if (!g->IsRoot()) {
        mpu = remapinfo.curmap[g->GetParent()];
    }



    /*
    * this below is misplaced and was causing bad remaps, and infinite loops
    if (!g->IsRoot() && mu == mpu)
        remapinfo.chain_lengths[g][mu] = remapinfo.chain_lengths[g->GetParent()][mu] + 1;
    else
        remapinfo.chain_lengths[g][mu] = 1;
    */

    if (g->IsRoot() || sindex < sindex_parent)
        remapinfo.chain_lengths[g][s] = 1;
    else
        remapinfo.chain_lengths[g][s] = remapinfo.chain_lengths[g->GetParent()][s] + 1;

    // dup changes, implements recurrences from paper
    if (g->IsRoot() || sindex < sindex_parent) {

        remapinfo.delta_gst[g][s][s] = std::max(remapinfo.chain_lengths[g][s], remapinfo.hashtable[sindex].size()) - remapinfo.hashtable[sindex].size();
        if (remapinfo.hashtable[mu_speciesindex].is_unique_max(g))
            remapinfo.delta_gst[g][s][mu] = -1;
        else
            remapinfo.delta_gst[g][s][mu] = 0;
        remapinfo.delta_gs[g][s] = remapinfo.delta_gst[g][s][s] + remapinfo.delta_gst[g][s][mu];
    }
    else {

        remapinfo.delta_gst[g][s][s] = std::max(remapinfo.chain_lengths[g][s], remapinfo.hashtable[sindex].size()) - remapinfo.hashtable[sindex].size();


        if (mu->GetIndex() < mpu->GetIndex())
            remapinfo.delta_gst[g->GetParent()][s][mu] = 0;			//WHY?  This sets something of the parent
        else	//WHY?  mu = mpu here, no?
            remapinfo.delta_gst[g->GetParent()][s][mu] = remapinfo.delta_gst[g->GetParent()][s][mpu];

        if (mu == mpu) {
            int hprim = remapinfo.hashtable[mu_speciesindex].size() + remapinfo.delta_gst[g->GetParent()][s][mu];
            int hmu = remapinfo.hashtable[mu_speciesindex].dupheight_at_subtree(g);
            int tmp = 0;
            if (hprim == hmu && remapinfo.hashtable[mu_speciesindex].is_unique_element(g, hprim))
                tmp = -1;
            remapinfo.delta_gst[g][s][mu] = remapinfo.delta_gst[g->GetParent()][s][mu] + tmp;
        }
        else {
            if (remapinfo.hashtable[mu_speciesindex].is_unique_max(g))
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
    else if (sindex <= sindex_parent) {
        remapinfo.lambda_gs[g][s] = lambda + GetSpeciesTreeDistance(s, mpu) - im_us_pu - GetSpeciesTreeDistance(mu, mpu) + impu;
    }
    else if (sindex > sindex_parent) {	//WHY not just else?
        remapinfo.lambda_gs[g][s] = remapinfo.lambda_gs[g->GetParent()][s] + lambda - GetSpeciesTreeDistance(s, mu);
    }


    double costchange = remapinfo.lambda_gs[g][s] * remapinfo.losscost + remapinfo.delta_gs[g][s] * remapinfo.dupcost;

    remapinfo.moves.AddUpMove(remapinfo.delta_gs[g][s], remapinfo.lambda_gs[g][s], costchange, g, s);



    if (debug_mode) {

        //do the remapping on a copy and check that computed delta dup is correct
        GSMap mapcopy = remapinfo.curmap;

        GNode* g_cur = g;
        bool done = false;
        while (!done) {
            mapcopy[g_cur] = s;
            if (g_cur->IsRoot()) {
                done = true;
            }
            else {
                g_cur = g_cur->GetParent();

                //we stop when g_cur is mapped to s or above, since the chain reaction of remaps can stop
                if (remapinfo.curmap[g_cur]->GetIndex() > s->GetIndex())
                    done = true;
            }
        }

        int dh_before = GetExhaustiveDupHeight(remapinfo.curmap);
        int dh_after = GetExhaustiveDupHeight(mapcopy);

        int delta_dh = dh_after - dh_before;


        if (delta_dh != remapinfo.delta_gs[g][s]) {
            cout << "ERROR: remapping g=" << g->GetLabel() << " to s=" << s->GetLabel() << " has different delta" << endl
                << "exhaustive delta_dh=" << delta_dh << "   delta_gs[g][s]=" << remapinfo.delta_gs[g][s] << endl;

        }
    }



}



bool SegmentalReconciler::ComputeBulkUpMove(SNode* s_src, SNode* s_dest, RemapperInfo& remapinfo) {

    unordered_set<GNode*> gnodes = remapinfo.hashtable[s_src->GetIndex()].return_max_heights();

    if (gnodes.size() == 0) {
        cout << "No gene mapped to src" << endl;
        return false;
    }

    if (gnodes.size() == 1) {   //bulk move must include two genes
        return false;
    }

    int nb_losses = 0;
    int max_dup_height_increase = -999999;

    for (GNode* g : gnodes) {
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

    remapinfo.moves.AddBulkUpMove(max_dup_height_increase, nb_losses, costchange, s_src, s_dest);

    return true;
}





void SegmentalReconciler::ApplyBulkUpMove(SNode* s_src, SNode* s_dest, RemapperInfo& remapinfo) {

    unordered_set<GNode*> gnodes = remapinfo.hashtable[s_src->GetIndex()].return_max_heights();

    if (gnodes.size() == 0) {
        cout << "Error: no gene is a dup in " << s_src->GetLabel() << endl;
    }

    for (Node* g : gnodes) {
        ApplyUpMove(g, s_dest, remapinfo);
    }
}






void SegmentalReconciler::ApplyUpMove(GNode* g, SNode* s, RemapperInfo& remapinfo) {

    GNode* g_cur = g;
    vector<GNode*> nodes_to_remap;

    bool done = false;
    while (!done) {

        nodes_to_remap.push_back(g_cur);

        if (g_cur->IsRoot())
            done = true;
        else {
            g_cur = g_cur->GetParent();

            //we stop when g_cur is mapped to s or above, since the chain reaction of remaps can stop
            if (remapinfo.curmap[g_cur]->GetIndex() > s->GetIndex())
                done = true;
        }
    }


    //important for hashtable data structure to do it in reverse
    for (int j = nodes_to_remap.size() - 1; j >= 0; j--) {
        GNode* g_cur = nodes_to_remap[j];

        if (IsDuplication(g_cur, remapinfo.curmap)) {
            SNode* s_cur = remapinfo.curmap[g_cur];
            int sindex = s_cur->GetIndex();
            int dupheight = GetDuplicationHeightUnder(g_cur, s_cur, remapinfo.curmap);
            remapinfo.hashtable[sindex].remove(dupheight, g_cur);
        }

    }


    //now we do the actual remap, making sure that we update the hash table
    for (GNode* g_cur : nodes_to_remap) {
        remapinfo.curmap[g_cur] = s;

        int dupheight = GetDuplicationHeightUnder(g_cur, s, remapinfo.curmap);
        if (dupheight > 0) {
            int sindex = s->GetIndex();
            remapinfo.hashtable[sindex].add_cell(dupheight, g_cur);
        }
    }
}





RemapperInfo SegmentalReconciler::BuildRemapperInfo() {

    RemapperInfo remapinfo;
    remapinfo.moves.stochastic_temperature = stochastic_temperature;
    remapinfo.curmap = GetLCAMapping();	//makes a copy, maybe slow
    this->lcamap = remapinfo.curmap;

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

                int sindex = s->GetIndex();

                remapinfo.hashtable[sindex].add_cell(dupheight, g);
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
                vector<SNode*> possibleremapping = GetPossibleUpRemaps(mu, this->lcamap[g]);


                if (!g->IsRoot() && remapinfo.curmap[g->GetParent()] == mu)
                    remapinfo.chain_lengths[g][mu] = remapinfo.chain_lengths[g->GetParent()][mu] + 1;
                else
                    remapinfo.chain_lengths[g][mu] = 1;



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




void SegmentalReconciler::ComputeAllBulkUpMoves(RemapperInfo& remapinfo) {
    //compute bulk moves for every species 
    TreeIterator* it_s = speciesTree->GetPreOrderIterator();
    while (SNode* s_src = it_s->next()) {
        if (remapinfo.hashtable[s_src->GetIndex()].size() >= 1) {
            vector<SNode*> possibleremapping;

            //******************
            //compute possible remaps, making sure if needed that no gene gets remapped higher than allowed
            if (max_remap_distance >= 99999) {
                possibleremapping = GetPossibleUpRemaps(s_src, nullptr);
            }
            else {
                unordered_set<GNode*> gnodes = remapinfo.hashtable[s_src->GetIndex()].return_max_heights();

                if (gnodes.size() >= 2) {
                    vector<SNode*> unfiltered_possibleremapping = GetPossibleUpRemaps(s_src, this->lcamap[*gnodes.begin()]);

                    //suboptimal
                    for (SNode* s : unfiltered_possibleremapping) {
                        bool all_ok = true;
                        for (GNode* g : gnodes) {
                            if (GetSpeciesTreeDistance(s, this->lcamap[g]) > max_remap_distance) {
                                all_ok = false;
                                break;
                            }
                        }
                        if (all_ok) {
                            possibleremapping.push_back(s);
                        }
                    }
                }
            }
            //******************

            for (SNode* s_dest : possibleremapping) {
                bool added_move = ComputeBulkUpMove(s_src, s_dest, remapinfo);
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
    int stochasticLoops = 0;

    ofstream costs;
    costs.open("stochastic_cost_changes.txt");

    while (!done) {

        RemapMove bestmove;

        if (stochastic) {

            stochasticLoops++;

            ComputeAllUpMoves(remapinfo);
            ComputeAllBulkUpMoves(remapinfo);
            ComputeAllDownMoves(remapinfo);
            ComputeAllBulkDownMoves(remapinfo);

            //At this point, remapinfo.moves has all the moves.  If nothing can be done, we exit
            if (remapinfo.moves.GetNbMoves() == 0) {
                done = true;
                break;	//evil break
            }

            bestmove = remapinfo.GetRandomMove();

            stochasticLoops++;
            if (stochasticLoops % 50 == 0) {
                cout << "Loops: " << stochasticLoops << endl;
            }

            //if (bestmove.delta_cost >= 0){

            if (stochasticLoops >= nbstochasticLoops) { // If nothing improves after 100 moves, we stop
                done = true;
                break;	//evil break
            }
            //}

            costs << bestmove.delta_cost << endl;
        }


        else {

            ComputeAllUpMoves(remapinfo);

            //At this point, remapinfo.moves has all the moves.  If nothing can be done, we exit
            //ML Jul31: Changed this, in case we need to do a down move
            /*if (remapinfo.moves.GetNbMoves() == 0){
                done = true;
                break;	//evil break
            }*/

            if (remapinfo.moves.GetNbMoves() > 0) {
                bestmove = remapinfo.GetBestMove();
            }

            if (remapinfo.moves.GetNbMoves() == 0 || bestmove.delta_cost > 0) {

                //if no move is good, try bulk/down moves
                ComputeAllDownMoves(remapinfo);
                ComputeAllBulkDownMoves(remapinfo);
                ComputeAllBulkUpMoves(remapinfo);
                bestmove = remapinfo.GetBestMove();

                //if still nothing, give up
                if (remapinfo.moves.GetNbMoves() == 0 || bestmove.delta_cost >= 0) {
                    done = true;
                    break;	//evil break
                }
            }
        }


        bestmove.debug_print();


        if (bestmove.type == UpMove)
            ApplyUpMove(bestmove.g, bestmove.s, remapinfo);
        else if (bestmove.type == DownMove)
            ApplyDownMove(bestmove.g, bestmove.s, remapinfo);
        else if (bestmove.type == BulkUpMove)
            ApplyBulkUpMove(bestmove.s_src, bestmove.s_dest, remapinfo);
        else if (bestmove.type == BulkDownMove)
            ApplyBulkDownMove(bestmove.s_src, bestmove.s_dest, remapinfo);

        if (stochastic) {
            dupHeightSum = GetDupHeightSum(remapinfo.hashtable);
            nbLosses = GetNbLosses(remapinfo.curmap);
            cost = dupHeightSum * remapinfo.dupcost + nbLosses * remapinfo.losscost;
            if (cost < remapinfo.bestcost)
                remapinfo.UpdateBestMapping(cost, dupHeightSum, nbLosses);
        }


        remapinfo.Reset();
    }

    SegmentalReconcilerInfo return_info;

    if (stochastic) {

        return_info.dupHeightSum = remapinfo.bestdupHeightSum;
        return_info.nbLosses = remapinfo.bestNblosses;
        return_info.curmap = remapinfo.bestmap;

    }

    else {

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







bool SegmentalReconciler::IsMapped(Node* g, GSMap& curmap) {
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
            if (n->IsLeaf()) {
                lcamap[n] = this->leaf_species_map[n];
            }
            else {
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

    for (int i = 0; i < geneTrees.size(); i++) {
        Node* genetree = geneTrees[i];

        TreeIterator* it = genetree->GetPostOrderIterator();
        while (GNode* g = it->next()) {
            if (!g->IsLeaf()) {
                bool isdup = this->IsDuplication(g, mapping);

                int d1 = GetSpeciesTreeDistance(mapping[g], mapping[g->GetChild(0)]);
                int d2 = GetSpeciesTreeDistance(mapping[g], mapping[g->GetChild(1)]);

                int losses_tmp = (double)(d1 + d2);
                if (!isdup) {
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







int SegmentalReconciler::GetExhaustiveDupHeight(GSMap& curmap) {
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
    while (Node* s = it_s->next()) {
        duplicationHeights[s] = 0;
    }
    speciesTree->CloseIterator(it_s);


    int nblosses = 0;

    for (int i = 0; i < geneTrees.size(); i++) {
        Node* genetree = geneTrees[i];

        TreeIterator* it = genetree->GetPostOrderIterator();
        while (GNode* g = it->next()) {
            if (!g->IsLeaf()) {

                bool isdup = this->IsDuplication(g, lcamap);

                int d1 = GetSpeciesTreeDistance(lcamap[g], lcamap[g->GetChild(0)]);
                int d2 = GetSpeciesTreeDistance(lcamap[g], lcamap[g->GetChild(1)]);

                int losses_tmp = (double)(d1 + d2);
                if (!isdup) {
                    losses_tmp -= 2;
                }

                nblosses += losses_tmp;


                if (isdup) {
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
    while (Node* s = it_s2->next()) {
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
    while (Node* s = it_s->next()) {
        duplicationHeights[s] = 0;
    }
    speciesTree->CloseIterator(it_s);

    int nblosses = 0;

    for (int i = 0; i < geneTrees.size(); i++) {
        Node* genetree = geneTrees[i];

        TreeIterator* it = genetree->GetPostOrderIterator();
        while (GNode* g = it->next()) {
            if (!g->IsLeaf()) {

                Node* s = GetSimphyMapping(g, curmap);
                curmap[g] = s;


                bool isdup = this->IsDuplication(g, curmap);

                int d1 = GetSpeciesTreeDistance(curmap[g], curmap[g->GetChild(0)]);
                int d2 = GetSpeciesTreeDistance(curmap[g], curmap[g->GetChild(1)]);

                int losses_tmp = (double)(d1 + d2);
                if (!isdup) {
                    losses_tmp -= 2;
                }

                nblosses += losses_tmp;


                if (isdup) {
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
    SNode* look_species_ = curmap[g->GetChild(0)]->FindLCAWith(curmap[g->GetChild(1)]);

    while (!look_species_->IsRoot()) {

        string look_species = look_species_->GetLabel();

        if (species == look_species) {
            return look_species_;
        }
        if (!look_species_->IsRoot())
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





void SegmentalReconciler::ComputeAllDownMoves(RemapperInfo& remapinfo) {
    //TODO: this is copy-pasted from ComputeAllUpMoves, consider merging

    //compute down moves for every gene tree node
    for (int i = 0; i < geneTrees.size(); i++) {
        GNode* g_root = geneTrees[i];

        TreeIterator* it = g_root->GetPreOrderIterator();	//needs to be preorder for DP to work
        while (GNode* g = it->next()) {
            if (!g->IsLeaf()) {
                SNode* mu = remapinfo.curmap[g];
                vector<SNode*> possibleremapping = GetPossibleDownRemaps(g, remapinfo.curmap);

                for (SNode* s : possibleremapping) {
                    ComputeDownMove(g, s, remapinfo);
                }
            }
        }
        g_root->CloseIterator(it);
    }

}



bool SegmentalReconciler::ComputeDownMove(GNode* g, SNode* s, RemapperInfo& remapinfo) {

    SNode* mu_orig = remapinfo.curmap[g];

    //get ancestors of g also mapped to mu, last element is the root of the dup subtree in mu_orig
    vector<GNode*> ancestors_with_same_mu = { g };
    while (!ancestors_with_same_mu.back()->IsRoot() && remapinfo.curmap[ancestors_with_same_mu.back()->GetParent()] == mu_orig) {
        ancestors_with_same_mu.push_back(ancestors_with_same_mu.back()->GetParent());
    }

    int dh_s_pre = GetDupHeight(s, remapinfo);
    int dh_mu_pre = GetDupHeight(mu_orig, remapinfo);

    int nb_losses_pre = 0;

    nb_losses_pre += this->GetNbLossesOnParentBranch(g, remapinfo.curmap);
    nb_losses_pre += this->GetNbLossesOnParentBranch(g->GetChild(0), remapinfo.curmap);
    nb_losses_pre += this->GetNbLossesOnParentBranch(g->GetChild(1), remapinfo.curmap);




    //idea is to temporarily remap g to s, and put it back after.
    //if you add an return "return" to this function for some reason, make sure you put mu_orig back
    int nb_losses_post = 0;
    remapinfo.curmap[g] = s;

    nb_losses_post += this->GetNbLossesOnParentBranch(g, remapinfo.curmap);
    nb_losses_post += this->GetNbLossesOnParentBranch(g->GetChild(0), remapinfo.curmap);
    nb_losses_post += this->GetNbLossesOnParentBranch(g->GetChild(1), remapinfo.curmap);

    //check changes in dup heights for s and mu_orig
    int local_dh_s_post = GetDuplicationHeightUnder(g, s, remapinfo.curmap);
    int local_dh_mu_post = GetDuplicationHeightUnder(ancestors_with_same_mu.back(), mu_orig, remapinfo.curmap);

    remapinfo.curmap[g] = mu_orig;

    //dup height in s has changed only if remapping g to s created a higher dup subtree
    int dh_s_post = dh_s_pre;
    if (local_dh_s_post > remapinfo.hashtable[s->GetIndex()].size())
        dh_s_post = local_dh_s_post;

    //dup height in mu changed only if highest ancestor is now the root of a smaller subtree, and no other large subtree exists
    int dh_mu_post = dh_mu_pre;
    int muindex = mu_orig->GetIndex();

    //case 1: dup subtree in mu height reduced by 1, and was unique max
    if (local_dh_mu_post == dh_mu_pre - 1 && remapinfo.hashtable[muindex].is_unique_max(ancestors_with_same_mu.back())) {
        dh_mu_post = local_dh_mu_post;
    }
    //case 2: dup subtree height reduced by 2 because parent of g became a spec
    else if (local_dh_mu_post == dh_mu_pre - 2) {
        if (remapinfo.hashtable[muindex].is_unique_max(ancestors_with_same_mu.back())) {
            dh_mu_post = dh_mu_pre - 1;

            if (remapinfo.hashtable[muindex].is_unique_after_max(ancestors_with_same_mu.end()[-2])) {
                dh_mu_post = dh_mu_pre - 2;
            }
        }
    }
    else if (local_dh_mu_post > dh_mu_pre) {
        std::cout << "Remapping g to s made local_dh_mu_post=" << local_dh_mu_post << " and dh_mu_pre=" << dh_mu_pre << endl;
        throw "local_dh_mu_post makes no sense";
    }




    int diff_loss = nb_losses_post - nb_losses_pre;
    int diff_dh = dh_s_post - dh_s_pre + dh_mu_post - dh_mu_pre;


    remapinfo.lambda_gs[g][s] = diff_loss;
    remapinfo.delta_gs[g][s] = diff_dh;
    remapinfo.delta_gst[g][s][s] = dh_s_post - dh_s_pre;
    remapinfo.delta_gst[g][s][mu_orig] = dh_mu_post - dh_mu_pre;

    double costchange = diff_loss * remapinfo.losscost + diff_dh * remapinfo.dupcost;

    remapinfo.moves.AddDownMove(diff_dh, diff_loss, costchange, g, s);

    return true;
}


int SegmentalReconciler::GetNbLossesOnParentBranch(GNode* g_child, GSMap& curmap) {
    if (g_child->IsRoot())
        return 0;

    bool parent_isdup = this->IsDuplication(g_child->GetParent(), curmap);

    int losses = GetSpeciesTreeDistance(curmap[g_child->GetParent()], curmap[g_child]);

    if (!parent_isdup) {
        losses -= 1;
    }

    return losses;

}


int SegmentalReconciler::GetDupHeight(SNode* s, RemapperInfo& remapinfo) {
    int sindex = s->GetIndex();
    return remapinfo.hashtable[sindex].size();
}



vector<SNode*> SegmentalReconciler::GetPossibleDownRemaps(GNode* gnode, GSMap& gsmap)
{
    vector<SNode*> sp;

    if (gnode->IsLeaf())
        return sp;

    SNode* curnode = gsmap[gnode->GetChild(0)]->FindLCAWith(gsmap[gnode->GetChild(1)]);

    while (curnode != gsmap[gnode])
    {
        sp.push_back(curnode);
        curnode = curnode->GetParent();
    }

    return sp;
}




void SegmentalReconciler::ApplyDownMove(GNode* g, SNode* s, RemapperInfo& remapinfo) {

    SNode* mu = remapinfo.curmap[g];

    //get ancestors of g also mapped to mu, last element is the root of the dup subtree in mu_orig
    //TODO: duped code in ComputeDownMove
    vector<GNode*> ancestors_with_same_mu = { g };
    while (!ancestors_with_same_mu.back()->IsRoot() && remapinfo.curmap[ancestors_with_same_mu.back()->GetParent()] == mu) {
        ancestors_with_same_mu.push_back(ancestors_with_same_mu.back()->GetParent());
    }




    //important for hashtable data structure to do it in reverse
    for (int j = ancestors_with_same_mu.size() - 1; j >= 0; j--) {
        GNode* g_cur = ancestors_with_same_mu[j];

        if (IsDuplication(g_cur, remapinfo.curmap)) {
            int muindex = mu->GetIndex();
            int dupheight = GetDuplicationHeightUnder(g_cur, mu, remapinfo.curmap);
            remapinfo.hashtable[muindex].remove(dupheight, g_cur);
        }
    }


    //now we do the actual remap, making sure that we update the hash table
    remapinfo.curmap[g] = s;

    int dupheight = GetDuplicationHeightUnder(g, s, remapinfo.curmap);
    if (dupheight > 0) {
        int sindex = s->GetIndex();
        remapinfo.hashtable[sindex].add_cell(dupheight, g);
    }

    //start at 1 to skip g
    for (int j = 1; j < ancestors_with_same_mu.size(); ++j) {
        GNode* g_cur = ancestors_with_same_mu[j];

        if (IsDuplication(g_cur, remapinfo.curmap)) {
            int muindex = mu->GetIndex();
            int dupheight = GetDuplicationHeightUnder(g_cur, mu, remapinfo.curmap);
            remapinfo.hashtable[muindex].add_cell(dupheight, g_cur);
        }
    }
}




void SegmentalReconciler::ComputeAllBulkDownMoves(RemapperInfo& remapinfo) {

    //compute bulk moves for every species 
    TreeIterator* it_s = speciesTree->GetPreOrderIterator();
    while (SNode* s_src = it_s->next()) {


        int max_height = remapinfo.hashtable[s_src->GetIndex()].size();
        unordered_set<GNode*> gnodes = remapinfo.hashtable[s_src->GetIndex()].return_max_heights();

        //gnodes contains the topmost dups, but here we want those at depth max_height
        vector<GNode*> deepest_nodes;
        for (GNode* g_cur : gnodes) {
            GetDupsAtDepth(g_cur, s_src, remapinfo.curmap, max_height, 1, deepest_nodes);
        }

        set<SNode*> possibleremaps;
        bool firstfill = true;

        //omg this is a slooow way of finding the species that all deepestnodes can map to
        for (GNode* g : deepest_nodes) {
            vector<SNode*> vec = GetPossibleDownRemaps(g, remapinfo.curmap);
            set<SNode*> setvec = set<SNode*>(vec.begin(), vec.end());


            if (firstfill) {
                possibleremaps = setvec;
                firstfill = false;
            }
            else {

                set<SNode*> inter;
                std::set_intersection(possibleremaps.begin(), possibleremaps.end(), setvec.begin(), setvec.end(),
                    std::inserter(inter, inter.begin()));
                possibleremaps = inter;
            }
        }


        for (SNode* s_dest : possibleremaps) {
            ComputeBulkDownMove(s_src, s_dest, deepest_nodes, remapinfo);
        }

    }

    speciesTree->CloseIterator(it_s);

}



bool SegmentalReconciler::ComputeBulkDownMove(SNode* s_src, SNode* s_dest, vector<GNode*>& deepest_nodes, RemapperInfo& remapinfo) {

    if (deepest_nodes.size() <= 1) {    //bulk move must include two genes
        return false;
    }


    int nb_losses = 0;
    int max_delta_gss = 0;
    bool all_reduce_by_two = true;

    for (GNode* g : deepest_nodes) {
        nb_losses += remapinfo.lambda_gs[g][s_dest];

        /*
        Idea of this below: each deepest may increase the height in s_dest by 1, and if one of them does the max will record a +1
        As for s_src, either we reduce the height by 1 by remapping every deepest node, or we reduce by 2 if all
        remaps reduce by 2.
        */
        max_delta_gss = max(max_delta_gss, remapinfo.delta_gst[g][s_dest][s_dest]);
        if (remapinfo.delta_gst[g][s_dest][s_src] > -2)
            all_reduce_by_two = false;
    }

    int delta_gs_mu = -1;
    if (all_reduce_by_two)
        delta_gs_mu = -2;

    int delta_dup = delta_gs_mu + max_delta_gss;

    double costchange = nb_losses * remapinfo.losscost + delta_dup * remapinfo.dupcost;

    remapinfo.moves.AddBulkDownMove(delta_dup, nb_losses, costchange, s_src, s_dest);

    return true;



}



void SegmentalReconciler::GetDupsAtDepth(GNode* curnode, SNode* s, GSMap& curmap, int target_depth, int cur_depth, vector<GNode*>& vec_to_fill) {
    if (curmap[curnode] != s || !this->IsDuplication(curnode, curmap))
        return;
    if (cur_depth > target_depth)
        return;

    //at this point, g is a dup mapped to s
    if (cur_depth == target_depth) {
        vec_to_fill.push_back(curnode);
        return; //we won't find anything deeper, so stop it
    }
    else {
        GetDupsAtDepth(curnode->GetChild(0), s, curmap, target_depth, cur_depth + 1, vec_to_fill);
        GetDupsAtDepth(curnode->GetChild(1), s, curmap, target_depth, cur_depth + 1, vec_to_fill);
    }
}




void SegmentalReconciler::ApplyBulkDownMove(SNode* s_src, SNode* s_dest, RemapperInfo& remapinfo) {
    //TODO: code is duped with ComputeAllBulkDownMoves
    int max_height = remapinfo.hashtable[s_src->GetIndex()].size();
    unordered_set<GNode*> gnodes = remapinfo.hashtable[s_src->GetIndex()].return_max_heights();

    vector<GNode*> deepest_nodes;
    for (GNode* g_cur : gnodes) {
        GetDupsAtDepth(g_cur, s_src, remapinfo.curmap, max_height, 1, deepest_nodes);
    }

    for (GNode* g : deepest_nodes) {
        ApplyDownMove(g, s_dest, remapinfo);
    }
}