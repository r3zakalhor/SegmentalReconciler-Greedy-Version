#include "SegmentalReconcile.h"

SegmentalReconcile::SegmentalReconcile(vector<Node*>& geneTrees, Node* speciesTree, unordered_map<Node*, Node*>& geneSpeciesMapping, double dupcost, double losscost, int maxDupHeight, int nbspecies)
{
    this->geneTrees = geneTrees;
    this->speciesTree = speciesTree;
    this->geneSpeciesMapping = geneSpeciesMapping;
    this->dupcost = dupcost;
    this->losscost = losscost;
    this->maxDupHeight = maxDupHeight;
    this->nbspecies = nbspecies;
    hashtable.resize(nbspecies);
}



SegmentalReconcileInfo SegmentalReconcile::Reconcile()
{
    ComputeLCAMapping();

    unordered_map<Node*, Node*> partialMapping(this->geneSpeciesMapping);

    unordered_map<Node*, int> duplicationHeights;
    TreeIterator* it = speciesTree->GetPostOrderIterator();
    while (Node* s = it->next())
    {
        duplicationHeights[s] = 0;
    }
    speciesTree->CloseIterator(it);

    vector<Node*> minimalNodes = GetMinimalUnmappedNodes(partialMapping);
    int dupheight = 0;

    //int added_losses = CleanupPartialMapping(partialMapping, duplicationHeights, minimalNodes);

    while (minimalNodes.size() > 0)
    {
        for (int j = minimalNodes.size() - 1; j >= 0; j--)
        {
            Node* g = minimalNodes[j];
            g->SetDup(false);
            //bool canBeSpec = !IsRequiredDuplication(g, partialMapping);
            //bool isEasyDup = IsEasyDuplication(g, partialMapping, duplicationheights);
            //if (canBeSpec || isEasyDup)
            //{
            Node* s = GetLowestPossibleMapping(g, partialMapping);
            partialMapping[g] = s;
            dupheight = GetDuplicationHeightUnder(g, s, partialMapping);
            //cout << " slbl " << s->GetLabel();
            if (dupheight > 0) {
                int slbl = Util::ToInt(s->GetLabel());
                int glbl = Util::ToInt(g->GetLabel());
                //cout << " s " << slbl << " g " << g->GetLabel() << " dupheight " << dupheight << endl;
                g->SetDup(true);
                hashtable[slbl].add_cell(dupheight, g);
            }
                //nblosses += GetSpeciesTreeDistance(s, partialMapping[g->GetChild(0)]);
                //nblosses += GetSpeciesTreeDistance(s, partialMapping[g->GetChild(1)]);
                //if (canBeSpec)
                    //nblosses -= 2;
                //the parent of g might become minimal - we'll add it in this case.
                //since we are iterating over new_minimals in the reverse order, this below works
            if (!g->IsRoot() && IsMinimalUnmapped(g->GetParent(), partialMapping))
            {
                minimalNodes.push_back(g->GetParent());
            }
            //}
            //no point in considering g from now on - we remove it from further consideration.
            minimalNodes.erase(minimalNodes.begin() + j);
        }
    }


    for (int i = 0; i < hashtable.size(); i++) {
        cout << "Species " << i << " ";
        hashtable[i].print();
        //cout << "size " << hashtable[i].size() << endl;
    }

    int nblosses = GetnbLosses(partialMapping);
    int dupheightsum = GetdupHeightSum(hashtable);
    cout << "/////////////////////////////////////////////////////////////////////////////////////" << endl;
    cout << " nb Losses of LCA : " << nblosses << " nb Dupheightsum of LCA : " << dupheightsum << endl;


    currentBestInfo.dupHeightSum = dupheightsum;
    currentBestInfo.nbLosses = nblosses;
    currentBestInfo.isBad = true;
    double cost = currentBestInfo.GetCost(dupcost, losscost);

    cout << " cost of LCA : " << cost << endl;

    //hashtable[13].remove(2, "160");
    //GreedyRemapping(partialMapping, hashtable, geneTrees, speciesTree, cost, 10, 0.1);
    SegmentalReconcileInfo info;
    info.dupHeightSum = dupheightsum;
    info.nbLosses = nblosses;
    info.isBad = false;
    info = GreedyRemapping(partialMapping, hashtable, geneTrees, speciesTree, cost, dupcost, losscost);
    //info.dupHeightSum = 0;
    //info.nbLosses = added_losses;
    info.partialMapping = partialMapping;

    SegmentalReconcileInfo retinfo = ReconcileRecursive(info, duplicationHeights);

    return retinfo;
}


SegmentalReconcileInfo SegmentalReconcile::GreedyRemapping(unordered_map<Node*, Node*>& partialMapping, vector<hashlist> hashtable, vector<Node*>& geneTrees, Node* speciesTree, double LCAcost, double dupcost, double losscost)
{   
    vector<hashlist> backuphash = hashtable;
    unordered_map<Node*, Node*> backuppartialMapping = partialMapping;
    vector<Node*> minimalNodes;
    SegmentalReconcileInfo greedyinfo;
    greedyinfo.dupHeightSum = GetdupHeightSum(hashtable);
    greedyinfo.nbLosses = GetnbLosses(partialMapping);
    SegmentalReconcileInfo greedyinfotmp;
    double currentbestcost = LCAcost;

    for (int i = 0; i < geneTrees.size(); i++)
    {
        Node* g = geneTrees[i];
        Node* currents;
        vector<Node*> possibleremapping;

        TreeIterator* it = g->GetPostOrderIterator();
        while (Node* n = it->next())
        {
            if (!n->IsLeaf())
            {
                currents = partialMapping[n];
                Node* s = currents;
                if (!s->IsRoot()) {
                    s = s->GetParent();
                    while (!s->IsRoot()) {
                        possibleremapping.push_back(s);
                        s = s->GetParent();
                    }
                    possibleremapping.push_back(s);
                }
                while (possibleremapping.size() > 0) {
                    s = possibleremapping.back();
                    possibleremapping.pop_back();
                    //cout << "remap " << n->GetLabel() << " from " << currents->GetLabel() << " to " << s->GetLabel() << endl;
                    if (IsDuplication(n, partialMapping)) {
                        int slbl = Util::ToInt(currents->GetLabel());
                        int dupheight = GetDuplicationHeightUnder(n, currents, partialMapping);
                        //cout << "is dup" << endl;
                        hashtable[slbl].remove(dupheight, n->GetLabel());
                    }
                    //cout << "before remapping" << endl;
                    partialMapping[n] = s;
                    //cout << "after remapping" << endl;
                    if (IsDuplication(n, partialMapping)) {
                        int dupheight = GetDuplicationHeightUnder(n, s, partialMapping);
                        //cout << "is new dup" << endl;
                        if (dupheight > 0) {
                            int slbl = Util::ToInt(s->GetLabel());
                            hashtable[slbl].add_cell(dupheight, n);
                        }
                    }

                    if (!n->IsRoot()) {
                        int slbl = Util::ToInt(s->GetLabel());
                        int slbl_parent = Util::ToInt(partialMapping[n->GetParent()]->GetLabel());
                        //cout << "slbl: " << slbl << " slbl parent " << slbl_parent << endl;
                        if (slbl >= slbl_parent)
                        {
                            minimalNodes.push_back(n->GetParent());
                        }
                    }

                    while (minimalNodes.size() > 0) {
                        for (int j = minimalNodes.size() - 1; j >= 0; j--)
                        {
                            Node* m = minimalNodes[j];
                            Node* currents1 = partialMapping[m];
                            if (IsDuplication(m, partialMapping)) {
                                int slbl = Util::ToInt(currents1->GetLabel());
                                int dupheight = GetDuplicationHeightUnder(m, currents1, partialMapping);
                                //cout << "remove " << dupheight << ", " << m->GetLabel() << " from " << slbl << endl;
                                hashtable[slbl].remove(dupheight, m->GetLabel());
                                //cout << "removed " << dupheight << ", " << m->GetLabel() << " from " << slbl << endl;
                            }
                            partialMapping.erase(m);
                            Node* s1 = GetLowestPossibleMapping(m, partialMapping); // should be same as s
                            //cout << "effected remap " << m->GetLabel() << " from " << currents1->GetLabel() << " to " << s1->GetLabel() << endl;
                            partialMapping[m] = s1;
                            int dupheight = GetDuplicationHeightUnder(m, s1, partialMapping);
                            if (dupheight > 0) {
                                int slbl = Util::ToInt(s1->GetLabel());
                                hashtable[slbl].add_cell(dupheight, m);
                            }
                            minimalNodes.erase(minimalNodes.begin() + j);
                            if (!m->IsRoot()) {
                                int slbl = Util::ToInt(s1->GetLabel());
                                int slbl_parent = Util::ToInt(partialMapping[m->GetParent()]->GetLabel());
                                if (slbl > slbl_parent)
                                {
                                    //cout << slbl << ">" << slbl_parent << endl;
                                    //cout << " next is " << m->GetParent()->GetLabel() << endl;
                                    minimalNodes.push_back(m->GetParent());
                                }
                            }
                        }
                    }
                    greedyinfotmp.nbLosses = GetnbLosses(partialMapping);
                    greedyinfotmp.dupHeightSum = GetdupHeightSum(hashtable);
                    double tmpcost = greedyinfotmp.GetCost(dupcost, losscost);
                    if (currentbestcost < tmpcost) {
                        hashtable = backuphash;
                        partialMapping = backuppartialMapping;
                    }
                    else {
                        backuphash = hashtable;
                        backuppartialMapping = partialMapping;
                        currentbestcost = tmpcost;
                        greedyinfo.nbLosses = greedyinfotmp.nbLosses;
                        greedyinfo.dupHeightSum = greedyinfotmp.dupHeightSum;
                        cout << "Good Remap!" << endl;
                    }

                }
            }

        }
        g->CloseIterator(it);
    }
    cout << "After run of Greedy remapping : " << endl;
    for (int i = 0; i < hashtable.size(); i++) {
        cout << "Species " << i << " ";
        hashtable[i].print();
        //cout << "size " << hashtable[i].size() << endl;
    }
    cout << "/////////////////////////////////////////////////////////////////////////////////////" << endl;
    cout << " nb Losses of Greedy: " << greedyinfo.nbLosses << " nb Dupheightsum of Greedy : " << greedyinfo.dupHeightSum << endl;
    cout << " cost of Greedy : " << currentbestcost << endl;

    return greedyinfo;
}


SegmentalReconcileInfo SegmentalReconcile::ReconcileRecursive(SegmentalReconcileInfo& info, unordered_map<Node*, int>& duplicationHeights)
{
    //IMPORTANT ASSERTION: partialMapping is clean

    //ASSERTION 2 : dupheights is smaller than maxDupheight
    if (info.dupHeightSum > maxDupHeight)
    {
        SegmentalReconcileInfo retinfo;
        retinfo.isBad = true;
        return retinfo;
    }

    //this makes this more of a branch-and-bound algorithm now...
    if (!currentBestInfo.isBad && currentBestInfo.GetCost(dupcost, losscost) < info.GetCost(dupcost, losscost))
    {
        info.isBad = true;
        return info;
    }


    //TODO: we could be more clever and avoid recomputing this at every recursion
    unordered_map<Node*, Node*> partialMapping = info.partialMapping;
    vector<Node*> minimalNodes = GetMinimalUnmappedNodes(partialMapping);

    if (minimalNodes.size() == 0) //normally, this means the mapping is complete
    {

        if (currentBestInfo.isBad || info.GetCost(dupcost, losscost) < currentBestInfo.GetCost(dupcost, losscost))
            currentBestInfo = info;

        return info;
    }
    else
    {
        Node* lowest = GetLowestMinimalNode(minimalNodes, partialMapping);

        vector<Node*> sps = GetPossibleSpeciesMapping(lowest, partialMapping);

        //we'll try mapping lowest to every possible species, and keep the best
        SegmentalReconcileInfo bestInfo;
        bestInfo.nbLosses = 999999;
        bestInfo.dupHeightSum = 999999;
        bestInfo.isBad = true;

        for (int i = 0; i < sps.size(); i++)
        {
            int local_nblosses = info.nbLosses;
            Node* s = sps[i];
            unordered_map<Node*, Node*> local_partialMapping(partialMapping);   //copy constructor called here
            unordered_map<Node*, int> local_duplicationHeights(duplicationHeights);   //copy constructor called here

            local_duplicationHeights[s] = duplicationHeights[s] + 1;    //requires proof, see paper

            local_partialMapping[lowest] = s;

            local_nblosses += GetSpeciesTreeDistance(s, local_partialMapping[lowest->GetChild(0)]);
            local_nblosses += GetSpeciesTreeDistance(s, local_partialMapping[lowest->GetChild(1)]);

            vector<Node*> new_minimals;
            //if parent has become minimal, we'll have to add it
            if (!lowest->IsRoot() && IsMinimalUnmapped(lowest->GetParent(), local_partialMapping))
            {
                new_minimals.push_back(lowest->GetParent());
            }

            //Map every minimal node that can be mapped to s
            for (int j = 0; j < minimalNodes.size(); j++)
            {
                Node* g = minimalNodes[j];

                if (g != lowest)
                {
                    Node* sg = GetLowestPossibleMapping(g, local_partialMapping);

                    if (sg->HasAncestor(s))
                    {
                        local_partialMapping[g] = s;

                        local_nblosses += GetSpeciesTreeDistance(s, local_partialMapping[g->GetChild(0)]);
                        local_nblosses += GetSpeciesTreeDistance(s, local_partialMapping[g->GetChild(1)]);

                        if (!g->IsRoot() && IsMinimalUnmapped(g->GetParent(), local_partialMapping))
                        {
                            new_minimals.push_back(g->GetParent());
                        }

                    }
                }
            }

            //CLEANUP PHASE
            int added_losses = CleanupPartialMapping(local_partialMapping, local_duplicationHeights, new_minimals);
            local_nblosses += added_losses;

            SegmentalReconcileInfo recursiveCallInfo;
            recursiveCallInfo.dupHeightSum = info.dupHeightSum + 1;
            recursiveCallInfo.nbLosses = local_nblosses;
            recursiveCallInfo.partialMapping = local_partialMapping;

            SegmentalReconcileInfo recursiveRetinfo = ReconcileRecursive(recursiveCallInfo, local_duplicationHeights);

            if (!recursiveRetinfo.isBad)
            {
                if (recursiveCallInfo.GetCost(dupcost, losscost) < bestInfo.GetCost(dupcost, losscost))
                    bestInfo = recursiveRetinfo;
            }
        }

        return bestInfo;
    }

}




int SegmentalReconcile::CleanupPartialMapping(unordered_map<Node*, Node*>& partialMapping, unordered_map<Node*, int>& duplicationheights, vector<Node*>& minimalNodes)
{
    int nblosses = 0;
    //CLEANUP PHASE
    //we want to do: while there is an easy node, map it
    //at this point, we know that minimals not in new_minimals cannot be speciations, nor mapped to s
    //so there is no point in checking them.
    while (minimalNodes.size() > 0)
    {
        for (int j = minimalNodes.size() - 1; j >= 0; j--)
        {
            Node* g = minimalNodes[j];
            bool canBeSpec = !IsRequiredDuplication(g, partialMapping);
            bool isEasyDup = IsEasyDuplication(g, partialMapping, duplicationheights);

            if (canBeSpec || isEasyDup)
            {
                Node* s = GetLowestPossibleMapping(g, partialMapping);

                partialMapping[g] = s;

                nblosses += GetSpeciesTreeDistance(s, partialMapping[g->GetChild(0)]);
                nblosses += GetSpeciesTreeDistance(s, partialMapping[g->GetChild(1)]);

                if (canBeSpec)
                    nblosses -= 2;

                //the parent of g might become minimal - we'll add it in this case.
                //since we are iterating over new_minimals in the reverse order, this below works
                if (!g->IsRoot() && IsMinimalUnmapped(g->GetParent(), partialMapping))
                {
                    minimalNodes.push_back(g->GetParent());
                }
            }

            //no point in considering g from now on - we remove it from further consideration.
            minimalNodes.erase(minimalNodes.begin() + j);
        }
    }

    return nblosses;
}



bool SegmentalReconcile::IsEasyDuplication(Node* g, unordered_map<Node*, Node*>& partialMapping, unordered_map<Node*, int>& duplicationHeights)
{
    if (!IsMinimalUnmapped(g, partialMapping))
    {
        cout << "Error in IsEasyDuplication: g is not minimal." << endl;
        throw "Error in IsEasyDuplication: g is not minimal.";
    }

    Node* lca = GetLowestPossibleMapping(g, partialMapping);

    //get dup height of lca under g
    int d1 = GetDuplicationHeightUnder(g->GetChild(0), lca, partialMapping);
    int d2 = GetDuplicationHeightUnder(g->GetChild(1), lca, partialMapping);

    int h = 1 + max(d1, d2);

    return (h <= duplicationHeights[lca]);

}


int SegmentalReconcile::GetDuplicationHeightUnder(Node* g, Node* species, unordered_map<Node*, Node*>& partialMapping)
{
    if (!IsDuplication(g, partialMapping) || partialMapping[g] != species)
        return 0;

    int d1 = GetDuplicationHeightUnder(g->GetChild(0), species, partialMapping);
    int d2 = GetDuplicationHeightUnder(g->GetChild(1), species, partialMapping);

    return 1 + max(d1, d2);
}




bool SegmentalReconcile::IsDuplication(Node* g, unordered_map<Node*, Node*>& partialMapping)
{
    if (g->IsLeaf())
        return false;

    Node* s = partialMapping[g];
    Node* s1 = partialMapping[g->GetChild(0)];
    Node* s2 = partialMapping[g->GetChild(1)];

    if (s1->HasAncestor(s2) || s2->HasAncestor(s1))
        return true;

    if (s != s1->FindLCAWith(s2))
        return true;

    return false;

}



vector<Node*> SegmentalReconcile::GetPossibleSpeciesMapping(Node* minimalNode, unordered_map<Node*, Node*>& partialMapping)
{
    Node* s = GetLowestPossibleMapping(minimalNode, partialMapping);

    vector<Node*> sps;

    bool done = false;
    while (!done)
    {
        sps.push_back(s);

        if (s->IsRoot() || sps.size() >= (int)(dupcost / losscost))
        {
            done = true;
        }
        else
        {
            s = s->GetParent();
        }
    }

    return sps;
}



Node* SegmentalReconcile::GetLowestMinimalNode(vector<Node*>& minimalNodes, unordered_map<Node*, Node*>& partialMapping)
{
    Node* curmin = minimalNodes[0];
    Node* curlca = GetLowestPossibleMapping(curmin, partialMapping);

    for (int i = 1; i < minimalNodes.size(); i++)
    {
        Node* lca = GetLowestPossibleMapping(minimalNodes[i], partialMapping);

        //if lowest possible mapping of i-th node is strictly below current, it becomes current
        if (lca->HasAncestor(curlca) && curlca != lca)
        {
            curmin = minimalNodes[i];
            curlca = lca;
        }
    }

    return curmin;
}



void SegmentalReconcile::ComputeLCAMapping()
{
    lcaMapping.clear();
    for (int i = 0; i < geneTrees.size(); i++)
    {
        Node* g = geneTrees[i];

        TreeIterator* it = g->GetPostOrderIterator();
        while (Node* n = it->next())
        {
            if (n->IsLeaf())
            {
                lcaMapping[n] = geneSpeciesMapping[n];
            }
            else
            {
                lcaMapping[n] = lcaMapping[n->GetChild(0)]->FindLCAWith(lcaMapping[n->GetChild(1)]);
            }
        }
        g->CloseIterator(it);
    }
}



bool SegmentalReconcile::IsMapped(Node* g, unordered_map<Node*, Node*>& partialMapping)
{
    return (partialMapping.find(g) != partialMapping.end());
}


vector<Node*> SegmentalReconcile::GetMinimalUnmappedNodes(unordered_map<Node*, Node*>& partialMapping)
{
    vector<Node*> minimalNodes;


    for (int i = 0; i < geneTrees.size(); i++)
    {
        Node* genetree = geneTrees[i];
        TreeIterator* it = genetree->GetPostOrderIterator();
        while (Node* g = it->next())
        {
            if (IsMinimalUnmapped(g, partialMapping))
            {
                minimalNodes.push_back(g);
            }
        }
        genetree->CloseIterator(it);
    }


    return minimalNodes;
}



bool SegmentalReconcile::IsMinimalUnmapped(Node* g, unordered_map<Node*, Node*>& partialMapping)
{
    return (!IsMapped(g, partialMapping) &&
        IsMapped(g->GetChild(0), partialMapping) &&
        IsMapped(g->GetChild(1), partialMapping));
}



Node* SegmentalReconcile::GetLowestPossibleMapping(Node* g, unordered_map<Node*, Node*>& partialMapping)
{
    string err = "";
    if (g->IsLeaf())
    {
        err = "g is a leaf.";
    }
    if (IsMapped(g, partialMapping))
    {
        err = "g is already mapped.";
    }
    if (!IsMapped(g->GetChild(0), partialMapping))
    {
        err = "g's child 0 is not mapped.";
    }
    if (!IsMapped(g->GetChild(1), partialMapping))
    {
        err = "g's child 1 is not mapped.";
    }

    if (err != "")
    {
        cout << "Error in GetLowestPossibleMapping: " << err << endl;
        throw "Error in GetLowestPossibleMapping: " + err;
    }

    return partialMapping[g->GetChild(0)]->FindLCAWith(partialMapping[g->GetChild(1)]);

}


bool SegmentalReconcile::IsRequiredDuplication(Node* g, unordered_map<Node*, Node*>& partialMapping)
{
    //See the required duplication Lemma in the paper to see that this works

    if (g->IsLeaf())
        return false;

    Node* lca = lcaMapping[g];
    Node* s1 = partialMapping[g->GetChild(0)];
    Node* s2 = partialMapping[g->GetChild(1)];

    return (lca->HasAncestor(s1) || lca->HasAncestor(s2));
}


int SegmentalReconcile::GetnbLosses(unordered_map<Node*, Node*>& fullMapping)
{
    double cost = 0;
    int nbdups = 0;
    int nblosses = 0;

    for (int i = 0; i < geneTrees.size(); i++)
    {
        Node* genetree = geneTrees[i];

        TreeIterator* it = genetree->GetPostOrderIterator();

        while (Node* g = it->next())
        {
            if (!g->IsLeaf())
            {
                bool isdup = this->IsDuplication(g, fullMapping);

                int d1 = GetSpeciesTreeDistance(fullMapping[g], fullMapping[g->GetChild(0)]);
                int d2 = GetSpeciesTreeDistance(fullMapping[g], fullMapping[g->GetChild(1)]);

                int losses_tmp = (double)(d1 + d2);
                if (!isdup)
                {
                    losses_tmp -= 2;
                }
                else
                {
                    //nbdups++;
                    //cost += this->dupcost;
                }

                nblosses += losses_tmp;
                cost += losses_tmp * this->losscost;
            }
        }

        genetree->CloseIterator(it);

    }

    /*int dupheight = 0;

    //COMPUTE DUP HEIGHTS THE HARD WAY...
    TreeIterator* its = speciesTree->GetPostOrderIterator();
    while (Node* s = its->next())
    {
        int maxh = 0;
        for (int i = 0; i < geneTrees.size(); i++)
        {
            Node* genetree = geneTrees[i];

            TreeIterator* it = genetree->GetPostOrderIterator();
            while (Node* g = it->next())
            {
                //I know I know this is suboptimal.
                int h = this->GetDuplicationHeightUnder(g, s, fullMapping);
                if (h > maxh)
                    maxh = h;
            }
            genetree->CloseIterator(it);
        }

        dupheight += maxh;

    }
    speciesTree->CloseIterator(its);

    cost += dupheight * this->dupcost;*/

    return nblosses;
}


int SegmentalReconcile::GetdupHeightSum(vector<hashlist>& hashtable)
{
    int dupHeightSum = 0;
    for (int i = 0; i < hashtable.size(); i++) {
        
        dupHeightSum += hashtable[i].size();
    }
    return dupHeightSum;
}



double SegmentalReconcile::GetMappingCost(unordered_map<Node*, Node*>& fullMapping)
{
    double cost = 0;
    int nbdups = 0;
    int nblosses = 0;

    for (int i = 0; i < geneTrees.size(); i++)
    {
        Node* genetree = geneTrees[i];

        TreeIterator* it = genetree->GetPostOrderIterator();

        while (Node* g = it->next())
        {
            if (!g->IsLeaf())
            {
                bool isdup = this->IsDuplication(g, fullMapping);

                int d1 = GetSpeciesTreeDistance(fullMapping[g], fullMapping[g->GetChild(0)]);
                int d2 = GetSpeciesTreeDistance(fullMapping[g], fullMapping[g->GetChild(1)]);

                int losses_tmp = (double)(d1 + d2);
                if (!isdup)
                {
                    losses_tmp -= 2;
                }
                else
                {
                    //nbdups++;
                    //cost += this->dupcost;
                }

                nblosses += losses_tmp;
                cost += losses_tmp * this->losscost;
            }
        }

        genetree->CloseIterator(it);

    }

    int dupheight = 0;

    //COMPUTE DUP HEIGHTS THE HARD WAY...
    TreeIterator* its = speciesTree->GetPostOrderIterator();
    while (Node* s = its->next())
    {
        int maxh = 0;
        for (int i = 0; i < geneTrees.size(); i++)
        {
            Node* genetree = geneTrees[i];

            TreeIterator* it = genetree->GetPostOrderIterator();
            while (Node* g = it->next())
            {
                //I know I know this is suboptimal.
                int h = this->GetDuplicationHeightUnder(g, s, fullMapping);
                if (h > maxh)
                    maxh = h;
            }
            genetree->CloseIterator(it);
        }

        dupheight += maxh;

    }
    speciesTree->CloseIterator(its);

    cost += dupheight * this->dupcost;

    return cost;
}




int SegmentalReconcile::GetSpeciesTreeDistance(Node* x, Node* y)
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
