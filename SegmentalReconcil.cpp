#include "SegmentalReconcile.h"


SegmentalReconcile::SegmentalReconcile(vector<Node*>& geneTrees, Node* speciesTree, unordered_map<Node*, Node*>& geneSpeciesMapping, double dupcost, double losscost, int maxDupHeight, int nbspecies, int nbgenes, string algorithm)
{
    this->geneTrees = geneTrees;
    this->speciesTree = speciesTree;
    this->geneSpeciesMapping = geneSpeciesMapping;
    this->dupcost = dupcost;
    this->losscost = losscost;
    this->maxDupHeight = maxDupHeight;
    this->nbspecies = nbspecies;
    this->nbgenes = nbgenes;
    this->algorithm = algorithm;
    hashtable.resize(nbspecies);
    DupChanges = new int*[nbgenes];
    LossChanges = new int*[nbgenes];
    Chain = new int* [nbgenes];
    for (int i = 0; i < nbgenes; ++i) {
        DupChanges[i] = new int[nbspecies];
        LossChanges[i] = new int[nbspecies];
        Chain[i] = new int[nbspecies];
    }
    
    TDupChanges = new int** [nbgenes];
    for (int i = 0; i < nbgenes; i++) {
        TDupChanges[i] = new int* [nbspecies];
        for (int j = 0; j < nbspecies; j++) {
            TDupChanges[i][j] = new int[nbspecies];
        }
    }
}



SegmentalReconcileInfo SegmentalReconcile::Reconcile()
{
    SegmentalReconcileInfo info;
    
    if (algorithm == "lca") {
        info = lca_algorithm();
        return info;
    }
    else if (algorithm == "simphy") {
        info = simphy_algorithm();
        return info;
    }
    else if (algorithm == "greedy") {
        info = greedy_algorithm();
        return info;
    }
    else if (algorithm == "ultragreedy") {
        info = ultragreedy_algorithm();
        return info;
    }

    ComputeLCAMapping();

    unordered_map<Node*, Node*> partialMapping(this->geneSpeciesMapping);
    Node* gg;
    unordered_map<Node*, int> duplicationHeights;
    TreeIterator* it = speciesTree->GetPostOrderIterator();
    while (Node* s = it->next())
    {
        duplicationHeights[s] = 0;
    }
    speciesTree->CloseIterator(it);

    vector<Node*> minimalNodes = GetMinimalUnmappedNodes(partialMapping);
    int dupheight = 0;
    int numofnodes = 0;
    //int added_losses = CleanupPartialMapping(partialMapping, duplicationHeights, minimalNodes);

    while (minimalNodes.size() > 0)
    {
        for (int j = minimalNodes.size() - 1; j >= 0; j--)
        {
            Node* g = minimalNodes[j];
            //g->SetDup(false);
            //bool canBeSpec = !IsRequiredDuplication(g, partialMapping);
            //bool isEasyDup = IsEasyDuplication(g, partialMapping, duplicationheights);
            //if (canBeSpec || isEasyDup)
            //{
            
            Node* s = GetLowestPossibleMapping(g, partialMapping);
            
            //Node* s = GetSimphyMapping(g, partialMapping);
            //cout << " s " << s->GetLabel() << " g " << g->GetLabel() << " dupheight " << dupheight << endl;
            numofnodes++;
            partialMapping[g] = s;
            dupheight = GetDuplicationHeightUnder(g, s, partialMapping);
            //cout << " slbl " << s->GetLabel();
            if (dupheight > 0) {
                int slbl = Util::ToInt(s->GetLabel());
                int glbl = Util::ToInt(g->GetLabel());
                //cout << " s " << slbl << " g " << g->GetLabel() << " dupheight " << dupheight << endl;
                //g->SetDup(true);
                hashtable[slbl].add_cell(dupheight, g);
                gg = g;
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

    cout << " Number of Nodes : " << numofnodes << endl;


    
    info.dupHeightSum = dupheightsum;
    info.nbLosses = nblosses;
    info.isBad = false;
    bool improve = true;
    int numofruns = 0;
    /*while (improve)
    {
        cost = info.GetCost(dupcost, losscost);
        info = UltraGreedyRemapping(partialMapping, hashtable, geneTrees, speciesTree, cost, dupcost, losscost, &improve);
        numofruns++;
    }*/

    // Here we are able to run ultra greedy or greedy algorithm
    //info = GreedyRemapping(partialMapping, hashtable, geneTrees, speciesTree, cost, dupcost, losscost, &improve);
    //info = UltraGreedyRemapping(partialMapping, hashtable, geneTrees, speciesTree, cost, dupcost, losscost, &improve);
    // 
    // 
    //info.dupHeightSum = 0;
    //info.nbLosses = added_losses;
    cout << "Number of Runs of Greedy: " << numofruns << endl;
    info.partialMapping = partialMapping;

    //SegmentalReconcileInfo retinfo = ReconcileRecursive(info, duplicationHeights);

    return info;
}

SegmentalReconcileInfo SegmentalReconcile::lca_algorithm() {

    ComputeLCAMapping();

    unordered_map<Node*, Node*> partialMapping(this->geneSpeciesMapping);
    Node* gg;
    unordered_map<Node*, int> duplicationHeights;
    TreeIterator* it = speciesTree->GetPostOrderIterator();
    while (Node* s = it->next())
    {
        duplicationHeights[s] = 0;
    }
    speciesTree->CloseIterator(it);

    vector<Node*> minimalNodes = GetMinimalUnmappedNodes(partialMapping);
    int dupheight = 0;
    int numofnodes = 0;
    //int added_losses = CleanupPartialMapping(partialMapping, duplicationHeights, minimalNodes);

    while (minimalNodes.size() > 0)
    {
        for (int j = minimalNodes.size() - 1; j >= 0; j--)
        {
            Node* g = minimalNodes[j];
            //g->SetDup(false);
            //bool canBeSpec = !IsRequiredDuplication(g, partialMapping);
            //bool isEasyDup = IsEasyDuplication(g, partialMapping, duplicationheights);
            //if (canBeSpec || isEasyDup)
            //{

            Node* s = GetLowestPossibleMapping(g, partialMapping);

            //Node* s = GetSimphyMapping(g, partialMapping);
            //cout << " s " << s->GetLabel() << " g " << g->GetLabel() << " dupheight " << dupheight << endl;
            numofnodes++;
            partialMapping[g] = s;
            dupheight = GetDuplicationHeightUnder(g, s, partialMapping);
            //cout << " slbl " << s->GetLabel();
            if (dupheight > 0) {
                int slbl = Util::ToInt(s->GetLabel());
                int glbl = Util::ToInt(g->GetLabel());
                //cout << " s " << slbl << " g " << g->GetLabel() << " dupheight " << dupheight << endl;
                //g->SetDup(true);
                hashtable[slbl].add_cell(dupheight, g);
                gg = g;
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


    /*for (int i = 0; i < hashtable.size(); i++) {
        cout << "Species " << i << " ";
        hashtable[i].print();
        //cout << "size " << hashtable[i].size() << endl;
    }*/

    int nblosses = GetnbLosses(partialMapping);
    int dupheightsum = GetdupHeightSum(hashtable);
    //cout << "/////////////////////////////////////////////////////////////////////////////////////" << endl;
    //cout << " nb Losses of LCA : " << nblosses << " nb Dupheightsum of LCA : " << dupheightsum << endl;


    currentBestInfo.dupHeightSum = dupheightsum;
    currentBestInfo.nbLosses = nblosses;
    currentBestInfo.isBad = true;
    double cost = currentBestInfo.GetCost(dupcost, losscost);

    //cout << " cost of LCA : " << cost << endl;

    //cout << " Number of Nodes : " << numofnodes << endl;


    SegmentalReconcileInfo info;
    info.dupHeightSum = dupheightsum;
    info.nbLosses = nblosses;
    info.isBad = false;
    bool improve = true;
    int numofruns = 0;
    /*while (improve)
    {
        cost = info.GetCost(dupcost, losscost);
        info = UltraGreedyRemapping(partialMapping, hashtable, geneTrees, speciesTree, cost, dupcost, losscost, &improve);
        numofruns++;
    }*/

    // Here we are able to run ultra greedy or greedy algorithm
    //info = GreedyRemapping(partialMapping, hashtable, geneTrees, speciesTree, cost, dupcost, losscost, &improve);
    //info = UltraGreedyRemapping(partialMapping, hashtable, geneTrees, speciesTree, cost, dupcost, losscost, &improve);
    // 
    // 
    //info.dupHeightSum = 0;
    //info.nbLosses = added_losses;
    //cout << "Number of Runs of Greedy: " << numofruns << endl;
    info.partialMapping = partialMapping;

    //SegmentalReconcileInfo retinfo = ReconcileRecursive(info, duplicationHeights);
    cout << "LCA reconciliation is finished!" << endl;
    return info;
}

SegmentalReconcileInfo SegmentalReconcile::simphy_algorithm() {

    ComputeLCAMapping();

    unordered_map<Node*, Node*> partialMapping(this->geneSpeciesMapping);
    Node* gg;
    unordered_map<Node*, int> duplicationHeights;
    TreeIterator* it = speciesTree->GetPostOrderIterator();
    while (Node* s = it->next())
    {
        duplicationHeights[s] = 0;
    }
    speciesTree->CloseIterator(it);

    vector<Node*> minimalNodes = GetMinimalUnmappedNodes(partialMapping);
    int dupheight = 0;
    int numofnodes = 0;
    //int added_losses = CleanupPartialMapping(partialMapping, duplicationHeights, minimalNodes);

    while (minimalNodes.size() > 0)
    {
        for (int j = minimalNodes.size() - 1; j >= 0; j--)
        {
            Node* g = minimalNodes[j];
            //g->SetDup(false);
            //bool canBeSpec = !IsRequiredDuplication(g, partialMapping);
            //bool isEasyDup = IsEasyDuplication(g, partialMapping, duplicationheights);
            //if (canBeSpec || isEasyDup)
            //{

            //Node* s = GetLowestPossibleMapping(g, partialMapping);

            Node* s = GetSimphyMapping(g, partialMapping);
            //cout << " s " << s->GetLabel() << " g " << g->GetLabel() << " dupheight " << dupheight << endl;
            numofnodes++;
            partialMapping[g] = s;
            dupheight = GetDuplicationHeightUnder(g, s, partialMapping);
            //cout << " slbl " << s->GetLabel();
            if (dupheight > 0) {
                int slbl = Util::ToInt(s->GetLabel());
                int glbl = Util::ToInt(g->GetLabel());
                //cout << " s " << slbl << " g " << g->GetLabel() << " dupheight " << dupheight << endl;
                //g->SetDup(true);
                hashtable[slbl].add_cell(dupheight, g);
                gg = g;
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


    /*for (int i = 0; i < hashtable.size(); i++) {
        cout << "Species " << i << " ";
        hashtable[i].print();
        //cout << "size " << hashtable[i].size() << endl;
    }*/

    int nblosses = GetnbLosses(partialMapping);
    int dupheightsum = GetdupHeightSum(hashtable);
    //cout << "/////////////////////////////////////////////////////////////////////////////////////" << endl;
    //cout << " nb Losses of LCA : " << nblosses << " nb Dupheightsum of LCA : " << dupheightsum << endl;


    currentBestInfo.dupHeightSum = dupheightsum;
    currentBestInfo.nbLosses = nblosses;
    currentBestInfo.isBad = true;
    double cost = currentBestInfo.GetCost(dupcost, losscost);

    //cout << " cost of LCA : " << cost << endl;

    //cout << " Number of Nodes : " << numofnodes << endl;


    SegmentalReconcileInfo info;
    info.dupHeightSum = dupheightsum;
    info.nbLosses = nblosses;
    info.isBad = false;
    bool improve = true;
    int numofruns = 0;
    /*while (improve)
    {
        cost = info.GetCost(dupcost, losscost);
        info = UltraGreedyRemapping(partialMapping, hashtable, geneTrees, speciesTree, cost, dupcost, losscost, &improve);
        numofruns++;
    }*/

    // Here we are able to run ultra greedy or greedy algorithm
    //info = GreedyRemapping(partialMapping, hashtable, geneTrees, speciesTree, cost, dupcost, losscost, &improve);
    //info = UltraGreedyRemapping(partialMapping, hashtable, geneTrees, speciesTree, cost, dupcost, losscost, &improve);
    // 
    // 
    //info.dupHeightSum = 0;
    //info.nbLosses = added_losses;
    //cout << "Number of Runs of Greedy: " << numofruns << endl;
    info.partialMapping = partialMapping;

    //SegmentalReconcileInfo retinfo = ReconcileRecursive(info, duplicationHeights);
    cout << "simphy reconciliation is finished!" << endl;
    return info;
}

SegmentalReconcileInfo SegmentalReconcile::greedy_algorithm() {

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
    int numofnodes = 0;
    //int added_losses = CleanupPartialMapping(partialMapping, duplicationHeights, minimalNodes);

    while (minimalNodes.size() > 0)
    {
        for (int j = minimalNodes.size() - 1; j >= 0; j--)
        {
            Node* g = minimalNodes[j];
            //g->SetDup(false);
            //bool canBeSpec = !IsRequiredDuplication(g, partialMapping);
            //bool isEasyDup = IsEasyDuplication(g, partialMapping, duplicationheights);
            //if (canBeSpec || isEasyDup)
            //{

            Node* s = GetLowestPossibleMapping(g, partialMapping);

            //Node* s = GetSimphyMapping(g, partialMapping);
            //cout << " s " << s->GetLabel() << " g " << g->GetLabel() << " dupheight " << dupheight << endl;
            numofnodes++;
            partialMapping[g] = s;
            dupheight = GetDuplicationHeightUnder(g, s, partialMapping);
            //cout << " slbl " << s->GetLabel();
            if (dupheight > 0) {
                int slbl = Util::ToInt(s->GetLabel());
                int glbl = Util::ToInt(g->GetLabel());
                //cout << " s " << slbl << " g " << g->GetLabel() << " dupheight " << dupheight << endl;
                //g->SetDup(true);
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


    /*for (int i = 0; i < hashtable.size(); i++) {
        cout << "Species " << i << " ";
        hashtable[i].print();
        //cout << "size " << hashtable[i].size() << endl;
    }*/

    int nblosses = GetnbLosses(partialMapping);
    int dupheightsum = GetdupHeightSum(hashtable);
    //cout << "/////////////////////////////////////////////////////////////////////////////////////" << endl;
    //cout << " nb Losses of LCA : " << nblosses << " nb Dupheightsum of LCA : " << dupheightsum << endl;


    currentBestInfo.dupHeightSum = dupheightsum;
    currentBestInfo.nbLosses = nblosses;
    currentBestInfo.isBad = true;
    double cost = currentBestInfo.GetCost(dupcost, losscost);

    //cout << " cost of LCA : " << cost << endl;

    //cout << " Number of Nodes : " << numofnodes << endl;


    SegmentalReconcileInfo info;
    info.dupHeightSum = dupheightsum;
    info.nbLosses = nblosses;
    info.isBad = false;
    bool improve = true;
    int numofruns = 0;
    /*while (improve)
    {
        cost = info.GetCost(dupcost, losscost);
        info = UltraGreedyRemapping(partialMapping, hashtable, geneTrees, speciesTree, cost, dupcost, losscost, &improve);
        numofruns++;
    }*/

    // Here we are able to run ultra greedy or greedy algorithm
    //info = GreedyRemapping(partialMapping, hashtable, geneTrees, speciesTree, cost, dupcost, losscost, &improve);
    info = FastGreedyRemapping(partialMapping, hashtable, Chain, TDupChanges, DupChanges, LossChanges, geneTrees, speciesTree, cost, dupcost, losscost, &improve);
    //info = UltraGreedyRemapping(partialMapping, hashtable, geneTrees, speciesTree, cost, dupcost, losscost, &improve);
    // 
    // 
    //info.dupHeightSum = 0;
    //info.nbLosses = added_losses;
    //cout << "Number of Runs of Greedy: " << numofruns << endl;
    info.partialMapping = partialMapping;

    //SegmentalReconcileInfo retinfo = ReconcileRecursive(info, duplicationHeights);
    cout << "greedy reconciliation is finished!" << endl;
    return info;
}

SegmentalReconcileInfo SegmentalReconcile::ultragreedy_algorithm() {

    ComputeLCAMapping();

    unordered_map<Node*, Node*> partialMapping(this->geneSpeciesMapping);
    Node* gg;
    unordered_map<Node*, int> duplicationHeights;
    TreeIterator* it = speciesTree->GetPostOrderIterator();
    while (Node* s = it->next())
    {
        duplicationHeights[s] = 0;
    }
    speciesTree->CloseIterator(it);

    vector<Node*> minimalNodes = GetMinimalUnmappedNodes(partialMapping);
    int dupheight = 0;
    int numofnodes = 0;
    //int added_losses = CleanupPartialMapping(partialMapping, duplicationHeights, minimalNodes);

    while (minimalNodes.size() > 0)
    {
        for (int j = minimalNodes.size() - 1; j >= 0; j--)
        {
            Node* g = minimalNodes[j];
            //g->SetDup(false);
            //bool canBeSpec = !IsRequiredDuplication(g, partialMapping);
            //bool isEasyDup = IsEasyDuplication(g, partialMapping, duplicationheights);
            //if (canBeSpec || isEasyDup)
            //{

            Node* s = GetLowestPossibleMapping(g, partialMapping);

            //Node* s = GetSimphyMapping(g, partialMapping);
            //cout << " s " << s->GetLabel() << " g " << g->GetLabel() << " dupheight " << dupheight << endl;
            numofnodes++;
            partialMapping[g] = s;
            dupheight = GetDuplicationHeightUnder(g, s, partialMapping);
            //cout << " slbl " << s->GetLabel();
            if (dupheight > 0) {
                int slbl = Util::ToInt(s->GetLabel());
                int glbl = Util::ToInt(g->GetLabel());
                //cout << " s " << slbl << " g " << g->GetLabel() << " dupheight " << dupheight << endl;
                //g->SetDup(true);
                hashtable[slbl].add_cell(dupheight, g);
                gg = g;
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


    /*for (int i = 0; i < hashtable.size(); i++) {
        cout << "Species " << i << " ";
        hashtable[i].print();
        //cout << "size " << hashtable[i].size() << endl;
    }*/

    int nblosses = GetnbLosses(partialMapping);
    int dupheightsum = GetdupHeightSum(hashtable);
    //cout << "/////////////////////////////////////////////////////////////////////////////////////" << endl;
    //cout << " nb Losses of LCA : " << nblosses << " nb Dupheightsum of LCA : " << dupheightsum << endl;


    currentBestInfo.dupHeightSum = dupheightsum;
    currentBestInfo.nbLosses = nblosses;
    currentBestInfo.isBad = true;
    double cost = currentBestInfo.GetCost(dupcost, losscost);

    //cout << " cost of LCA : " << cost << endl;

    //cout << " Number of Nodes : " << numofnodes << endl;


    SegmentalReconcileInfo info;
    info.dupHeightSum = dupheightsum;
    info.nbLosses = nblosses;
    info.isBad = false;
    bool improve = true;
    int numofruns = 0;
    /*while (improve)
    {
        cost = info.GetCost(dupcost, losscost);
        info = UltraGreedyRemapping(partialMapping, hashtable, geneTrees, speciesTree, cost, dupcost, losscost, &improve);
        numofruns++;
    }*/

    // Here we are able to run ultra greedy or greedy algorithm
    //info = GreedyRemapping(partialMapping, hashtable, geneTrees, speciesTree, cost, dupcost, losscost, &improve);
    info = UltraGreedyRemapping(partialMapping, hashtable, geneTrees, speciesTree, cost, dupcost, losscost, &improve);
    // 
    // 
    //info.dupHeightSum = 0;
    //info.nbLosses = added_losses;
    //cout << "Number of Runs of Greedy: " << numofruns << endl;
    info.partialMapping = partialMapping;

    //SegmentalReconcileInfo retinfo = ReconcileRecursive(info, duplicationHeights);
    cout << "ultragreedy reconciliation is finished!" << endl;
    return info;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int SegmentalReconcile::CalculateCostChange(Node* n, Node* s, unordered_map<Node*, Node*>& partialMapping, vector<hashlist>& hashtable, int** Chain, int*** TDupChanges, int** DupChanges, int** LossChanges, vector<Node*>& geneTrees, Node* speciesTree, double dupcost, double losscost) {

    int slbl = Util::ToInt(s->GetLabel());
    Node* mu = partialMapping[n];
    int mu_lbl = Util::ToInt(mu->GetLabel());
    int mu_index = mu->GetIndex();
    int g_index = n->GetIndex();
    int s_index = s->GetIndex();

    int delta_dup = 0;

    int slbl_parent = slbl;
    Node* mpu = mu;
    int mpu_index = mu_index;
    int g_parent_index = g_index;

    if (!n->IsRoot()) {
        slbl_parent = Util::ToInt(partialMapping[n->GetParent()]->GetLabel());
        mpu = partialMapping[n->GetParent()];
        mpu_index = mpu->GetIndex();
        g_parent_index = (n->GetParent())->GetIndex();
    }

    if (!n->IsRoot() && mu_index == mpu_index)
        Chain[g_index][mu_index] = Chain[g_parent_index][mu_index] + 1;
    else
        Chain[g_index][mu_index] = 1;

    // dup changes
    if (n->IsRoot() || slbl < slbl_parent) {
        Chain[g_index][s_index] = 1;
        TDupChanges[g_index][s_index][s_index] = std::max(Chain[g_index][s_index], hashtable[slbl].size()) - hashtable[slbl].size();
        if (hashtable[slbl].is_unique_max(n))
            TDupChanges[g_index][s_index][mu_index] = -1;
        else 
            TDupChanges[g_index][s_index][mu_index] = 0;
        DupChanges[g_index][s_index] = TDupChanges[g_index][s_index][s_index] + TDupChanges[g_index][s_index][mu_index];
    }
    else {
        /*if (slbl == slbl_parent) {
            Chain[g_parent_index][s_index] = 1;
        }
        else {
            Chain[g_index][s_index] = Chain[g_parent_index][s_index] + 1;
        }*/
        Chain[g_index][s_index] = Chain[g_parent_index][s_index] + 1;
        TDupChanges[g_index][s_index][s_index] = std::max(Chain[g_index][s_index], hashtable[slbl].size()) - hashtable[slbl].size();
        if (mu->GetLabel() < mpu->GetLabel())
            TDupChanges[g_parent_index][s_index][mu_index] = 0;
        else
            TDupChanges[g_parent_index][s_index][mu_index] = TDupChanges[g_parent_index][s_index][mpu_index];
        if (mu->GetLabel() == mpu->GetLabel()) {
            int hprim = hashtable[mu_lbl].size() + TDupChanges[g_parent_index][s_index][mu_index];
            int hmu = hashtable[mu_lbl].dupheight_at_subtree(n);
            int tmp = 0;
            if (hprim == hmu && hashtable[mu_lbl].is_unique_element(n,hprim))
                tmp = -1;
            TDupChanges[g_index][s_index][mu_index] = TDupChanges[g_parent_index][s_index][mu_index] + tmp;
        }
        else {
            if(hashtable[slbl].is_unique_max(n))
                TDupChanges[g_index][s_index][mu_index] = -1;
            else 
                TDupChanges[g_index][s_index][mu_index] = 0;
        }
        if (mpu_index == s_index) {
            DupChanges[g_parent_index][s_index] = 0;
            TDupChanges[g_parent_index][s_index][s_index] = 0;
            TDupChanges[g_parent_index][s_index][mu_index] = 0;
        }
        DupChanges[g_index][s_index] = DupChanges[g_parent_index][s_index] - TDupChanges[g_parent_index][s_index][s_index] - TDupChanges[g_parent_index][s_index][mu_index] +
            TDupChanges[g_index][s_index][s_index] + TDupChanges[g_index][s_index][mu_index];

        if (DupChanges[g_index][s_index] > 100) {
            cout << "TDupChanges[g_parent_index][s_index][mu_index] " << TDupChanges[g_parent_index][s_index][mu_index] << endl;
            cout << "g_parent_index: " << g_parent_index << "s_index: " << s_index << endl;
        }
    }

    // loss changes
    int imu = 0;
    int impu = 0;
    int im_us_pu = 0;

    if (IsDuplication(n, partialMapping))
        imu = 0;
    else
        imu = 2;

    if (!n->IsRoot()) {
        if (IsDuplication(n->GetParent(), partialMapping))
            impu = 0;
        else
            impu = 2;

        partialMapping[n] = s;
        if (IsDuplication(n->GetParent(), partialMapping))
            im_us_pu = 0;
        else
            im_us_pu = 2;
        partialMapping[n] = mu;
    }

    int lambda = GetSpeciesTreeDistance(s, partialMapping[n->GetChild(0)]) + GetSpeciesTreeDistance(s, partialMapping[n->GetChild(1)]) -
        GetSpeciesTreeDistance(mu, partialMapping[n->GetChild(0)]) - GetSpeciesTreeDistance(mu, partialMapping[n->GetChild(1)]) + imu;
    if (n->IsRoot()) {
        LossChanges[g_index][s_index] = lambda;
    }
    else if (slbl <= slbl_parent) {
        LossChanges[g_index][s_index] = lambda + GetSpeciesTreeDistance(s, mpu) - im_us_pu - GetSpeciesTreeDistance(mu, mpu) + impu;
    }
    else if (slbl > slbl_parent) {
        LossChanges[g_index][s_index] = LossChanges[g_parent_index][s_index] + lambda - GetSpeciesTreeDistance(s, mu);
    }



    int costchange = LossChanges[g_index][s_index] * losscost + DupChanges[g_index][s_index] * dupcost;

    return costchange;
}

int SegmentalReconcile::ApplyChange(Node* n, Node* s, unordered_map<Node*, Node*>& partialMapping, vector<hashlist>& hashtable) {

    vector<Node*> minimalNodes;
    minimalNodes.push_back(n);

    while (minimalNodes.size() > 0) {
        for (int j = minimalNodes.size() - 1; j >= 0; j--)
        {
            Node* m = minimalNodes[j];
            Node* currents = partialMapping[m];

            if (IsDuplication(m, partialMapping)) {
                int slbl = Util::ToInt(currents->GetLabel());
                int dupheight = GetDuplicationHeightUnder(m, currents, partialMapping);
                //cout << "remove " << dupheight << ", " << m->GetLabel() << " from " << slbl << endl;
                bool fl = hashtable[slbl].remove(dupheight, m);
            }


            partialMapping.erase(m);

            //cout << "effected remap " << m->GetLabel() << " from " << currents1->GetLabel() << " to " << s1->GetLabel() << endl;
            partialMapping[m] = s;

            int dupheight = GetDuplicationHeightUnder(m, s, partialMapping);
            if (dupheight > 0) {
                int slbl = Util::ToInt(s->GetLabel());
                hashtable[slbl].add_cell(dupheight, m);
            }

            minimalNodes.erase(minimalNodes.begin() + j);

            if (!m->IsRoot()) {
                int slbl = Util::ToInt(s->GetLabel());
                int slbl_parent = Util::ToInt(partialMapping[m->GetParent()]->GetLabel());
                if (slbl >= slbl_parent)
                {
                    //cout << slbl << ">" << slbl_parent << endl;
                    //cout << " next is " << m->GetParent()->GetLabel() << endl;
                    minimalNodes.push_back(m->GetParent());
                }
            }
        }
    }
    return -1;
}

SegmentalReconcileInfo SegmentalReconcile::FastGreedyRemapping(unordered_map<Node*, Node*>& partialMapping, vector<hashlist>& hashtable, int** Chain, int*** TDupChanges, int** DupChanges, int** LossChanges, vector<Node*>& geneTrees, Node* speciesTree, double LCAcost, double dupcost, double losscost, bool* improve)
{
    int BestCostChange, CurrentCostChange;
    BestCostChange = 0;
    vector<Node*> minimalNodes;
    SegmentalReconcileInfo greedyinfo;
    bool remap_find = false;
    Node* BestChange[2];


    ofstream myfile;
    myfile.open("example.txt");

    int cntn = 0;
    int loop = 1;
    *improve = true;
    bool newremmap = true;
    myfile << "start!" << endl;
    while (newremmap) {
        BestCostChange = 0;
        cntn = 0;
        newremmap = false;
        for (int i = 0; i < geneTrees.size(); i++)
        {
            Node* g = geneTrees[i];
            Node* currents;
            vector<Node*> possibleremapping;

            TreeIterator* it = g->GetPreOrderIterator();
            while (Node* n = it->next())
            {
                myfile << "index: " << Util::ToString(n->GetIndex()) << endl;
                if (!n->IsLeaf())
                {
                    currents = partialMapping[n];
                    Node* s = currents;

                    Node* mu = partialMapping[n];
                    int mu_index = mu->GetIndex();
                    int g_index = n->GetIndex();
                    Node* mpu = mu;
                    int mpu_index = mu_index;
                    int g_parent_index = g_index;

                    if (!n->IsRoot()) {
                        mpu = partialMapping[n->GetParent()];
                        mpu_index = mpu->GetIndex();
                        g_parent_index = (n->GetParent())->GetIndex();
                    }
                    if (!n->IsRoot() && mu_index == mpu_index)
                        Chain[g_index][mu_index] = Chain[g_parent_index][mu_index] + 1;
                    else
                        Chain[g_index][mu_index] = 1;


                    cntn++;
                    if (cntn % 2000 == 0) {
                        cout << cntn << endl;
                    }
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
                        
                        CurrentCostChange = CalculateCostChange(n, s, partialMapping, hashtable, Chain, TDupChanges, DupChanges, LossChanges, geneTrees, speciesTree, dupcost, losscost);
                        myfile << "gene index: " << Util::ToString(n->GetIndex()) << " from " << Util::ToString(currents->GetIndex()) << " -> " << Util::ToString(s->GetIndex()) 
                            << " Dup change: " << Util::ToString(DupChanges[n->GetIndex()][s->GetIndex()]) << " loss change: " << Util::ToString(LossChanges[n->GetIndex()][s->GetIndex()]) 
                            << " cost change: "<< CurrentCostChange << endl;
                        
                        if (CurrentCostChange <= BestCostChange) {
                            BestCostChange = CurrentCostChange;
                            BestChange[0] = n;  
                            BestChange[1] = s;
                            newremmap = true;
                            remap_find = true;
                            cout << "Good Remap!" << endl;
                            cout << "gene index: " << Util::ToString(n->GetIndex()) << " from " << Util::ToString(currents->GetIndex()) << " -> " << Util::ToString(s->GetIndex()) << endl;
                        }

                    }

                }

            }
            g->CloseIterator(it);
        }
        if (remap_find) {
            ApplyChange(BestChange[0], BestChange[1], partialMapping, hashtable);
            myfile << "at iteration: " << loop << "gene index: " << Util::ToString(BestChange[0]->GetIndex()) << " -> " << Util::ToString(BestChange[1]->GetIndex()) << " Cost change: " << BestCostChange << endl;
            loop++;
            remap_find = false;
        }
    }
    /*cout << "After run of Greedy remapping : " << endl;
    for (int i = 0; i < hashtable.size(); i++) {
        cout << "Species " << i << " ";
        hashtable[i].print();
        //cout << "size " << hashtable[i].size() << endl;
    }
    cout << "/////////////////////////////////////////////////////////////////////////////////////" << endl;
    cout << " nb Losses of Greedy: " << bestgreedyinfo.nbLosses << " nb Dupheightsum of Greedy : " << bestgreedyinfo.dupHeightSum << endl;
    cout << " cost of Greedy : " << currentbestcost << endl;*/
    greedyinfo.dupHeightSum = GetdupHeightSum(hashtable);
    greedyinfo.nbLosses = GetnbLosses(partialMapping);
    double Cost = greedyinfo.GetCost(dupcost, losscost);

    return greedyinfo;
}


/// <summary>
/// ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// </summary>
/// <param name="partialMapping"></param>
/// <param name="hashtable"></param>
/// <param name="geneTrees"></param>
/// <param name="speciesTree"></param>
/// <param name="LCAcost"></param>
/// <param name="dupcost"></param>
/// <param name="losscost"></param>
/// <param name="improve"></param>
/// <returns></returns>

SegmentalReconcileInfo SegmentalReconcile::GreedyRemapping(unordered_map<Node*, Node*>& partialMapping, vector<hashlist>& hashtable, vector<Node*>& geneTrees, Node* speciesTree, double LCAcost, double dupcost, double losscost, bool* improve)
{
    vector<hashlist> backuphash = hashtable;
    unordered_map<Node*, Node*> backuppartialMapping = partialMapping;
    vector<hashlist> besthash = hashtable;
    unordered_map<Node*, Node*> bestpartialMapping = partialMapping;

    vector<Node*> minimalNodes;
    SegmentalReconcileInfo backupgreedyinfo, bestgreedyinfo;
    backupgreedyinfo.dupHeightSum = GetdupHeightSum(hashtable);
    backupgreedyinfo.nbLosses = GetnbLosses(partialMapping);
    bestgreedyinfo.dupHeightSum = GetdupHeightSum(hashtable);
    bestgreedyinfo.nbLosses = GetnbLosses(partialMapping);
    int currenttotalnblosses = backupgreedyinfo.nbLosses;
    int currentbestnblosses = backupgreedyinfo.nbLosses;
    SegmentalReconcileInfo greedyinfotmp;

    Triple T;
    vector<Triple> modify_add;
    vector<Triple> modify_remove;
    unordered_map<Node*, Node*> firstmap;
    unordered_map<Node*, Node*> currentmap;
    unordered_map<Node*, Node*> remap;
    unordered_map<Node*, Node*> remap1;
    ofstream myfile;
    myfile.open("example.txt");

    double currentbestcost = LCAcost;
    int cntn = 0;
    bool isdup;
    int d1, d2, losses_tmp;
    *improve = true;
    bool newremmap = true;
    int num = 0;
    while (newremmap) {
        cntn = 0;
        newremmap = false;
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
                    firstmap[n] = currents;
                    myfile << "1.firstmap[n]: " << n->GetLabel() << " to " << currents->GetLabel() << endl;
                    Node* s = currents;
                    cntn++;
                    if (cntn % 2000 == 0) {
                        cout << cntn << endl;
                    }
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
                        isdup = IsDuplication(n, backuppartialMapping);

                        //calulate number of losses for n 
                        d1 = GetSpeciesTreeDistance(backuppartialMapping[n], backuppartialMapping[n->GetChild(0)]);
                        d2 = GetSpeciesTreeDistance(backuppartialMapping[n], backuppartialMapping[n->GetChild(1)]);
                        losses_tmp = (double)(d1 + d2);

                        if (IsDuplication(n, partialMapping)) {
                            int slbl = Util::ToInt(currents->GetLabel());
                            int dupheight = GetDuplicationHeightUnder(n, currents, partialMapping);
                            //cout << "is dup" << endl;
                            bool fl = hashtable[slbl].remove(dupheight, n);
                            /*if (fl) {
                                T.slbl = slbl;
                                T.dupHeight = dupheight;
                                T.n = n;
                                modify_add.push_back(T);
                            }*/
                        }
                        if(!isdup)
                        {
                            losses_tmp -= 2;
                        }

                        currenttotalnblosses -= losses_tmp;
                        //cout << "before remapping" << endl;
                        partialMapping[n] = s;
                        remap[n] = s;
                        myfile << "1.remap[n]: " << n->GetLabel() << " to " << s->GetLabel() << endl;
                        //cout << "after remapping" << endl;
                        isdup = IsDuplication(n, partialMapping);
                        //calulate new number of losses for n 
                        d1 = GetSpeciesTreeDistance(partialMapping[n], partialMapping[n->GetChild(0)]);
                        d2 = GetSpeciesTreeDistance(partialMapping[n], partialMapping[n->GetChild(1)]);
                        losses_tmp = (double)(d1 + d2);

                        if (isdup) {
                            int dupheight = GetDuplicationHeightUnder(n, s, partialMapping);
                            //cout << "is new dup" << endl;
                            if (dupheight > 0) {
                                int slbl = Util::ToInt(s->GetLabel());
                                hashtable[slbl].add_cell(dupheight, n);
                                /*T.slbl = slbl;
                                T.dupHeight = dupheight;
                                T.n = n;
                                modify_remove.push_back(T);*/
                            }
                        }
                        else
                        {
                            losses_tmp -= 2;
                        }
                        currenttotalnblosses += losses_tmp;
                        bool check_losses = true;
                        Node* last_ancestor = n;

                        if (!n->IsRoot()) {
                            int slbl = Util::ToInt(s->GetLabel());
                            int slbl_parent = Util::ToInt(partialMapping[n->GetParent()]->GetLabel());
                            //cout << "slbl: " << slbl << " slbl parent " << slbl_parent << endl;
                            last_ancestor = n->GetParent();
                            if (slbl >= slbl_parent)
                            {
                                minimalNodes.push_back(n->GetParent());
                            }
                            if (slbl == slbl_parent) {
                                check_losses = false;
                            }
                        }

                        while (minimalNodes.size() > 0) {
                            for (int j = minimalNodes.size() - 1; j >= 0; j--)
                            {
                                Node* m = minimalNodes[j];
                                Node* currents1 = partialMapping[m];
                                currentmap[m] = currents1;
                                myfile << "2.currentmap[m]: " << m->GetLabel() << " to " << currents1->GetLabel() << endl;
                                isdup = IsDuplication(m, backuppartialMapping);

                                //calulate number of losses for n 
                                d1 = GetSpeciesTreeDistance(backuppartialMapping[m], backuppartialMapping[m->GetChild(0)]);
                                d2 = GetSpeciesTreeDistance(backuppartialMapping[m], backuppartialMapping[m->GetChild(1)]);
                                losses_tmp = (double)(d1 + d2);
                                if (!isdup) {
                                    losses_tmp -= 2;
                                }
                                currenttotalnblosses -= losses_tmp;

                                if (IsDuplication(m, partialMapping)) {
                                    int slbl = Util::ToInt(currents1->GetLabel());
                                    int dupheight = GetDuplicationHeightUnder(m, currents1, partialMapping);
                                    //cout << "remove " << dupheight << ", " << m->GetLabel() << " from " << slbl << endl;
                                    bool fl = hashtable[slbl].remove(dupheight, m);
                                    /*if (fl) {
                                        T.slbl = slbl;
                                        T.dupHeight = dupheight;
                                        T.n = n;
                                        modify_add.push_back(T);
                                    }*/
                                    //cout << "removed " << dupheight << ", " << m->GetLabel() << " from " << slbl << endl;
                                }


                                partialMapping.erase(m);
                                Node* s1 = GetLowestPossibleMapping(m, partialMapping); // should be same as s
                                //cout << "effected remap " << m->GetLabel() << " from " << currents1->GetLabel() << " to " << s1->GetLabel() << endl;
                                partialMapping[m] = s1;
                                remap[m] = s1;
                                myfile << "2.remap[m]: " << m->GetLabel() << " to " << s1->GetLabel() << endl;
                                isdup = IsDuplication(m, partialMapping);
                                //calulate new number of losses for n 
                                d1 = GetSpeciesTreeDistance(partialMapping[m], partialMapping[m->GetChild(0)]);
                                d2 = GetSpeciesTreeDistance(partialMapping[m], partialMapping[m->GetChild(1)]);
                                losses_tmp = (double)(d1 + d2);
                                if (!isdup) {
                                    losses_tmp -= 2;
                                }
                                currenttotalnblosses += losses_tmp;

                                int dupheight = GetDuplicationHeightUnder(m, s1, partialMapping);
                                if (dupheight > 0) {
                                    int slbl = Util::ToInt(s1->GetLabel());
                                    hashtable[slbl].add_cell(dupheight, m);
                                    /*T.slbl = slbl;
                                    T.dupHeight = dupheight;
                                    T.n = n;
                                    modify_remove.push_back(T);*/
                                }
                                minimalNodes.erase(minimalNodes.begin() + j);
                                if (!m->IsRoot()) {
                                    int slbl = Util::ToInt(s1->GetLabel());
                                    int slbl_parent = Util::ToInt(partialMapping[m->GetParent()]->GetLabel());
                                    last_ancestor = m->GetParent();
                                    if (slbl >= slbl_parent)
                                    {
                                        //cout << slbl << ">" << slbl_parent << endl;
                                        //cout << " next is " << m->GetParent()->GetLabel() << endl;
                                        minimalNodes.push_back(m->GetParent());
                                    }
                                    if (slbl == slbl_parent) {
                                        check_losses = false;
                                    }
                                }
                            }
                        }

                        if (check_losses) {

                            isdup = IsDuplication(last_ancestor, backuppartialMapping);
                            //calulate number of losses for n 
                            d1 = GetSpeciesTreeDistance(backuppartialMapping[last_ancestor], backuppartialMapping[last_ancestor->GetChild(0)]);
                            d2 = GetSpeciesTreeDistance(backuppartialMapping[last_ancestor], backuppartialMapping[last_ancestor->GetChild(1)]);
                            losses_tmp = (double)(d1 + d2);
                            if (!isdup)
                            {
                                losses_tmp -= 2;
                            }
                            currenttotalnblosses -= losses_tmp;

                            isdup = IsDuplication(last_ancestor, partialMapping);
                            //calulate new number of losses for n 
                            d1 = GetSpeciesTreeDistance(partialMapping[last_ancestor], partialMapping[last_ancestor->GetChild(0)]);
                            d2 = GetSpeciesTreeDistance(partialMapping[last_ancestor], partialMapping[last_ancestor->GetChild(1)]);
                            losses_tmp = (double)(d1 + d2);
                            if (!isdup) {
                                losses_tmp -= 2;
                            }
                            currenttotalnblosses += losses_tmp;
                        }

                        greedyinfotmp.nbLosses = currenttotalnblosses;
                        greedyinfotmp.dupHeightSum = GetdupHeightSum(hashtable);
                        double tmpcost = greedyinfotmp.GetCost(dupcost, losscost);
                        if (currentbestcost <= tmpcost) {
                            hashtable = backuphash; // X
                            /*while (modify_remove.size() > 0) {
                                T = modify_remove.back();
                                modify_remove.pop_back();
                                hashtable[T.slbl].remove(T.dupHeight, T.n);
                            }
                            while (modify_add.size() > 0) {
                                T = modify_add.back();
                                modify_add.pop_back();
                                hashtable[T.slbl].add_cell(T.dupHeight, T.n);
                            }*/

                            //partialMapping = backuppartialMapping;
                            for (auto x : firstmap) {
                                partialMapping[x.first] = backuppartialMapping[x.first];
                                myfile << "3.firstmap: " << x.first->GetLabel() << " to " << x.second->GetLabel() << endl;
                            }
                            for (auto x : currentmap) {
                                partialMapping[x.first] = backuppartialMapping[x.first];
                                myfile << "3.currentmap: " << x.first->GetLabel() << " to " << x.second->GetLabel() << endl;
                            }
                            currentmap.clear();
                            remap.clear();
                            currenttotalnblosses = backupgreedyinfo.nbLosses;
                            *improve = false;
                            //newremmap = false;
                        }
                        else {
                            besthash = hashtable; // X
                            //bestpartialMapping = partialMapping;
                            for (auto x : remap) {
                                bestpartialMapping[x.first] = partialMapping[x.first];
                                remap1[x.first] = partialMapping[x.first];
                                myfile << "4.remap: " << x.first->GetLabel() << " to " << x.second->GetLabel() << endl;
                            }
                           for (auto x : partialMapping) {
                                if (partialMapping[x.first]->GetLabel() != bestpartialMapping[x.first]->GetLabel()) {
                                    bestpartialMapping[x.first] = partialMapping[x.first];
                                    myfile << "Error!" << endl;
                                    myfile << "partial: " << x.first->GetLabel() << "->" << partialMapping[x.first]->GetLabel() << " best: " << x.first->GetLabel() << "->" << bestpartialMapping[x.first]->GetLabel() << endl;
                                }
                            }
                            currentbestcost = tmpcost;
                            bestgreedyinfo.nbLosses = greedyinfotmp.nbLosses;
                            bestgreedyinfo.dupHeightSum = greedyinfotmp.dupHeightSum;
                            cout << "Good Remap!" << endl;
                            *improve = true;
                            newremmap = true;
                            hashtable = backuphash;
                            //partialMapping = backuppartialMapping;
                            for (auto x : firstmap) {
                                partialMapping[x.first] = backuppartialMapping[x.first];
                                myfile << "4.firstmap: " << x.first->GetLabel() << " to " << x.second->GetLabel() << endl;
                            }
                            for (auto x : currentmap) {
                                partialMapping[x.first] = backuppartialMapping[x.first];
                                myfile << "4.currentmap: " << x.first->GetLabel() << " to " << x.second->GetLabel() << endl;
                            }

                            currentmap.clear();
                            remap.clear();
                            currenttotalnblosses = backupgreedyinfo.nbLosses;
                        }

                    }
                    firstmap.clear();
                    currentmap.clear();
                    remap.clear();
                }
                firstmap.clear();
                currentmap.clear();
                remap.clear();
            }
            g->CloseIterator(it);
        }
        num++;
        backuphash = besthash;
        backuppartialMapping = bestpartialMapping;
        hashtable = besthash;
        partialMapping = bestpartialMapping;
        backupgreedyinfo.nbLosses = bestgreedyinfo.nbLosses;
        backupgreedyinfo.dupHeightSum = bestgreedyinfo.dupHeightSum;
        currentbestcost = bestgreedyinfo.GetCost(dupcost, losscost);
        cout << " number of remapping: " << num << endl;
    }
    /*cout << "After run of Greedy remapping : " << endl;
    for (int i = 0; i < hashtable.size(); i++) {
        cout << "Species " << i << " ";
        hashtable[i].print();
        //cout << "size " << hashtable[i].size() << endl;
    }
    cout << "/////////////////////////////////////////////////////////////////////////////////////" << endl;
    cout << " nb Losses of Greedy: " << bestgreedyinfo.nbLosses << " nb Dupheightsum of Greedy : " << bestgreedyinfo.dupHeightSum << endl;
    cout << " cost of Greedy : " << currentbestcost << endl;*/

    return bestgreedyinfo;
}


SegmentalReconcileInfo SegmentalReconcile::UltraGreedyRemapping(unordered_map<Node*, Node*>& partialMapping, vector<hashlist>& hashtable, vector<Node*>& geneTrees, Node* speciesTree, double LCAcost, double dupcost, double losscost, bool *improve)
{   
    vector<hashlist> backuphash = hashtable;
    unordered_map<Node*, Node*> backuppartialMapping = partialMapping;
    vector<Node*> minimalNodes;
    SegmentalReconcileInfo greedyinfo;
    greedyinfo.dupHeightSum = GetdupHeightSum(hashtable);
    greedyinfo.nbLosses = GetnbLosses(partialMapping);
    int currenttotalnblosses = greedyinfo.nbLosses;
    SegmentalReconcileInfo greedyinfotmp;
    double currentbestcost = LCAcost;
    int cntn = 0;
    bool isdup;
    int d1, d2, losses_tmp;
    *improve = false;

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
                cntn++;
                if (cntn % 500 == 0) {
                    cout << cntn << endl;
                }
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
                    isdup = IsDuplication(n, partialMapping);

                    //calulate number of losses for n 
                    d1 = GetSpeciesTreeDistance(backuppartialMapping[n], backuppartialMapping[n->GetChild(0)]);
                    d2 = GetSpeciesTreeDistance(backuppartialMapping[n], backuppartialMapping[n->GetChild(1)]);
                    losses_tmp = (double)(d1 + d2);

                    if (isdup) {
                        int slbl = Util::ToInt(currents->GetLabel());
                        int dupheight = GetDuplicationHeightUnder(n, currents, partialMapping);
                        //cout << "is dup" << endl;
                        hashtable[slbl].remove(dupheight, n);
                    }
                    else
                    {
                        losses_tmp -= 2;
                    }
                    currenttotalnblosses -= losses_tmp;
                    //cout << "before remapping" << endl;
                    partialMapping[n] = s;
                    //cout << "after remapping" << endl;
                    isdup = IsDuplication(n, partialMapping);
                    //calulate new number of losses for n 
                    d1 = GetSpeciesTreeDistance(partialMapping[n], partialMapping[n->GetChild(0)]);
                    d2 = GetSpeciesTreeDistance(partialMapping[n], partialMapping[n->GetChild(1)]);
                    losses_tmp = (double)(d1 + d2);

                    if (isdup) {
                        int dupheight = GetDuplicationHeightUnder(n, s, partialMapping);
                        //cout << "is new dup" << endl;
                        if (dupheight > 0) {
                            int slbl = Util::ToInt(s->GetLabel());
                            hashtable[slbl].add_cell(dupheight, n);
                        }
                    }
                    else
                    {
                        losses_tmp -= 2;
                    }
                    currenttotalnblosses += losses_tmp;
                    bool check_losses = true;
                    Node* last_ancestor = n;

                    if (!n->IsRoot()) {
                        int slbl = Util::ToInt(s->GetLabel());
                        int slbl_parent = Util::ToInt(partialMapping[n->GetParent()]->GetLabel());
                        //cout << "slbl: " << slbl << " slbl parent " << slbl_parent << endl;
                        last_ancestor = n->GetParent();
                        if (slbl >= slbl_parent)
                        {
                            minimalNodes.push_back(n->GetParent());
                        }
                        if (slbl == slbl_parent) {
                            check_losses = false;
                        }
                    }

                    while (minimalNodes.size() > 0) {
                        for (int j = minimalNodes.size() - 1; j >= 0; j--)
                        {
                            Node* m = minimalNodes[j];
                            Node* currents1 = partialMapping[m];

                            isdup = IsDuplication(m, backuppartialMapping);

                            //calulate number of losses for n 
                            d1 = GetSpeciesTreeDistance(backuppartialMapping[m], backuppartialMapping[m->GetChild(0)]);
                            d2 = GetSpeciesTreeDistance(backuppartialMapping[m], backuppartialMapping[m->GetChild(1)]);
                            losses_tmp = (double)(d1 + d2);
                            if (!isdup) {
                                losses_tmp -= 2;
                            }
                            currenttotalnblosses -= losses_tmp;

                            if (IsDuplication(m, partialMapping)) {
                                int slbl = Util::ToInt(currents1->GetLabel());
                                int dupheight = GetDuplicationHeightUnder(m, currents1, partialMapping);
                                //cout << "remove " << dupheight << ", " << m->GetLabel() << " from " << slbl << endl;
                                hashtable[slbl].remove(dupheight, m);
                                //cout << "removed " << dupheight << ", " << m->GetLabel() << " from " << slbl << endl;
                            }

                            
                            partialMapping.erase(m);
                            Node* s1 = GetLowestPossibleMapping(m, partialMapping); // should be same as s
                            //cout << "effected remap " << m->GetLabel() << " from " << currents1->GetLabel() << " to " << s1->GetLabel() << endl;
                            partialMapping[m] = s1;

                            isdup = IsDuplication(m, partialMapping);
                            //calulate new number of losses for n 
                            d1 = GetSpeciesTreeDistance(partialMapping[m], partialMapping[m->GetChild(0)]);
                            d2 = GetSpeciesTreeDistance(partialMapping[m], partialMapping[m->GetChild(1)]);
                            losses_tmp = (double)(d1 + d2);
                            if (!isdup) {
                                losses_tmp -= 2;
                            }
                            currenttotalnblosses += losses_tmp;

                            int dupheight = GetDuplicationHeightUnder(m, s1, partialMapping);
                            if (dupheight > 0) {
                                int slbl = Util::ToInt(s1->GetLabel());
                                hashtable[slbl].add_cell(dupheight, m);
                            }
                            minimalNodes.erase(minimalNodes.begin() + j);
                            if (!m->IsRoot()) {
                                int slbl = Util::ToInt(s1->GetLabel());
                                int slbl_parent = Util::ToInt(partialMapping[m->GetParent()]->GetLabel());
                                last_ancestor = m->GetParent();
                                if (slbl >= slbl_parent)
                                {
                                    //cout << slbl << ">" << slbl_parent << endl;
                                    //cout << " next is " << m->GetParent()->GetLabel() << endl;
                                    minimalNodes.push_back(m->GetParent());
                                }
                                if (slbl == slbl_parent) {
                                    check_losses = false;
                                }
                            }
                        }
                    }

                    if (check_losses) {

                        isdup = IsDuplication(last_ancestor, backuppartialMapping);
                        //calulate number of losses for n 
                        d1 = GetSpeciesTreeDistance(backuppartialMapping[last_ancestor], backuppartialMapping[last_ancestor->GetChild(0)]);
                        d2 = GetSpeciesTreeDistance(backuppartialMapping[last_ancestor], backuppartialMapping[last_ancestor->GetChild(1)]);
                        losses_tmp = (double)(d1 + d2);
                        if (!isdup)
                        {
                            losses_tmp -= 2;
                        }
                        currenttotalnblosses -= losses_tmp;

                        isdup = IsDuplication(last_ancestor, partialMapping);
                        //calulate new number of losses for n 
                        d1 = GetSpeciesTreeDistance(partialMapping[last_ancestor], partialMapping[last_ancestor->GetChild(0)]);
                        d2 = GetSpeciesTreeDistance(partialMapping[last_ancestor], partialMapping[last_ancestor->GetChild(1)]);
                        losses_tmp = (double)(d1 + d2);
                        if (!isdup) {
                            losses_tmp -= 2;
                        }
                        currenttotalnblosses += losses_tmp;
                    }

                    greedyinfotmp.nbLosses = currenttotalnblosses;
                    greedyinfotmp.dupHeightSum = GetdupHeightSum(hashtable);
                    double tmpcost = greedyinfotmp.GetCost(dupcost, losscost);
                    if (currentbestcost < tmpcost) {
                        hashtable = backuphash;
                        partialMapping = backuppartialMapping;
                        currenttotalnblosses = greedyinfo.nbLosses;
                    }
                    else {
                        backuphash = hashtable;
                        backuppartialMapping = partialMapping;
                        currentbestcost = tmpcost;
                        greedyinfo.nbLosses = greedyinfotmp.nbLosses;
                        greedyinfo.dupHeightSum = greedyinfotmp.dupHeightSum;
                        cout << "Good Remap!" << endl;
                        cout << "gene index: " << Util::ToString(n->GetIndex()) << " from " << Util::ToString(currents->GetIndex()) << " -> " << Util::ToString(s->GetIndex()) << endl;
                        *improve = true;
                    }

                }
            }

        }
        g->CloseIterator(it);
    }
    /*cout << "After run of Greedy remapping : " << endl;
    for (int i = 0; i < hashtable.size(); i++) {
        cout << "Species " << i << " ";
        hashtable[i].print();
        //cout << "size " << hashtable[i].size() << endl;
    }
    cout << "/////////////////////////////////////////////////////////////////////////////////////" << endl;
    cout << " nb Losses of Greedy: " << greedyinfo.nbLosses << " nb Dupheightsum of Greedy : " << greedyinfo.dupHeightSum << endl;
    cout << " cost of Greedy : " << currentbestcost << endl;*/

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

Node* SegmentalReconcile::GetSimphyMapping(Node* g, unordered_map<Node*, Node*>& partialMapping)
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
    char separator = '-';
    int i = 0;
    string glbl = g->GetLabel();
    bool flag = false;
    bool flag2 = false;
    string s = glbl;
    //cout << glbl << endl;
    /*while (glbl[i] != '\0') {
        if (glbl[i] != separator && !flag) {
            // Append the char to the temp string.
            s += glbl[i];
        }
        else {
            flag = true;
        }
        i++;
    }*/
    /*if (!flag) {
        return partialMapping[g->GetChild(0)]->FindLCAWith(partialMapping[g->GetChild(1)]);
    }*/
    //else {

    char delimiter = '_';

    // Create a stringstream from the input string
    std::stringstream ss(s);

    std::vector<std::string> tokens;
    std::string token;

    // Tokenize the input string by the delimiter and store each token in the vector
    while (std::getline(ss, token, delimiter)) {
        tokens.push_back(token);
    }
    string species;
    //cout << "mapped to " << endl;
    for (const auto& t : tokens) {
        //std::cout << t << std::endl;
        species = t;
    }
    //cout << species << endl;
    Node * look_species_ = partialMapping[g->GetChild(0)]->FindLCAWith(partialMapping[g->GetChild(1)]);

    while (!look_species_->IsRoot()) {

        string look_species = look_species_->GetLabel();
        // Find the position of the first occurrence of the delimiter
        size_t pos = look_species.find(separator);
        if (pos != std::string::npos) {
            // Extract the substring just before the delimiter
            look_species = look_species.substr(0, pos);
            //std::cout << "Result: " << look_species << std::endl;
        }

        if (species == look_species) {
            //cout << "simphy mapped to: "<< look_species << endl;
            return look_species_;
        }
        if(!look_species_->IsRoot())
            look_species_ = look_species_->GetParent();
    }
    //cout << "simphy mapped to: " << look_species_->GetLabel() << endl;
    return look_species_;
    //}

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

int SegmentalReconcile::GetnbLossesexceptgenetree(unordered_map<Node*, Node*>& fullMapping, int genetreenum)
{
    double cost = 0;
    int nbdups = 0;
    int nblosses = 0;

    for (int i = 0; i < geneTrees.size(); i++)
    {
        if (i != genetreenum) {
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

    }



    return nblosses;
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
