#include "ReconciliationTester.h"
#include "newicklex.h"





int ReconciliationTester::GetDupHeightUnder(GNode* g, SNode* s, GSMap& map) {

    if (!IsDup(g, map))
        return 0;

    if (map[g] != s)
        return 0;

    return max(GetDupHeightUnder(g->GetChild(0), s, map), GetDupHeightUnder(g->GetChild(1), s, map)) + 1;
}




int ReconciliationTester::GetSpeciesDistance(SNode* s, SNode* t) {
    int dist = 0;

    while (s != t) {
        if (s->IsRoot())
            throw "t is not an ancestor of s";
        s = s->GetParent();
        dist++;
    }

    return dist;
}


int ReconciliationTester::GetNbLossesOnParentBranch(GNode* g, GSMap& map) {
    if (g->IsRoot())
        return 0;

    int dist = GetSpeciesDistance(map[g], map[g->GetParent()]);

    if (!IsDup(g->GetParent(), map))
        dist -= 1;

    return dist;
}





map<SNode*, int> ReconciliationTester::GetDupHeights(SNode* species_tree, vector<GNode*>& gene_trees, GSMap& curmap) {

    map<SNode*, int> dhs;

    for (int i = 0; i < gene_trees.size(); ++i) {
        GNode* g_root = gene_trees[i];

        TreeIterator* it = g_root->GetPostOrderIterator();
        while (GNode* g = it->next()) {
            if (!g->IsLeaf()) {
                SNode* mu = curmap[g];
                int dh = GetDupHeightUnder(g, mu, curmap);

                if (!dhs.count(mu))
                    dhs[mu] = dh;
                else
                    dhs[mu] = max(dhs[mu], dh);
            }
        }
        g_root->CloseIterator(it);
    }

    

    return dhs;
}




bool ReconciliationTester::IsDup(GNode* g, GSMap& map) {
    if (g->IsLeaf())
        return false;

    SNode* sl = map[g->GetChild(0)];
    SNode* sr = map[g->GetChild(1)];
    SNode* lca = sl->FindLCAWith(sr);

    if (lca == sl || lca == sr)
        return true;

    if (lca == map[g])
        return false;

    return true;
}




int ReconciliationTester::GetNbLosses(SNode* species_tree, vector<GNode*>& gene_trees, GSMap& curmap) {
    
    int nblosses = 0;
    for (int i = 0; i < gene_trees.size(); ++i) {
        GNode* g_root = gene_trees[i];

        TreeIterator* it = g_root->GetPostOrderIterator();
        while (GNode* g = it->next()) {
            nblosses += GetNbLossesOnParentBranch(g, curmap);
        }
        g_root->CloseIterator(it);
    }
    return nblosses;
}



bool ReconciliationTester::IsValidMap(SNode* species_tree, vector<GNode*> &gene_trees, GSMap &map, GSMap &lcamap) {
    for (int i = 0; i < gene_trees.size(); ++i) {
        GNode* g_root = gene_trees[i];

        TreeIterator* it = g_root->GetPostOrderIterator();
        while (GNode* g = it->next()) {
            if (g->IsLeaf()) {
                if (!map[g]->IsLeaf()) {
                    cout << g->GetLabel() + " not mapped to a leaf" << endl;
                    g_root->CloseIterator(it);
                    return false;
                }
            }
            else {
                
                SNode* sl = map[g->GetChild(0)];
                SNode* sr = map[g->GetChild(1)];
                SNode* lca = sl->FindLCAWith(sr);
                SNode* mu = map[g];

                if (mu != lca && !lca->HasAncestor(mu)) {
                    cout << g->GetLabel() + " mapped to " + mu->GetLabel() + " which is below lca" << endl;
                    g_root->CloseIterator(it);
                    return false;
                }

            }
        }
        g_root->CloseIterator(it);
    }
    return true;
}



GNode* ReconciliationTester::FindGeneTreeNode(GNode* g_root, string name) {
    GNode* found_node = nullptr;
    TreeIterator* it = g_root->GetPostOrderIterator();
    while (GNode* g = it->next()) {
        string lbl = g->GetLabel();
        int pos = lbl.find_first_of('_');
        string lbl2 = lbl.substr(pos + 1);

        if (lbl2 == name){
            found_node = g;
            break;
        }
    }
    g_root->CloseIterator(it);

    return found_node;
}




GSMap ReconciliationTester::GetLcaMap(vector<GNode*> gene_trees, SNode* species_tree) {
    GSMap map;
    
    for (int i = 0; i < gene_trees.size(); ++i) {
        GNode* g_root = gene_trees[i];

        TreeIterator* it = g_root->GetPostOrderIterator();
        while (GNode* g = it->next()){
            if (g->IsLeaf()) {
                vector<string> pz = Util::Split(g->GetLabel(), "_");
                string slbl = "'" + pz[0] + "'";
                SNode* s = species_tree->GetLeafByLabel(slbl);

                if (!s) {
                    throw "Species " + slbl + " not found";
                }

                map[g] = s;
            }
            else {
                map[g] = map[g->GetChild(0)]->FindLCAWith(map[g->GetChild(1)]);
            }
        }
        g_root->CloseIterator(it);
    }

    return map;
}


GSMap ReconciliationTester::ReadSegrecFile(string filename)
{
    

    double cost = 0, dupheight = 0, nblosses = 0;
    SNode* species_tree = nullptr;
    vector<GNode*> gene_trees;
    vector<string> lines = Util::GetFileLines(filename);

    GSMap lcamap;


    vector<string> dups_per_species;

    string mode = "";

    for (string line : lines) {

        if (line[0] == '<' && line[1] == '/')
            mode = "";
        else if (line[0] == '<')
            mode = line;
        else if (line != ""){
            if (mode == "<COST>") {
                cost = Util::ToDouble(line);
            }
            else if (mode == "<DUPHEIGHT>") {
                dupheight = Util::ToDouble(line);
            }
            else if (mode == "<NBLOSSES>") {
                nblosses = Util::ToDouble(line);
            }
            else if (mode == "<SPECIESTREE>") {
                species_tree = NewickLex::ParseNewickString(line);
            }
            else if (mode == "<GENETREES>") {
                GNode* g = NewickLex::ParseNewickString(line);
                gene_trees.push_back(g);
                
            }
            else if (mode == "<DUPS_PER_SPECIES>") {
                dups_per_species.push_back(line);
            }
        }
        
    }



    if (!species_tree) {
        cout << "No species tree found    *FAIL*" << endl;
        return lcamap;
    }


    lcamap = GetLcaMap(gene_trees, species_tree);

    GSMap newmap = lcamap;

    for (string dps : dups_per_species) {
        vector<string> pz = Util::Split(dps, " ", false);

        string sp = Util::ReplaceAll(pz[0], "[", "");
        sp = Util::ReplaceAll(sp, "]", "");

        SNode* spnode = species_tree->GetNodeWithLabel(sp);

        for (int i = 1; i < pz.size(); i += 2) {
            string gname = pz[i];
            string gstrindex = Util::ReplaceAll(Util::ReplaceAll(pz[i + 1], "(", ""), ")", "");
            gstrindex = Util::ReplaceAll(gstrindex, "G", "");
            int gindex = Util::ToInt(gstrindex);

            GNode* g = FindGeneTreeNode(gene_trees[gindex - 1], gname);

            if (!g) {
                throw "Label " + gname + " not found in " + pz[i + 1];
            }

            newmap[g] = spnode;
        }
    }



    bool validmap = IsValidMap(species_tree, gene_trees, newmap, lcamap);
    cout << "validmap = " << (validmap ? 1 : 0) << " ";


    map<SNode*, int> dhs = GetDupHeights(species_tree, gene_trees, newmap);
    int dhsum = 0;

    for (auto& keyval : dhs) {
        dhsum += keyval.second;
    }

    cout << "mydh=" << dhsum << "  filedh=" << dupheight << " ";


    int my_nblosses = GetNbLosses(species_tree, gene_trees, newmap);

    cout << "mylosses=" << my_nblosses << "  filedh=" << nblosses << " ";

    if (validmap && dhsum == dupheight && my_nblosses == nblosses) {
        cout << "   *PASS*" << endl;
    }
    else
    {
        cout << "   *FAIL*" << endl;
    }




    return newmap;
}
