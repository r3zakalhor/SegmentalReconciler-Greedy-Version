#include <iostream>

#include <map>
#include "util.h"
#include "newicklex.h"
#include "node.h"
#include "genespeciestreeutil.h"
#include "treeiterator.h"
#include "SegmentalReconciler.h"


#include "ReconciliationTester.h"

using namespace std;





//sets the preorderid as index for every node of SpeciesTree, then returns the number of nodes
int SetSpeciesIndex(Node* speciesTree) {
    int cpt = 0;
    TreeIterator* it = speciesTree->GetPostOrderIterator();
    while (Node* g = it->next())
    {
        g->SetIndex(cpt);
        cpt++;
    }
    return cpt;
}

int GetNbLeaves(Node* speciesTree) {
    int cpt = 0;
    TreeIterator* it = speciesTree->GetPostOrderIterator();
    while (Node* g = it->next())
    {
        if(g->IsLeaf())
            cpt++;
    }
    return cpt;
}

int SetGenesIndex(vector<Node*> geneTrees) {
    int cpt = 0;
    for (int i = 0; i < geneTrees.size(); i++)
    {
        Node* genetree = geneTrees[i];
        TreeIterator* it = genetree->GetPreOrderIterator();
        while (Node* g = it->next())
        {
            g->SetIndex(cpt);
            /*if (!g->IsLeaf()) {
                g->SetLabel(Util::ToString(cpt));
            }*/
            cpt++;
        }
    }
    return cpt;
}




//labels the gene trees to prepare them for output.  Adds the species mapping to the out,
//plus _Spec or _Dup_nbx, where x is a dup id.  Also returns a map of dups per species,
//since we're computing it in this function anyway.  The value is a vector of int/Node pairs,
//where the int is the gene tree index and the node is a dup node in this tree.
//I agree that this function does more than just labelling with species mapping, but hey, life is tough.
map<Node*, vector< pair<int, Node*> > > LabelGeneTreesWithSpeciesMapping(vector<Node*> geneTrees, Node* speciesTree, SegmentalReconciler& reconciler, SegmentalReconcilerInfo& info, bool resetLabels = true)
{
    map<Node*, vector< pair<int, Node*> > > dups_per_species;
    int dup_counter = 1;
    for (int i = 0; i < geneTrees.size(); i++)
    {
        Node* genetree = geneTrees[i];

        TreeIterator* it = genetree->GetPostOrderIterator();
        while (Node* g = it->next())
        {
            if (!g->IsLeaf())
            {
                string lbl = g->GetLabel();
                if (lbl != "")
                    lbl += "_";

                if (resetLabels)
                    lbl = "";

                lbl += info.curmap[g]->GetLabel();

                if (reconciler.IsDuplication(g, info.curmap))
                {
                    lbl += "_Dup_nb" + Util::ToString(dup_counter);
                    dup_counter++;

                    vector< pair<int, Node*> > dups_for_s;
                    if (dups_per_species.find(info.curmap[g]) != dups_per_species.end())
                        dups_for_s = dups_per_species[info.curmap[g]];
                    pair<int, Node*> p = make_pair(i + 1, g);
                    dups_for_s.push_back(p);
                    dups_per_species[info.curmap[g]] = dups_for_s;
                }
                else
                    lbl += "_Spec";
                g->SetLabel(lbl);
            }
        }
        genetree->CloseIterator(it);
    }

    return dups_per_species;
}



GSMap GetGeneSpeciesMapping(vector<GNode*> geneTrees, SNode* speciesTree, string species_separator, int species_index)
{
    //compute gene leaf to species leaf mapping.  Could be faster by using a map for the species leaves by name...
    GSMap geneSpeciesMapping;

    for (int i = 0; i < geneTrees.size(); i++)
    {
        GSMap tmpmap = GeneSpeciesTreeUtil::Instance()->GetGeneSpeciesMappingByLabel(geneTrees[i], speciesTree, species_separator, species_index);

        for (GSMap::iterator it = tmpmap.begin(); it != tmpmap.end(); ++it)
        {
            geneSpeciesMapping[it->first] = it->second;
        }
    }

    return geneSpeciesMapping;
}



void PrintHelp()
{
	cout<<"Please update me"<<endl;
	return;
}





SegmentalReconcilerInfo Execute(map<string, string> &args) {
    SegmentalReconcilerInfo info;
    
    vector<GNode*> geneTrees;
    SNode* speciesTree = NULL;


    string algorithm;
    string species_separator = "_";
    int species_index = 0;
    double dupcost = 2;
    double losscost = 1;
    double stochastic_temperature = 1;
	int nbstochasticLoops = 2000;

    int max_remap_distance = 999999;
	
    //parse dup loss cost and max dup height
    if (args.find("d") != args.end()){
        dupcost = Util::ToDouble(args["d"]);
    }
    if (args.find("l") != args.end()){
        losscost = Util::ToDouble(args["l"]);
    }
    if (args.find("al") != args.end()){
        algorithm = args["al"];
    }
	if (args.find("tmp") != args.end()){
        stochastic_temperature = Util::ToDouble(args["tmp"]);
    }
	if (args.find("stoloops") != args.end()){
        nbstochasticLoops = Util::ToInt(args["stoloops"]);
    }
    if (args.find("maxremap") != args.end()) {
        max_remap_distance = Util::ToInt(args["maxremap"]);
        cout << "Max remap distance set to " << max_remap_distance << endl;
    }
    string outfile = "";
    if (args.find("o") != args.end()){
        outfile = args["o"];
    }
	
	for (auto const& x : args)
    {
		std::cout << x.first  // string (key)
              << ':' 
              << x.second // string's value 
              << std::endl;
	}
	
    //parse gene trees, either from command line or from file
    if (args.find("g") != args.end()){
        vector<string> gstrs = Util::Split(Util::ReplaceAll(args["g"], "\n", ""), ";", false);

        for (int i = 0; i < gstrs.size(); i++){
            string str = gstrs[i];
            GNode* tree = NewickLex::ParseNewickString(str, false);

            if (!tree){
                cout << "Error: there is a problem with input gene tree " << str << endl;
                return info;
            }

            geneTrees.push_back(tree);
        }
    }
    else if (args.find("gf") != args.end()){
        string gcontent = Util::GetFileContent(args["gf"]);

        vector<string> lines = Util::Split(gcontent, "\n");
        vector<string> gstrs;
        for (int l = 0; l < lines.size(); l++){
            vector<string> trees_on_line = Util::Split(lines[l], ";", false);
            for (int t = 0; t < trees_on_line.size(); t++){
                if (trees_on_line[t] != "")
                    gstrs.push_back(trees_on_line[t]);
            }
        }


        for (int i = 0; i < gstrs.size(); i++){
            string str = gstrs[i];
            GNode* tree = NewickLex::ParseNewickString(str, false);

            if (!tree)
            {
                cout << "Error: there is a problem with input gene tree " << str << endl;
                return info;
            }

            geneTrees.push_back(tree);
        }
    }


    //parse species trees, either from command line or from file
    if (args.find("s") != args.end()){
        speciesTree = NewickLex::ParseNewickString(args["s"], false);

        if (!speciesTree){
            cout << "Error: there is a problem with the species tree." << endl;
            return info;
        }
    }
    else if (args.find("sf") != args.end()){
        string scontent = Util::GetFileContent(args["sf"]);
        speciesTree = NewickLex::ParseNewickString(scontent, false);

        if (!speciesTree)
        {
            cout << "Error: there is a problem with the species tree." << endl;
            return info;
        }
    }




    //parse species separator and index
    if (args.find("spsep") != args.end()){
        species_separator = args["spsep"];
    }

    if (args.find("spindex") != args.end()){
        species_index = Util::ToInt(args["spindex"]);
    }




    if (geneTrees.size() == 0){
        cout << "No gene tree given.  Program will exit." << endl;
        return info;
    }
    if (!speciesTree){
        cout << "No species tree given.  Program will exit." << endl;
        return info;
    }

    if (dupcost < 0 || losscost <= 0){
        cout << "dupcost < 0 or losscost <= 0 are prohibited.  Program will exit." << endl;
        return info;
    }
    //cout << algorithm << endl;
    string all[6] = { "lca", "greedy", "ultragreedy", "simphy", "fastgreedy", "stochastic"};
    bool flagx = false;
    for (int i = 0; i < 6; i++) {
        if (algorithm == all[i]) {
            flagx = true;
        }
    }
    if (!flagx) {
        cout << "Algorithm not set to an existing type, will use default greedy." << endl;
    }


    unordered_map<Node*, Node*> geneSpeciesMapping = GetGeneSpeciesMapping(geneTrees, speciesTree, species_separator, species_index);
    int nbspecies = SetSpeciesIndex(speciesTree);
    int numleaves = GetNbLeaves(speciesTree);
    int nbgenes = SetGenesIndex(geneTrees);
    
    int numintnodes = nbspecies;
    nbspecies = nbspecies + numleaves + 1;
	bool stochastic = false;
	if (algorithm == "stochastic"){
		stochastic = true;
	}
    SegmentalReconciler reconciler(geneTrees, speciesTree, geneSpeciesMapping, dupcost, losscost, nbspecies, nbgenes, numintnodes, stochastic, stochastic_temperature, nbstochasticLoops);
    cout << "dupcost: " << dupcost << " losscost: " << losscost << endl;
    
    reconciler.SetMaxRemapDistance(max_remap_distance);

	if (algorithm == "lca"){
		info = reconciler.ReconcileWithLCAMap();
	}
	else if (algorithm == "simphy"){
		info = reconciler.ReconcileWithSimphy();
	}
	else{
		info = reconciler.Reconcile();
	}


    string output = "";
    if (info.isBad)
    {
        output = "NO SOLUTION FOUND";
    }
    else
    {
        output += "<COST>\n" + Util::ToString(info.GetCost(dupcost, losscost)) + "\n</COST>\n";
        output += "<DUPHEIGHT>\n" + Util::ToString(info.dupHeightSum) + "\n</DUPHEIGHT>\n";
        output += "<NBLOSSES>\n" + Util::ToString(info.nbLosses) + "\n</NBLOSSES>\n";
        output += "<SPECIESTREE>\n" + NewickLex::ToNewickString(speciesTree) + "\n</SPECIESTREE>\n";

        map<Node*, vector< pair<int, Node*> > > dups_per_species = LabelGeneTreesWithSpeciesMapping(geneTrees, speciesTree, reconciler, info, false);

        output += "<GENETREES>\n";
        for (int t = 0; t < geneTrees.size(); t++)
        {
            output += NewickLex::ToNewickString(geneTrees[t]) + "\n";
        }
        output += "</GENETREES>\n";

        output += "<DUPS_PER_SPECIES>\n";
        TreeIterator* itsp = speciesTree->GetPostOrderIterator();
        while (Node* s = itsp->next())
        {
            if (dups_per_species.find(s) != dups_per_species.end())
            {
                output += "[" + s->GetLabel() + "] ";
                vector< pair<int, Node*> > dups_for_s = dups_per_species[s];

                for (int d = 0; d < dups_for_s.size(); d++)
                {
                    pair<int, Node*> p = dups_for_s[d];
                    string lbl = p.second->GetLabel();
                    lbl = Util::GetSubstringAfter(lbl, "_");

                    output += lbl + " (G" + Util::ToString(p.first) + ") ";

                }
                output += "\n";
            }
        }
        speciesTree->CloseIterator(itsp);
        output += "</DUPS_PER_SPECIES>\n";
        
        
    }

    if (outfile == "")
        cout << output;
    else
        Util::WriteFileContent(outfile, output);




    for (int i = 0; i < geneTrees.size(); i++)
    {
        delete geneTrees[i];
    }

    delete speciesTree;

    return info;
}
















int main(int argc, char* argv[])
{

    map<string, string> args;

    bool hasHelp = false;
    bool hasTest = false;

    //BUILD DICTIONARY OF ARGS
    string prevArg = "";
    for (int i = 0; i < argc; i++){
        if (string(argv[i]) == "-v"){
            prevArg = "";
        }
        else if (string(argv[i]) == "--help"){
            hasHelp = true;
        }
        else{
            if (prevArg != "" && prevArg[0] == '-'){
                args[Util::ReplaceAll(prevArg, "-", "")] = string(argv[i]);
            }

            prevArg = string(argv[i]);
        }
    }

    
    if (args.find("help") != args.end() || hasHelp)
    {
        PrintHelp();
    }
    else if (args.count("testout")) {
        string filename = args["testout"];
        //filename = "C:\\cygwin64\\home\\Manuel\\git\\SegmentalReconciler-Greedy-Version\\data\\sim_4\\out_stochastic100.txt";
        ReconciliationTester tester;
        tester.ReadSegrecFile(filename);
    }
    else
    {
        //ML's ad-hoc testing stuff for Windows.
        /*args["d"] = "10";
        args["l"] = "1";
        args["gf"] = "sample_data/all_genetrees_edited.txt";
        args["sf"] = "sample_data/s_tree.newick";
        args["o"] = "sample_data/out_greedy.txt";
        args["al"] = "stochastic";*/

        /*
        args["d"] = "5";
        args["l"] = "1";
        args["gf"] = "data/gene_trees.txt";
        args["sf"] = "data/s_tree.newick";
        //args["o"] = "sample_data/out_greedy.txt";
        args["al"] = "fastgreedy";

        Example call:
        ./segrec -d 5 -l 1 -gf ../data/gene_trees.txt -sf ../data/s_tree.newick -spsep "_" -spindex 0 -al fastgreedy

        from build dir:
        ./segrec -sf "../data/sim_4/s_tree.newick" -gf "../data/sim_4/applied_loss_fix_all_genetrees_edited.txt" -d 10 -l 1 -al fastgreedy -maxremap 5
        */

        /*args["d"] = "100";
        args["l"] = "1";
        args["gf"] = "data/sim_1/all_genetrees.txt";
        args["gf"] = "data/sim_1/applied_loss_fix_all_genetrees_edited.txt";
        args["sf"] = "data/sim_1/s_tree.newick";
        args["al"] = "greedy";
        args["maxremap"] = "2";*/
        


        time_t start, end;
        time(&start);


        
        Execute(args);

        

        // Recording end time.
        time(&end);

        // Calculating total time taken by the program.
        double time_taken = double(end - start);
        cout << "Time taken by program is : " << fixed
            << time_taken ;
        cout << " sec " << endl;
        
        return 0;
    }
}

