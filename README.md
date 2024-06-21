# SegmentalReconciler-Greedy-Version

#########
COMPILING
#########

You can compile with the classic cmake - make approach, e.g. to build the exe in a build directory
mkdir build
cd build
cmake ..
make


This will produce the 
segrec 
executable



* If you'd really like to compile manually, you could do:
g++ -std=c++11 main.cpp define.h genespeciestreeutil.h hashtable.h newicklex.h node.h SegmentalReconciler.h treeinfo.h treeiterator.h util.h genespeciestreeutil.cpp newicklex.cpp node.cpp SegmentalReconciler.cpp treeinfo.cpp treeiterator.cpp  -o segrec
(possibly add optimization flags)

#########
USING
#########

Here is an example call 

./segrec -d 5 -l 1 -gf ../data/gene_trees.txt -sf ../data/s_tree.newick -spsep "_" -spindex 0 -o output.txt -al greedy 


The argument -al specifies the type of reconciliation which has 3 options:

- simphy: Consider the original mapping from Simphy simulator.
- LCA: Considers the LCA mapping.
- anything else: Does greedy remapping


### Simphy

To run a simphy simulation: https://github.com/adamallo/SimPhy/wiki/Manual

#### Postorder labeling

After simulating Simphy, we need to do a postorder labeling of the species tree with the following command: (where you can find the Python code in the supplementaries)

 - $python post-order-labeling.py $SpeciesTree $Output_SpeciesTree

Then, we need to map gene trees nodes based on the Simphy mapping by following command: (where you can find the Python code in the supplementaries)

- $python map_gene_trees.py $GeneTrees $Output_GeneTrees $Simphy_simulation_directory

#### Now, we can run the segmental reconciler for four different algorithms and produce the appropriate outputs.

### Calculate distance between mappings

Finally, by using the segmental reconciler outputs we can calculate the distance between any two mapping by fallowing command: (where you can find the Python code in the supplementaries)

- $python compare_mapping.py $firstmapping $secondmapping $SpeciesTree $Output

### All in one script:

You can find all these steps in a bash script that might be useful. (supplementaries/simphy_simulation.sh)

### Generate number of duplication per species figure

To generate plots of nb duplication and species, you first need to use a script to calculate number of duplications per species and save the results at CSV files using output_nbdup_csv.sh bash script. (where you can find it in the supplementaries) 
Then use the pyhton code to generate plots for all of the simulations. (create_figures.py)
