# SegmentalReconciler-Greedy-Version

##### First of all, to compile segmental reconciler we can use the following command:

- $g++ -std=c++11 main.cpp define.h genespeciestreeutil.h hashtable.h multigenereconciler.h newicklex.h node.h SegmentalReconcile.h treeinfo.h treeiterator.h util.h genespeciestreeutil.cpp multigenereconciler.cpp newicklex.cpp node.cpp SegmentalReconcil.cpp treeinfo.cpp treeiterator.cpp  -o segmentalreconciler

##### Then, we would be able to use ./segmentalreconciler by following command:

- $./segmentalreconciler -d $dupcost -l $losscost -gf $GeneTreeFile -sf $SpeciesTreeFile -o $OutputFile -al $algorithm

Where $algorithm specify the type of reconciliation which has 5 options:

- simphy: Consider the original mapping from Simphy simulator.
- LCA: Considers the LCA mapping.
- ultragreedy: It goes through all the nodes and as soon as it finds a remap that reduces the cost, it applies it.
- fastgreedy: It goes through all nodes and finds the best remap (has the lowest cost between all possible remapping for all nodes), then applies it and again iterates over all nodes until it finds a remap that improves the cost. (this algorithm is based on the changes which fast)
- greedy: It works like fastgreedy but is slower.

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

  
