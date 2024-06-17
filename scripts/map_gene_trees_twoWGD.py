import re
from ete3 import Tree
import os,sys

gene_trees_file = sys.argv[1]
edited_gene_trees_file = sys.argv[2]
sim_out_directory = sys.argv[3]

def remove_branch_lengths(newick_file, output_file, genetree_id):
    # Read the Newick tree from the file
    tree = Tree(newick_file, format=1)
    id_counter = 0
    #print(tree)
    num_nodes = len(tree.search_nodes())
    num_nodes -= 1
    #print(num_nodes)
    root = tree.get_tree_root()
    root.name = str(num_nodes)
    # Clear branch lengths from the tree
    for node in tree.traverse():

        #print(node.name)
        section = re.split('_', node.name)
        if section.__len__() > 1:
            # print(section)
            flag1 = False
        else:
            species_name, event = look_sim_species(node.name, genetree_id)
            node.name = str(node.name) + "_" + species_name + "_" + event
            #print(node.name)
        node.dist = 0

    try:
        # Open the file in append mode
        #with open(output_file, 'a') as file:
        content = tree.write(format=1,  format_root_node=lambda node: f"[&name={node.name}]")
        modified_content = content[:-3] + ";"
        output_file.write(modified_content + "\n")
        #print("Lines appended successfully.")
    except Exception as e:
        print("Error appending lines to the file:", str(e))

    #file.close()
    # Write the modified tree to the output file
    # tree.write(outfile=output_file, format=1)

def setsimphyMapping():
    counter = 0
    gene_trees = open(gene_trees_file, 'r')
    edited_gene_trees = open(edited_gene_trees_file, 'w')
    lines = gene_trees.readlines()
    for line2 in lines:
        counter = counter + 1
        remove_branch_lengths(line2, edited_gene_trees, counter)

def look_sim_species(gene_id, genetree_id):

    if genetree_id < 10:
        genetree_id = "0" + str(genetree_id)
    elif genetree_id < 100:
        genetree_id = "" + str(genetree_id)

    name = str(genetree_id) + "l1g.maplg.modded.modded"
    #print(name)
    gene_tree_mapping = open(sim_out_directory + "/" + name, 'r')
    lines = gene_tree_mapping.readlines()
    locus_id = "-1"
    flag = False
    for line in lines:
        #print(line)
        words = line.split()
        if words and words[0] == gene_id:
            #print(words[0])
            locus_id = words[1]
            flag = True
    gene_tree_mapping.close()
    if flag == False:
        locus_id = gene_id
    # print(gene_id + " -> " + locus_id)

    name2 = str(genetree_id) + ".mapsl.modded.modded"
    # print(name2)
    gene_tree_mapping = open(sim_out_directory + "/" + name2, 'r')
    lines = gene_tree_mapping.readlines()
    species_id = "-1"
    flag = False
    event = "ND"
    for line in lines:
        words = line.split()
        if words and words[0] == locus_id:
            species_id = words[3]
            event = words[2]
            flag = True
    gene_tree_mapping.close()
    if flag == False:
        species_id = locus_id

    #translation_table = str.maketrans("", "", "'")
    #species_id = species_id.translate(translation_table)
    #print(locus_id + " -> " + species_id)
    return species_id, event


if __name__ == "__main__":

    setsimphyMapping()
    print("simphy mapping is computed!")