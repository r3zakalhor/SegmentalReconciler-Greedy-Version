import re
from ete3 import Tree
import os,sys


output = sys.argv[4]
species_tree_file = sys.argv[3]
filename = sys.argv[1]
filename2 = sys.argv[2]

#output = "output.txt"
#species_tree_file = "s_tree.newick"
#filename = "out_ultragreedy.txt"
#filename2 = "out_simphy.txt"

f = open(output, 'w')
sys.stdout = f


if __name__ == "__main__":

    list_of_distances = []
    total_distance = 0

    file1 = open(species_tree_file, 'r')
    lines = file1.readlines()
    species_tree = lines[0]
    # print(species_tree)
    s_tree = Tree(species_tree, format=1)
    for node in s_tree.traverse():
        node.dist = 1
        if node.is_leaf():
            node.name = node.name + "_0"
    output_file = 'tree_without_branch_lengths.newick'
    s_tree.write(format=1, outfile=output_file, format_root_node=True)
    # B = s_tree.search_nodes(name="51")[0]
    root_node = s_tree.get_tree_root()
    # print(root_node)
    # s_tree.show()
    # s_tree.render("tree.png", w=600, units="px")
    file1.close()
    # filename2 = 'out_simphy.txt'
    file2 = open(filename2, 'r')
    lines = file2.readlines()
    species_tree = ""

    print("Cost and Number of segmental duplication of " + filename2 + ": ")
    for i in range(9):
        print(lines[i], end="")
    file2.close()

    file = open(filename, 'r')
    lines = file.readlines()
    flag = False
    s_flag = False
    species_tree = ""
    count = 0
    print("Cost and Number of segmental duplication of " + filename + ": ")
    for i in range(9):
        print(lines[i], end="")
    file.close()
    for line in lines:
        if "</GENETREES>" in line:
            flag = False
            # print(line)

        if flag:
            # reading each word
            count = count + 1
            for word in re.split('[( )]', line):
                # displaying the words
                # print(word)
                count1 = 0
                section = re.split('_', word)
                if section.__len__() > 3:
                    # print(section)
                    sim_map = section[1]
                    algorithm_map = section[3]
                    # print(section)
                    # print(sim_map[0])
                    # print(greedy_map[0])
                    # sim_species = look_sim_species(sim_map[0], count)
                    # print(sim_species + ", " + greedy_map[0])
                    distance = 0
                    if sim_map != algorithm_map:
                        # print (sim_map)
                        # print (algorithm_map)
                        translation_table = str.maketrans("", "", "'")
                        #algorithm_map = algorithm_map.translate(translation_table)
                        #sim_map = sim_map.translate(translation_table)
                        # A = s_tree.search_nodes(name=sim_map)
                        #if len(A) > 0:
                            # A = A[0]
                        # B = s_tree.search_nodes(name=algorithm_map)
                        # if len(B) > 0:
                           # B = B[0]
                        # A = s_tree.get_descendants_by_name(sim_species)[0]
                        # B = s_tree.get_descendants_by_name(greedy_map[0])[0]
                        nodes_with_algorihtm_name = s_tree.search_nodes(name=algorithm_map)
                        if nodes_with_algorihtm_name:
                            algorithm_map = algorithm_map
                        else:
                            algorithm_map = algorithm_map + "_0"
                            # print(algorithm_map)
                        nodes_with_algorihtm_name.clear()
                        nodes_with_sim_name = s_tree.search_nodes(name=sim_map)
                        if nodes_with_sim_name:
                            sim_map = sim_map
                        else:
                            sim_map = sim_map + "_0"
                            #print(sim_map)
                        nodes_with_sim_name.clear()
                        distance = s_tree.get_distance(sim_map, algorithm_map)
                        # if len(A) > 0 and len(B) > 0:
                            # dist_A = s_tree.get_distance(A, s_tree.get_tree_root())
                            # dist_B = s_tree.get_distance(B, s_tree.get_tree_root())
                            # distance = abs(dist_A - dist_B)
                            # distance = s_tree.get_distance(A, B)
                    # genetree number, gene tree node, simphy mapping, greedy mapping, distance
                    element = [count, section[0], sim_map, algorithm_map, distance]
                    list_of_distances.append(element)
                    total_distance = total_distance + distance
                    # print("distance: " + str(distance))
                    # compute_distance(sim_species, greedy_map[0])
                # for sec in section:
                # count1 = count1 + 1
                # if count1 > 3:
                # print(word)

        if "<GENETREES>" in line:
            flag = True
            # print(line)

        # print(line)

    print("total distance between mappings: ")
    print(str(total_distance))
    print("number of gene trees: " + str(count))

    for element in list_of_distances:
        if element[4] != 0:
            print("Element:", element)
            print("gene tree number: ", element[0])
            print("gene tree node: ", element[1])
            print("simphy mapping", element[2])
            print("greedy mapping", element[3])
            print("distance", element[4])
            print("------------------------")

    print("total distance between mappings: " + str(total_distance))
