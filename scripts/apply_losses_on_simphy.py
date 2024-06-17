from ete3 import Tree
import random
import os, sys

out_name = sys.argv[1]
gene_trees = sys.argv[2]
numerator = int(sys.argv[3])
denominator = int(sys.argv[4])
numerator_fix = int(sys.argv[5])
denominator_fix = int(sys.argv[6])


# species_tree_file = sys.argv[1]
# num_gene_trees = int(sys.argv[2])
# duprate = int(sys.argv[3])
# numerator = int(sys.argv[4])
# denominator = int(sys.argv[5])
# numerator_fix = int(sys.argv[6])
# denominator_fix = int(sys.argv[7])

def rename_leaves(tree):
    for leaf in tree.iter_leaves():
        #print(leaf.name)
        node_name_parts = leaf.name.split('_')
        leaf.name = f"{node_name_parts[0]}_0_0"  # You can modify the naming convention here
        #print(leaf.name)


def count_leaves(tree):
    leaf_count = {}
    for leaf in tree.iter_leaves():
        leaf_name = leaf.name
        leaf_count[leaf_name] = leaf_count.get(leaf_name, 0) + 1
    return leaf_count


def do_loss(gene_tree, leaf_name, pr, numerator, numerator_fix, denominator, denominator_fix):
    current_pr = pr - numerator
    leaves = list(gene_tree.iter_leaves())
    random.shuffle(leaves)
    for leaf in leaves:
        if leaf_name == leaf.name:
            random_pr = random.random()
            if random_pr <= current_pr / pr:
                leaf.detach()
                if current_pr != 0 and numerator_fix == 0:
                    current_pr = current_pr - numerator
                if pr != 0 and denominator_fix == 0:
                    pr = pr - denominator




def remove_single_child_internal_nodes(tree):
    for node in tree.traverse("postorder"):
        if node.is_leaf():
            node_name_parts = node.name.split('_')
            if node_name_parts[2]!="0":
                parent = node.up
                if parent:
                    if len(parent.children) == 2:
                        node.detach()
                    else:
                        node.detach()
                        leng = len(parent.children)
                        while leng == 1:
                            if parent:
                                node = parent
                                parent = node.up
                                leng = len(parent.children)
                                node.detach()
                            else:
                                node.detach()
                                leng = 0
                else:
                    node.detach()
    for node in tree.traverse("postorder"):
        if not node.is_leaf() and len(node.children) == 1:
            child = node.children[0]
            parent = node.up
            # Remove the internal node
            node.detach()
            # Link the parent to the child
            if parent:
                parent.add_child(child)
            else:
                for node1 in tree.traverse("postorder"):
                    node1.detach()
                tree.name = ""
                #print(tree.write())

def generate_new_gene_trees(set_of_gene_trees):
    gene_trees = []
    for tree in set_of_gene_trees:

        #root = gene_tree.get_tree_root()
        rename_leaves(tree)
        leaf_counts = count_leaves(tree)
        for leaf_name, count in leaf_counts.items():
            do_loss(tree, leaf_name, count, numerator, numerator_fix, denominator, denominator_fix)
            # print(f"{leaf_name}: {count}")
        remove_single_child_internal_nodes(tree)
        gene_trees.append(tree)
    return gene_trees


def main():
    # species_tree = Tree(species_tree_file, format=1)
    # gene_trees = "all_genetrees_edited.txt"
    Set_of_gene_trees = []
    with open(gene_trees, 'r') as file:
        for line in file:
            tree_str = line.strip()  # Remove leading/trailing whitespace
            tree = Tree(line, format=1)
            Set_of_gene_trees.append(tree)
            #print(tree.write())

    gene_trees_out = generate_new_gene_trees(Set_of_gene_trees)

    with open(out_name, "w") as output_file:
        for i, gene_tree in enumerate(gene_trees_out, start=1):
            num_nodes = len(gene_tree.get_tree_root().get_descendants()) + 1
            if num_nodes > 1:
                output_file.write(gene_tree.write(format=1, format_root_node=lambda node: f"[&name={node.name}]") + "\n")
            else:
                output_file.write("\n")

if __name__ == "__main__":
    main()
