from ete3 import Tree

def generate_recphyloxml(species_tree, gene_trees_file_path, output_filename):
    def recursive_write_SP(node, indent, output_file):
        output_file.write(f"{indent}<clade>\n")
        output_file.write(f"{indent}<name>{node.name}</name>\n")
        for child in node.children:
            recursive_write_SP(child, indent + '  ', output_file)
        output_file.write(f"{indent}</clade>\n")

    def loop_gene_trees(gene_trees_file_path, xml_file):
        with open(gene_trees_file_path, 'r') as file:
            for line in file:
                gene_tree_string = line.strip()
                gene_tree = Tree(gene_tree_string, format=1)
                xml_file.write("<recGeneTree>\n")
                xml_file.write("<phylogeny rooted=\"true\">\n")
                recursive_write_GT(gene_tree, '  ', xml_file)
                xml_file.write("</phylogeny>\n")
                xml_file.write("</recGeneTree>\n")

    def recursive_write_GT(node, indent, output_file):
        output_file.write(f"{indent}<clade>\n")
        output_file.write(f"{indent}<name>{node.name}</name>\n")
        output_file.write(f"{indent}<eventsRec>\n")
        components = node.name.split('_')
        if node.is_leaf():
            event = "leaf"
            output_file.write(f"{indent}<{event} speciesLocation=\"\'{components[0]}\'\"/>\n")
        else:
            if components[4] == "Spec":
                event = "speciation"
            elif components[4] == "Dup":
                event = "duplication"
            output_file.write(f"{indent}<{event} speciesLocation=\"{components[3]}\"/>\n")

        output_file.write(f"{indent}</eventsRec>\n")
        for child in node.children:
            recursive_write_GT(child, indent + '  ', output_file)
        output_file.write(f"{indent}</clade>\n")

    with open(output_filename, 'w') as xml_file:
        xml_file.write("<recPhylo>\n")
        xml_file.write("<spTree>\n")
        xml_file.write("<phylogeny>\n")
        recursive_write_SP(species_tree, '  ', xml_file)
        xml_file.write("</phylogeny>\n")
        xml_file.write("</spTree>\n")
        loop_gene_trees(gene_trees_file_path, xml_file)
        xml_file.write("</recPhylo>\n")

if __name__ == "__main__":
    species_tree_path = "species_tree.txt"
    gene_trees_path = "gene_trees.txt"
    species_tree = Tree(species_tree_path, format=1)
    generate_recphyloxml(species_tree, gene_trees_path, "output.xml")
