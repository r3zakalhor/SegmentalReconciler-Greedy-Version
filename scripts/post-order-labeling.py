import sys
from ete3 import Tree

# Step 1: Read the Newick tree from the file
#file_path = "s_tree.trees"
# print(sys.argv[1])
# print(sys.argv[2])
file_path = sys.argv[1]
output = sys.argv[2]
tree = Tree(file_path)

root = tree.get_tree_root()
# print(root.name)
# Step 2: Perform postorder traversal and assign new IDs
def postorder_traversal(node, id_counter):
    if node.is_leaf():
        # node.name = str(id_counter[0])
        node.name = "'" + node.name + "'"
        id_counter[0] += 1
    else:
        for child in node.children:
            postorder_traversal(child, id_counter)
        node.name = str(id_counter[0])
        id_counter[0] += 1

id_counter = [0]  # Create a list to hold the counter as a mutable object
postorder_traversal(tree, id_counter)

id_counter[0] -= 1

# Step 3: Print the modified Newick tree with updated IDs
# print(tree.write(format=1, format_root_node=lambda node: f"[&name={node.name}]"))
tree.write(format=1, outfile=output, format_root_node=lambda node: f"[&name={node.name}]")

file_path = output
with open(file_path, "r") as file:
    content = file.read()

# Step 2: Modify the content to remove the last two characters and replace them with one character
modified_content = content[:-3] + ";"  # Replace "X" with your desired single character
# print(modified_content)

# Step 3: Write the modified content back to the file
with open(file_path, "w") as file:
    file.write(modified_content)

print("species tree is on post-order!")