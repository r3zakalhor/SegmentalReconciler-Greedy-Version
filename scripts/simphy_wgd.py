from ete3 import Tree
import random
import os,sys

species_tree_file = "s_tree.newick"


max_leaf_id = 0

def copy_species_subtree( node ):
    global max_leaf_id
    
    ncopy = Tree()
    
    ncopy.dist = node.dist
    
    if node.is_leaf():
        ncopy.add_feature("origin_leaf_id", node.name.replace("'", ""))
        #ncopy.name = "'" + str(max_leaf_id + 1) + "'"
        ncopy.name = str(max_leaf_id + 1)
        max_leaf_id += 1
    else:
        
        orig_preid = getattr(node, "preorderid")
        ncopy.add_feature("origin_preorderid", orig_preid)
        
        orig_postid = getattr(node, "postorderid")
        ncopy.add_feature("origin_postorderid", orig_postid)
        
        for c in node.children:
            child_copy = copy_species_subtree( c )
            ncopy.add_child(child_copy)
            
    
    return ncopy 




#input = list of lines of a mapsl file 
#output = dict, where key = spid, value = max val in a line of the form 'spid_val' 
def get_max_loci_per_species(lines):
    
    retdict = {}
    
    for line in lines:
        pz = line.split("\t")
        
        if len(pz) >= 4 and pz[0].startswith("'") and not pz[0].startswith("'Lost"):
            pk = pz[0].replace("'", "").split("_")
            spid = pk[0]
            val = int(pk[1])
            
            if spid in retdict:
                retdict[spid] = max( retdict[spid], val )
            else:
                retdict[spid] = val
    
    return retdict



def main(args):
    global max_leaf_id
    
    print("Reading " + args["-spfile"])
    species_tree = Tree(args["-spfile"], format=1)
    
    #format 5 = leaf names + all branch lengths, no internal node name (required for simphy)
    #dist_formatter also required
    oldnewick = species_tree.write(format=5, dist_formatter="%0.10f")
    
    print(oldnewick)
    
    cptnode = 0
    
    num_leaves = len([leaf for leaf in species_tree.iter_leaves()])
    node_to_dup = None
    candidates = []
    
    for node in species_tree.traverse("preorder"):

        if not node.is_leaf():
            node.add_feature("preorderid", cptnode)
            

            #TODO: subtree to dup is first we see with good number of leaves.  Maybe pick randomly instead
            if node_to_dup is None and len(node) <= num_leaves-2 and len(node) >= num_leaves/2:
                candidates.append(node)
                #node_to_dup = node
            
        
        else: 
            leafid = (int(node.name.replace("'", "")))
            max_leaf_id = max(max_leaf_id, leafid)
        cptnode += 1
    
    
    if not candidates:
        print("No node to duplicate")
        sys.exit()
        
    node_to_dup = random.choice(candidates)
    cptnode = 0
    for node in species_tree.traverse("postorder"):
        if not node.is_leaf():
            node.add_feature("postorderid", cptnode)
        cptnode += 1

    if node_to_dup is None:
        print("No node to duplicate")
        sys.exit()
        
    if node_to_dup.is_root():
        print("Node to dup is root, not supported, too many special cases to handle")
        sys.exit()
            
            
    node_copy = copy_species_subtree( node_to_dup )
    
    #node_to_dup is removed from its parent
    prev_parent = node_to_dup.up
    
    d = node_to_dup.dist 
    
    
    node_to_dup.detach()
    
    node_new_parent = Tree()
    
    node_new_parent.add_feature("isextra", True)
    
    
    node_new_parent.add_child(node_to_dup)
    node_new_parent.add_child(node_copy)
    prev_parent.add_child(node_new_parent)
    
    #branch above node is split in two halves
    node_to_dup.dist = float(d)/2.0
    node_new_parent.dist = float(d)/2.0
    node_copy.dist = float(d)/2.0
    
    
    newnewick = species_tree.write(format=5, dist_formatter="%0.10f")
    
    
    #find new preorder ids
    cptnode = 1
    
    #list of (p, q) where p=newpreorderid q=prevpreorderid
    preorder_maps_str = ""
    
    #list of (p, q) where p=newpostorderid q=prevpostorderid
    postorder_maps_str = ""
    
    #(p, q) where p=newpreorderid q=old preorderid of node it corresponds to
    newinternal_preordermaps_str = ""
    
    #(p, q) where p=newpostorderid q=old postorderid of node it corresponds to
    newinternal_postordermaps_str = ""
    
    extranode_preorder_id = ""
    
    extranode_postorder_id = ""
    
    
    duped_node_preorder_id = str(getattr(node_to_dup, "preorderid"))
    
    duped_node_postorder_id = str(getattr(node_to_dup, "postorderid"))
    
    
    
    
    
    #(p, q) where p=new leaf id   q = corresponding leaf id.  The ids contain apostrophes, e.g. '20'
    newleaf_maps_str = ""
    
    cptnode = 0
    for node in species_tree.traverse("preorder"):
        if not node.is_leaf():
            node.add_feature("newpreorderid", str(cptnode))
             
            if hasattr(node, "preorderid"):
                preorder_maps_str += " (" + getattr(node, "newpreorderid") + "," + str(getattr(node, "preorderid")) + ")"
            
            if hasattr(node, "isextra"):
                extranode_preorder_id = getattr(node, "newpreorderid")
            
            if hasattr(node, "origin_preorderid"):
                newinternal_preordermaps_str += " (" + getattr(node, "newpreorderid") + "," + str(getattr(node, "origin_preorderid")) + ")"
            
        
        else:
            if hasattr(node, "origin_leaf_id"):
                newleaf_maps_str += " (" + node.name + "," + str(getattr(node, "origin_leaf_id")) + ")"
        cptnode += 1
    cptnode = 0
    for node in species_tree.traverse("postorder"):
        if not node.is_leaf():
            node.add_feature("newpostorderid", str(cptnode))
            
            if hasattr(node, "postorderid"):
                postorder_maps_str += " (" + getattr(node, "newpostorderid") + "," + str(getattr(node, "postorderid")) + ")"
            
            if hasattr(node, "isextra"):
                extranode_postorder_id = getattr(node, "newpostorderid")
            
            if hasattr(node, "origin_postorderid"):
                newinternal_postordermaps_str += " (" + getattr(node, "newpostorderid") + "," + str(getattr(node, "origin_postorderid")) + ")"
            
            
        cptnode += 1
    
    #to output
    #oldnewick
    #newpreorderid -> oldpreorderid 
    #extranode preorderid 
    #newpreorderid -> original node preorder id 
    #newleaf -> oldleaf 
    strout = ""
    strout += "oldnewick=" + oldnewick + "\n"
    strout += "newnewick=" + newnewick + "\n"
    strout += "preorderid_map=" + preorder_maps_str + "\n"
    strout += "postorderid_map=" + postorder_maps_str + "\n"
    strout += "extranode_preorderid=" + extranode_preorder_id + "\n"
    strout += "extranode_postorderid=" + extranode_postorder_id + "\n"
    strout += "duped_node_preorder_id=" + duped_node_preorder_id + "\n"
    strout += "duped_node_postorder_id=" + duped_node_postorder_id + "\n"
    strout += "newinternal_preordermaps_str=" + newinternal_preordermaps_str + "\n"
    strout += "newinternal_postordermaps_str=" + newinternal_postordermaps_str + "\n"
    strout += "newleaf_map=" + newleaf_maps_str + "\n"
    
    
    if "-sout" in args:
        f = open(args["-sout"], "w")
        f.write(newnewick)
        f.close()
        
    if "-auxfile" in args:
        f = open(args["-auxfile"], "w")
        f.write(strout)
        f.close()
    else:
        print(strout)
    



#parses a string of the form (a, b) (c, d)
#couples separated by space
#returns { "a" => "b", "c" => "d" }
def get_couples_dict( str_with_couples ):
    
    ret = dict()
    couples = str_with_couples.split(" ")
    
    
    for c in couples:
        if "," in c:
            pz = c.split(",")
            ret[ pz[0].replace("(", "").replace(")", "").strip() ] = pz[1].replace("(", "").replace(")", "").strip()
    
    return ret 



def remap_simphy_file( auxfilename, mapsl_filename, maplg_filename, genetree_filename, out_prefix = 'wgd_' ):
    
    f_aux = open(auxfilename, 'r')
    lines = f_aux.readlines()
    f_aux.close()
    
    params = dict()
    for line in lines:
        line = line.replace("\n", "").strip()

        if "=" in line:
            pz = line.split("=", 1)
            params[ pz[0] ] = pz[1]
    
    
    output_str = ""
    
    
    #todo: check for errors, eg missing lines 
    
    new_to_old_leafids = get_couples_dict( params["newleaf_map"] )
    
    #old nodes: newpostorderid => orig 
    new_to_old_postorderids = get_couples_dict( params["postorderid_map"] )
    
    #newnodes: newpostorderid => orig 
    newnode_corresp_postorders = get_couples_dict( params["newinternal_postordermaps_str"] )
    
    
    extranode_postorderid = params["extranode_postorderid"]
    
    #params["duped_node_postorder_id"] = "10"
    
    duped_node_postorder_id = ""
    if "duped_node_postorder_id" in params:
        duped_node_postorder_id = params["duped_node_postorder_id"]
        #little trick: tell remapper to remap extranode_postorderid to duped_node_postorder_id
        newnode_corresp_postorders[ extranode_postorderid ] = duped_node_postorder_id
    else:
        print("Error: duped_node_postorder_id not set in aux file")
    
    #print(params)
    #print(new_to_old_leafids)
    #print(new_to_old_postorderids)
    #print(newnode_corresp_postorders)
    
    f_sl = open(mapsl_filename, 'r')
    lines = f_sl.readlines()
    f_sl.close()
    
    
    loci_per_species = get_max_loci_per_species(lines)
    loci_to_new_loci = {}   #eg loci_to_new_loci['23_1'] = '15_2'
    
    '''
    This modifes each line by remapping whatever needs to be remapped
    Leaf case:
    if leaf '39' is new and corresponds to original leaf '10', lines of the form 
    '39_0'    0    Leaf    '39'
    become 
    '10_0'    0    Leaf    '10'
    
    A leaf could also be a loss, form 
    'Lost-45_0'    0    Loss    '20'
    here only the last '20' needs replacement
    OR
    'Lost-45_0'    0    Loss    20
    where now 20 is a postorder map
    
    Internal node case:
    40  0   Sp  26
    Here 26 is the species postorder id.  This one needs to be replaced with corresponding sp id.
    
    '''
    output_str2 = []
    for line in lines:
        
        newline = line.replace("\n", "").strip()
        
        changed = False
        
        pz = newline.split("\t")
        
        
        event = pz[2]
        newevent = event
        
        
        sp = pz[3]  #species that the locus node is mapped to 
        newsp = sp 
        spid = sp.replace("'", "")
        
        #first check pz[3] = species -> updates newsp
        #species leaf case
        if sp.startswith("'"):
            if spid in new_to_old_leafids:
                newsp = "'" + new_to_old_leafids[spid] + "'"
        #species internal node -> postorderid case
        else:
            if sp in newnode_corresp_postorders:
                newsp = newnode_corresp_postorders[sp]
            elif sp in new_to_old_postorderids:
                newsp = new_to_old_postorderids[sp]
            else:
                print("ERROR: species " + sp + " not a postorderid of newnodes nor of oldnodes")
                
                
            #special case: things mapped to extra node must be Dups
            if sp == extranode_postorderid:
                newevent = "Dup"
        
        #then check pz[0] = locus -> updates newloc
        loc = pz[0]
        newloc = loc
        
        para = pz[1]
        newpara = para
        #one case left: locus leaf of the form '39_0'   0   Sp  '39'
        if loc.startswith("'") and not loc.startswith("'Lost"):
            if not sp.startswith("'"):
                print(f"ERROR: leaf {pz[0]} mapped to species internal node {sp}")
            elif not loc.startswith("'" + sp.replace("'", "") + "_"):
                print(f"ERROR: leaf {pz[0]} mapped to incorrect species leaf {sp}")
            else:
                if sp != newsp:
                    spname = newsp.replace("'", "")
                    if spname not in loci_per_species:
                        loci_per_species[spname] = -1	#will become 0
					
                    maxval_loci = loci_per_species[spname]
                    newloc = "'" + spname + "_" + str(maxval_loci + 1) + "'"
                    
                    
                    
                    newpara = str(maxval_loci + 1)
                    loci_per_species[spname] += 1
                    
                    #newloc = loc.replace("'" + sp.replace("'", "") + "_", "'" + newsp.replace("'", "") + "_")

        loci_to_new_loci[loc] = newloc
        
        if sp != newsp or loc != newloc or event != newevent:
            #kids today would use pz.join("\t") or something, whatever
            newline = newloc + "\t" + newpara + "\t" + newevent + "\t" + newsp 
            changed = True
        
        if len(newline) > 0:
            output_str2.append(newline)
        
        
        if changed:
            print(line.replace("\n", "") + "   -->")
            print(newline)
    
    #print(output_str)
    
    output_str = "\n".join(output_str2)
    with open(mapsl_filename + ".modded", 'w') as f:
        print(output_str, file=f)




    print("Now modding lg file.  Note: will not work if using ILS!")

    #and finally the maplg file
    #lines of the form 
    #'39_0_0'    '39_0'    0
    #must be changed to 
    #'10_0_0'    '10_0'    0
    #Also we need to use updated loci vals
    output_str_lg = ""
    f_lg = open(maplg_filename, 'r')
    lines = f_lg.readlines()
    f_lg.close()

    for line in lines:
        
        newline = line.replace("\n", "").strip()
        changed = False
        
        if newline.startswith("'"):
            pz = newline.split("\t")
            
            if pz[1] in loci_to_new_loci and pz[1] != loci_to_new_loci[pz[1]]:
                newloc = loci_to_new_loci[pz[1]]
                newline = "'" + newloc.replace("'", "") + "_0'" + "\t" + newloc + "\t0"
                
                changed = True

    
        if changed:
            print(line.replace("\n", "") + "   -->")
            print(newline)
    
        output_str_lg += newline + "\n"

    with open(maplg_filename + ".modded", 'w') as f:
        print(output_str_lg, file=f)

    
    
    # Open the Newick tree file
    with open(genetree_filename, "r") as file:
        newick_string = file.read()

    # Create a Tree object from the Newick string
    tree = Tree(newick_string)
    #print (new_to_old_leafids)
    # Traverse the tree and visit each leaf node
    for leaf in tree.iter_leaves():

        words = leaf.name.split('_')
        flag=False
        if words[0] in new_to_old_leafids:
            flag=True
            #print(new_to_old_leafids[words[0]])
            words[0] = new_to_old_leafids[words[0]]
        updated_string = "_".join(words)
        if flag:
            print("Leaf name: ", leaf.name)
            print("updated Leaf name: ", updated_string)
            leaf.name = updated_string
        
        
    tree.write(outfile=genetree_filename + ".modded")

    # Open the file in append mode with newline parameter set to an empty string
    with open(genetree_filename + ".modded", 'a', newline='') as file:
        # Add a newline character to the end of the file
        file.write('\n')



if __name__ == "__main__":
    
    
    
    
    if "-h" in sys.argv or "--help" in sys.argv:
        print("This program has two modes.  The default is take an input species tree file, in newick as generated by simphy, and the program duplicates a subtree of it and outputs a modified species tree (usable with simphy) and an auxfile.")
        print("With mode=remap, the program takes an auxfile and a mapsl file, and modifies it to revert to original species.")
        print("Default mode: arguments are\n -sin=[path to input sp tree] -sout=[path to output modded sp tree] -auxfile=[path to output aux data]")
        print("If -spfile is not specified, modded sp tree is not output.  If -auxfile not specified, we just output it to stdout.")
        print("Example: python simphy_wgd.py -spfile=mysptree.txt -auxfile=myaux.txt")
        print("mode=remap: outputs modded mapsl file to same mapsl.modded")
        print("python simphy_wgd.py -auxfile=myaux.txt -slfile=1.mapsl -lgfile=1.maplg")
        sys.exit()
    
    args = dict(arg.split('=') for arg in sys.argv[1:])
    
    if "-mode" in args and args["-mode"] == "remap":
        aux_filename = args["-auxfile"]
        mapsl_filename = args["-slfile"]
        maplg_filename = args["-lgfile"]
        genetree_filename = args["-genetreefile"]
        remap_simphy_file( aux_filename, mapsl_filename, maplg_filename, genetree_filename )
    else:
        if not "-spfile" in args:
            print("You need to give an input sp tree with -spfile")
            sys.exit()
        main(args)
        
        
        
        
        
        
        
        
        

