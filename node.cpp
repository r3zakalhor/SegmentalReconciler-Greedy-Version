#include "node.h"



//This is used for pruning species trees.  The Node::Restrict function takes a fct pointer as argument.
bool Node__RestrictToLeafsetFunction(Node* n, void* arg)
{
    set<Node*>* leavesToKeep = (set<Node*>*)arg;
    if (n->IsLeaf())
    {
        return (leavesToKeep->find(n) != leavesToKeep->end());
    }
    else
    {
        return (n->IsRoot() || n->GetNbChildren() > 1);
    }
}




Node::Node(bool maintainTreeInfo)
{
    this->parent = NULL;
    this->nodeInfo = NULL;
    this->pathBits = 0;
    this->branchLength = 0.0;
    this->depth = -1;

    if (maintainTreeInfo)
    {
        this->treeInfo = new TreeInfo(this);
        this->treeInfo->OnNodeInserted(this, 0);
    }
    else
    {
        this->treeInfo = NULL;
    }

    state = 0;
}


Node::Node(TreeInfo* treeInfo)
{
    this->depth = -1;
    this->parent = NULL;
    this->treeInfo = treeInfo;
    this->nodeInfo = NULL;
    this->state = 0;
    this->pathBits = 0;
    this->branchLength = 0.0;
}



Node::~Node()
{
    for (int i = 0; i < children.size(); i++)
    {
        delete children[i];
    }
    children.clear();

    if (nodeInfo)
        delete nodeInfo;

    if (this->parent == NULL && treeInfo)
        delete treeInfo;

}

Node* Node::AddChild()
{
    Node* n = this->CreateNode(treeInfo);
    this->AddSubTree(n);

    if (treeInfo)
        treeInfo->OnNodeInserted(n, this->children.size() - 1);

    return n;
}

Node* Node::CreateNode(TreeInfo* treeInfo)
{
    return new Node(treeInfo);
}

Node* Node::GetParent()
{
    return parent;
}

Node* Node::InsertChild(int pos)
{
    vector<Node*>::iterator it = children.begin();
    Node* n = this->CreateNode(treeInfo);
    children.insert(it + pos, n);
    n->SetParent(this);


    if (treeInfo)
        treeInfo->OnNodeInserted(n, pos);

    return n;
}

int Node::GetNbChildren()
{
    return children.size();
}

Node* Node::GetChild(int pos)
{
    if (pos < 0 || pos >= children.size())
        throw "Bad child position";

    return children[pos];
}

void Node::SetLabel(string lbl)
{
    string prevLabel = this->label;
    this->label = lbl;

    if (treeInfo)
        treeInfo->OnLabelChanged(this, prevLabel, lbl);
}

void Node::SetIndex(int index)
{
    this->index = index;
}

int Node::GetIndex()
{
    return this->index;
}

void Node::SetDup(bool dup)
{
    this->isdup = dup;
}

bool Node::IsDup()
{
    return isdup;
}


string Node::GetLabel()
{
    return label;
}



Node* Node::GetLeftSibling()
{
    if (this->IsRoot())
        return NULL;
    int pindex = -1;

    //this pretty much sucks, having to iterate through all parent's children to find "this"
    for (int i = 0; i < this->parent->GetNbChildren(); i++)
    {
        if (parent->GetChild(i) == this)
            pindex = i;
    }

    if (pindex == 0)
        return NULL;

    return parent->GetChild(pindex - 1);
}

Node* Node::GetRightSibling()
{
    if (this->IsRoot())
        return NULL;
    int pindex = -1;

    //this pretty much sucks, having to iterate through all parent's children to find "this"
    for (int i = 0; i < this->parent->GetNbChildren(); i++)
    {
        if (parent->GetChild(i) == this)
            pindex = i;
    }

    if (pindex == parent->GetNbChildren() - 1)
        return NULL;

    return parent->GetChild(pindex + 1);
}


TreeIterator* Node::GetPostOrderIterator(bool leavesOnly)
{
    TreeIterator* it = new PostOrderTreeIterator(this);
    it->SetLeavesOnly(leavesOnly);
    return it;
}


TreeIterator* Node::GetPreOrderIterator(bool leavesOnly)
{
    TreeIterator* it = new PreOrderTreeIterator(this);
    it->SetLeavesOnly(leavesOnly);
    return it;
}

void Node::CloseIterator(TreeIterator* it)
{
    delete it;
}

TreeInfo* Node::GetTreeInfo()
{
    return treeInfo;
}


void Node::SetDepth(int depth)
{
    this->depth = depth;
}

int Node::GetDepth()
{
    if (this->depth < 0)
    {
        int cpt = 0;
        Node* n = this;
        while (n->parent)
        {
            cpt++;
            n = n->parent;
        }
        return cpt;
    }
    else
        return this->depth;
}

void Node::SetPathBits(uint64 pathBits)
{
    this->pathBits = pathBits;
}

uint64 Node::GetPathBits()
{
    return pathBits;
}


bool Node::IsRoot()
{
    return (this->parent == NULL);
}

bool Node::IsLeaf()
{
    return (this->children.size() == 0);
}


void Node::SetParent(Node* parent)
{
    //TODO : treeInfo might need to be updated
    //if (!this->parent)
    //{
    //    this->parent = parent;
    //}

    //TODO : parent must have added node to its children beforehand

    this->parent = parent;
}


void Node::AddSubTree(Node* n)
{
    //if (!n->IsRoot())
    //    return;

    if (n->treeInfo || this->treeInfo)
    {
        throw "AddSubTree not supported with treeInfo";
    }

    this->children.push_back(n);
    n->SetParent(this);

}


void Node::DeleteTreeInfo()
{
    for (int i = 0; i < this->GetNbChildren(); i++)
    {
        this->GetChild(i)->DeleteTreeInfo();
    }

    if (this->treeInfo)
    {
        if (this->IsRoot())
            delete treeInfo;

        this->treeInfo = NULL;
    }




}


void Node::RemoveChild(Node* node)
{
    vector<Node*>::iterator it = this->children.begin();

    while (it != this->children.end())
    {
        if ((*it) == node)
        {
            it = children.erase(it);

            if (treeInfo)
                treeInfo->OnNodeDeleted();
        }
        else
        {
            it++;
        }
    }

}


void Node::RemoveChildren(bool setParentToNullOnChildren)
{
    if (setParentToNullOnChildren)
    {
        for (int i = 0; i < children.size(); i++)
        {
            this->children[i]->SetParent(NULL);
        }
    }
    children.clear();


    if (treeInfo)
        treeInfo->OnNodeDeleted();
}


int Node::GetState()
{
    return state;
}

void Node::SetState(int state)
{
    this->state = state;
}

void Node::CopyFrom(Node* n, set<Node*> ignoreNodes)
{
    this->depth = n->depth;
    this->label = n->label;
    this->pathBits = n->pathBits;
    this->state = n->state;
    this->branchLength = n->GetBranchLength();

    if (n->nodeInfo)
    {
        this->nodeInfo = n->nodeInfo->GetClone();
    }

    for (int i = 0; i < n->children.size(); i++)
    {
        if (ignoreNodes.find(n->children[i]) == ignoreNodes.end())
        {
            Node* nc = this->AddChild();
            nc->CopyFrom(n->children[i], ignoreNodes);
        }
    }
}



void Node::Restrict(bool (*fncDelete)(Node*, void*), void* arg)
{
    TreeIterator* it = this->GetPostOrderIterator();

    Node* n = it->next();
    while (n)
    {
        bool ok = (*fncDelete)(n, arg);

        if (!ok)
        {
            n = it->DeleteCurrent();
        }
        else
        {
            n = it->next();
        }
    }

    this->CloseIterator(it);

}


void Node::DeleteSingleChildDescendants()
{
    TreeIterator* it = this->GetPostOrderIterator();

    Node* n = it->next();
    while (n)
    {
        if (n->GetNbChildren() == 1)
        {
            n = it->DeleteCurrent();
        }
        else
        {
            n = it->next();
        }
    }
    this->CloseIterator(it);
}


Node* Node::FindLCAWith(Node* n)
{
    //NAIVE WAY: list all nodes from this to root
    //then list all from n to root, return first met that was visited previously
    unordered_set<Node*> visited;

    Node* cur = this;
    visited.insert(cur);

    while (cur)
    {
        cur = cur->parent;
        visited.insert(cur);
    }

    cur = n;
    while (cur)
    {
        if (visited.find(cur) != visited.end())
            return cur;
        cur = cur->parent;
    }

    //this return should never happen in theory
    return NULL;
}


Node* Node::FindLCA(vector<Node*> nodes)
{
    if (nodes.size() == 0)
        return NULL;

    Node* lca = nodes[0];
    for (int i = 1; i < nodes.size(); i++)
    {
        lca = lca->FindLCAWith(nodes[i]);
    }

    return lca;
}


/*void Node::SetMappingLabel(string lbl)
{
    this->mappingLabel = lbl;
}

string Node::GetMappingLabel()
{
    return this->mappingLabel;
}*/


bool Node::HasAncestor(Node* ancestor)
{

    Node* cur = this;
    while (cur)
    {
        if (cur == ancestor)
            return true;
        else
            cur = cur->GetParent();
    }

    return false;
}




void Node::SetBranchLength(double length)
{
    branchLength = length;
}

double Node::GetBranchLength()
{
    return branchLength;
}




Node* Node::GetNodeWithLabel(string lbl, bool ignoreCase)
{
    TreeIterator* it = this->GetPostOrderIterator();

    Node* found = NULL;
    while (Node* n = it->next())
    {
        if (ignoreCase)
        {
            if (Util::ToUpper(lbl) == Util::ToUpper(n->GetLabel()))
            {
                found = n;
                break;
            }
        }
        else
        {
            if (lbl == n->GetLabel())
            {
                found = n;
                break;
            }
        }
    }
    this->CloseIterator(it);

    return found;
}


set<Node*> Node::GetLeafSet()
{
    set<Node*> leaves;
    TreeIterator* it = this->GetPostOrderIterator(true);

    while (Node* n = it->next())
    {
        leaves.insert(n);
    }
    this->CloseIterator(it);

    return leaves;

}


vector<Node*> Node::GetLeafVector()
{
    vector<Node*> leaves;
    TreeIterator* it = this->GetPostOrderIterator(true);

    while (Node* n = it->next())
    {
        leaves.push_back(n);
    }
    this->CloseIterator(it);

    return leaves;

}



Node* Node::InsertParentWith(Node* sibling)
{
    Node* prevParent = this->parent;

    if (!prevParent)
        return NULL;

    int sindex = -1;
    for (int i = 0; i < prevParent->GetNbChildren(); i++)
    {
        if (prevParent->GetChild(i) == sibling)
        {
            sindex = i;
            break;
        }
    }

    if (sindex == -1)
        return NULL;

    Node* newParent = prevParent->AddChild();

    prevParent->RemoveChild(this);
    prevParent->RemoveChild(sibling);

    newParent->AddSubTree(this);
    newParent->AddSubTree(sibling);

    return newParent;

}



void Node::BinarizeRandomly()
{
    TreeIterator* it = this->GetPostOrderIterator();
    while (Node* n = it->next())
    {
        if (n != this)
        {
            n->BinarizeRandomly();
        }
    }
    this->CloseIterator(it);

    while (this->GetNbChildren() > 2)
    {

        int ic1 = rand() % this->GetNbChildren();
        int ic2 = rand() % this->GetNbChildren();

        if (ic1 != ic2)
        {
            Node* c1 = this->GetChild(ic1);
            Node* c2 = this->GetChild(ic2);

            c1->InsertParentWith(c2);
        }
    }
}

vector<Node*> Node::GetChildrenVector()
{
    vector<Node*> v(this->children);

    return v;
}

vector<Node*> Node::GetPostOrderedNodes()
{
    vector<Node*> v;
    TreeIterator* it = this->GetPostOrderIterator();
    while (Node* n = it->next())
    {
        v.push_back(n);
    }
    this->CloseIterator(it);

    return v;
}


Node* Node::SetAsRootInCopy()
{
    Node* res = this->SetAsRootInCopy(NULL);
    res->DeleteSingleChildDescendants();
    return res;
}



Node* Node::SetAsRootInCopy(Node* ignore)
{
    Node* copy = new Node(false);
    set<Node*> ignoreSet;
    if (ignore)
        ignoreSet.insert(ignore);
    copy->CopyFrom(this, ignoreSet);

    if (this->parent)
    {
        Node* parentCopy = this->parent->SetAsRootInCopy(this);
        copy->AddSubTree(parentCopy);
    }

    return copy;
}


Node* Node::SetRootOnParentEdgeInCopy()
{
    if (!this->parent)
        throw "Node has no parent edge.";

    Node* newRoot = new Node(false);

    Node* thisCopy = new Node(false);
    thisCopy->CopyFrom(this);

    Node* parentCopy = this->parent->SetAsRootInCopy(this);

    newRoot->AddSubTree(thisCopy);
    newRoot->AddSubTree(parentCopy);

    return newRoot;
}


void Node::SetCustomField(string name, string val)
{
    this->customFields[name] = val;
}

string Node::GetCustomField(string name)
{
    if (this->customFields.find(name) == this->customFields.end())
        return "";

    return this->customFields[name];
}


Node* Node::GraftOnParentEdge(Node* nodeToGraft)
{
    if (this->IsRoot())
        throw "Can graft on root";
    Node* parent = this->GetParent();

    parent->RemoveChild(this);
    Node* n = parent->AddChild();
    n->AddSubTree(this);
    n->AddSubTree(nodeToGraft);

    return n;
}


int Node::GetNbLeaves()
{
    vector<Node*> v = GetLeafVector();

    return v.size();
}

Node* Node::GetLeafByLabel(string label)
{

    vector<Node*> v = this->GetLeafVector();
    for (int i = 0; i < v.size(); i++)
    {
        if (v[i]->GetLabel() == label)
            return v[i];
    }

    return NULL;
}


set<string> Node::GetLeafLabels()
{
    set<string> labels;
    TreeIterator* it = this->GetPostOrderIterator(true);

    while (Node* n = it->next())
    {
        labels.insert(n->GetLabel());
    }
    this->CloseIterator(it);

    return labels;
}


//static
void Node::RestrictToLeafset(Node* root, set<Node*> leavesToKeep)
{
    root->Restrict(&Node__RestrictToLeafsetFunction, (void*)&leavesToKeep);
}
