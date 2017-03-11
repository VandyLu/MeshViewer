#ifndef DISJOINTSET_H_INCLUDED
#define DISJOINTSET_H_INCLUDED
#include <vector>
class DisjointSet{
private:
    int size;
    int*parent;
public:
    DisjointSet(int s);
    ~DisjointSet(){delete[]parent;}
    void Union(int root1,int root2);
    int Find(int x);
    int Find2(int x);
    void Display()const;
    std::vector<int> Types(int &n);
};


#endif // DISJOINTSET_H_INCLUDED
