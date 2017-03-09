#include "../include/DisjointSet.h"
#include <stack>
#include <iostream>


using namespace std;

DisjointSet::DisjointSet(int s)
{
    size =s;
    parent=new int[s];
    for(int i=0;i<size;++i)parent[i]=-1;
}
int DisjointSet::Find(int x)
{
    if(parent[x]<0)return x;            //根则直接返回
    return parent[x]=Find(parent[x]);   //非根则找父节点，同时路径压缩
}

void DisjointSet::Union(int root1,int root2)
{
    if(root1==root2)return;
    if(parent[root1]>parent[root2]) //root1规模小（因为是负数） 把root1作为子树
    {
        parent[root2] += parent[root1];
        parent[root1] = root2;
    }else{
        parent[root1] += parent[root2];
        parent[root2] = root1;
    }
}

int DisjointSet::Find2(int x)
{
    stack<int> s;
    while(parent[x]>=0)
    {
        s.push(x);
        x = parent[x];
    }
    // x 是根节点
    while(!s.empty())
    {
        int tmp = s.top();s.pop();
        parent[tmp] = x;
    }
    return x;
}
void DisjointSet::Display()const
{
    for(int i=0;i<size;i++)
    {
        cout << parent[i]<<' ';
    }
    cout << endl;
}
int DisjointSet::RootSize()const
{
	int cnt=0;
	for(int i=0;i<size;i++)
		if(parent[i]<0)cnt++;
	return cnt;
}
