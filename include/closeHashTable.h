#ifndef CLOSEHASHTABLE_H_INCLUDED
#define CLOSEHASHTABLE_H_INCLUDED

#include <iostream>
#include <vector>

template<class Type>
class hashTable{
public:
    virtual bool find(const Type&x)const=0;
    virtual bool insert(const Type&x)=0;
    virtual bool remove(const Type&x)=0;
protected:
    int (*key)(const Type &x);
    static int defaultKey(const int&k){return k;}
};

template<class Type>
class closeHashTable:public hashTable<Type>{
private:
    struct node{
        Type data;
        int state;
        node(){state =0;}
    };
    node*array;
    int size;
    int dataSize;
    int deleteSize;
    int maxDeleteSize;

    void autoRehash();
    void doubleSpace();
public:
    closeHashTable(int length=101,int _maxDeleteSize=10,int(*f)(const Type&x)=hashTable<Type>::defaultKey);
    ~closeHashTable(){ delete []array;}
    bool find(const Type&x)const;
    bool insert(const Type&x);
    bool remove(const Type&x);

    void rehash();
    std::vector<Type> list()const;
};

template<class Type>
closeHashTable<Type>::closeHashTable(int length,int _maxDeleteSize,int(*f)(const Type&x))
{
    maxDeleteSize = _maxDeleteSize;
    dataSize=0;
    deleteSize=0;
    size = length;
    hashTable<Type>::key = f;
    array = new node[size];
}
template<class Type>
bool closeHashTable<Type>::insert(const Type&x)
{
    int initPos,pos;
    if(dataSize>size/2) doubleSpace();
    initPos = pos = hashTable<Type>::key(x) % size;
    dataSize++;
    do{
        if(array[pos].state!=1)
        {
            array[pos].data = x;
            array[pos].state = 1;
            return true;
        }
        if(array[pos].state == 1&&array[pos].data ==x){
            return true;
        }
        pos = (pos+1)%size;
    }while(pos!=initPos);
    dataSize--;
    return false;
}
template<class Type>
bool closeHashTable<Type>::remove(const Type&x)
{
    int initPos,pos;
    autoRehash();
    initPos = pos = hashTable<Type>::key(x) % size;
    do{
        if(array[pos].state ==0)return false;
        if(array[pos].state == 1&&array[pos].data ==x)
        {
            array[pos].state = 2;
            dataSize--;
            deleteSize++;
            return true;
        }
        pos = (pos+1)%size;
    }while(pos!=initPos);
    return false;
}
template<class Type>
bool closeHashTable<Type>::find(const Type&x)const
{
    int initPos,pos;
    initPos = pos= hashTable<Type>::key(x)%size;

    do{
        if(array[initPos].state==0)return false;
        if(array[initPos].state==1 && array[initPos].data ==x)return true;
        pos = (pos+1)%size;
    }while(initPos!=pos);
    return false;
}
template<class Type>
void closeHashTable<Type>::rehash()
{
    deleteSize = 0;
    dataSize = 0;
    node* tmp = array;
    array = new node[size];
    for(int i=0;i<size;i++)
        if(tmp[i].state==1) insert(tmp[i].data);
    delete [] tmp;
}
template<class Type>
void closeHashTable<Type>::autoRehash()
{
    if(deleteSize>maxDeleteSize) {
        //std::cout << "rehash()被调用，此时 有效数据："<<dataSize<<",无效数据："<<deleteSize<<std::endl;
        rehash();
    }
}
template<class Type>
void closeHashTable<Type>::doubleSpace()
{
    //std::cout << "doubleSpace()被调用，此时 有效数据："<<dataSize<<",容量："<<size<<std::endl;
    node* tmp = array;
    array = new node[2*size];
    dataSize =0;
    deleteSize =0;
    for(int i=0;i<size;i++)
        if(tmp[i].state==1) insert(tmp[i].data);
    size*=2;
    delete []tmp;
}
template<class Type>
std::vector<Type> closeHashTable<Type>::list()const
{
    std::vector<Type> l;
    for(int i=0;i<size;i++)
	if(array[i].state==1) l.push_back(array[i].data);
    return l;
}
/*
template<class Type>
int closeHashTable<Type>::currentSize()const
{
    int cnt=0;
    for(int i=0;i<size;i++)
        if(array[i].state==1)
            cnt++;
    return cnt;
}*/
#endif // CLOSEHASHTABLE_H_INCLUDED
