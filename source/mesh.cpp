#include "../include/mesh.h"
#include "../include/meshmatrix.h"
#include "../include/closeHashTable.h"
#include "../include/DisjointSet.h"
#include <cstring>
#include <iostream>
#include <strstream>
#include <fstream>
#include <cmath>
#include <float.h>
#include <queue>
#include <stack>
using namespace std;

/////////////////////////////////////////
// helping inline functions
inline double Cot(const Vector3d & p1, const Vector3d & p2, const Vector3d & p3) {
	Vector3d v1 = p1 - p2;
	Vector3d v2 = p3 - p2;

	v1 /= v1.L2Norm();
	v2 /= v2.L2Norm();
	double tmp = v1.Dot(v2);
	return 1.0 / tan(acos(tmp));
}

inline double Area(const Vector3d & p1, const Vector3d & p2, const Vector3d & p3) {
	Vector3d v1 = p2 - p1;
	Vector3d v2 = p3 - p1;
	return v1.Cross(v2).L2Norm() / 2.0;
}


/////////////////////////////////////////
// implementation of OneRingHEdge class
OneRingHEdge::OneRingHEdge(const Vertex * v) {
	if (v == NULL) start = next = NULL;
	else start = next = v->HalfEdge();
}

HEdge * OneRingHEdge::NextHEdge() {
	HEdge *ret = next;
	if (next && next->Prev()->Twin() != start)
		next = next->Prev()->Twin();
	else
		next = NULL;
	return ret;
}

/////////////////////////////////////////
// implementation of Mesh class
//
// function AddFace
// it's only for loading obj model, you do not need to understand it
void Mesh::AddFace(int v1, int v2, int v3) {
	int i;
	HEdge *he[3], *bhe[3];
	Vertex *v[3];
	Face *f;

	// obtain objects
	for (i=0; i<3; i++) he[i] = new HEdge();
	for (i=0; i<3; i++) bhe[i] = new HEdge(true);
	v[0] = vList[v1];
	v[1] = vList[v2];
	v[2] = vList[v3];
	f = new Face();

	// connect prev-next pointers
	SetPrevNext(he[0], he[1]);
	SetPrevNext(he[1], he[2]);
	SetPrevNext(he[2], he[0]);
	SetPrevNext(bhe[0], bhe[1]);
	SetPrevNext(bhe[1], bhe[2]);
	SetPrevNext(bhe[2], bhe[0]);

	// connect twin pointers
	SetTwin(he[0], bhe[0]);
	SetTwin(he[1], bhe[2]);
	SetTwin(he[2], bhe[1]);

	// connect start pointers for bhe
	bhe[0]->SetStart(v[1]);
	bhe[1]->SetStart(v[0]);
	bhe[2]->SetStart(v[2]);
	for (i=0; i<3; i++) he[i]->SetStart(v[i]);

	// connect start pointers
	// connect face-hedge pointers
	for (i=0; i<3; i++) {
		v[i]->SetHalfEdge(he[i]);
		v[i]->adjHEdges.push_back(he[i]);
		SetFace(f, he[i]);
	}
	v[0]->adjHEdges.push_back(bhe[1]);
	v[1]->adjHEdges.push_back(bhe[0]);
	v[2]->adjHEdges.push_back(bhe[2]);

	// mearge boundary if in need
	for (i=0; i<3; i++) {
		Vertex *start = bhe[i]->Start();
		Vertex *end   = bhe[i]->End();
		for (size_t j=0; j<end->adjHEdges.size(); j++) {
			HEdge *curr = end->adjHEdges[j];
			if (curr->IsBoundary() && curr->End()==start) {
				SetPrevNext(bhe[i]->Prev(), curr->Next());
				SetPrevNext(curr->Prev(), bhe[i]->Next());
				SetTwin(bhe[i]->Twin(), curr->Twin());
				bhe[i]->SetStart(NULL);	// mark as unused
				curr->SetStart(NULL);	// mark as unused
				break;
			}
		}
	}

	// finally add hedges and faces to list
	for (i=0; i<3; i++) heList.push_back(he[i]);
	for (i=0; i<3; i++) bheList.push_back(bhe[i]);
	fList.push_back(f);
}

// function LoadObjFile
// it's only for loading obj model, you do not need to understand it
bool Mesh::LoadObjFile(const char *filename) {
	if (filename==NULL || strlen(filename)==0) return false;
	ifstream ifs(filename);
	if (ifs.fail()) return false;

	Clear();

	char buf[1024], type[1024];
	do {
		ifs.getline(buf, 1024);
		istrstream iss(buf);
		iss >> type;

		// vertex
		if (strcmp(type, "v") == 0) {
			double x, y, z;
			iss >> x >> y >> z;	
			Vertex* vtmp = new Vertex(x,y,z);
            AddVertex(vtmp);
		}
		// face
		else if (strcmp(type, "f") == 0) {
			int index[3];
			iss >> index[0] >> index[1] >> index[2];
			AddFace(index[0]-1, index[1]-1, index[2]-1);
		}
	} while (!ifs.eof());
	ifs.close();

	size_t i;
	Vector3d box = this->MaxCoord() - this->MinCoord();
	for (i=0; i<vList.size(); i++) vList[i]->SetPosition(vList[i]->Position() / box.X());

	Vector3d tot;
	for (i=0; i<vList.size(); i++) tot += vList[i]->Position();
	Vector3d avg = tot / vList.size();
	for (i=0; i<vList.size(); i++) vList[i]->SetPosition(vList[i]->Position() - avg);

	HEdgeList list;
	for (i=0; i<bheList.size(); i++)
		if (bheList[i]->Start()) list.push_back(bheList[i]);
	bheList = list;

	for (i=0; i<vList.size(); i++) 
	{
		vList[i]->adjHEdges.clear(); 
		vList[i]->SetIndex((int)i);
		vList[i]->SetFlag(0);
	}

	// SphereClassifier
	//SC.showStatus();
	return true;
}

void Mesh::DisplayMeshInfo()
{
	int NO_VERTICES = (int)vList.size();
	int NO_FACES = (int)fList.size();
	int NO_HEDGES = (int)heList.size()+(int)bheList.size();
	int NO_B_LOOPS = CountBoundaryLoops();
	//int NO_COMPONENTS = CountConnectedComponents();
	int NO_COMPONENTS = 0;
	int NO_GENUS = NO_COMPONENTS - (NO_VERTICES - NO_HEDGES/2 +  NO_FACES + NO_B_LOOPS)/2;

	
	cout << "the number of vertices is " <<NO_VERTICES << endl;
	cout << "the number of faces is " <<NO_FACES << endl;
    cout << "the number of half edges is "<< NO_HEDGES << endl;
	cout << "the number of boundary loops is "<< NO_B_LOOPS << endl;
    cout << "the number of connected components is "<< NO_COMPONENTS << endl;
	cout << "the number of genus is "<< NO_GENUS << endl;
    cout << "the size of the boundary loop is  "<< bheList.size() << endl;
	
}

int Mesh::CountBoundaryLoops()
{
    int no_loop =0;//count the number of boundary loops
	size_t i;
    for (i=0; i< bheList.size(); i++)
	{
	   HEdge *cur=bheList[i];
       HEdge *nex=cur;
	   while(nex->Start()->visit!=1)
	   {
	     nex->Start()->visit=1;
         nex=nex->Next();
		 if (nex==cur)
		 {no_loop++;break;} 
	   }
	}
	return no_loop;

}

int Mesh::CountConnectedComponents()
{

    int no_component =0;//count the number no_component
	size_t i;
    for (i=0; i< vList.size(); i++) // initialize
	   vList[i]->visit=2;
	size_t j;
    for (j=0; j< vList.size(); j++) // for each iteration, vertices with visit=2 forms a connected component and will not be visited anymore.
	   if(vList[j]->visit==2)
	   {no_component++; DFSVisit(vList[j]);}
	return no_component;

}

void Mesh::DFSVisit(Vertex * v)
{
	if(v->HalfEdge()==NULL)return; // to avoid NULL pointer
	v->visit=3;
	HEdge *cur = v->HalfEdge()->Twin();
	HEdge *nex = cur;
	while (nex && nex->Twin()->Prev()!=cur) //assign 'visit' for the neighbors of v
	{
        if(nex->Start()->visit==2)
		DFSVisit(nex->Start());
		nex = nex->Twin()->Prev();	
	}
    v->visit=4;
}


// -------------------------------------------------------
// DO NOT TOUCH THE FOLLOWING FOR NOW
// -------------------------------------------------------
void Mesh::ComputeVertexNormals()
{
	//Area-weighted scheme
	size_t i;
    for (i=0; i< vList.size(); i++) 
	{
	   Vector3d  v_normal(0.0,0.0,0.0);
	   double  v_A_1=0.0; //initialize the area of one face incident to each vertex v

       HEdge *cur=vList[i]->HalfEdge();
       HEdge *nex=cur;
	   while(nex && nex->Prev()->Twin()!=cur) 
	   {
		 if(nex && !nex->IsBoundary())//exclude the holes formed by boundary loops
		 {		 
	     
		   const Vector3d & pos1 = vList[i]->Position();//p1
           const Vector3d & pos2 = nex->End()->Position();//p2
           const Vector3d & pos3 =  nex->Prev()->Start()->Position();//p3
		   v_A_1 = Area(pos1, pos2, pos3);
           v_normal = v_normal + v_A_1*nex->LeftFace()->Normal_f();//SUM
		 }
		 nex = nex->Prev()->Twin();	  
	   }
       
       //process the last nex
		 if(nex && !nex->IsBoundary())
		 {		 
	     
		   const Vector3d & pos1 = vList[i]->Position();//p1
           const Vector3d & pos2 = nex->End()->Position();//p2
           const Vector3d & pos3 =  nex->Prev()->Start()->Position();//p3
		   v_A_1 = Area(pos1, pos2, pos3);
           v_normal = v_normal + v_A_1*nex->LeftFace()->Normal_f();//SUM
		 }

       v_normal /= v_normal.L2Norm(); 
	   vList[i]->SetNormal(v_normal);
	} 
}

void Mesh::UmbrellaSmooth() 
{
	cout<<"UmbrellaSmooth starts"<<endl;

	//construct matrix L
    double cot_sum=0.0;//initialize the total weight of neighboring vertices of the current vertex
    double cot_nex=0.0;//weight due to the current neighboring vertex

	double landa = 0.8;//set landa = 0.8 for 'Xt1 = Xt + landa*L*Xt'
	double cot_next =0.0;//reset cot_nex

    Matrix L((int)vList.size(), (int)vList.size());

    int i;
    for (i=0; i< (int)vList.size(); i++)
	{

	  L.AddElement(i, i, -landa);//add element to L
	  //compute 'cot_sum': the sum of weight incident to the current vertex
	  
	  HEdge *cur=vList[i]->HalfEdge();
      HEdge *nex=cur;
	  while(nex && nex->Prev()->Twin()!=cur) 
	  {
	    if( !vList[i]->IsBoundary())//interior vertices
	    {		 
		  const Vector3d & pos1 = vList[i]->Position();//p1
          const Vector3d & pos2 = nex->End()->Position();//p2
          const Vector3d & pos3 =  nex->Prev()->Start()->Position();//p3
		  const Vector3d & pos4 =  nex->Prev()->Twin()->Prev()->Start()->Position();//p4
          cot_nex = Cot(pos1, pos2, pos3) + Cot(pos1, pos4, pos3);
	      cot_sum = cot_sum + cot_nex; 
	   
		}//end if
        nex = nex->Prev()->Twin();
	  }//end while
      //process the last 'nex'
	  if( !vList[i]->IsBoundary())//interior vertices
	  {		 
		const Vector3d & pos1 = vList[i]->Position();//p1
        const Vector3d & pos2 = nex->End()->Position();//p2
        const Vector3d & pos3 =  nex->Prev()->Start()->Position();//p3
		const Vector3d & pos4 =  nex->Prev()->Twin()->Prev()->Start()->Position();//p4
        cot_nex = Cot(pos1, pos2, pos3) + Cot(pos1, pos4, pos3);
	    cot_sum = cot_sum + cot_nex;
	  }//end if
      cot_nex=0.0;//reset cot_nex

      // compute matrix 'L = landa*L'
      HEdge *curr=vList[i]->HalfEdge();
      HEdge *next=curr;
	  while(next && next->Prev()->Twin()!=curr) 
	  {
	    if( !vList[i]->IsBoundary())//interior vertices
	    {		 
		  const Vector3d & pos1 = vList[i]->Position();//p1
          const Vector3d & pos2 = next->End()->Position();//p2
          const Vector3d & pos3 =  next->Prev()->Start()->Position();//p3
		  const Vector3d & pos4 =  next->Prev()->Twin()->Prev()->Start()->Position();//p4
          cot_next = Cot(pos1, pos2, pos3) + Cot(pos1, pos4, pos3);




          //add element to matrix L
          L.AddElement(i, next->Prev()->Start()->Index(), landa*cot_next/cot_sum);
	   
		}//end if
        next = next->Prev()->Twin();
	  }//end while
      //process the last 'next'
	  if( !vList[i]->IsBoundary())//interior vertices
	  {		 
		const Vector3d & pos1 = vList[i]->Position();//p1
        const Vector3d & pos2 = next->End()->Position();//p2
        const Vector3d & pos3 =  next->Prev()->Start()->Position();//p3
		const Vector3d & pos4 =  next->Prev()->Twin()->Prev()->Start()->Position();//p4
        cot_next = Cot(pos1, pos2, pos3) + Cot(pos1, pos4, pos3);




        //add element to matrix L
        L.AddElement(i, next->Prev()->Start()->Index(), landa*cot_next/cot_sum);

        
	  }//end if

      cot_next =0.0;//reset cot_next
	  cot_sum = 0.0;//reset cot_sum

	}//end the first 'for' loop in this function

    //sort matrix L
    L.SortMatrix();

    //define arrays
	double **Xt_array = new double*[3];
	double **Xt1_array= new double*[3];
    
	for (i=0; i< 3; i++)
	Xt_array[i] = new double[vList.size()];
	for (i=0; i< 3; i++) 
	Xt1_array[i] = new double[vList.size()];

	//compute  Xt_array
    for (i=0; i< (int)vList.size(); i++)
    {
      Xt_array[0][i] =  vList[i]->Position().X();
      Xt_array[1][i] =  vList[i]->Position().Y();
      Xt_array[2][i] =  vList[i]->Position().Z();
    }	

    //compute  matrix 'landa*L*Xt'
    L.Multiply(Xt_array[0], Xt1_array[0]);
	L.Multiply(Xt_array[1], Xt1_array[1]);
	L.Multiply(Xt_array[2], Xt1_array[2]);
    
	//compute  matrix 'Xt1 = Xt + landa*L*Xt'
    for (i=0; i< (int)vList.size(); i++)
	{
	   Xt1_array[0][i] = Xt1_array[0][i] + Xt_array[0][i];
	   Xt1_array[1][i] = Xt1_array[1][i] + Xt_array[1][i];
	   Xt1_array[2][i] = Xt1_array[2][i] + Xt_array[2][i];	
	}
    
    //release space for matrix Xt
	for (i=0; i< 3; i++) 
	delete[] Xt_array[i];

	delete[] Xt_array;


ofstream outfile("mannequin.txt"); //20150531


    //reset the poistions of the interior vertices
	for (i=0; i< (int)vList.size(); i++)
	{
outfile << "v" <<" "<<Xt1_array[0][i]<<" "<<Xt1_array[1][i]<<" "<<Xt1_array[2][i]<<"\n";//20150531
	  if(!vList[i]->IsBoundary())//interior vertices
	  {
	    double x1 = Xt1_array[0][i];
        double y1 = Xt1_array[1][i];
	    double z1 = Xt1_array[2][i];  

        //ofstream outfile("bu_head_test.txt"); //20150531
//outfile << "v" <<" "<<Xt1_array[0][i]<<" "<<Xt1_array[1][i]<<" "<<Xt1_array[2][i]<<"\n";//20150531
        //outfile.close();



        Vector3d new_position(x1,y1,z1);
        vList[i]->SetPosition(new_position);
	  }//end if

	}

outfile.close();//20150531



	//release space for matrix Xt1
    for (i=0; i< 3; i++) 
	delete[] Xt1_array[i];

	delete[] Xt1_array;
	
    cout<<"done now "<<endl;
}

void Mesh::ImplicitUmbrellaSmooth()
{
    cout<< "ImplicitUmbrellaSmooth starts..."<<endl;
	//construct matrix A for AX=b, A=I-landa*L
    double cot_sum=0.0;//initialize the total weight of neighboring vertices of the current vertex
    double cot_nex=0.0;//weight due to the current neighboring vertex

	double landa = 0.8;//set landa = 0.8 for 'Xt1 = Xt + landa*L*Xt'
	double cot_next =0.0;//reset cot_nex

    Matrix A((int)vList.size(), (int)vList.size());

    int i;
    for (i=0; i< (int)vList.size(); i++)
	{

	  A.AddElement(i, i, 1+landa);//add element to L
	  //compute 'cot_sum': the sum of weight incident to the current vertex
	  
	  HEdge *cur=vList[i]->HalfEdge();
      HEdge *nex=cur;
	  while(nex && nex->Prev()->Twin()!=cur) 
	  {
	    if( !vList[i]->IsBoundary())//interior vertices
	    {		 
		  const Vector3d & pos1 = vList[i]->Position();//p1
          const Vector3d & pos2 = nex->End()->Position();//p2
          const Vector3d & pos3 =  nex->Prev()->Start()->Position();//p3
		  const Vector3d & pos4 =  nex->Prev()->Twin()->Prev()->Start()->Position();//p4
          cot_nex = Cot(pos1, pos2, pos3) + Cot(pos1, pos4, pos3);
	      cot_sum = cot_sum + cot_nex; 
	   
		}//end if
        nex = nex->Prev()->Twin();
	  }//end while
      //process the last 'nex'
	  if( !vList[i]->IsBoundary())//interior vertices
	  {		 
		const Vector3d & pos1 = vList[i]->Position();//p1
        const Vector3d & pos2 = nex->End()->Position();//p2
        const Vector3d & pos3 =  nex->Prev()->Start()->Position();//p3
		const Vector3d & pos4 =  nex->Prev()->Twin()->Prev()->Start()->Position();//p4
        cot_nex = Cot(pos1, pos2, pos3) + Cot(pos1, pos4, pos3);
	    cot_sum = cot_sum + cot_nex;
	  }//end if
      cot_nex=0.0;//reset cot_nex

      // compute matrix 'L = landa*L'
      HEdge *curr=vList[i]->HalfEdge();
      HEdge *next=curr;
	  while(next && next->Prev()->Twin()!=curr) 
	  {
	    if( !vList[i]->IsBoundary())//interior vertices
	    {		 
		  const Vector3d & pos1 = vList[i]->Position();//p1
          const Vector3d & pos2 = next->End()->Position();//p2
          const Vector3d & pos3 =  next->Prev()->Start()->Position();//p3
		  const Vector3d & pos4 =  next->Prev()->Twin()->Prev()->Start()->Position();//p4
          cot_next = Cot(pos1, pos2, pos3) + Cot(pos1, pos4, pos3);

          //add element to matrix L
          A.AddElement(i, next->Prev()->Start()->Index(), 1-landa*cot_next/cot_sum);
	   
		}//end if
        next = next->Prev()->Twin(); 
	  }//end while
      //process the last 'next'
	  if( !vList[i]->IsBoundary())//interior vertices
	  {		 
		const Vector3d & pos1 = vList[i]->Position();//p1
        const Vector3d & pos2 = next->End()->Position();//p2
        const Vector3d & pos3 =  next->Prev()->Start()->Position();//p3
		const Vector3d & pos4 =  next->Prev()->Twin()->Prev()->Start()->Position();//p4
        cot_next = Cot(pos1, pos2, pos3) + Cot(pos1, pos4, pos3);

        //add element to matrix L
        A.AddElement(i, next->Prev()->Start()->Index(), 1-landa*cot_next/cot_sum);

        
	  }//end if

      cot_next =0.0;//reset cot_next
	  cot_sum = 0.0;//reset cot_sum

	}//end the first 'for' loop in this function

    //sort matrix L
    A.SortMatrix();

	//construct arrays
	double **Xt_array = new double*[3];
	double **Xt1_array= new double*[3];
    
	for (i=0; i< 3; i++)
	Xt_array[i] = new double[vList.size()];
	for (i=0; i< 3; i++) 
	Xt1_array[i] = new double[vList.size()];


	//compute  Xt_array
    for (i=0; i< (int)vList.size(); i++)
    {
      Xt_array[0][i] =  vList[i]->Position().X();
      Xt_array[1][i] =  vList[i]->Position().Y();
      Xt_array[2][i] =  vList[i]->Position().Z(); 
    }


   	//initialize Xt1_array 
    for (i=0; i< (int)vList.size(); i++)
    {
      Xt1_array[0][i] =  0.001;
      Xt1_array[1][i] =  0.001;
      Xt1_array[2][i] =  0.001; 
    }


    //call the BiConjugate Gradient solver to construct Xt1_array
	  A.BCG(Xt_array[0], Xt1_array[0], 1, 2e-2);
	  A.BCG(Xt_array[1], Xt1_array[1], 1, 2e-2);
	  A.BCG(Xt_array[2], Xt1_array[2], 1, 2e-2);	

	//release space for matrix Xt
	for (i=0; i< 3; i++) 
	delete[] Xt_array[i];

	delete[] Xt_array;

    //reset the poistions of the interior vertices
	for (i=0; i< (int)vList.size(); i++)
	{
	  if(!vList[i]->IsBoundary())//interior vertices
	  {
	    double x1 = Xt1_array[0][i];
        double y1 = Xt1_array[1][i];
	    double z1 = Xt1_array[2][i];
		//cout<<"x1, y1, z1 is "<<x1<<y1<<z1<<endl;
        Vector3d new_position(x1,y1,z1);
        vList[i]->SetPosition(new_position);
	  }//end if
	}

	//release space for matrix Xt1
    for (i=0; i< 3; i++) 
	delete[] Xt1_array[i];

	delete[] Xt1_array;
	
	cout<<"done now "<<endl;
	/*************************/
	/* insert your code here */
	/*************************/
}
void Mesh::ComputeVertexCurvatures(int Mode)
{
	//Gaussian curvature
	size_t i;
	double min_c=200000;//initialize the maximum value of curvature
	double max_c=-200000;//initialize the minimum value of curvature

    for (i=0; i< vList.size(); i++) 
	{
    if(!vList[i]->IsBoundary()) //for interior vertices only
	{
    double  v_A_1 = 0.0; //initialize the area of any face incident to each vertex v
    double  v_A = 0.0;//the area of all face incident to vertex v
	double  Angle_v_1 = 0.0; //initialize any angle incident on each vertex v
    double  Angle_v = 0.0;//the angles incident on vertex v

	double cot_sum=0.0;//coefficient of (vj - vi) in the mean curvature function
    Vector3d color_v(0.0,0.0,0.0);

    HEdge *cur=vList[i]->HalfEdge();
    HEdge *nex=cur;
	while(nex && nex->Prev()->Twin()!=cur) 
	{
	  if(nex && !nex->IsBoundary())//exclude the holes formed by boundary loops
	    {		 
		   const Vector3d & pos1 = vList[i]->Position();//p1
           const Vector3d & pos2 = nex->End()->Position();//p2
           const Vector3d & pos3 =  nex->Prev()->Start()->Position();//p3
           if(!vList[i]->IsBoundary())//if the current vertex is an interior vertex
		   {
		     const Vector3d & pos4 =  nex->Prev()->Twin()->Prev()->Start()->Position();//p4, for computing the mean curvature
             cot_sum = Cot(pos1, pos2, pos3) + Cot(pos1, pos4, pos3);
             color_v = color_v + cot_sum*(pos3 - pos1);  

		   }
		   v_A_1 = Area(pos1, pos2, pos3);
           v_A = v_A+ v_A_1;//SUM

		   Vector3d v1 = pos2 - pos1;
	       Vector3d v2 = pos3 - pos1;
		   v1 /= v1.L2Norm();
	       v2 /= v2.L2Norm();
	       double temp = v1.Dot(v2);
           Angle_v_1 = acos(temp);
           Angle_v = Angle_v +  Angle_v_1;
		 }
         
		 nex = nex->Prev()->Twin();	  
         
	}
       
    //process the last nex
    if(nex && !nex->IsBoundary())
	  {		    
		const Vector3d & pos1 = vList[i]->Position();//p1
        const Vector3d & pos2 = nex->End()->Position();//p2
        const Vector3d & pos3 =  nex->Prev()->Start()->Position();//p3
		if(!vList[i]->IsBoundary())//if the current vertex is an interior vertex
		{
		  const Vector3d & pos4 =  nex->Prev()->Twin()->Prev()->Start()->Position();//p4, for computing the mean curvature
          cot_sum = Cot(pos1, pos2, pos3) + Cot(pos1, pos4, pos3);
          color_v = color_v + cot_sum*(pos3 - pos1);  

		}

		v_A_1 = Area(pos1, pos2, pos3);
        v_A =  v_A + v_A_1;//SUM
		
        //computing the sum of angles for Gaussian curvature
		Vector3d v1 = pos2 - pos1;
	    Vector3d v2 = pos3 - pos1;
		v1 /= v1.L2Norm();
	    v2 /= v2.L2Norm();
	    double temp = v1.Dot(v2);
        Angle_v_1 = acos(temp);
        Angle_v = Angle_v +  Angle_v_1;
	}//end if

    double PI = 3.14159265358979323846; 
	double G = (2.0/v_A)*(2*PI - Angle_v);//Gaussian curvature
	double H = -color_v.L2Norm()/(2*v_A); //mean curvature   

	//compute principle curvatures
	double K1 = H + sqrt(fabs(H*H - G));//max principle curvature
	
	double K2 = H - sqrt(fabs(H*H - G));//min principle curvature
    //remember different types of curvatures for each 
    vList[i]->G = G; //record gaussian curvature at vList[i]
    vList[i]->H = H; //record mean curvature at vList[i]
	vList[i]->K1 = K1; //record max principle curvature at vList[i]
	vList[i]->K2 = K2; //record min principle curvature at vList[i]
	
	//"max principal curvature" scheme
    if(Mode == MAX_CURVATURE)
	{
		if(K1 < min_c) // reset min curvatures
		min_c = K1;
		if(K1 > max_c) //reset max curvatures
		max_c = K1;	
	}
	if(Mode == MIN_CURVATURE)
	{
		//"min principal curvature" scheme
		if(K2 < min_c) // reset min curvatures
		min_c = K2;
		if(K2 > max_c) //reset max curvatures
		max_c = K2;	
	}

	
    //"mean curvature" scheme
    if(Mode == MEAN_CURVATURE)
	{
		if(H < min_c) // reset min curvatures
		min_c = H;
		if(H > max_c) //reset max curvatures
		max_c = H;	
	}
  

	if(Mode == GAUSSIAN_CURVATURE)
	{
		//"Gaussian curvature" scheme
		if(G < min_c) // reset min curvatures
		min_c = G;
		if(G > max_c) //reset max curvatures
		max_c = G;	
	}

	}//end 'if' condition for interior vertices

	}//end the first 'for' loop
	
	/**********interpolation based on the scheme of [0->1] equals [blue->red]**********/
	//"max principal curvature" scheme

	if(Mode == MAX_CURVATURE)
	{
		for (i=0; i< vList.size(); i++) 
		{ 
		  if(!vList[i]->IsBoundary())//process interior vertices
		  {
			Vector3d color_max((vList[i]->K1 - min_c)/(max_c - min_c),0.0,1-(vList[i]->K1 - min_c)/(max_c - min_c));
			vList[i]->SetColor(color_max); 
		  }
		  else //process boundary vertices
		  {
			Vector3d color_max(0.0,0.0,0.0);
			vList[i]->SetColor(color_max); 
		  }
		}//end this for
	}
	if(Mode == MIN_CURVATURE)
	{
		//"min principal curvature" scheme
		for (i=0; i< vList.size(); i++) 
		{ 
		  if(!vList[i]->IsBoundary())//process interior vertices
		  {
			Vector3d color_min((vList[i]->K2 - min_c)/(max_c - min_c),0.0,1-(vList[i]->K2 - min_c)/(max_c - min_c));
			vList[i]->SetColor(color_min); 
		  }
		  else //process boundary vertices
		  {
			Vector3d color_max(0.0,0.0,0.0);
			vList[i]->SetColor(color_max); 
		  }
		}//end this for
	}

	if(Mode==MEAN_CURVATURE)
	{
		//"mean curvature" scheme
		for (i=0; i< vList.size(); i++) 
		{ 
		  if(!vList[i]->IsBoundary())//process interior vertices
		  {
			Vector3d color_mean((vList[i]->H - min_c)/(max_c - min_c),0.0,1-(vList[i]->H - min_c)/(max_c - min_c));
			vList[i]->SetColor(color_mean); 
		  }
		  else //process boundary vertices
		  {
			Vector3d color_mean(0.0,0.0,0.0);
			vList[i]->SetColor(color_mean); 
		  }
		}//end this for
	}

	//"Gaussian curvature" scheme
	if(Mode == GAUSSIAN_CURVATURE)
	{
		/*
		for (i=0; i< vList.size(); i++) 
		{ 
		  if(!vList[i]->IsBoundary())//process interior vertices
		  {
			Vector3d color_Gau((vList[i]->G - min_c)/(max_c - min_c),0.0,1-(vList[i]->G - min_c)/(max_c - min_c));
			//if(vList[i]->G >2000)
			//	 color_Gau = Vector3d(0,0.9,0);
			vList[i]->SetColor(color_Gau); 
		  }
		  else //process boundary vertices
		  {
			Vector3d color_max(0.0,0.0,0.0);
			vList[i]->SetColor(color_max); 
		  }
		  
		
		}//end this for
		*/
		cout << "Gaussian Curvature Computed"<<endl;
	}

}

extern int currSelectedVertex;
//一.切割法
// 取高斯曲率为负的点 放到“沟”里 漫水填充
// 1.找出所有曲率满足条件的点
// 2.对这些点拟合出适当的平面
// 3.用这些平面分割模型 分别储存
// 4.判断currSelectedPoint在哪块模型中
// 5.取出模型块
//二.高斯曲率块
//三.高斯曲率为种子 选择周边块
void Mesh::Separate()
{
	// calculate normals and curvatures first
	this->ComputeVertexCurvatures();
	this->ComputeVertexNormals();
	
	
}
/*
void Mesh::MarkCurvPoints(int threshold)
{
	vGList.clear();
	queue<Vertex*> s;
	//找到所有高斯曲率大点
	for(int i=0;i<vList.size();i++)
	{
		if(vList[i]->G < threshold) continue;
		else{
			vList[i]->GVisited(true);
			vGList.push_back(vList[i]);
		}
	}

	cout << "vGList size "<<vGList.size()<<endl;
	if(vGList.empty())return;

	for(int i=0;i<vList.size();i++)
		vList[i]->SetColor(Vector3d(0,0,1));
	for(int i=0;i<fList.size();i++)
		fList[i]->setInBlock(false);

	feat.clear();
	int selectMax = 100; //选择高斯曲率大的点周围100个点 
	for(int i=0;i<vGList.size();i++)
	{
		if(!vGList[i]->GVisited())continue; // 高斯曲率不符合 或者已经被Block包含
		// create new block
		// 这里是创建block并往里面塞点 塞点的规则可以在这里修改！！
		int selected=0;
		feat.start();
		vGList[i]->GVisited(false);
		s.push(vGList[i]);
		feat.push(vGList[i]);
		selected++;
		while(!s.empty() && selected <selectMax)
		{
			Vertex* curVertex = s.front();s.pop();
			HEdge * curEdge = curVertex->HalfEdge();
			Vertex* endV = curEdge->End();
			do
			{
				if(curEdge->End()->GVisited())curEdge->End()->GVisited(false); // 此处改变了GVisited不是很好 但所有满足的点都在vGList中了 GVisited失效了 
				s.push(curEdge->End());
				feat.push(curEdge->End());
				if(!curEdge->LeftFace()->isInBlock())
				{
					curEdge->LeftFace()->setInBlock(true);
					feat.pushface(curEdge->LeftFace());
				}
				selected++;
				curEdge = curEdge->Twin()->Next();
			}while(curEdge->End()!=endV);
		}
		feat.end();
	}
	
	//for(int i=0;i<vGList.size();i++)
	//{
	//	if(!vGList[i]->GVisited())continue;
	//	// create new block
	//	// 这里是创建block并往里面塞点 塞点的规则可以在这里修改！！
	//	feat.start();
	//	vGList[i]->GVisited(false);
	//	s.push(vGList[i]);
	//	feat.push(vGList[i]);
	//	while(!s.empty())
	//	{
	//		Vertex* curVertex = s.top();
	//		s.pop();
	//		HEdge * curEdge = curVertex->HalfEdge();
	//		Vertex* endV = curEdge->End();
	//		do
	//		{
	//			if(curEdge->End()->GVisited()){	
	//				curEdge->End()->GVisited(false); // 此处改变了GVisited不是很好 但所有满足的点都在vGList中了 GVisited失效了 
	//				s.push(curEdge->End());
	//				feat.push(curEdge->End());
	//			}
	//			curEdge = curEdge->Twin()->Next();
	//		}while(curEdge->End()!=endV);
	//	}
	//	// 调试时使用的 
	//	if(feat[feat.size()-1].size()==1)
	//	{
	//		cout << i<<"---"<<vGList[i]->G << endl;
	//		HEdge * curEdge = vGList[i]->HalfEdge();
	//		Vertex* endV = curEdge->End();
	//		do
	//		{
	//			cout << curEdge->End()->G <<" "<<endl;
	//			curEdge = curEdge->Twin()->Next();
	//		}while(curEdge->End()!=endV);
	//		vGList[i]->SetColor(Vector3d(1,0,0));
	//	}
	//	feat.end();
	//}
	// 对特征区块进行筛选
	feat.Select();
	
	cout << "feat size " << feat.size()<< endl;
	
	for(int i=0;i<feat.size();i++)cout << feat.blockSize(i)<< " ";
	
	// 可视化特征区块
	feat.SetColor();

	// 计算直方图
	
}*/

/*
void Mesh::SeparateLines()
{
	for(int i=0;i<vGList.size();i++)
		vGList[i]->SetGDepth(0);
	for(int i=0;i<vGList.size();i++)
	{
		if(vGList[i]->GinLine())continue;
		// 否则创建一条新路
		queue<Vertex*> vQue;
		vGList[i]->GinLine(true);
		vQue.push(vGList[i]);
		// 从选择点向外遍历 
		while(!vQue.empty())
		{
			Vertex* curV = vQue.front();
			vQue.pop();
			HEdge * curEdge = curV->HalfEdge();
			Vertex* endV = curEdge->End();
		
			if(!curEdge->End()->GinLine() && curEdge->End()->GVisited())
			{
				curEdge->End()->GinLine(true);
				if(curEdge->End()->GetGDepth()<curV->GetGDepth()+1)
					curEdge->End()->SetGDepth(curV->GetGDepth()+1);
				vQue.push(curEdge->End());
			}
			curEdge = curEdge->Twin()->Next();
			while( curEdge->End()!= endV)
			{	
				if(!curEdge->End()->GVisited() && curEdge->End()->GVisited())
				{
					curEdge->End()->GinLine(true);
					if(curEdge->End()->GetGDepth()<curV->GetGDepth()+1)
						curEdge->End()->SetGDepth(curV->GetGDepth()+1);
					vQue.push(curEdge->End());
				}
				curEdge = curEdge->Twin()->Next();
			}
		}
	}
	//算完后找到depth No1和No2 的点 反向找到路径
}*/
/*
// 分水岭算法
void Watershed()
{
	
}
*/
// 局部特征

/*
void Mesh::MarkCurvPoints(int threshold)
{
	// select a vertex
	if(currSelectedVertex==-1 || vList[currSelectedVertex]->IsBoundary())
	{
		cout << "Select a vertex" << endl;
		return;
	}
	queue<Vertex*> vGList;
	vList[currSelectedVertex]->GVisited(true);
	vGList.push(vList[currSelectedVertex]);

	// 从选择点向外遍历 

	while(!vGList.empty())
	{
		Vertex* curV = vGList.front();
		vGList.pop();
		HEdge * curEdge = curV->HalfEdge();
		Vertex* endV = curEdge->End();
		
		if(!curEdge->End()->GVisited())
		{
			curEdge->End()->GVisited(true);
			vGList.push(curEdge->End());
		}
		curEdge = curEdge->Twin()->Next();
		while( curEdge->End()!= endV)
		{	
			if(!curEdge->End()->GVisited())
			{
				curEdge->End()->GVisited(true);
				vGList.push(curEdge->End());
			}
			curEdge = curEdge->Twin()->Next();
		}
	}
}*/

void Mesh::Move(const Vector3d p)
{
	for(int i=0;i<vList.size();i++)
		vList[i]->SetPosition(vList[i]->Position()+p);
}

void Mesh::Rotate(Quaternion q)
{
	for(int i=0;i<vList.size();i++)
	{
		Quaternion pos(vList[i]->Position());
		Quaternion res = q*pos*q.inv();
		Vector3d curPos = res.getPos();
		vList[i]->SetPosition(curPos);
	}
}
void Mesh::Noise(double max)
{
	for(int i=0;i<vList.size();i++)
	{
		Vector3d pos = vList[i]->Position();
		pos+= Vector3d(	(double)rand()/(RAND_MAX+1)*max,
						(double)rand()/(RAND_MAX+1)*max,
						(double)rand()/(RAND_MAX+1)*max);
		vList[i]->SetPosition(pos);
	}
}
void Mesh::ComputeFaceNormals()
{
	for (size_t i=0; i<fList.size(); i++) {
			if(fList[i]!=NULL && fList[i]->HalfEdge()->LeftFace()!=NULL )
			{
			Face *f = fList[i];
			const Vector3d & pos1 = f->HalfEdge()->Start()->Position();
			const Vector3d & pos2 = f->HalfEdge()->End()->Position();
			const Vector3d & pos3 = f->HalfEdge()->Next()->End()->Position();
			Vector3d normal = (pos2-pos1).Cross(pos3-pos1);
			normal /= normal.L2Norm();
        
			f->SetNormal_f(normal);//1007
			}
	}
}
void Mesh::Flip()
{
	for (size_t i = 0; i < fList.size(); i++)
	{
		HEdge *he = fList[i]->HalfEdge();
		Vertex *p1 = he->Start();
		Vertex *p2 = he->End();
		Vertex *p3 = he->Next()->End();

	}
}
bool Mesh::SaveObjFile(const char * filename)const
{
	std::ofstream fout(filename);
	fout << "g object\n";
	fout.precision(16);
	//output coordinates of each vertex
	for (int i=0; i<vList.size(); i++)
	{
		fout << "v " << std::scientific << vList[i]->Position().X()
			<< " " << vList[i]->Position().Y() << " " << vList[i]->Position().Z() << "\n";
	}

	// 		for (viter = pvertices_list_->begin();viter!=pvertices_list_->end(); viter++) 
	// 		{
	// 			fout<<"vn "<< std::scientific <<(*viter)->normal_.x() 
	// 				<<" "<<(*viter)->normal_.y() <<" "<<(*viter)->normal_.z() <<"\n";
	// 		}
	//output the valence of each face and its vertices_list' id
	for (int i=0; i < fList.size(); i++)
	{
		fout << "f";
		HEdge * he = fList[i]->HalfEdge();
		HEdge *start = he;
		do {
			fout << " " << he->Twin()->End()->Index()+1;
			he = he->Next();
		} while (he!=start);
		fout << "\n";
	}

	fout.close();
	return true;
}
int hashKey(HEdge* const&he)
{
	int * p = (int*)he;
	return *p;	
}
void Mesh::AutoSupports()
{
	FaceList hface;
	ComputeFaceNormals();
	for(int i=0;i<fList.size();i++)
	{
		Vector3d n = fList[i]->Normal_f();
		n /= n.L2Norm();
		double cosa = n.Dot(Vector3d(0.0,0.0,-1.0));
		if(cosa > 0.866)hface.push_back(fList[i]);
	}
	// merge areas
	closeHashTable<HEdge*> table(hface.size()*3,hface.size(),hashKey);
	for(int i=0;i<hface.size();i++)
	{
		HEdge *he = hface[i]->HalfEdge();
			table.insert(he);
			table.insert(he->Next());
			table.insert(he->Next()->Next());	
	}

	HEdgeList halfEdges = table.list();
	bool *boundary = new bool[halfEdges.size()];
	for(int i=0;i<halfEdges.size();i++) boundary[i] = false;
	DisjointSet ds(hface.size());
	for(int i=0;i<halfEdges.size();i++)
	{
		if(!table.find(halfEdges[i]->Twin())) boundary[i] = true;
		for(int j=0;j<halfEdges.size();j++)
		{
			if(i==j) continue;
			if(halfEdges[i]->Twin()==halfEdges[j]) ds.Union(ds.Find(i),ds.Find(j));
		}
	}	
	std::vector<int> types = ds.Types();
	for(int i=0;i<types.size();i++) std::cout << types[i] << ' ';
	// split hanging vertex
	
}


double Supports::dist(const Vector3d &v)const
{
	return (center-v).L2Norm();
}
void Supports::addVertex(Vertex* v)
{
	vList.push_back(v);
}
void Supports::addFace(Face* f)
{
	
}
