
#include "GL/glut.h"
#include "../include/mesh.h"
#include "../include/meshmatrix.h"
#include "../include/OpenGLProjector.h"
#include <fstream>//20110819


#define TEST_NORMAL
//#define DEMO
// 计算两个向量夹角 弧度制 带符号
#define sgn(a) (a>=0)
double getalpha(const Vector3d& p1, const Vector3d& p2, const Vector3d&n)
{
	return acos(p1.Dot(p2) / p1.L2Norm() / p2.L2Norm()) * (sgn(p1.Cross(p2).Z()) == sgn(n.Z()) ? 1 : -1);
}

bool selectDone =false;
Vertex* currentVertex;

// Enumeration
enum VIEWMESHNO{ No1=0,No2=1,ALL=2};
enum EnumDisplayMode { HIDDENLINE, FLATSHADED, SMOOTHSHADED, COLORSMOOTHSHADED,DELETESELECTEDVERTEX,MATCHBLOCKS};//WIREFRAME, 

VIEWMESHNO ViewMeshNo=No1;



enum Mode 
{ 
	Viewing, 
	Selection,
};
Mode currentMode = Viewing;


// variables
int displayMode = FLATSHADED;	// current display mode
int mainMenu, displayMenu;		// glut menu handlers
int winWidth, winHeight;		// window width and height
double winAspect;				// winWidth / winHeight;
int lastX, lastY;				// last mouse motion position
int currSelectedVertex = -1;         // current selected vertex
bool leftDown,leftUp, rightUp, rightDown, middleDown, middleUp, shiftDown;		// mouse down and shift down flags
double sphi = 90.0, stheta = 45.0, sdepth = 10;	// for simple trackball
double xpan = 0.0, ypan = 0.0;				// for simple trackball
double zNear = 1.0, zFar = 100.0;
double g_fov = 45.0;
Vector3d g_center;
double g_sdepth;
Mesh mesh1;	// our mesh1
//Mesh mesh2;
Mesh mesh2;

// functions
void SetBoundaryBox(const Vector3d & bmin, const Vector3d & bmax);
void InitGL();
void InitMenu();
void InitGeometry();
void MenuCallback(int value);
void ReshapeFunc(int width, int height);
void DisplayFunc();
void DrawWireframe();
void DrawHiddenLine();
void DrawFlatShaded();
void DrawSmoothShaded();
void DrawColorSmoothShaded();
void DrawSelectedVertices();
//void Partition();

void KeyboardFunc(unsigned char ch, int x, int y);
void MouseFunc(int button, int state, int x, int y);
void MotionFunc(int x, int y);
void SelectVertexByPoint(int mode=0);
void DeleteSelectedVertex(int vertex);


void SetBoundaryBox(const Vector3d & bmin, const Vector3d & bmax) {
	double PI = 3.14159265358979323846;
	double radius = bmax.Distance(bmin);
	g_center = 0.5 * (bmin+bmax);
	zNear    = 0.2 * radius / sin(0.5 * g_fov * PI / 180.0);
	zFar     = zNear + 2.0 * radius;
	g_sdepth = zNear + radius;
	zNear *= 0.1;
	zFar *= 10;
	sdepth = g_sdepth;
}

// init openGL environment
void InitGL() {
	GLfloat light0Position[] = { 0, 1, 0, 1.0 }; 

	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
	glutInitWindowSize(500, 500);
	glutCreateWindow("Comp541 Mesh Viewer");
	glClearColor(1.0, 1.0, 1.0, 1.0);
	glPolygonOffset(1.0, 1.0);
	glDepthFunc(GL_LEQUAL);
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_CULL_FACE);
	glEnable(GL_COLOR_MATERIAL);
	glColorMaterial(GL_FRONT, GL_DIFFUSE);
	glLightfv (GL_LIGHT0, GL_POSITION, light0Position);
	glEnable(GL_LIGHT0);
	// 接口函数 line 95 
	glutReshapeFunc(ReshapeFunc);
	glutDisplayFunc(DisplayFunc);
	glutKeyboardFunc(KeyboardFunc);
	glutMouseFunc(MouseFunc);
	glutMotionFunc(MotionFunc);
}

// init right-click menu
void InitMenu() {
	displayMenu = glutCreateMenu(MenuCallback);
	//glutAddMenuEntry("Wireframe", WIREFRAME);
	glutAddMenuEntry("Hidden Line", HIDDENLINE);
	glutAddMenuEntry("Flat Shaded", FLATSHADED);
	glutAddMenuEntry("Smooth Shaded", SMOOTHSHADED);
	glutAddMenuEntry("Color Smooth Shaded", COLORSMOOTHSHADED);
    //glutAddMenuEntry("Delete Selected Vertex", DELETESELECTEDVERTEX);
	mainMenu = glutCreateMenu(MenuCallback);
	glutAddSubMenu("Display", displayMenu);
	glutAddMenuEntry("Exit", 99);
	glutAttachMenu(GLUT_RIGHT_BUTTON);//glutAttachMenu(GLUT_MIDDLE_BUTTON);
}

// init geometry (if no input argument is provided)
void InitGeometry() {
	const int VSIZE = 4;
	const int HESIZE = 12;
	const int FSIZE = 4;
	int i;
	Vertex *v[VSIZE];
	HEdge *he[HESIZE];
	Face *f[FSIZE];
	
	for (i=0; i<VSIZE; i++) {
		v[i] = new Vertex();
		mesh1.vList.push_back(v[i]);
	}
	v[0]->SetPosition(Vector3d(0.0, 0.0, 0.0));
	v[1]->SetPosition(Vector3d(10.0, 0.0, 0.0));
	v[2]->SetPosition(Vector3d(0.0, 10.0, 0.0));
	v[3]->SetPosition(Vector3d(0.0, 0.0, 10.0));

	v[0]->SetNormal(Vector3d(-0.577, -0.577, -0.577));
	v[1]->SetNormal(Vector3d(0.0, -0.7, -0.7));
	v[2]->SetNormal(Vector3d(-0.7, 0.0, -0.7));
	v[3]->SetNormal(Vector3d(-0.7, -0.7, 0.0));

	for (i=0; i<FSIZE; i++) {
		f[i] = new Face();
		mesh1.fList.push_back(f[i]);
	}

	for (i=0; i<HESIZE; i++) {
		he[i] = new HEdge();
		mesh1.heList.push_back(he[i]);
	}
	for (i=0; i<FSIZE; i++) {
		int base = i*3;
		SetPrevNext(he[base], he[base+1]);
		SetPrevNext(he[base+1], he[base+2]);
		SetPrevNext(he[base+2], he[base]);
		SetFace(f[i], he[base]);
	}
	SetTwin(he[0], he[4]);
	SetTwin(he[1], he[7]);
	SetTwin(he[2], he[10]);
	SetTwin(he[3], he[8]);
	SetTwin(he[5], he[9]);
	SetTwin(he[6], he[11]);
	he[0]->SetStart(v[1]); he[1]->SetStart(v[2]); he[2]->SetStart(v[3]);
	he[3]->SetStart(v[0]); he[4]->SetStart(v[2]); he[5]->SetStart(v[1]);
	he[6]->SetStart(v[0]); he[7]->SetStart(v[3]); he[8]->SetStart(v[2]);
	he[9]->SetStart(v[0]); he[10]->SetStart(v[1]); he[11]->SetStart(v[3]);
	v[0]->SetHalfEdge(he[3]);
	v[1]->SetHalfEdge(he[0]);
	v[2]->SetHalfEdge(he[1]);
	v[3]->SetHalfEdge(he[2]);
}

// GLUT menu callback function
void MenuCallback(int value) {
	switch (value) {
	case 99: exit(0); break;
	default: 
		displayMode = value;
		glutPostRedisplay();
		break;
	}
}

// GLUT reshape callback function
void ReshapeFunc(int width, int height) {
	winWidth = width;
	winHeight = height;
	winAspect = (double)width/(double)height;
	glViewport(0, 0, width, height);
	glutPostRedisplay();
}

// GLUT display callback function
void DisplayFunc() {
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(g_fov, winAspect, zNear, zFar);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity(); 
	glTranslatef(xpan, ypan, -sdepth);
	glRotatef(-stheta, 1.0, 0.0, 0.0);
	glRotatef(sphi, 0.0, 1.0, 0.0);
	glTranslatef(-g_center[0], -g_center[1], -g_center[2]);

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	switch (displayMode) {
	//case WIREFRAME: DrawWireframe(); break;
	case HIDDENLINE: DrawHiddenLine(); break;
	case FLATSHADED: DrawFlatShaded(); break;
	case SMOOTHSHADED: DrawSmoothShaded(); break;
	case COLORSMOOTHSHADED: DrawColorSmoothShaded(); break;
	case DELETESELECTEDVERTEX: DeleteSelectedVertex(currSelectedVertex); break;
	
	}
	
	DrawSelectedVertices();

	//ComputeSS();//20110817

	glutSwapBuffers();
}

// Wireframe render function
void DrawWireframe() {
	HEdgeList heList = mesh1.Edges();
	HEdgeList bheList = mesh1.BoundaryEdges();
	glColor3f(0.3, 0.3, 1.0);
	glBegin(GL_LINES);
	size_t i;
	for (i=0; i<heList.size(); i++) {
		glVertex3dv(heList[i]->Start()->Position().ToArray());
		glVertex3dv(heList[i]->End()->Position().ToArray());
	}
    
	glColor3f(1, 0, 0);
	for (i=0; i<bheList.size(); i++) {
		glVertex3dv(bheList[i]->Start()->Position().ToArray());
		glVertex3dv(bheList[i]->End()->Position().ToArray());
	}
	
	glEnd();
}

// Hidden Line render function
void DrawHiddenLine() {
	FaceList fList = mesh1.Faces();
	glShadeModel(GL_FLAT); 
	glEnable(GL_POLYGON_OFFSET_FILL);
	glColor3f(0, 0, 0);
	glBegin(GL_TRIANGLES);
	for (size_t i=0; i<fList.size(); i++) {
		Face *f = fList[i];
		const Vector3d & pos1 = f->HalfEdge()->Start()->Position();
		const Vector3d & pos2 = f->HalfEdge()->End()->Position();
		const Vector3d & pos3 = f->HalfEdge()->Next()->End()->Position();
		glVertex3dv(pos1.ToArray());
		glVertex3dv(pos2.ToArray());
		glVertex3dv(pos3.ToArray());
	}
	glEnd();
	glDisable(GL_POLYGON_OFFSET_FILL);

	DrawWireframe();
}

// Flat Shaded render function
extern double plane[4];
void DrawFlatShaded() {
	FaceList fList;
	
	if(ViewMeshNo == No1 || ViewMeshNo==ALL)
	{
		fList = mesh1.Faces();
		glShadeModel(GL_FLAT); 
		glEnable(GL_LIGHTING);
		glColor3f(0.4f, 0.4f, 1.0f);
		glBegin(GL_TRIANGLES);
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
        
			glNormal3dv(normal.ToArray());
			glVertex3dv(pos1.ToArray());
			glVertex3dv(pos2.ToArray());
			glVertex3dv(pos3.ToArray());

			}
		}
		glEnd();
		glDisable(GL_LIGHTING);
	}
	if(ViewMeshNo == No2 || ViewMeshNo == ALL)
	{
		fList = mesh2.Faces();
		glShadeModel(GL_FLAT); 
		glEnable(GL_LIGHTING);
		glColor3f(0.4f, 0.4f, 1.0f);
		glBegin(GL_TRIANGLES);
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
        
			glNormal3dv(normal.ToArray());
			glVertex3dv(pos1.ToArray());
			glVertex3dv(pos2.ToArray());
			glVertex3dv(pos3.ToArray());
			}
		}
		glEnd();
		glDisable(GL_LIGHTING);
	}
	
}

// Smooth Shaded render function
void DrawSmoothShaded() { 
	FaceList fList;
	if(ViewMeshNo==No1 || ViewMeshNo==ALL )
	{
		fList = mesh1.Faces();
		glShadeModel(GL_SMOOTH); 
		glEnable(GL_LIGHTING);
		glColor3f(0.4f, 0.4f, 1.0f);
		glBegin(GL_TRIANGLES) ;
		for (size_t i=0; i<fList.size(); i++) {
			Face *f = fList[i];
			Vertex * v1 = f->HalfEdge()->Start();
			Vertex * v2 = f->HalfEdge()->End();
			Vertex * v3 = f->HalfEdge()->Next()->End();
			glNormal3dv(v1->Normal().ToArray());
			glVertex3dv(v1->Position().ToArray());
			glNormal3dv(v2->Normal().ToArray());
			glVertex3dv(v2->Position().ToArray());
			glNormal3dv(v3->Normal().ToArray());
			glVertex3dv(v3->Position().ToArray());
		}
		glEnd();
		glDisable(GL_LIGHTING);
	}
	if(ViewMeshNo==No2 || ViewMeshNo==ALL){
		fList = mesh2.Faces();
		glShadeModel(GL_SMOOTH); 
		glEnable(GL_LIGHTING);
		glColor3f(0.4f, 0.4f, 1.0f);
		glBegin(GL_TRIANGLES) ;
		for (size_t i=0; i<fList.size(); i++) {
			Face *f = fList[i];
			Vertex * v1 = f->HalfEdge()->Start();
			Vertex * v2 = f->HalfEdge()->End();
			Vertex * v3 = f->HalfEdge()->Next()->End();
			glNormal3dv(v1->Normal().ToArray());
			glVertex3dv(v1->Position().ToArray());
			glNormal3dv(v2->Normal().ToArray());
			glVertex3dv(v2->Position().ToArray());
			glNormal3dv(v3->Normal().ToArray());
			glVertex3dv(v3->Position().ToArray());
		}
		glEnd();
		glDisable(GL_LIGHTING);
	}
}

void DrawColorSmoothShaded() {
	FaceList fList;
	if(ViewMeshNo==No1 ||ViewMeshNo==ALL ){
		fList = mesh1.Faces();
		glShadeModel(GL_SMOOTH); 
		glEnable(GL_LIGHTING);
		glColor4f(0.4f, 0.4f, 1.0f,0.1);
		glBegin(GL_TRIANGLES) ;
		for (size_t i=0; i<fList.size(); i++) {
			Face *f = fList[i];
			Vertex * v1 = f->HalfEdge()->Start();
			Vertex * v2 = f->HalfEdge()->End();
			Vertex * v3 = f->HalfEdge()->Next()->End();
			glNormal3dv(v1->Normal().ToArray());
			glColor3dv(v1->Color().ToArray());
			glVertex3dv(v1->Position().ToArray());
			glNormal3dv(v2->Normal().ToArray());
			glColor3dv(v2->Color().ToArray());
			glVertex3dv(v2->Position().ToArray());
			glNormal3dv(v3->Normal().ToArray());
			glColor3dv(v3->Color().ToArray());
			glVertex3dv(v3->Position().ToArray());
		}
		glEnd();
		glDisable(GL_LIGHTING);

	}
	if(ViewMeshNo==No2 || ViewMeshNo==ALL)
	{
		fList = mesh2.Faces();
		glShadeModel(GL_SMOOTH);
		glEnable(GL_LIGHTING);
		glColor4f(0.4f, 0.4f, 1.0f, 0.1);
		glBegin(GL_TRIANGLES);
		for (size_t i = 0; i<fList.size(); i++) {
			Face *f = fList[i];
			Vertex * v1 = f->HalfEdge()->Start();
			Vertex * v2 = f->HalfEdge()->End();
			Vertex * v3 = f->HalfEdge()->Next()->End();
			glNormal3dv(v1->Normal().ToArray());
			glColor3dv(v1->Color().ToArray());
			glVertex3dv(v1->Position().ToArray());
			glNormal3dv(v2->Normal().ToArray());
			glColor3dv(v2->Color().ToArray());
			glVertex3dv(v2->Position().ToArray());
			glNormal3dv(v3->Normal().ToArray());
			glColor3dv(v3->Color().ToArray());
			glVertex3dv(v3->Position().ToArray());
		}
		glEnd();
		glDisable(GL_LIGHTING);
	}
}



//partition the model
/*void Partition();
{

	VertexList vList = mesh1.Vertices();
	glColor3f(1.0, 0.0, 0.0);
	glPointSize(10.0);
	glBegin(GL_POINTS);
	size_t i;

 
}
*/
// draw the selected ROI vertices on the mesh
void DrawSelectedVertices()
{
	VertexList vList;
	if(ViewMeshNo==0)vList= mesh1.Vertices();
	else if(ViewMeshNo==1)vList = mesh2.Vertices();
	glColor3f(1.0, 0.0, 0.0);
	glPointSize(10.0);
	glBegin(GL_POINTS);
	size_t i;

    for (i=0; i<vList.size(); i++) 
	{
		
		if (vList[i]->Flag())
		{
			switch (vList[i]->Flag() % 3)
			{
			case 0: // handle vertices
				glColor3f(1.0, 0.3, 0.3);
				break;
			case 1: 
				glColor3f(1.0, 0.0, 0.0); 
				break;
			case 2: 
				glColor3f(0.3, 0.3, 1.0); 
				break;
			}
			glVertex3dv(vList[i]->Position().ToArray());
            //cout<< vList[i]->Position().X() <<"  "<< vList[i]->Position().Y() <<"  "<<vList[i]->Position().Z() <<"  "<<endl; //20110816
		}
	}
	glEnd();
}

//delete selected vertex and its incident faces and half-edges

void DeleteSelectedVertex(int vertex)
{
     
	//case 1: vList[vertex] is not on any boundary triangle
	//case 2: vList[vertex] is on some boundary triangle, and it is not on any bounary loop
	//case 3: vList[vertex] is on some boundary triangle, and it is on some boundary loop

	VertexList vList = mesh1.Vertices();
	HEdge *cur = vList[vertex]->HalfEdge();
	HEdge *nex = cur;
    //delete the neighboring faces of v
	while (nex && nex->Prev()->Twin()!=cur) 
	{
	  if(nex->LeftFace()!=NULL)//take care of 3 cases
	  nex->SetFace(NULL);
      if(nex->Next()->LeftFace()!=NULL)//take care of 3 cases
	  nex->Next()->SetFace(NULL);
      if(nex->Prev()->LeftFace()!=NULL)//take care of 3 cases
	  nex->Prev()->SetFace(NULL);
	  nex = nex->Prev()->Twin();
	}
	  //process the last "nex"
	  if(nex->LeftFace()!=NULL)//take care of 3 cases
	  nex->SetFace(NULL);
      if(nex->Next()->LeftFace()!=NULL)//take care of 3 cases
	  nex->Next()->SetFace(NULL);
      if(nex->Prev()->LeftFace()!=NULL)//take care of 3 cases
	  nex->Prev()->SetFace(NULL);
    //cur->Twin()->SetFace(NULL);
	//cur->Twin()->Next()->SetFace(NULL);
	//cur->Twin()->Prev()->SetFace(NULL);
	//if(cur->Twin()->Prev()->LeftFace()==NULL) //test whether the left is really null
    cout<< "adjacent faces are deleted and will not be rendered"<<endl;
	
	//for processing case 2, mark the nodes to be deleted on the surrounding loop of the selected vertex
    HEdge *next_0 = cur;
	while(next_0 && next_0->Prev()->Twin()!=cur) 
	{
	  next_0=next_0->Prev()->Twin();
	  if(next_0->Next()->Twin()->IsBoundary() && next_0->Twin()->Prev()->Twin()->IsBoundary())
	  {
	    next_0->End()->mark=1;//mark this vertex as to be deleted later in the next phase
	  }
	  if(!next_0->Next()->Twin()->IsBoundary() && !next_0->Twin()->Prev()->Twin()->IsBoundary())
	  {next_0->End()->mark=0;}//mark this vertex as not on any boundary loop
	
	}
	//process the last next_0
    if(next_0->Next()->Twin()->IsBoundary() && next_0->Twin()->Prev()->Twin()->IsBoundary())
	 {
	   next_0->End()->mark=1;//mark this vertex as to be deleted later in the next phase
	 }
	 if(!next_0->Next()->Twin()->IsBoundary() && !next_0->Twin()->Prev()->Twin()->IsBoundary())
	 {next_0->End()->mark=0;}//mark this vertex as not on any boundary loop
    //process cur
	 if(cur->Next()->Twin()->IsBoundary() && cur->Twin()->Prev()->Twin()->IsBoundary())
	 {
	   cur->End()->mark=1;//mark this vertex as to be deleted later in the next phase
	 }
	 if(!cur->Next()->Twin()->IsBoundary() && !cur->Twin()->Prev()->Twin()->IsBoundary())
	 {cur->End()->mark=0;}//mark this vertex as not on any boundary loop
    

    //initialize
	int no_case =0;
    HEdge *next_3 = cur;

	while (next_3 && next_3->Prev()->Twin()!=cur)//determine case 3 or "1 and 2"
	{
	   if (next_3->IsBoundary() || next_3->Twin()->IsBoundary())
	   {no_case =3;break;}//mark the half edge emanating out of the selcted vertex, such that wither itself or its twin is a boundary edge
       next_3=next_3->Prev()->Twin();
	}
	//process the last "next_3"
    if (next_3->IsBoundary() || next_3->Twin()->IsBoundary())
	{no_case =3;}//mark the half edge emanating out of the selected vertex, such that wither itself or its twin is a boundary edge
	
	if(no_case ==3) //case 3
	{
	  HEdge *current = cur;
      while (current && current->Prev()->Twin()!= cur)
	  {    
		 SetPrevNext(current->Twin()->Prev(), current->Next());
		 current=current->Prev()->Twin();
		 if(current->Twin()->Next()!=cur)
		 {current->Twin()->Next()->SetTwin(NULL);current->Twin()->SetNext(NULL);}
		 //current->SetTwin(NULL); //conflict 111

	  }
	  //process the last "current" 
	  SetPrevNext(current->Twin()->Prev(), current->Next());
      current->SetTwin(NULL);
	  current=NULL;
	  //process "cur"
	  cur->SetTwin(NULL); 
	  cur=NULL;	
	}//end if

	else//case 1 and 2 
	{
    //delete half edges opposite to the selected vertex
	HEdge *next_1 = cur;
	while(next_1 && next_1->Prev()->Twin()!=cur) 
	{
      if(next_1->End()->mark==1 ||  next_1->Twin()->Prev()->Twin()->IsBoundary())//next_1->End() is on some boundary loop
	  {

		if(next_1->End()->mark==1)
        next_1->End()->mark=2;//means one side of the vertex has been cleared

		if(!next_1->Next()->Twin()->IsBoundary() &&  next_1->Twin()->Prev()->Twin()->IsBoundary())//this vertex is not deleted, but the prev and Next relation is rebuilt
		{
		  SetPrevNext(next_1->Twin()->Prev()->Twin()->Prev(),next_1->Next());
		}
        
		//delete the marked vertex
		if(next_1->Twin()->Prev()->Start()->mark==2)
        next_1->Twin()->Prev()->SetStart(NULL);
		//delete boundary the half edge and its twin
		next_1->Twin()->Prev()->SetTwin(NULL);
		next_1->Twin()->SetPrev(NULL);
	  }//end if 
	  if(next_1->End()->mark==0)//next_1->End() is not on any boundary loop
	  {
	    SetPrevNext(next_1->Twin()->Prev(),next_1->Next());
	  }
      //delete half edges incident to the selected vertex
	  //next_1->SetTwin(NULL);
	  next_1=next_1->Prev()->Twin();

	  if(next_1->Twin()->Next()!=cur)
	  {
		next_1->Twin()->Next()->SetTwin(NULL);next_1->Twin()->SetNext(NULL);

	  } 
      
	}//end while 
	
    //process the last next_1 below
  
     if(next_1->End()->mark==1 ||  next_1->Twin()->Prev()->Twin()->IsBoundary())//next_1->End() is on some boundary loop
	  {

		if(next_1->End()->mark==1)
        next_1->End()->mark=2;//means one side of the vertex has been cleared

		if(!next_1->Next()->Twin()->IsBoundary() &&  next_1->Twin()->Prev()->Twin()->IsBoundary())//this vertex is not deleted, but the prev and Next relation is rebuilt
		{
		  SetPrevNext(next_1->Twin()->Prev()->Twin()->Prev(),next_1->Next());
		}
        
		//delete the marked vertex
		if(next_1->Twin()->Prev()->Start()->mark==2)
        next_1->Twin()->Prev()->SetStart(NULL);
		//delete boundary the half edge and its twin
		next_1->Twin()->Prev()->SetTwin(NULL);
		next_1->Twin()->SetPrev(NULL);
	 }//end if 
	  if(next_1->End()->mark==0)//next_1->End() is not on any boundary loop
	  {
	    SetPrevNext(next_1->Twin()->Prev(),next_1->Next());
	  }
      //delete half edges incident to the selected vertex
	  //next_1->SetTwin(NULL);
	  next_1=next_1->Prev()->Twin();

	  if(next_1->Twin()->Next()!=cur)
	  {
		next_1->Twin()->Next()->SetTwin(NULL);next_1->Twin()->SetNext(NULL);

	  }  

      //process "next_1" and "cur"
		next_1->SetTwin(NULL);
		next_1=NULL;
		cur->SetTwin(NULL);
		cur=NULL;

    }//end else
    vList[vertex]=NULL; //delete the selected vertex
	if(vList[vertex]==NULL)
	cout<< "the selected vertex is gone, its index is "<<vertex<<endl;
} 

void testRotate()
{
	Vertex pt;
	pt.SetPosition(Vector3d(1,1,1));
	Quaternion rotate(3.14/2,0,0,1,1);
	cout << pt.Position()<<endl;
	Quaternion abc = rotate.inv();

	Quaternion pos(pt.Position());
	Quaternion res =rotate*pos;
	res = res*rotate.inv();
	Vector3d curPos = res.getPos();
	pt.SetPosition(curPos);
	cout << pt.Position().X()<<","<< pt.Position().Y()<<","<< pt.Position().Z() <<endl;
}

bool getTwoBlock = false;
 //DWORD WINAPI QThread(void* pParam)

int Query(void)
{
	currentMode = Selection;
	ViewMeshNo=No1;
	cout << "Choose Block 1"<< endl;
	selectDone = false;
	while(!selectDone);
	ViewMeshNo=No2;
	cout << "Choose Block 2"<< endl;
	selectDone = false;
	while(!selectDone);
	selectDone = false;
	getTwoBlock = true;
	return 1;
}

// GLUT keyboard callback function
void KeyboardFunc(unsigned char ch, int x, int y) { 
	switch (ch) { 
	//case 'u':
	case '3':
		/************************************************************************/
		/* activate the following code if you finish the corresponding functions*/
// 		mesh1.UmbrellaSmooth();
 		mesh1.ComputeVertexNormals();
		cout<<"the normals of vertices are computed already"<<endl;
// 		mesh1.ComputeVertexCurvatures();
		/************************************************************************/
		break; 
	
	case '4':
		/************************************************************************/
		/* activate the following code if you finish the corresponding functions*/
// 		mesh1.ImplicitUmbrellaSmooth();
// 		mesh1.ComputeVertexNormals();
 		mesh1.ComputeVertexCurvatures(); 
		cout<<"the Curvatures of vertices are computed already"<<endl;
		/************************************************************************/
		break;
	case '5':
        mesh1.UmbrellaSmooth();
        break;
	case '6':
        mesh1.ImplicitUmbrellaSmooth();
        break;

	// edit by lfb
	case 's':
		if(ViewMeshNo==No1)
		{
			ViewMeshNo=No2;
			cout <<"Current Mesh is Mesh2"<<endl;
		}else if(ViewMeshNo==No2)
		{
			ViewMeshNo=ALL;
			cout <<"Current Mesh is ALL"<<endl;
		}else if(ViewMeshNo==ALL)
		{
			ViewMeshNo=No1;
			cout << "Current Mesh is Mesh1"<<endl;
		}
		break;
	case 'r':
		// rotate the model
		mesh1.Rotate(Quaternion(3.14/4,Vector3d(0,0,1)));
		break;
	
	case 't':
		mesh1.Move(Vector3d(0.1,0.1,0.1));
		break;
	case 'p':

		mesh1.ComputeFaceNormals();
		mesh1.ComputeVertexNormals();
		mesh1.AutoSupports();
		break;
	case 'j':
		
		currentMode = Selection;

		mesh1.ComputeVertexCurvatures();
		mesh1.ComputeFaceNormals();
		mesh1.ComputeVertexNormals();

		mesh2.ComputeVertexCurvatures();
		mesh2.ComputeFaceNormals();
		mesh2.ComputeVertexNormals();

		//上色
		{	
		int threshold = 300;
		Vector3d keyColor(0.0,1.0,0.0);
		Vector3d defaultColor(0.0, 0.0, 1.0);
		for (int i = 0; i < mesh1.Vertices().size(); i++)
			if (mesh1.Vertices()[i]->G > threshold)
				mesh1.Vertices()[i]->SetColor(keyColor);
			else mesh1.Vertices()[i]->SetColor(defaultColor);

		for (int i = 0; i < mesh2.Vertices().size(); i++)
			if (mesh2.Vertices()[i]->G > threshold)
				mesh2.Vertices()[i]->SetColor(keyColor);
			else mesh2.Vertices()[i]->SetColor(defaultColor);
			//cout << "mesh2 color:\n\t" << mesh2.Vertices()[0]->Color();
			//cout << "mesh2 position:\n\t" << mesh2.Vertices()[0]->Position();
		}
		cout << "preprocess done!" << endl;
		break;

	case 'm':	
	{

		// 选择分别三个点
		int pt1, pt2, pt3, pta, ptb, ptc;
		cout << "Mesh1 three vertices" << endl;
		cin >> pt1 >> pt2 >> pt3;
		cout << "Mesh2 three vertices" << endl;
		cin >> pta >> ptb >> ptc;

		Vector3d p1 = mesh1.vList[pt1]->Position();
		Vector3d p2 = mesh1.vList[pt2]->Position();
		Vector3d p3 = mesh1.vList[pt3]->Position();
		Vector3d pa = mesh2.vList[pta]->Position();
		Vector3d pb = mesh2.vList[ptb]->Position();
		Vector3d pc = mesh2.vList[ptc]->Position();

		// 计算重心
		Vector3d g1 = (p1 + p2 + p3) / 3;
		Vector3d g2 = (pa + pb + pc) / 3;
		//cout << g1 << endl << g2 << endl;

		// 旋转至法向重合
		Vector3d p12 = p2 - p1;
		Vector3d p13 = p3 - p1;
		Vector3d pab = pb - pa;
		Vector3d pac = pc - pa;
		Vector3d n1 = p12.Cross(p13);
		Vector3d n2 = pab.Cross(pac);
		Vector3d axis = n1.Cross(n2); axis /= axis.L2Norm(); //单位向量 旋转轴 为 n1 x n2
		double angle = getalpha(n1, n2, n1.Cross(n2));
		cout << "axis:" << endl << "\t" << axis << endl;
		cout << "angle:" << endl << "\t" << angle << endl;
		mesh1.Rotate(Quaternion(-angle, axis));
		// update position
		p1 = mesh1.vList[pt1]->Position();
		p2 = mesh1.vList[pt2]->Position();
		p3 = mesh1.vList[pt3]->Position();
		pa = mesh2.vList[pta]->Position();
		pb = mesh2.vList[ptb]->Position();
		pc = mesh2.vList[ptc]->Position();
		p12 = p2 - p1;
		p13 = p3 - p1;
	    pab = pb - pa;
		pac = pc - pa;
		n1 = p12.Cross(p13);
		n2 = pab.Cross(pac);
		g1 = (p1 + p2 + p3) / 3;
		g2 = (pa + pb + pc) / 3;
		cout << n1 - n2 << endl;  //should be all zero

		Vector3d axis2 = n2; // 未移动过 因此直接选择n2为转轴

		// method 1.simple

		double b1 = getalpha(p1 - g1, pa - g2, n2);
		double b2 = getalpha(p2 - g1, pb - g2, n2);
		double b3 = getalpha(p3 - g1, pc - g2, n2);
		double bavg = -(b1 + b2 + b3)/3;
		mesh1.Rotate(Quaternion(bavg, axis2));
		cout << "\n\naxis2:\n\t" << axis2 << endl;
		cout << "beta angles:\n\t";
		cout << b1 << ',' << b2 << ',' << b3 << endl;
		cout << "beta avg(angle2):\n\t" << bavg<<endl;
		
		
		// 平移至重心重合
		p1 = mesh1.vList[pt1]->Position();
		p2 = mesh1.vList[pt2]->Position();
		p3 = mesh1.vList[pt3]->Position();
		pa = mesh2.vList[pta]->Position();
		pb = mesh2.vList[ptb]->Position();
		pc = mesh2.vList[ptc]->Position();
		p12 = p2 - p1;
		p13 = p3 - p1;
		pab = pb - pa;
		pac = pc - pa;
		n1 = p12.Cross(p13);
		n2 = pab.Cross(pac);
		g1 = (p1 + p2 + p3) / 3;
		g2 = (pa + pb + pc) / 3;
		cout << "mesh1 move:" << g2 - g1 << endl;
		mesh1.Move(g2 - g1);

		char choice='n';
		cout << "save models?(y/n)" << endl;
		cin >> choice;
		if (choice == 'y' || choice == 'Y')
		{
			mesh1.SaveObjFile("mesh1.obj");
			mesh2.SaveObjFile("mesh2.obj");
		}
		break;
	}
	case 'n':
		mesh1.Noise(0.00005);
		break;
	//case '7':
		//ComputeSS();//20110918
        //DrawSS();
        //break;

	case '1':	// key '1'
		currentMode = Viewing;
		cout << "Viewing mode" << endl;
		break;
	case '2':	// key '2'
		currentMode = Selection;
		cout << "Selection mode" << endl;
		break;

	case '9': 
		DeleteSelectedVertex(currSelectedVertex);
		break;
	case 27:
		exit(0);
		break;
	}
	glutPostRedisplay();
}

// GLUT mouse callback function
void MouseFunc(int button, int state, int x, int y) {
	
	lastX = x;
	lastY = y;
	leftDown = (button == GLUT_LEFT_BUTTON) && (state == GLUT_DOWN);
	leftUp = (button == GLUT_LEFT_BUTTON) && (state == GLUT_UP);
    rightDown = (button == GLUT_RIGHT_BUTTON) && (state == GLUT_DOWN);
	rightUp = (button == GLUT_RIGHT_BUTTON) && (state == GLUT_UP);
	middleDown = (button == GLUT_MIDDLE_BUTTON) && (state == GLUT_DOWN);
	middleUp = (button == GLUT_MIDDLE_BUTTON) && (state == GLUT_UP);
	shiftDown = (glutGetModifiers() & GLUT_ACTIVE_SHIFT);

	if (currentMode == Selection && state == GLUT_UP)
	{
	    if (middleUp)
		{
		    if (currSelectedVertex != -1)
		    {
		        mesh1.Vertices()[currSelectedVertex]->SetFlag(0);
		        currSelectedVertex = -1;
		    }
		}
		else 
		{
			SelectVertexByPoint();
			if(currSelectedVertex && ViewMeshNo!=ALL)
			{
				const char *strout = (ViewMeshNo == No1 ? "Mesh No1:" : "Mesh No2:");
				cout << strout << "\n\t Vertex ID:"<<currSelectedVertex << endl;
			}
			//cout << mesh1.Vertices()[currSelectedVertex]->G<<endl;
		}
		//if (leftUp)
	    //{
		  //  SelectVertexByPoint(0); 
			
		//}
        //if(middleUp)
		//{DeleteSelectedVertex(currSelectedVertex);}
        
		lastX = lastY = 0;
		glutPostRedisplay();
	}
}

// GLUT mouse motion callback function
void MotionFunc(int x, int y) {
	if (leftDown)
		if(!shiftDown) { // rotate
			sphi += (double)(x - lastX) / 4.0;
			stheta += (double)(lastY - y) / 4.0;
		} else { // pan
			xpan += (double)(x - lastX)*sdepth/zNear/winWidth;
			ypan += (double)(lastY - y)*sdepth/zNear/winHeight;
		}
	// scale
	if (middleDown) sdepth += (double)(lastY - y) / 10.0;

	lastX = x;
	lastY = y;
	glutPostRedisplay();
}


// select a mesh point
void SelectVertexByPoint(int mode)
{

	// get the selection position
	int x = lastX, y = winHeight - lastY;
	Vector3d u(x,y,0);
	std::cout << x << ',' << y << std::endl;
	OpenGLProjector projector; 

	VertexList vList;
	if(mode==0 && ViewMeshNo==No1)
		vList = mesh1.Vertices();
	else if(mode==0 && ViewMeshNo==No2)
		vList = mesh2.Vertices();
	//else vList = SC.getVertexList();

    double mindis = 1e6; int selectedIndex = -1;
	for (size_t i=0; i<vList.size(); i++) 
	{
		Vector3d v = projector.Project(vList[i]->Position());
		
		vList[i]->SetFlag(0); // set back to unselected
		
		if (projector.GetDepthValue((int)v.X(), (int)v.Y()) - v.Z() <= -1e-8)
		{
		    continue;
		}
		v.Z() = 0;
		
		double dist = (u - v).L2Norm();
		if (dist < mindis)
		{
		    mindis = dist;
		    selectedIndex = (int)i;
		}
	}
	
	if (selectedIndex != -1)
	{
	    currSelectedVertex = selectedIndex;
	    vList[selectedIndex]->SetFlag(1);
		Vector3d aaa = (projector.Project(vList[selectedIndex]->Position()));
		aaa.Z() = 0;
		std::cout << mindis << std::endl;
		std::cout << '\n' << aaa << '\n' << (u - aaa).L2Norm() << std::endl;
		std::cout << "the selected vertex index is " << currSelectedVertex << std::endl;
        //cout <<"the slected vertex index is " <<  currSelectedVertex <<endl;
	}
	/*if(ViewMeshNo==0)
	{
		currentVertex = vList[selectedIndex];
	}else if(ViewMeshNo==1)
	{
	}*/
	selectDone = true;

	DrawSelectedVertices();
}

#ifndef DEMO
// main function
int main(int argc, char **argv) {

	glutInit(&argc, argv);
	InitGL();
	InitMenu();
#ifdef TEST_NORMAL	
	if (argc==2)
	{
		mesh1.LoadObjFile(argv[1]);
	}
	else if (argc==3)
	{
		mesh1.LoadObjFile(argv[1]);
		mesh2.LoadObjFile(argv[2]);
	}
	else InitGeometry();
#endif
#ifdef TEST_MYOUT
	mesh1.LoadObjFile("mesh1");
	mesh2.LoadObjFile("mesh2");
#endif
#ifdef TEST_BU
	mesh1.LoadObjFile("/home/lfb/Documents/models/bu.obj");
	mesh2.LoadObjFile("/home/lfb/Documents/models/bu.obj");
#endif
	SetBoundaryBox(mesh1.MinCoord(), mesh1.MaxCoord());
	
	/************************************************************************/
	/* activate the following code if you finish the corresponding functions*/
 	mesh1.DisplayMeshInfo();
	/************************************************************************/

	/* for arch 
	Cube c(Vector3d(-0.25,-0.25,-0.5),Vector3d(0.5,0.5,1));
	s.setCube(c);
	s.generate();
	*/
	glutMainLoop();
	return 0;
}
#endif
