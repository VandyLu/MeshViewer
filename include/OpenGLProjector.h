#ifndef __OPENGLPROJECTOR__
#define __OPENGLPROJECTOR__

#include "GL/glut.h"
#include <cstdlib>
#include "vector3.h"

class OpenGLProjector 
{
private:
	double modelView[16];
	double projection[16];
	int viewport[4];

    float* depthBuffer;
    
public:
	double* ModelViewMatrix() { return modelView; }
	double* ProjectionMatrix() { return projection; }
	int* Viewport() { return viewport; }

	OpenGLProjector()
	{
		glGetDoublev(GL_MODELVIEW_MATRIX, modelView);
		glGetDoublev(GL_PROJECTION_MATRIX, projection);
		glGetIntegerv(GL_VIEWPORT, viewport);
		
		depthBuffer = new float[viewport[2] * viewport[3]];

		glReadPixels(viewport[0], viewport[1], viewport[2], viewport[3], GL_DEPTH_COMPONENT, GL_FLOAT, depthBuffer);

	}
	Vector3d UnProject(double inX, double inY, double inZ)
	{
		double x,y,z;
		gluUnProject(inX, inY, inZ, modelView, projection, viewport, &x, &y, &z);
		return Vector3d(x,y,z);
	}
	Vector3d UnProject(Vector3d p)
	{
		double x,y,z;
		gluUnProject(p.X(), p.Y(), p.Z(), modelView, projection, viewport, &x, &y, &z);
		return Vector3d(x,y,z);
	}
	Vector3d Project(double inX, double inY, double inZ)
	{
		double x,y,z;
		gluProject(inX, inY, inZ, modelView, projection, viewport, &x, &y, &z);
		return Vector3d(x,y,z);
	}
	Vector3d Project(Vector3d p)
	{
		double x,y,z;
		gluProject(p.X(), p.Y(), p.Z(), modelView, projection, viewport, &x, &y, &z);
		return Vector3d(x,y,z);
	}
	double GetDepthValue(int x, int y) 
	{
		return depthBuffer[(y-viewport[1])*viewport[2] + (x-viewport[0])];
	}
};

#endif // __OPENGLPROJECTOR__
