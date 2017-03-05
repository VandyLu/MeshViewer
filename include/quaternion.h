#ifndef __QUATERNION_H
#define __QUATERNION_H

#include <iostream>
#include "vector3.h"

using namespace std;

class Quaternion{
private:

public:
	double q0,q1,q2,q3;
	Quaternion(){}
	Quaternion(double _q0,double _q1,double _q2,double _q3):q0(_q0),q1(_q1),q2(_q2),q3(_q3){}
	Quaternion(double alpha,double x,double y,double z,int mode){
		Enter( alpha, x, y, z);
	}
	Quaternion(Vector3d v){ q0=0; q1=v.X();q2=v.Y();q3=v.Z();}                       //用于生成待旋转坐标
	Quaternion(double alpha,Vector3d v){ Enter(alpha,v.X(),v.Y(),v.Z());}//用于生成旋转四元数
	~Quaternion(){}
	void Enter(double alpha,double x,double y,double z)
	{
		double r = sqrt(x*x+y*y+z*z);
		x /=r; y/=r; z/=r;
		q0 = cos(alpha/2);
		double k = sin(alpha/2);
		q1 = k*x;
		q2 = k*y;
		q3 = k*z;		
	}
	void Normalize()
	{
		double len = Length();
		q0/=len; q1/=len; q2/=len;q3/=len;
	}
	Vector3d getPos()const {return Vector3d(q1,q2,q3);}
	double Length()const{return sqrt(q0*q0+q1*q1+q2*q2+q3*q3);}
	Quaternion inv()const{ return conj()/Length();}
	Quaternion conj()const{return Quaternion(q0,-q1,-q2,-q3);}
	Quaternion& operator=(const Quaternion& q){ q0 = q.q0; q1=q.q1; q2=q.q2;q3=q.q3;return *this;}
	Quaternion operator*(double a)const{ return Quaternion(q0*a,q1*a,q2*a,q3*a);}
	Quaternion& operator*=(double a){q0*=a;q1*=a;q2*=a;q3*=a; return *this;}
	Quaternion operator*(const Quaternion& q)const{
		return Quaternion(	q0*q.q0 - q1*q.q1 -q2*q.q2 -q3*q.q3,
							q0*q.q1 +q1*q.q0 + q3*q.q2 -q2*q.q3,
							q0*q.q2 +q2*q.q0 + q1*q.q3 -q3*q.q1,
							q0*q.q3 +q3*q.q0 + q2*q.q1 -q1*q.q2);
	}
	Quaternion& operator*=(const Quaternion& q){
		*this = Quaternion(	q0*q.q0 - q1*q.q1 -q2*q.q2 -q3*q.q3,
							q0*q.q1 +q1*q.q0 + q3*q.q2 -q2*q.q3,
							q0*q.q2 +q2*q.q0 + q1*q.q3 -q3*q.q1,
							q0*q.q3 +q3*q.q0 + q2*q.q1 -q1*q.q2);
		return *this;
	}
	Quaternion operator/(double a){return (*this)*(1.0/a);}
	Quaternion& operator/=(double a){ *this = (*this)/a;return *this;}
	void display()const
	{
		cout<<'['<<q0<<','<<q1<<','<<q2<<','<<q3<<']'<<endl; 
	}
};

#endif