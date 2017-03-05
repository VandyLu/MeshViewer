/**
* Sparse Matrix
*/

#ifndef __MESHMATRIX_H__
#define __MESHMATRIX_H__

#include <cstdlib>
#include <vector>
#include <algorithm>



// entry in a matrix
struct MatrixElement
{
	int row, col;
	double value;

	MatrixElement(int r, int c, double v)
		: row(r), col(c), value(v) 
	{}

	// upper left elementes are "smaller"
	static bool order (MatrixElement e1, MatrixElement e2)
	{
		if (e1.row < e2.row) return true;
		if (e1.row == e2.row) return (e1.col < e2.col);
		return false;
	}
};

typedef std::vector<MatrixElement> MatrixElementList;

// Sparse Matrix
class Matrix 
{
private:
	int m, n;
	MatrixElementList elements;


	// stores the index of each element in a row
	int * rowIndex; 



	// fields for CG method
	double* diagInv;
	double* r;
	double* r2;
	double* d;
	double* d2;
	double* q;
	double* s;
	double* s2;
	double* tmp;

public:
	Matrix(int m, int n) : m(m), n(n) 
	{
		rowIndex = new int[m+1]; 
		diagInv = new double[m];
		r = new double[m];
		r2 = new double[m];
		d = new double[m];
		d2 = new double[m];
		q = new double[m];
		s = new double[m];
		s2 = new double[m];
		tmp = new double[m];
	}
	~Matrix() 
	{
		delete[] rowIndex; 
		delete[] r;
		delete[] r2;
		delete[] d;
		delete[] d2;
		delete[] q;
		delete[] s;
		delete[] s2;
	}

	void AddElement(int row, int col, double value)
	{
		elements.push_back(MatrixElement(row, col, value));
	}
	int RowSize() const { return m; }
	int ColSize() const { return n; }

	const MatrixElementList & Elements() const { return elements;}

	int * RowIndex() const { return rowIndex; }

	void SortMatrix()
	{
		sort(elements.begin( ), elements.end( ), MatrixElement::order);

		for (int i=0; i<m+1; i++)
			rowIndex[i] = 0;
		for (int i=0; i<(int)elements.size(); i++)
			rowIndex[elements[i].row + 1] = i + 1;

		for (int i=0; i<m; i++)
			diagInv[i] = 0;
		for (int i=0; i<(int)elements.size(); i++)
			if (elements[i].row == elements[i].col)
				diagInv[elements[i].row] = 1.0 / elements[i].value;
	}

/*
//fetch element value
double GetElement(int row_1, int col_1)//1013
{
return elements[row_1*n + col_1].value;
}
*/
	void Multiply(double* xIn, double* xOut)//Suppose the Matrix is M, then xOut = M*xIn;
	{
		for (int i=0; i<m; i++)//IN OUR SETTING, XOUT HAS 3*M ELEMENTS
		{
			double sum = 0;
			for (int j=rowIndex[i]; j<rowIndex[i+1]; j++)
				sum += elements[j].value * xIn[elements[j].col];
			xOut[i] = sum;
		}
	}

	void PreMultiply(double* xIn, double* xOut)//Suppose the Matrix is M, then xOut = xIn*M
	{
		for (int i=0; i<n; i++) xOut[i] = 0;

		for (int i=0; i<m; i++)
		{
			for (int j=rowIndex[i]; j<rowIndex[i+1]; j++)
				xOut[elements[j].col] += elements[j].value * xIn[i];
		}
	}


	// solve the linear system Ax = b by using Bi-Conjugete Gradient Method
	void BCG(double* b, double* x, int maxIter, double tolerance)
	{
		int iter;
		double errNew;
		double errOld;
		double err;

		iter = 0;
		Multiply(x, r);
		// initialize the variables
		for (int i=0; i<m; i++) 
			d2[i] = d[i] = r2[i] = r[i] = b[i] - r[i];
		errNew = 0;
		for (int i=0; i<m; i++) errNew += r2[i] * r[i];
		err = abs(errNew);

		while (iter<maxIter && abs(errNew) > tolerance*err)
		{
			cout << errNew << " " << tolerance*err << endl;

			Multiply(d, q);

			double alpha = 0;
			for (int i=0; i<m; i++) alpha += d2[i] * q[i];
			alpha = errNew / alpha;
			for (int i=0; i<m; i++) x[i] += alpha * d[i];

			if (iter % 50 == 0)
			{
				Multiply(x, r);
				for (int i=0; i<m; i++) r[i] = b[i] - r[i];
				PreMultiply(d2, q);
				for (int i=0; i<m; i++) r2[i] -= alpha * q[i];
			}
			else
			{
				for (int i=0; i<m; i++) r[i] -= alpha * q[i];
				PreMultiply(d2, q);
				for (int i=0; i<m; i++) r2[i] -= alpha * q[i];
			}

			errOld = errNew;
			errNew = 0;
			for (int i=0; i<m; i++) errNew += r2[i] * r[i];
			double beta = errNew / errOld;
			for (int i=0; i<m; i++) d[i] = r[i] + beta * d[i];
			for (int i=0; i<m; i++) d2[i] = r2[i] + beta * d2[i];
			iter++;
		}

		cout << errNew << " " << tolerance*err << endl;
		cout << "iter: " << iter << endl;
	}


   void VecSubtract(double* lhs, double* rhs, double* result, int n) 
   {
	 for (int i = 0; i < n; i++) 
	  {
		result[i] = lhs[i] - rhs[i];
	  }
   }
	
   double VecDot(double* rhs, double* lhs, int n) 
   {
	 double r = 0;
	 for (int i = 0; i < n; i++) 
	   {
		r += rhs[i] * lhs[i];
	   }
	 return r;
   }
 

	friend ostream & operator<< (ostream & out, const Matrix & r) 
	{
		for (int i=0; i<r.m; i++)
		{
			for(int j=r.rowIndex[i]; j<r.rowIndex[i+1]; j++)
				out << r.elements[j].value << " ";
			out << endl;
		}
		return out;
	}
};

#endif __MATRIX_H__