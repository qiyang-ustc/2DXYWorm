
//*******************************************************************
// Ising model on the square Lattice

// Error bars are calculated using the blocking technique.
// Let ' \chi^2 = \sum_{t=1}^{T} <(O_t -<O>)^2> ' be the squared chi for
// 'T' blocks of observable 'O'. Assuming each block of data is independent
// of each other, the error bar of 'O' is given by 'Dev = \sqrt{\chi^2/T(T-1)}.

// Reliabity of the obtained errors is monitored by t=1 correlation,
// for which tolerance is set by variable 'tol' (default: tol=0.20d0).

// Composite quantities like Binder ratios are calculated in each block, and
// the associated error bars are obtained from their fluctuations.

// Results are written into a special file 'dat.***' if the number of
// blocks is less than 125 or correlation is too big. Data in each
// block will be also printed out in this case.

// Default number of extensive simulation is 'NBlck=1024'.

// For test purpose, for which huge amount of information will be
// printed out, 'NBlck' should be set smaller but >2.

// Dynamical behavior is not studied.

// 'my_vrbls.f90', 'carlo.f90', 'monte.f90', 'measure.f90',
// 'write2file.f90' and etc need to be modified for new projects.

//  Look for 'PROJECT-DEPENDENT'.

//  Author: Yuan Huang
//  Date  : April 19th, 2012.
//*******************************************************************

// Look for 'PROJECT-DEPENDENT' for different projects
#include <vector>
#include<iostream>
#include"my_rng.h"
#include <math.h>
#include "bessi.h"
#ifndef _MY_VRBLS_
#define _MY_VRBLS_
//------------------Project-Independed---------
extern const double tm32;	 //1/2^(32)
extern const double eps;     //very small number
extern const double tol;	 //tolerance for correlation 0.2
extern const int MaxInt;	 //Maximun integer
extern const int MinInt;	 //Minimun integer
extern double Pi;      //PI
extern const int MaxBlock;	     //Maximun number of block
extern const int MinBlock;		 //Minimun number of block
extern int NBlock;               //N=number of
extern int NSample;				// samples in a block
extern int TotalSample;
extern int NToss;				//samples to be trown away
extern int Collect_Interval;
//----------------Project-depended-------------
extern const int Dimension;
extern int Lx, Ly;
extern double Jcp;
extern int Vol; extern double Vol2; extern double Vol4;
extern int ident;//used to identify cluster
extern const int NObs;
extern const int NQuan;
extern int steps;
template<class T>
/* Some algorithm which may be useful*/
void Delete_Num(typename std::vector<T>& v, T del)
//Can not Delete Points!!!
{
	for (typename std::vector<T>::iterator iter = v.begin();iter != v.end();)
	{
		if (del == *iter)
		{
			iter = v.erase(iter);
		}
		else
		{
			++iter;
		}
	}
}
//----------ALGORITHM------END!!!

//---------------Block-----------------------
class Block {
public:
	std::vector<double> Quan;
	Block()
	{
		this->p_block = new int*[Lx];
		for (int i = 0;i < Lx;i++)
		{
			(this->p_block)[i] = new int[Ly];
		}

		this->left_block = new int*[Lx];
		for (int i = 0;i < Lx;i++)
		{
			(this->left_block)[i] = new int[Ly];
		}

		this->down_block = new int*[Lx];
		for (int i = 0;i < Lx;i++)
		{
			(this->down_block)[i] = new int[Ly];
		}
		this->fresh();
		Quan.resize(NObs);
		//initialize temp  which is the hash function.
	}
	~Block() {}
	void Print_ele(std::ostream &fout)
	{
		for (int i = 0;i < Lx;i++)
		{
			for (int j = 0;j < Ly;j++)
			{
				fout << this->ele(i, j) << "	";
			}
			fout << std::endl;
		}
	}//out put block configuration.
	int& ele(int i, int j)
	{
		return (this->p_block)[(i+Lx)%Lx][(j+Ly)%Ly];
	}//return p_block(i,j)
	int& left_flow(int i,int j)
	{	
		return (this->left_block)[(i + Lx) % Lx][(j + Ly) % Ly];
	}
	int& down_flow(int i, int j)
	{
		return (this->down_block)[(i + Lx) % Lx][(j + Ly) % Ly];
	}
	void fresh()
	{
		for (int i = 0;i < Lx;i++)
			for (int j = 0;j < Ly;j++)
			{
				p_block[i][j] = 0;
				left_block[i][j] = 0;
				down_block[i][j] = 0;
			}
	}
	void Sample()
	{
		this->Quan[0] = -Jcp*(this->ele(1, 0) + this->ele(-1, 0) + this->ele(0, 1) + this->ele(0, -1)) / this->ele(0, 0)/2;
		int sum=0;
		for (int i = 0;i < Lx;i++)
		{
			for (int j = 0;j < Ly;j++)
			{
				sum += this->ele(i, j);
			}
		}
		this->Quan[1]= 1.0*sum/Vol/this->ele(0,0);
	}
	void Calculate_Quan()
	{
		this->Sample();
		/*
		0.<M>
		*/
	}
private:
	int **p_block;
	int **left_block;
	int **down_block;
};

class FlowWeight
{
public:
	FlowWeight()
	{
		this->length = 10;
		this->w.resize(this->length);
		for (int i = 0;i < this->length;i++)
		{
			w[i] = BESSI(i, Jcp);
		}
	}
	void Expand2(int n)
	{
		this->w.resize(n);
		for (int i = this->length;i < n;i++)
		{
			w[i] = BESSI(i, Jcp);
		}
		this->length = n;
	}
	double weight(int n)
	{
		if (abs(n) < this->length)
		{
			return w[abs(n)];
		}
		else
		{
			this->Expand2(2*abs(n));
			return w[abs(n)];
		}
	}
	int Print_ele(std::ostream &fout)
	{
		for (int i = 0;i < length;i++)
		{
			fout << w[i]<<"	";
		}
		fout << std::endl;
		return 0;
	}
private:
	int length;
	std::vector<double> w;

};

#endif
