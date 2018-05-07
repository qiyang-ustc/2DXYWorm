#include <vector>
#include<iostream>
#include"my_rng.h"
#include"my_vrbls.h"
#ifndef _MY_ALGO_
#define _MY_ALGO_
extern bool intersected(std::vector<int> up,std::vector<int> down);// return 1 if a ann b are intersected.

class Union_Find_int
{
public:
	Union_Find_int(int size)
	{
		labels.resize(size);
		for (int i = 0;i < size;i++)
		{
			labels[i] = i;
		}
	}
	Union_Find_int() {}
	void Refresh()
	{
		int size = labels.size();
		for (int i = 0;i < size;i++)
		{
			labels[i] = i;
		}
	}
	int getfather(int x)
	{
		int y=x;
		int z;
		while (y != labels[y])
		{
			y = labels[y];
		}
		while (x != labels[x])
		{
			z = labels[x];
			labels[x] = y;
			x = z;
		}//compress path;
		return y;
	}
	void Union(int x, int y)
	{
		if(x!=y)
		labels[getfather(x)] = getfather(y);
	}
	void Resize(int t)
	{
		labels.resize(t);
		this->Refresh();
	}
	void Print_labels(std::ostream& fout)
	{
		for (int i = 0;i < this->labels.size();i++)
		{
			fout << i<<"	";
		}
		fout << std::endl;
		for (int i = 0;i < this->labels.size();i++)
		{
			fout << labels[i]<< "	";
		}
		fout << std::endl;
	}
private:
	std::vector<int> labels;
};
class Cluster
{
public:
	Cluster() {
		this->p_block = new int*[Lx];
		for (int i = 0;i < Lx;i++)
		{
			(this->p_block)[i] = new int[Ly];
		}
		this->fresh();
		this->ident = 0;
	}
	~Cluster() {
		for (auto i = 0;i<Lx;++i) {
			delete[] p_block[i];
		}
		delete[] p_block;
	}
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
		return (this->p_block)[(i + Lx) % Lx][(j + Ly) % Ly];
	}//return p_block(i,j)
	void fresh()
	{
		for (int i = 0;i < Lx;i++)
		{
			for (int j = 0;j < Ly;j++)
			{
				p_block[i][j] = 0;
			}
		}
	}
	int ident;
private:
	int **p_block;
};
/*The following algorithms are used for breadth first search in Wolff algorithm*/
class Coordinates2D
{
public:
	Coordinates2D(int a, int b)
	{
		this->x = a;
		this->y = b;
	}
	int x;
	int y;
};

//Wolff End
#endif
