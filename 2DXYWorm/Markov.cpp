#include "my_vrbls.h"
#include "my_rng.h"
#include <math.h>
#include "my_algorithm.h"
#include <list>
void Metropolis(Block& b)
{
	double dE;
	double bij;
	double temp; //random new spin
	for (int i = 0;i < Lx;i++)
	{
		for (int j = 0;j < Ly; j++)
		{
			temp = rn();
			bij = b.ele(i, j);
			dE = cos(2 * Pi*(b.ele(i + 1, j) - bij)) + cos(2 * Pi*(b.ele(i - 1, j) - bij)) + cos(2 * Pi*(b.ele(i, j + 1) - bij)) + cos(2 * Pi*(b.ele(i, j - 1) - bij))
				- (cos(2 * Pi*(b.ele(i + 1, j) - temp)) + cos(2 * Pi*(b.ele(i - 1, j) - temp)) + cos(2 * Pi*(b.ele(i, j + 1) - temp)) + cos(2 * Pi*(b.ele(i, j - 1) - temp)));
			if (dE<0 || exp(-Jcp*dE)>rn())
			{
				b.ele(i, j) = temp;
			}
		}
	}
}
void FindNeighbour(std::list<Coordinates2D>& q,std::list<Coordinates2D>::iterator& iter, Cluster& c, Block& b, double n)
{
	int x;int y;double temp1,temp2,temp3,temp4,temp;
	x = iter->x;
	y = iter->y;
	temp = cos(2*Pi*(b.ele(x, y)-n));
	temp1 = temp*cos(2 * Pi*(b.ele(x+1, y) - n));
	temp2 = temp*cos(2 * Pi*(b.ele(x-1, y) - n));
	temp3 = temp*cos(2 * Pi*(b.ele(x, y+1) - n));
	temp4 = temp*cos(2 * Pi*(b.ele(x, y-1) - n));
	if (c.ele(x + 1, y) == 0 && temp1<0 &&  rn()<(1-exp(2*Jcp*temp1)) )
	{
		c.ele(x + 1, y) = 1;
		q.push_back(Coordinates2D(x + 1, y));
	}
	if (c.ele(x - 1, y) == 0 && temp2<0 && rn()<(1 - exp(2 * Jcp*temp2)))
	{
		c.ele(x - 1, y) = 1;
		q.push_back(Coordinates2D(x - 1, y));
	}
	if (c.ele(x, y+1) == 0 && temp3<0 && rn()<(1 - exp(2 * Jcp*temp3)))
	{
		c.ele(x, y+1) = 1;
		q.push_back(Coordinates2D(x, y+1));
	}
	if (c.ele(x, y-1) == 0 && temp4<0 && rn()<(1 - exp(2 * Jcp*temp4)))
	{
		c.ele(x , y-1) = 1;
		q.push_back(Coordinates2D(x, y-1));
	}
}
void Wolff(Block& b)
{
	int x;int y;double n;
	double temp;
	static std::list<Coordinates2D> q;
	static Cluster c;
	static std::list<Coordinates2D>::iterator iter;
	x = rn_i() % Lx;
	y = rn_i() % Ly;// (x,y) is a random start-point;
	n = b.ele(x,y); // n is a random direction respective to (x,y)
	q.push_back(Coordinates2D(x, y));
	iter = q.begin();
	while (iter==q.end())
	{
		FindNeighbour(q,iter,c,b,n);
		iter++;
	}
	while (!q.empty())
	{
		c.ele(q.front().x, q.front().y) = 0;
		b.ele(q.front().x, q.front().y) = fmod(2 * n - b.ele(q.front().x, q.front().y) + 1, 1);
		q.pop_front();
	}
}
void Markov(Block& b)
{
	Metropolis(b);
	Wolff(b);
}