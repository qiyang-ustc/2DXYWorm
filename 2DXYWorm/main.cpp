#include<iostream>
#include<fstream>
#include"my_rng.h"
#include"my_vrbls.h"
#include"my_statistics.h"
#include"Markov.h"
#include<vector>
#include<stdlib.h>
void Initialize()
{
	Vol = Lx*Ly;
	Vol2 = pow(Vol,2);
	Vol4 = pow(Vol2,2);
	Block_Obs.resize(NBlock);
	Sample_Obs.resize(NSample);
	TotalSample = NBlock*NSample;
}
time_t tm;
void RandomWalk(Block& b)//Walk "steps" steps.....I hope you know
{
	static FlowWeight fw;
	int x = 0;
	int y = 0;
	b.fresh();
	int temp;
	double tempr;
	for (int i = 0;i < steps;i++)
	{
		tempr = rn();
		temp = 1 * (tempr < 0.25) + 2 * (tempr > 0.25&&tempr < 0.5) + 3 * (tempr > 0.5&&tempr < 0.75) + 4 * (tempr > 0.75);
		switch (temp)
		{
		case 1:
			if (rn() < (fw.weight(1+b.down_flow(x - 1, y)) / fw.weight(b.down_flow(x - 1, y))))
			{
				x = x - 1;
				b.down_flow(x, y)++;
			}
			b.ele(x, y)++;
			break;
		case 2:
			if (rn() < (fw.weight(1 + b.left_flow(x, y+1)) / fw.weight(b.left_flow(x, y+1))))
			{
				y = y + 1;
				b.left_flow(x, y)++;
			}
			b.ele(x, y)++;
			break;
		case 3:
			if (rn() < (fw.weight(-1 +b.down_flow(x, y)) / fw.weight(b.down_flow(x, y))))
			{
				b.down_flow(x, y)--;
				x = x + 1;
			}
			b.ele(x, y)++;
			break;
		case 4:
			if (rn() < (fw.weight(-1+b.left_flow(x, y)) / fw.weight(b.left_flow(x, y))))
			{
				b.left_flow(x, y)--;
				y = y -1;
			} 
			b.ele(x, y)++;
			break;
		default:
			break;
		}
		x = (x + Lx) % Lx;
		y = (y + Ly) % Ly;
	}
}
int main(int argc,char *argv[])
{
	set_elapse_time();
	/*
	Lx = atoi(argv[1]);
	Ly = atoi(argv[2]);
	Jcp=atof(argv[3]);
	NBlock=atoi(argv[4]);
	NSample=atoi(argv[5]);
	std::ofstream fout(OutFile, std::ios::app);
	*/
	//We use System parameters to get parameters;


	std::ifstream fin(InFile, std::ios::in);
	if (fin.is_open())
	{
		Get_Parameters(fin);
	}
	else
	{
		Get_Parameters(std::cin);
	}
	Initialize();
	Block b;

	//!!!!!*********NSample must be 1, in odrer to inherit all function in following*********
	for (int i = 0;i < NBlock;i++)
	{
		for (int j = 0;j < NSample;j++)
		{
			RandomWalk(b);
			b.Calculate_Quan();
			Collect_data(b,j);
		}
		Normalize_data(i);
	}
	//-------Analysis----------
	Analyze_data();
	//-------write to file-----
	write2file(std::cout);
	elapse_time();
	system("pause");
}
