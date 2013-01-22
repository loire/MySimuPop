/*
 * migration.cpp
 *
 *  Created on: 7 déc. 2012
 *      Author: etienne
 */
#include "depression.h"
#include <vector>
#include <boost/dynamic_bitset.hpp>
#include <cmath>
#include <bitset>
#include <string>
using namespace std;

extern MTRand rnd;
extern void (*recfunc)(Parameter&, chr_diplo &, chr_diplo &, chr_diplo &);


void twoDemeMigration(Parameter &param, int N_1,  double ** &Wij,  double * &wmax, chr_diplo ** &temp, chr_diplo ** &pop)
{
	int i,j,ind;
	int par1,par2;
	// Deme 0:
//	recfunc(param.Get_N());
	for (ind = 0;  ind < param.Get_N() ; ind++)
	{
		// Phillopatric

		if (rnd.rand()>param.Get_m())
		{
			do
			{
				par1 = rnd.randInt(N_1);
			}
			while (rnd.rand() > Wij[0][par1] / wmax[0]);

			do
			{
				par2 = rnd.randInt(N_1);
			}
			while (rnd.rand() > Wij[0][par2] / wmax[0]);
			// recombination:
			if (param.Get_freerec()!=1)
				rec_L(param, temp[0][ind], pop[0][par1], pop[0][par2]);
			else
				freerec(param, temp[0][ind], pop[0][par1], pop[0][par2]);
		}
		else  // migrants:
		{
			do
			{
				par1 = rnd.randInt(N_1);
			}
			while (rnd.rand() > Wij[1][par1] / wmax[1]);

			do
			{
				par2 = rnd.randInt(N_1);
			}
			while (rnd.rand() > Wij[1][par2] / wmax[1]);
			// recombination:
			if (param.Get_freerec()!=1)
				rec_L(param, temp[0][ind], pop[1][par1], pop[1][par2]);
			else
				freerec(param, temp[0][ind], pop[1][par1], pop[1][par2]);
		}
	}
	// Deme 1:

	for (ind = 0;  ind < param.Get_N() ; ind++)
	{

		// Phillopatric:
		if (rnd.rand()>param.Get_m())
		{
			do
			{
				par1 = rnd.randInt(N_1);
			}
			while (rnd.rand() > Wij[1][par1] / wmax[1]);

			do
			{
				par2 = rnd.randInt(N_1);
			}
			while (rnd.rand() > Wij[1][par2] / wmax[1]);
			// recombination:
			if (param.Get_freerec()!=1)
				rec_L(param, temp[1][ind], pop[1][par1], pop[1][par2]);
			else
				freerec(param, temp[1][ind], pop[1][par1], pop[1][par2]);
		}
		else   // migrants:
		{
			do
			{
				par1 = rnd.randInt(N_1);
			}
			while (rnd.rand() > Wij[0][par1] / wmax[0]);

			do
			{
				par2 = rnd.randInt(N_1);
			}
			while (rnd.rand() > Wij[0][par2] / wmax[0]);
			// recombination:
			if (param.Get_freerec()!=1)
				rec_L(param, temp[1][ind], pop[0][par1], pop[0][par2]);
			else
				freerec(param, temp[1][ind], pop[0][par1], pop[0][par2]);
		}
	}
	for (i=0; i < param.Get_n() ; i++)
	{
		for (j=0; j< param.Get_N() ; j++)
		{
			pop[i][j]=temp[i][j];
		}
	}

}



void steppingStone1D(Parameter &param, int N_1, double ** &Wij,double * &wmax, chr_diplo ** &temp, chr_diplo ** &pop)
{
	int i,j,ind;
	int par1,par2;
	double migprob;

	// First deme

	for (ind = 0;  ind < param.Get_N() ; ind++)
	{
		if (rnd.rand()>param.Get_m())
		{
			do
			{
				par1 = rnd.randInt(N_1);
			}
			while (rnd.rand() > Wij[0][par1] / wmax[0]);

			do
			{
				par2 = rnd.randInt(N_1);
			}
			while (rnd.rand() > Wij[0][par2] / wmax[0]);
			// recombination:
			if (param.Get_freerec()!=1)
				rec_L(param, temp[0][ind], pop[0][par1], pop[0][par2]);
			else
				freerec(param, temp[0][ind], pop[0][par1], pop[0][par2]);
		}
		else
		{
			do
			{
				par1 = rnd.randInt(N_1);
			}
			while (rnd.rand() > Wij[1][par1] / wmax[1]);

			do
			{
				par2 = rnd.randInt(N_1);
			}
			while (rnd.rand() > Wij[1][par2] / wmax[1]);
			// recombination:
			if (param.Get_freerec()!=1)
				rec_L(param, temp[0][ind], pop[1][par1], pop[1][par2]);
			else
				freerec(param, temp[0][ind], pop[1][par1], pop[1][par2]);
		}

	}
	//for deme 1 to nv-2
	double mig_half= param.Get_m()/2;
	double mig_half_1= 1 - mig_half;

	for (i=1; i < param.Get_n()-1 ; i++)
	{
		int i_l = i-1; // deme of left
		int i_r = i+1;// deme of right

		for (ind = 0;  ind < param.Get_N() ; ind++)
		{
			migprob=rnd.rand();
			if (migprob < mig_half )			// Indiv come from left
			{
				do
				{
					par1 = rnd.randInt(N_1);
				}
				while (rnd.rand() > Wij[i_l][par1] / wmax[i_l]);

				do
				{
					par2 = rnd.randInt(N_1);
				}
				while (rnd.rand() > Wij[i_l][par2] / wmax[i_l]);
				// recombination:
				if (param.Get_freerec()!=1)
					rec_L(param, temp[i][ind], pop[i_l][par1], pop[i_l][par2]);
				else
					freerec(param, temp[i][ind], pop[i_l][par1], pop[i_l][par2]);
			}
			else
			{
				if (migprob > mig_half_1)			// Indiv come from right
				{
					do
					{
						par1 = rnd.randInt(N_1);
					}
					while (rnd.rand() > Wij[i_r][par1] / wmax[i_r]);

					do
					{
						par2 = rnd.randInt(N_1);
					}
					while (rnd.rand() > Wij[i_r][par2] / wmax[i_r]);
					// recombination:
					if (param.Get_freerec()!=1)
						rec_L(param, temp[i][ind], pop[i_r][par1], pop[i_r][par2]);
					else
						freerec(param, temp[i][ind], pop[i_r][par1], pop[i_r][par2]);
				}
				else				// Phillopatric indiv
				{
					do
					{
						par1 = rnd.randInt(N_1);
					}
					while (rnd.rand() > Wij[i][par1] / wmax[i]);

					do
					{
						par2 = rnd.randInt(N_1);
					}
					while (rnd.rand() > Wij[i][par2] / wmax[i]);
					// recombination:
					if (param.Get_freerec()!=1)
						rec_L(param, temp[i][ind], pop[i][par1], pop[i][par2]);
					else
						freerec(param, temp[i][ind], pop[i][par1], pop[i][par2]);
				}
			}
		}
	}
	// And the last deme (indiv from the left with a prob of mv)
	int nv_1 = param.Get_n()-1;
	int nv_2 = param.Get_n()-2;

	for (ind = 0;  ind < param.Get_N() ; ind++)
	{
		if (rnd.rand()>param.Get_m())
		{
			do
			{
				par1 = rnd.randInt(N_1);
			}
			while (rnd.rand() > Wij[nv_1][par1] / wmax[nv_1]);

			do
			{
				par2 = rnd.randInt(N_1);
			}
			while (rnd.rand() > Wij[nv_1][par2] / wmax[nv_1]);
			// recombination:
			if (param.Get_freerec()!=1)
				rec_L(param, temp[nv_1][ind], pop[nv_1][par1], pop[nv_1][par2]);
			else
				freerec(param, temp[nv_1][ind], pop[nv_1][par1], pop[nv_1][par2]);
		}
		else
		{
			do
			{
				par1 = rnd.randInt(N_1);
			}
			while (rnd.rand() > Wij[nv_2][par1] / wmax[nv_2]);

			do
			{
				par2 = rnd.randInt(N_1);
			}
			while (rnd.rand() > Wij[nv_2][par2] / wmax[nv_2]);
			// recombination:
			if (param.Get_freerec()!=1)
				rec_L(param, temp[nv_1][ind], pop[nv_2][par1], pop[nv_2][par2]);
			else
				freerec(param, temp[nv_1][ind], pop[nv_2][par1], pop[nv_2][par2]);
		}

	}
	// old pop is remplaced
	for (i=0; i < param.Get_n() ; i++)
	{
		for (j=0; j< param.Get_N() ; j++)
		{
			pop[i][j]=temp[i][j];
		}
	}
}
