/*
 * fitness.cpp
 *
 *  Created on: 7 déc. 2012
 *      Author: etienne
 */

#include "depression.h"
#include <vector>

#include <cmath>
#include <bitset>
#include <string>
using namespace std;

extern MTRand rnd;


void fitness(Parameter &param, double ** &Wij,double * &wbar,double * &wmax,chr_diplo ** &pop,selCoeffs * &Sc)
{
	int i,j,k;
	double w;
	int nbS_1=param.Get_nbS()-1;
	// loop for first habitat / deme -> locally adapted allele have positive selective values (using Sc[i].a1)
	for (i = 0; i < param.Get_loc(); i++)
	{
//		cout << "deme: " << i << " loc: " << param.Get_loc() << endl;

		wbar[i] = 0;
		wmax[i] = 0;
		for (j=0; j < param.Get_N() ; j++)
		{

			w = 1.0;
			for (k = nbS_1; k >= 0; k--)
			{
//				cout << "Locus:"  <<  k << " Fitness: " << Sc[k].a1 << endl;
				if (pop[i][j].chr1[k] != 0 or pop[i][j].chr2[k] != 0)		// if something to do
				{
					if (pop[i][j].chr1[k] == 1 and pop[i][j].chr2[k] == 1) // fitness of homozygotes
					{
						w *= Sc[k].a1;
					}
					else if (pop[i][j].chr1[k] == 1 or pop[i][j].chr2[k] == 1) // fitness of heterozygotes
					{
						w*= Sc[k].ha1;
					}
				}
			}
			Wij[i][j] = w;
			wbar[i] += w;
			if (wmax[i] < w)
				wmax[i] = w;
		}
		wbar[i]= wbar[i] / param.Get_N();
	}
	// Loop for habitat where locally adapted alleles have negative selected values (using Sc[i].a2)
	for (i = param.Get_loc(); i < param.Get_n(); i++)
	{
//		cout << "deme: " << i << " loc: " << param.Get_loc() << endl;

		wbar[i] = 0;
		wmax[i] = 0;
		for (j=0; j < param.Get_N() ; j++)
		{

			w = 1.0;
			for (k = nbS_1; k >= 0; k--)
			{
//				cout << "Locus:"  <<  k << " Fitness: " << Sc[k].a2 << endl;
				if (pop[i][j].chr1[k] != 0 or pop[i][j].chr2[k] != 0)		// if something to do
				{
					if (pop[i][j].chr1[k] == 1 and pop[i][j].chr2[k] == 1) // fitness of homozygotes
					{
						w *= Sc[k].a2;
					}
					else if (pop[i][j].chr1[k] == 1 or pop[i][j].chr2[k] == 1) // fitness of heterozygotes
					{
						w*= Sc[k].ha2;
					}
				}
			}
			Wij[i][j] = w;
			wbar[i] += w;
			if (wmax[i] < w)
				wmax[i] = w;
		}
		wbar[i]= wbar[i] / param.Get_N();
	}

}
