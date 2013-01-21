#include "depression.h"
#include <vector>
#include <boost/dynamic_bitset.hpp>
#include <cmath>
#include <bitset>
#include <string>
using namespace std;

extern MTRand rnd;


// mutation all site: from site mutation rate Uv, the number of mutation #mut is picked in poisson distribution.
// #mut Random site along genomes is mutated (1 -> 0 or 0 -> 1).
// each allele is equally subject to mutation.



void mutation_all_site(Parameter &param, chr_diplo ** &pop)
{
	int i,j,k,mut;
	int nbS_1=param.Get_nbS()-1;


	for (i = 0; i < param.Get_n(); i++)
	{
		for (j=0; j< param.Get_N() ; j++)
		{
			mut = int(poisdev(param.Get_U()));
			for (k = 0; k < mut; k++)
			{
				int l = rnd.randInt(nbS_1);
				if (rnd.rand()>0.5)		// mutation fall on chr1
				{
					if (pop[i][j].chr1[l]==0)
						pop[i][j].chr1[l]=1;
					else
						pop[i][j].chr1[l]=0;
				}
				else					// mutation fall on chr2
				{
					if (pop[i][j].chr2[l]==0)
						pop[i][j].chr2[l]=1;
					else
						pop[i][j].chr2[l]=0;
				}
			}
		}
	}
}
void mutation_muck(Parameter &param, chr_diplo ** &pop)
{
}

void initpop(Parameter &param, chr_diplo ** &pop, selCoeffs * &Sc)
{
	int i,j,k;
	int nbS_1=param.Get_nbS()-1;
	// loop for first habitat / deme -> locally adapted allele have positive selective values (using Sc[i].a1)
	for (i = 0; i < param.Get_loc(); i++)
//	for (i = param.Get_loc(); i < param.Get_n() ;i++)
	{
		cout << "deme: " << i << " loc: " << param.Get_loc() << endl;

		for (j=0; j < param.Get_N() ; j++)
		{
			for (k = nbS_1; k >= 0; k--)
			{

				if (Sc[k].a1 > 1)
				{
					pop[i][j].chr1[k]=1;
					pop[i][j].chr2[k]=1;
				}
			}

		}
	}
}



