#include "depression.h"
#include <vector>
#include <boost/dynamic_bitset.hpp>
#include <cmath>
#include <bitset>
#include <string>
using namespace std;

extern MTRand rnd;

boost::dynamic_bitset<> RandomMask(int N)

{	
	int num=floor(N/32)+1;
	string s_masque = "";
	bitset<32> tt;
	for (int i =0 ; i< num ; i++)
		{
	bitset<32> tt((unsigned long int ) rnd.randInt());
		s_masque+=tt.to_string();
	}
	boost::dynamic_bitset<> maskfromstring(s_masque);	
	maskfromstring.resize(N);
	return maskfromstring;

}

////////////////////////////////////////////////////////////////////
// rec_r choose position of crossing over according to r, recombination rate.
// It needs a random number for each position in decide for recombination events, thus it's more appropriate for few loci
// If more loci are needed, the program will use rec_L function
////////////////////////////////////////////////////////////////////

void rec_r(Parameter &param, chr_diplo &res, chr_diplo &c1, chr_diplo &c2)
{
	int nbCo;
	vector<int> pos;
	int j;
	int nS_1= param.Get_nbS() -1; // Last locus position (neutral one)
	boost::dynamic_bitset<> rec;
	boost::dynamic_bitset<> off1;
	boost::dynamic_bitset<> off2;

	res.chr1.clear();
	res.chr2.clear();

	// First parent (c1)

	// number and positions of cross-overs
	//int nbCo = int(poisdev(param.Get_L()));
	//for (j = 0; j < nbCo; j++)
	//	pos.push_back(rnd.randInt(nS_2));
	//sort(pos.begin(), pos.end());
	nbCo=0;
	for (j=0; j< nS_1; j++)
	{
		if ( rnd.rand() < param.Get_r() )
		{
			pos.push_back(j);
			nbCo++;
		}
	}
	if (rnd.rand()>param.Get_rn())
	{
		pos.push_back(nS_1);
		nbCo++;
	}
	// recombination mask:

	for (j = 0; j < nbCo; j++)
		rec.resize(pos[j], (j % 2) == 0 ? 0 : 1);
	rec.resize(param.Get_nbS(), (nbCo % 2) == 0 ? 0 : 1);
	off1 = (c1.chr1 & rec);
	rec.flip();
	off2 = (c1.chr2 & rec);
	res.chr1 = (off1 | off2);

	// Second Parent (c2)

	// number and positions of cross-overs
	pos.clear();
	rec.clear();
	off1.clear();
	off2.clear();
	nbCo=0;
	for (j=0; j< nS_1; j++)
	{
		if ( rnd.rand() < param.Get_r() )
		{
			pos.push_back(j);
			nbCo++;
		}
	}
	if (rnd.rand()>param.Get_rn())
	{
		pos.push_back(nS_1);
		nbCo++;
	}
	// recombination mask:
	for (j = 0; j < nbCo; j++)
		rec.resize(pos[j], (j % 2) == 0 ? 0 : 1);
	rec.resize(param.Get_nbS(), (nbCo % 2) == 0 ? 0 : 1);
	off1 = (c2.chr1 & rec);
	rec.flip();
	off2 = (c2.chr2 & rec);

	res.chr2 = (off1 | off2);
}
/////////////////////////////////////////////////////////////////////////////////////
// Rec_L use parameter L, ie genetic length of the chromosome.
// It takes a random number of position for crossing over in the chromosome according to its length
// Suitable for many loci
/////////////////////////////////////////////////////////////////////////////////////

void rec_L(Parameter &param, chr_diplo &res, chr_diplo &c1, chr_diplo &c2)
{

	vector<int> pos;
	int j;
	int nS_1= param.Get_nbS() -1; // Last locus position (neutral one)
	int nS_2 = nS_1 - 1;		  // The one before
	boost::dynamic_bitset<> rec;
	boost::dynamic_bitset<> off1;
	boost::dynamic_bitset<> off2;
	
	res.chr1.clear();
	res.chr2.clear();
	
	// First parent (c1)

	// number and positions of cross-overs
	int nbCo = int(poisdev(param.Get_L()));
	for (j = 0; j < nbCo; j++)
		pos.push_back(rnd.randInt(nS_2));
	sort(pos.begin(), pos.end());
	
	if (rnd.rand()>param.Get_rn())
		pos.push_back(nS_2);

	// recombination mask:

	for (j = 0; j < nbCo; j++)
		rec.resize(pos[j], (j % 2) == 0 ? 0 : 1);
	rec.resize(param.Get_nbS(), (nbCo % 2) == 0 ? 0 : 1);
	off1 = (c1.chr1 & rec);
	rec.flip();
	off2 = (c1.chr2 & rec);
	res.chr1 = (off1 | off2);

	// Second Parent (c2)

	// number and positions of cross-overs
	pos.clear();
	rec.clear();
	off1.clear();
	off2.clear();

	nbCo = int(poisdev(param.Get_L()));
	for (j = 0; j < nbCo; j++)
		pos.push_back(rnd.randInt(nS_2));
	sort(pos.begin(), pos.end());
	if (rnd.rand()>param.Get_rn())
		pos.push_back(nS_2);

	// recombination mask:

	for (j = 0; j < nbCo; j++)
		rec.resize(pos[j], (j % 2) == 0 ? 0 : 1);
	rec.resize(param.Get_nbS(), (nbCo % 2) == 0 ? 0 : 1);
	off1 = (c2.chr1 & rec);
	rec.flip();
	off2 = (c2.chr2 & rec);

	res.chr2 = (off1 | off2);

}


void freerec(Parameter &param, chr_diplo &res, chr_diplo &c1, chr_diplo &c2)
{
	

	boost::dynamic_bitset<> rec;
	boost::dynamic_bitset<> off1;
	boost::dynamic_bitset<> off2;
	
	res.chr1.clear();
	res.chr2.clear();
		// recombination mask:
		
	rec=RandomMask(param.Get_nbS());
	off1 = (c1.chr1 & rec);
	rec.flip();
	off2 = (c1.chr2 & rec);
	
	res.chr1 = (off1 | off2);

	rec.clear();
	off1.clear();
	off2.clear();

	rec=RandomMask(param.Get_nbS());
	off1 = (c2.chr1 & rec);
	rec.flip();
	off2 = (c2.chr2 & rec);

	res.chr2 = (off1 | off2);

}




