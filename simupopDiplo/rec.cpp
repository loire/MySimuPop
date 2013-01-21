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



void rec(chr_diplo &res, chr_diplo &c1, chr_diplo &c2, double Sz, int nS)
{

	cout << "change in rec function"
	cout << "And another change !!"
	cout << "And another !!"
	vector<int> pos;
	int j;
	int nS_1 = nS - 1;
	boost::dynamic_bitset<> rec;
	boost::dynamic_bitset<> off1;
	boost::dynamic_bitset<> off2;
	
	res.chr1.clear();
	res.chr2.clear();
	
	// First parent (c1)

	// number and positions of cross-overs
	int nbCo = int(poisdev(Sz));
	for (j = 0; j < nbCo; j++)
		pos.push_back(rnd.randInt(nS_1));
	sort(pos.begin(), pos.end());
	
	// recombination mask:

	for (j = 0; j < nbCo; j++)
		rec.resize(pos[j], (j % 2) == 0 ? 0 : 1);
	rec.resize(nS, (nbCo % 2) == 0 ? 0 : 1);
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

	nbCo = int(poisdev(Sz));
	for (j = 0; j < nbCo; j++)
		pos.push_back(rnd.randInt(nS_1));
	sort(pos.begin(), pos.end());

	// recombination mask:

	for (j = 0; j < nbCo; j++)
		rec.resize(pos[j], (j % 2) == 0 ? 0 : 1);
	rec.resize(nS, (nbCo % 2) == 0 ? 0 : 1);
	off1 = (c2.chr1 & rec);
	rec.flip();
	off2 = (c2.chr2 & rec);

	res.chr2 = (off1 | off2);

}


void freerec(chr_diplo &res, chr_diplo &c1, chr_diplo &c2, int Taille)
{
	

	boost::dynamic_bitset<> rec;
	boost::dynamic_bitset<> off1;
	boost::dynamic_bitset<> off2;
	
	res.chr1.clear();
	res.chr2.clear();
		// recombination mask:
		
	rec=RandomMask(Taille);
	off1 = (c1.chr1 & rec);
	rec.flip();
	off2 = (c1.chr2 & rec);
	
	res.chr1 = (off1 | off2);

	rec.clear();
	off1.clear();
	off2.clear();

	rec=RandomMask(Taille);
	off1 = (c2.chr1 & rec);
	rec.flip();
	off2 = (c2.chr2 & rec);

	res.chr2 = (off1 | off2);

}




