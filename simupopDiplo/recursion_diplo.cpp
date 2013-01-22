#include "depression.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <ctime>
#include <boost/dynamic_bitset.hpp>
#include <cmath>
#include <csignal>
using namespace std;

extern MTRand rnd;
extern ofstream fichierS;
bool cntl_c_bool = false;

void cntl_c_handler (int bidon) 
{
	cntl_c_bool = true;
}

void recursion(Parameter &param)
{


	// Useful stuff
	int i, j, k, l, gen, mut, chr1, chr2, ind;
	double w, x, inc;

	//double index_h, w_bet, w_int1, w_int2, w_bet_mean, w_int1_mean, w_int2_mean;	// to compute hybrid index
	int h,i1,i2,j1,j2;

	chr_diplo progeny;

	///////////////////////////////////////////////////////////////////////////
	// Output file name //
	//////////////////////////////////////////////////////////////////////////
	char nomFichier[256];
	stringstream nomF;
	nomF << "result_n" << param.Get_n() << "_N" << param.Get_N() << "_sA" << param.Get_a() << "_s" << param.Get_s() << "_sigs" << param.Get_sig_s() << "_c" << param.Get_i() << "_sigc" << param.Get_sig_i() << "_U" << param.Get_U() << "_nbS" << param.Get_nbS() << "_L" << param.Get_L() << "_loc" << param.Get_loc() << "_freerec" << param.Get_freerec() << ".txt";
	nomF >> nomFichier;

	ofstream fichierR;
	fichierR.open(nomFichier,ofstream::app);
	fichierR << "#New Sim\n";


	////////////////////////////////////////////////////////////////////////
	// Time of execution handling //
	////////////////////////////////////////////////////////////////////////
	time_t debut, fin;
	struct tm *ptr;
	debut = time(0);

	cntl_c_bool = false;
	signal(SIGINT, cntl_c_handler);


	///////////////////////////////////////////////////////////////////////
	// Allocate memory for individuals
	//////////////////////////////////////////////////////////////////////

	int nbFix = 0;
	int pas_c = 20;
	boost::dynamic_bitset<> mask;	
	boost::dynamic_bitset<> masktemp;

	double nN = double(param.Get_n() * param.Get_N());
	double not_mig  = 1 - param.Get_m(); // Probability of no migration
	int N_1 = param.Get_N() - 1;
	int nbS_1 = param.Get_nbS() - 1;

	// Bidimensional array : 	
	chr_diplo ** pop = new chr_diplo *[param.Get_n()];
	chr_diplo ** temp = new chr_diplo *[param.Get_n()];

	double ** Wij = new double *[param.Get_n()];
	for ( i=0; i < param.Get_n();i++)
	{
		pop[i] = new chr_diplo[param.Get_N()]; // ALlocate memory to store chr_diplo in each deme
		temp[i] = new chr_diplo[param.Get_N()]; // copy of the other one
		Wij[i] = new double[param.Get_N()]; // Allocate memory to store fitnesses in each deme
	}		

	double * wbar = new double [param.Get_n()];
	double * wmax = new double [param.Get_n()];
	double * le = new double [param.Get_nbS()]; // to store global allele frequency



	///////////////////////////////////////////////////////////////////////////
	// Allocate memory for fitness of locus
	///////////////////////////////////////////////////////////////////////////
	// TODO : Put in function
	selCoeffs * Sco = new selCoeffs [param.Get_nbS()];


	//sampling selection coefficients:
	vector<int> sites;
	for (i = 0; i < param.Get_nbS(); i++)
		sites.push_back(i);	

	int count_a = param.Get_a();
	double tmp_s;
	for (i = 0; i < param.Get_nbS(); i++)
	{

		if(count_a>0)
		{
			do
			{
				tmp_s = param.Get_s() + param.Get_sig_s() * gasdev();

				Sco[sites[i]].a1 = 1 + tmp_s;
				Sco[sites[i]].ha1= 1 + param.Get_h() * tmp_s;
			}
			while (Sco[sites[i]].a1 < 1);
			Sco[sites[i]].a2=2-Sco[sites[i]].a1;
			Sco[sites[i]].ha2=2-Sco[sites[i]].ha1;

			Sco[sites[i]].c=0;
			Sco[sites[i]].locus=0;
		}
		else // neutral mutation
		{
			Sco[sites[i]].a1 = 1;
			Sco[sites[i]].a2 = 1;
			Sco[sites[i]].ha1 = 1;
			Sco[sites[i]].ha2 = 1;
			Sco[sites[i]].c=0;	
			Sco[sites[i]].locus=0;
		}
		count_a--;
	}
	for (i = 0; i < param.Get_n(); i++) // for each deme
	{
		for (j=0; j < param.Get_N() ; j++)	// for each chrom
		{
			pop[i][j].chr1.resize(param.Get_nbS());			// create chrom (dynamic biteset) with nbSv sites set to 0
			pop[i][j].chr2.resize(param.Get_nbS());			// create chrom (dynamic biteset) with nbSv sites set to 0
		}
	}

	//////////////////////////////////////////////////////////////////////////////
	// Print parameters values on standard output
	/////////////////////////////////////////////////////////////////////////////


	cout << " nv "<< param.Get_n() << "\n";
	cout << " Nv "<< param.Get_N()<< "\n";
	cout << " mv "<<param.Get_m() << "\n";
	cout << " av "<< param.Get_a()<< "\n";
	cout << " sv "<<param.Get_s() << "\n";
	cout << " hv "<< param.Get_h()<< "\n";
	cout << " sig_sv "<<param.Get_sig_s() << "\n";
	cout << " iv "<<  param.Get_i()<< "\n";
	cout << " sig_iv "<< param.Get_sig_i() << "\n";
	cout << " Uv "<<param.Get_U() << "\n";
	cout << " nbSv "<< param.Get_nbS()<< "\n";
	cout << " Lv "<<param.Get_L() << "\n";
	cout << " locv "<< param.Get_loc()<< "\n";
	cout << " NbGenv "<< param.Get_NbGen()<< "\n";
	cout << " adapgenv "<<param.Get_adapgen() << "\n";
	cout << " pasv "<<param.Get_pas() << "\n";
	cout << " freerecv "<< param.Get_freerec() << "\n";
	cout << " no_mutv " << param.Get_no_mut() << "\n";
	cout <<  " nbSv  "<<param.Get_nbS() << "\n";

	////// Some printing (to be erased later)
	cout << "mutations\n";	
	for (i=0 ; i < param.Get_nbS() ; i ++)
	{
		cout << i << " " << Sco[i].a1 << " " << Sco[i].ha1 << " "<< Sco[i].a2 << " " << Sco[i].ha2 << " " << Sco[i].locus << " " << Sco[i].c << "\n";
	}

	////// 


	//initialization: no mutation


	//	progeny.chr1.resize(param.Get_nbS());
	//	progeny.chr2.resize(param.Get_nbS());

	////////////////////////////////////////////////////////////////////////////////
	// Addressing life cycle functions according to parameter
	////////////////////////////////////////////////////////////////////////////////
	void (*recfunc)(Parameter&, chr_diplo &, chr_diplo &, chr_diplo &) = NULL;

	if (param.Get_freerec()==1)
		recfunc= &freerec;
	else
	{
		if ( param.Get_nbS() < param.Get_L() )
			recfunc = &rec_r;
		else
			recfunc= &rec_L;
	}

	void (*migfunc)(Parameter&, int, double**&, double*&, chr_diplo**&, chr_diplo**&) = NULL;
	if (param.Get_n() > 2)
		migfunc = &steppingStone1D;
	else
		migfunc = &twoDemeMigration;



	void (*mutfunc)(Parameter &, chr_diplo ** &p) = NULL;
	if (param.Get_no_mut())
	{
		mutfunc = &mutation_muck;
		cout << "no mutations, pop will be initialized according to habitats (adaptative mutation fixed at the beginning of the simul\n";
		initpop(param, pop, Sco);
	}
	else
	{
		mutfunc = &mutation_all_site;
		cout << "locus will experience mutation at a rate " << param.Get_U() << " at each generations\n";
	}
	cout << param.Get_nbS() << endl;
	cout << param.Get_L() << endl;
	cout << param.Get_r() << endl;
	cout << param.Get_rn() << endl;
	cout << param.Get_freerec() << endl;



	////////////////////////////////////////////////////////
	///////////////// THE LOOP ////////////////////////////
	//////////////////////////////////////////////////////
	//generations
	//mask.resize(param.Get_nbS(),1); // fixed mask of size nbSv filled with 1
	for (gen = 0; gen < param.Get_NbGen(); gen++)				// for each generation
	{ 
		//// Affichage de la pop
		//
		//
		// 		for (i = 0; i < param.Get_n(); i++)
		// 			cout << "pop numero" << i << "\t" ;
		// 		cout << "\n";
		// 		for (j=0; j< param.Get_N() ; j++)
		// 			{
		// 					for (i = 0 ; i < param.Get_n() ; i++)
		// 						cout << pop[i][j].sel << "\t";
		// 					cout << "\n";
		// 			}

		//////

		// comptage des mutations:
		////// ici a modifier pour tester l'occurence d'une meme mutation chez tout les individus !
//		if (gen % pas_c == 0)			// at every pas_c, check for fixed mutation. If fixed, set to 0. Could be done more nicely using bitset operators.
//		{
//			masktemp=mask;
//
//			for (i = 0; i < param.Get_n() ; i++)
//			{
//
//				for (int j = 0 ; j < param.Get_N() ; j++)
//				{
//					masktemp = (masktemp & pop[i][j].sel);
//				}
//			}
//			nbFix+=masktemp.count();
//
//			for (i = 0; i < param.Get_n() ; i++)
//			{
//				for (int j = 0 ; j < param.Get_N() ; j++)
//				{
//					pop[i][j].sel = (pop[i][j].sel ^ masktemp);
//				}
//
//			}
//		}

		//mutation:
		(*mutfunc)(param,pop);
		//fitnesses:
		fitness(param, Wij,wbar,wmax,pop,Sco);  // Calculate fitness of individuals, mean and max fitness of each populations
		// Migration
		(*migfunc)(param, N_1, Wij, wmax, temp, pop);


		//		cout << "coucou fin migration\n";
		//		cout << "pasv :" << param.Get_pas() << "\n";
		//		cout << "gen: " << gen << "\n";
		//		cout << "gen % pasv: " << gen % param.Get_pas() << "\n";
		if (gen % param.Get_pas() == 0)
		{
			cout << gen << "\t";
			for (i = 0; i < param.Get_n() ; i ++)
			{
				for (j = 0 ; j < param.Get_nbS() ; j++)
				{
					double freq_allele=0.0;
					for (k=0 ; k < param.Get_N() ; k++)
					{
						freq_allele+=(int) pop[i][k].chr1[j];
						freq_allele+=(int) pop[i][k].chr2[j];
					}

					freq_allele/=(float) 2*param.Get_N();
					cout << freq_allele << "\t";
					fichierR << freq_allele << "\t";

				}
				cout << "|\t";
			}
			cout << "\n";
			fichierR << "\n";
		}
		//~ fprintf(fichierR,"%f\t",index_h);

		// to print global allele frequency

		// for (i=0;  i<param.Get_nbS()+1 ; i++)
		//  				le[i]=0;
		//
		//  			for (i = 0; i < param.Get_n(); i++)
		//  			{
		//  				for (j = 0 ; j < param.Get_N() ; j++)
		//  				{
		//  					for (k=0; k <param.Get_nbS()+1; k++)
		//  					{
		//  							le[k]+=pop[i][j].sel[k];
		//  					}
		// 				}
		//  			}
		//
		//  			//cout << le << "\t" << param.Get_N() <<"\t" << param.Get_n() << "\t"<< le/ (param.Get_N()*param.Get_n())<< "\n";
		//  			for (i=0;  i < param.Get_nbS()+1 ; i++)
		//  				cout << (le[i] / nN) << "\t";
		//  			cout << nbFix << endl;
		//		}

		if (cntl_c_bool)
			break;
	}

	fichierR << "\n";

	fichierR.close();
	fin = time(0);  

	if (cntl_c_bool)
		fichierS << "\n\nInterrupted by user at generation " << gen << endl;



	// Ecriture des resultats:
	fichierS << "\n\nResultats dans fichier ";
	fichierS << nomFichier << endl;


	// Duree de la simulation:
	int temps = int(difftime(fin, debut));
	fichierS << "\n" << param.Get_NbGen() << " generations ont pris "<< temps / 3600 << " heure(s) " << (temps % 3600) / 60 << " minute(s) " << (temps % 60) << " secondes\n";

	// Date et heure:
	ptr=localtime(&fin);
	fichierS <<  asctime(ptr);
	for (i=0;i< param.Get_n() ; i++)
	{
		delete [] pop[i];
		delete [] temp[i];
		delete [] Wij[i];
	}
	delete [] pop;
	delete [] temp;
	delete [] Wij;		
	delete [] Sco;
	delete [] wbar;
	delete [] wmax;
}
