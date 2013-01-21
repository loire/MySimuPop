#ifndef DEPRESSION_H
#define DEPRESSION_H

#include "parameter.h"
#include <vector>
#include <iostream>
#include <boost/dynamic_bitset.hpp>
#include "MersenneTwister.h"
using namespace std;


// Variables globales:

#define fichierLecture "parametres.txt"     // nom des fichiers 
#define fichierEcriture "resultats.txt"		// entree et sortie


struct chr_diplo
{
	boost::dynamic_bitset<> chr1; // selected loci
	boost::dynamic_bitset<> chr2; // selected loci
};



struct selCoeffs
{
	double a1;		// Fitness of homozygous adaptative mutation in environment 1
	double ha1;		// Fitness of heterozygous adapatative mutation in environment 1
	double a2;		// Fitness of homozygous adaptative mutation in environment 2
	double ha2;		// Fitness of heterozygous adapatative mutation in environment 2
	int locus;		// Incompatible locus
	double c;		// Cost of incompatibility
};

// Prototypes de fonctions

void ouvrirFichierE();
void ouvrirFichierS();
void ecrireParametres(Parameter &param);
bool lireFichier(Parameter &param);
void recursion(Parameter &param);
double gammln(const double xx);
double poisdev(const double xm);
double binldev(const double pp, const int n);
double gasdev();
void rec(chr_diplo &res, chr_diplo &c1, chr_diplo &c2, double sz, int nS);
void freerec(chr_diplo &res, chr_diplo &c1,chr_diplo &c2, int Taille);
void cntl_c_handler(int bidon);

void fitness(Parameter &param, double ** &Wij,double * &wbar,double * &wmax,chr_diplo ** &pop,selCoeffs * &Sc);
void steppingStone1D(Parameter &param, int N_1, double ** &Wij,double * &wmax, chr_diplo ** &temp, chr_diplo ** &pop);
void twoDemeMigration(Parameter &param, int N_1,  double ** &Wij,  double * &wmax, chr_diplo ** &temp, chr_diplo ** &pop);

void initpop(Parameter &, chr_diplo ** &p, selCoeffs * &Sc);
void mutation_muck(Parameter &param, chr_diplo ** &pop);
void mutation_all_site(Parameter &param, chr_diplo ** &pop);
#endif
