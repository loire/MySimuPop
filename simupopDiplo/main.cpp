// Fonction main(): demande a l'utilisateur s'il y a un fichier entree 
// contenant les valeurs des parametres, auquel cas le fichier est lu,
// sinon l'utilisateur doit entrer les differentes valeurs. Simulation,
// puis ecriture des resultats dans un fichier sortie.

#include "depression.h"
#include "parameter.h"
#include "MersenneTwister.h"
#include <iostream>
#include <fstream>
#include <boost/dynamic_bitset.hpp>

MTRand rnd;
using namespace std;

ofstream fichierS;
Parameter param;

int main()
{	
	cout << "Simulation start..." << endl;
	cout << "This is from recombination branch of the project" << endl;
	// Ouverture des fichiers entree et sortie:
	bool fin;

	fichierS.open(fichierEcriture,ofstream::app);
	
	fin = false;
	fin=lireFichier(param);
	if (fin)
	{
			ecrireParametres(param);
			recursion(param);
	}
	else
			cout << "Problem opening file" << endl;

//	int no = 1;
//	do
//	{
//		fin = lireFichier(n, Nt, m, a, s, h, sig_s, c, sig_c, U, nbS, L, loc, NbGen, adapgen, pas, freerec);
//
//		if (!fin)
//		{
//			ecrireParametres(n, Nt, m, a, s,h, sig_s, c, sig_c, U, nbS, L, loc, NbGen, adapgen, pas, freerec);
//						 // ecrit valeurs des parametres
//						 // dans fichier sortie
//			recursion(n, Nt, m, a, s, h, sig_s, c, sig_c, U, nbS, L, loc, NbGen, adapgen, pas, freerec);
//					// simulation, ecrit resultats dans fichier sortie
//			no++;
//		}
//	} while (!fin);
//
	// Fermeture des fichiers:
	fichierS.close();
	return 0 ;
}
