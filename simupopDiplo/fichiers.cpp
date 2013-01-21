#include "depression.h"
#include "parameter.h"
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <fstream>
#include <boost/tokenizer.hpp>
#include <libconfig.h++>
#include <string>
using namespace std;
using namespace boost;
using namespace libconfig;

extern ofstream fichierS;
//extern Parameter param;


void ouvrirFichierS()   
{
	fichierS.open(fichierEcriture,ifstream::app);
}

bool lireFichier(Parameter &param)
					 // lit les valeurs du fichier de config et les stocke dans l'objet param
{
		libconfig::Config cfg;		// cfg object
		// some test on config file
		try
		{
			cfg.readFile(fichierLecture);
		}
		catch(const FileIOException &fioex)
		{
			cout << "I/O error while reading file." << endl;
			return(0);
		}
		catch(const ParseException &pex)
		{
		std::cout << "Parse error at " << pex.getFile() << ":" << pex.getLine() << " - " << pex.getError() << std::endl;
		return(0);
		}

		// Trying to parse paramaters

		try
		{
		param.Set_n(cfg.lookup("metapop.n"));
        param.Set_N(cfg.lookup("metapop.N"));
        param.Set_loc(cfg.lookup("metapop.loc"));
        param.Set_m(cfg.lookup("metapop.m"));
		}
		catch(const SettingNotFoundException &nfex)
		{
			cerr << "No metapop setting in configuration file." << endl;
			return(0);
		}
		try{
		param.Set_a(cfg.lookup("chrom.a"));
        param.Set_s(cfg.lookup("chrom.s"));
        param.Set_h(cfg.lookup("chrom.h"));
        param.Set_sig_s(cfg.lookup("chrom.sig_s"));
		param.Set_i(cfg.lookup("chrom.i"));
        param.Set_sig_i(cfg.lookup("chrom.sig_i"));
        param.Set_U(cfg.lookup("chrom.U"));
        param.Set_nbS(cfg.lookup("chrom.nbS"));
        param.Set_L(cfg.lookup("chrom.L"));
		param.Set_freerec(cfg.lookup("chrom.freerec"));
		param.Set_no_mut(cfg.lookup("chrom.no_mut"));
		}
		catch(const SettingNotFoundException &nfex)
		{
			cout << "No chrom setting in configuration file." << endl;
			return(0);
		}
		try
		{
		param.Set_NbGen(cfg.lookup("simul.NbGen"));
        param.Set_adapgen(cfg.lookup("simul.adapgen"));
        param.Set_pas(cfg.lookup("simul.pas"));

		}
		catch(const SettingNotFoundException &nfex)
		{
			cout << "No simul setting in configuration file." << endl;
			return(0);
		}
		return 1;


}








void ecrireParametres(Parameter &param)
							// ecrit les valeurs des parametres
{
	// Ecrit dans le fichier sortie:
	
	fichierS << "\n_________________________________________\n";
	fichierS << "\nn = " << param.Get_n();
	fichierS << "\nN = " << param.Get_N();
	fichierS << "\nm = " << param.Get_m();
	fichierS << "\a = " << param.Get_a();
	fichierS << ", s = " << param.Get_s();
	fichierS << ",h = " << param.Get_h();
	fichierS << ", sig_s = " << param.Get_sig_s();
	fichierS << ", i = " << param.Get_i();
	fichierS << ", sig_i = " << param.Get_sig_i();
	fichierS << ", U = " << param.Get_U();
	fichierS << ", nbS = " << param.Get_nbS();
	fichierS << ", L = " << param.Get_L();
	fichierS << ", loc1 = " <<  param.Get_loc();
	fichierS << "\ngenerations = " << param.Get_NbGen();
	fichierS << "\nmutation adaptative a partir de = " << param.Get_adapgen();
	fichierS << ", pas = " << param.Get_pas();
	if (param.Get_freerec() == 1)
		fichierS << "\free recombination !\n";

}
