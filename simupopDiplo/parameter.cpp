#include "parameter.h"
#include <iostream>
#include <fstream>
#include <fstream>
#include <boost/tokenizer.hpp>
#include <string>
// Date member function
using namespace std;



void Parameter::Set_n(int n)
{
    nv = n;
}
void Parameter::Set_N(int N)
{
    Nv = N;
}
void Parameter::Set_m(double m)
{
	mv=m;
}
void Parameter::Set_a(int a)
{
	av=a;
}
void Parameter::Set_s(double s)
{
	sv=s;
}
void Parameter::Set_h(double h)
{
	hv=h;
}
void Parameter::Set_sig_s(double sig_s)
{
	sig_sv=sig_s;
}
void Parameter::Set_i(double i)
{
	iv=i;
}
void Parameter::Set_sig_i(double sig_i)
{ 
	sig_iv=sig_i; 
}
void Parameter::Set_U(double U)
{
	Uv=U;
}
void Parameter::Set_nbS(int nbS)
{
	nbSv=nbS;
}
void Parameter::Set_L(double L)
{
	Lv=L;
}
void Parameter::Set_r(double r)
{
	rv=r;
}
void Parameter::Set_rn(double rn)
{
	rnv=rn;
}

void Parameter::Set_loc(int loc)
{
	locv=loc;
}
void Parameter::Set_NbGen(int NbGen)
{
	NbGenv=NbGen;
}
void Parameter::Set_adapgen(int adapgen)
{
	adapgenv=adapgen;
}
void Parameter::Set_pas(int pas)
{
	pasv=pas;
}
void Parameter::Set_freerec(bool freerec)
{ 
	freerecv=freerec;
}
void Parameter::Set_no_mut(bool no_mut)
{
	no_mutv=no_mut;
}


