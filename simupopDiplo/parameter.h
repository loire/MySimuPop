#ifndef PARAMETER_H
#define PARAMETER_H
 
class Parameter
{
private:
	int nv;
	int Nv;
	double mv;
	int av;
	double sv;
	double hv;
	double sig_sv;
	double iv;
	double sig_iv;
	double Uv;
	int nbSv;
	double Lv;
	int locv;
	int NbGenv;
	int adapgenv;	
	int pasv;
	bool freerecv;
	bool no_mutv;
	
public:
    Parameter() {};    
    void Set_n(int n);
    void Set_N(int N);
    void Set_m(double m);
    void Set_a(int a);
    void Set_s(double s);
    void Set_h(double h);
	void Set_sig_s(double sig_s);
	void Set_i(double i);
	void Set_sig_i(double sig_i);
	void Set_U(double U);
	void Set_nbS(int nbS);
	void Set_L(double L);
	void Set_loc(int loc);
	void Set_NbGen(int NbGen);
	void Set_adapgen(int adapgen);
	void Set_pas(int pas);
	void Set_freerec(bool freerec);
	void Set_no_mut(bool no_mut);
    int Get_n() { return nv;}
    int Get_N() {return Nv;}
    double Get_m() {return mv;}
    int Get_a() {return av;}
    double Get_s() {return sv;}
    double Get_h() {return hv;}
	double Get_sig_s() {return sig_sv; }
	double Get_i() { return iv; }
	double Get_sig_i() {return sig_iv; }
	double Get_U() {return Uv;}
	int Get_nbS() {return nbSv;}
	double Get_L() {return Lv;}
	int Get_loc() {return locv;}
	int Get_NbGen() {return NbGenv;}
	int Get_adapgen() {return adapgenv;}
	int Get_pas() {return pasv;}
	bool Get_freerec(){return freerecv;} 
	bool Get_no_mut(){return no_mutv;}
};
#endif
