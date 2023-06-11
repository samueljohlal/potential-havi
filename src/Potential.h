#ifndef POTENTIAL_H
#define POTENTIAL_H

#include"Field.h"

#include<cmath>
#include<fstream>
#include<iostream>

//setor konstanta yang kira-kira akan dipakai di sini
namespace Const
{
    const double EPS_0 = 8.8541878e-12; //epsilon nol
    const double QE = 1.602176565e-19; //muatan elektron coulomb (tanda masih positif)
    const double AMU = 1.660538921e-27; //atomic mass unit kg
    const double ME = 9.10938215e-31; //massa elektron kg
    const double K = 1.380648e-23; //konstanta Boltzmann
    const double PI = 3.141592653; //pi
    const double EvToK = QE/K; //electron volt ke kelvin
    const double SPEEDC = 299792458; //kecepatan cahaya dalam meter per detik
}

//ini kita definisikan juga domainnya, jadi nggak perlu dipisah antara world dengan solver.
//sementara untuk koordinat silinder, dengan jejari kecil lebih besar dari nol (dirichlet)
class Solver
{
public: 
    Solver(int nz, int nr);
    ~Solver();
    //kayaknya nggak perlu destructor ya
    void setextents(double zmax, double zmin, double rmax, double rmin);
    void setParam(int iter, double tol);
    bool solveGS();
    void setRhoGaussian(double cz, double cr, double varz, double varr); 
    double gaussFunc(double z, double r, double cz, double cr, double varz, double varr);

    /**
     * @brief Mencari Hasilnya
     * 
     * @param sz Jumlah langkah perubahan pusat distribusi sumbu-z
     * @param sr Jumlah langkah perubahan pusat distribusi sumbu-r
     * @param vz menunjukkan langkah perubahan varian sumbu-z
     * @param vr menunjukkan langkah perubahan varian sumbu-r
     * 
     * @note kalo mau banyak medan sekaligus, kita perlu nuliskan di sini
     * @note sementara untuk distribusi gaussian aja dengan satu puncak
     * 
     */
    void writeMult(int sz, int sr, int vz, int vr);

public:
    Field phi;
    Field rho;
    FieldI object_id;
    const int nz, nr;

private:
    double x0[2];
    double xm[2];
    double dh[2];
    int max_solver_it = 0;
    double tolerance = 0;
};

#endif