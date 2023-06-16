#include "sinus.h"
#include "Field.h"

#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;

Solver::Solver(int nz, int nr)
    : nz{nz},
      nr{nr},
      phi(nz, nr),
      rho(nz, nr),
      object_id(nz, nr)
{
}

Solver::~Solver()
{
    // this->phi.~Field_();
    // this->rho.~Field_();
    // this->object_id.~Field_();
}

// kayaknya nggak perlu destructor ya
void Solver::setextents(double zmax, double zmin, double rmax, double rmin)
{
    this->x0[0] = zmin;
    this->x0[1] = rmin;

    this->xm[0] = zmax;
    this->xm[1] = rmax;

    this->dh[0] = (zmax - zmin) / (nz - 1);
    this->dh[1] = (rmax - rmin) / (nr - 1);

    this->L[0] = zmax - zmin;
    this->L[1] = rmax - rmin;
}

void Solver::setParam(int iter, double tol)
{
    this->max_solver_it = iter;
    this->tolerance = tol;
}

bool Solver::solveGS()
{
    // di Solver itu phi dan rho sudah didefinisikan, jadi nggak perlu dipanggil lagi
    // dh juga udah diset kan
    double idz = 1 / this->dh[0];
    double idr = 1 / this->dh[1];
    double idz2 = idz * idz;
    double idr2 = idr * idr;

    double L2 = 0;
    bool converged = false;

    double crz = 0.5 / (idz2 + idr2);
    // berikut ini iterasinya
    for (unsigned it = 0; it < this->max_solver_it; it++)
    {
        for (int i = 0; i < this->nz; i++)
        {
            for (int j = 0; j < this->nr; j++)
            {
                // ini asumsinya kasus konduktor di empat syarat batas ya
                if (i == 0) // paling kiri
                {
                    continue; // nggak diapa-apain
                }
                else if (i == (this->nz) - 1) // paling kanan
                {
                    continue; // nggak diapa-apain
                }
                else if (j == 0) // paling bawah, nggak boleh r=0
                {
                    continue; // nggak diapa-apain
                }
                else if (j == (this->nr) - 1)
                {
                    continue; // nggak diapa-apain
                }
                else // selain di syarat batas
                {
                    // ini nggak sesederhana itu ya, karena rumus aslinya referensinya (j=0) itu ya di sumbu z.
                    double crj = 0.5 / ((this->x0[1] + j * this->dh[0]) * this->dh[0]);
                    // double phi_baru = crz*( (rho[i][j]/EPS_0) + idz2*(phi[i+1][j] + phi[i-1][j]) + phi[i][j+1]*(idr2 + 0.5*idr2/j) + phi[i][j-1]*(idr2 - 0.5*idr2/j) );
                    double phi_baru = crz * ((this->rho[i][j] / Const::EPS_0) + (idz2 * (this->phi[i + 1][j] + this->phi[i - 1][j])) + (this->phi[i][j + 1] * (idr2 + crj)) + (this->phi[i][j - 1] * (idr2 - crj)));
                    // lanjutkan dengan SOR
                    phi[i][j] += 1.4 * (phi_baru - this->phi[i][j]);
                }
            }
        }
        if (it % 50 == 0)
        {
            double sum = 0;
            for (int i = 1; i < (this->nz) - 1; i++)
            {
                for (int j = 1; j < (this->nr) - 1; j++)
                {
                    double crj = 0.5 / ((this->x0[1] + j * this->dh[0]) * this->dh[0]);
                    double R = -2 * this->phi[i][j] * (idr2 + idz2) + ((this->rho[i][j] / Const::EPS_0)
                            + idz2 * (this->phi[i + 1][j]
                            + this->phi[i - 1][j]) + this->phi[i][j + 1] * (idr2 + crj)
                            + this->phi[i][j - 1] * (idr2 - crj));
                    sum += R * R;
                }
            }
            L2 = sqrt(sum / (nz * nr));
            if (L2 < tolerance)
            {
                converged = true;
                break;
            }
        }
    }
    if (!converged)
    {
        cerr << "Gauss seidel standar gagal konvergen, L2=" << L2 << endl;
    }
    return converged;
}

void Solver::setRhoGaussian(int nz, int nr)
{
    // definisikan fungsinya dulu
    for (int i = 0; i < nz; i++)
    {
        for (int j = 0; j < nr; j++)
        {
            double z = x0[0] + i * dh[0];
            double r = x0[1] + j * dh[1];
            rho[i][j] = gaussFunc(z, r, nz, nr);
        }
    }
}

double Solver::gaussFunc(double z, double r, int nz, int nr)
{
    // definisi di bawah belum tentu tepat juga sih
    double partz = sin((nz * Const::PI * (z-zmin))/this->L[0]);
    double partr = sin((nr * Const::PI * (r-rmin))/this->L[1]);
    double fung = partz * partr;
    return fung;
}

void Solver::writeMult(int sz, int sr, int vz, int vr)
{
    ofstream gout("rhomult.csv"); // anggap file input (rapat muatan)
    ofstream fout("outmult.csv"); // file output (distribusi tegangan)

    // tentukan lebar langkahnya dulu
    double dvz = (xm[0] - x0[0]) / (sz - 1);
    double dvr = (xm[1] - x0[0]) / (sr - 1);

    // varminnya sama dengan dh aja kali ya, sementara varmax nya sama dengan xm
    double vvz = (xm[0] - x0[0] - dh[0]) / (vz - 1);
    // double vvr = (xm[1] - xm[1] - dh[1]) / (vr - 1); // TANDAI
    double vvr = (xm[1] - x0[1] - dh[1]) / (vr - 1);

    // bagian iterasinya
    for (int l = 0; l < vr; l++)
    {
        double varianceR = dh[1] + l * vvr;
        for (int k = 0; k < vz; k++)
        {
            double varianceZ = dh[0] + k * vvz;
            for (int j = 0; j < sr; j++)
            {
                double centerR = x0[1] + j * dvr;
                for (int i = 0; i < sz; i++)
                {
                    double centerZ = x0[0] + i * dvz;
                    // hitung distribusi rho
                    setRhoGaussian(centerZ, centerR, varianceZ, varianceR);
                    // hitung distribusi potensial
                    phi.clear();
                    solveGS();
                    // selanjutnya iterasi output (satu medan satu baris)
                    for (int z = 0; z < nz; z++)
                    {
                        for (int r = 0; r < nr; r++)
                        {
                            if (z == nz - 1 && r == nr - 1)
                            {
                                fout << phi[z][r];
                                gout << rho[z][r];
                            }
                            else
                            {
                                fout << phi[z][r] << ",";
                                gout << rho[z][r] << ",";
                            }
                        }
                    }
                    if ((i == sz - 1) && (j == sr - 1) && (k == vz - 1) && (l == vr - 1))
                    {
                        continue;
                    }
                    else
                    {
                        fout << "\n";
                        gout << "\n";
                    }
                }
            }
        }
    }

    gout.close();
    fout.close();
}
