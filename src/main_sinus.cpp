#include"sinus.h"
#include<iostream>
#include<fstream>
#include<iomanip>

int main()
{
    //definisikan solvernya
    Solver solve(20,20);
    //set domainnya (ini dalam meter ya)
    solve.setextents(14.5,5,14.5,5);
    //set parameter 
    solve.setParam(1000, 0.01);

    //fungsi gaussian bisa dipanggil di sini, bisa dibuat iterasi juga kalau diperlukan 
    solve.writeMult(10,10,10,10);
    return 0;
}