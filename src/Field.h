#ifndef FIELD_H
#define FIELD_H

#include<iostream>
#include<fstream>
#include<vector> //kayaknya nggak butuh, kita nggak bahas medan listrik kan

//kita bikin field 2 dimensi. pakai template juga
template<typename T>
class Field_
{
public:
    /**
     * @brief Matriks nz x nr
     * 
     * @param nz Banyak row
     * @param nr Banyak column
     */
    Field_(int nz, int nr) : nz{nz}, nr{nr}
    {
        data = new T*[nz]; //ini dua alamat ya
        for (int z=0; z<nz; z++)
        {
            data[z] = new T[nr];
        }
        clear(); //isi semua jadi nol
    }

    // Copy Assignment (Special Member Function)
    Field_ (const Field_ &other) : Field_{other.nz, other.nr}
    {
        for (int z=0; z<nz; z++)
        {
            for (int r=0; r<nr; r++)
            {
                data[z][r] = other(z,r); //perlu overload
            }
        }
    }

    // Copy Assignment (Special Member Function)
    Field_& operator= (Field_ &&f)
    {
        return (*this);
    }

    // Destructor
    ~Field_()
    {
        if(data==nullptr)
        {
            return;
        }
        for(int z=0; z < nz; z++)
        {
            delete[] data[z];
        }
        delete[] data;
    }

    /**
     * @brief Menge-nolkan semua value dari alokasi heap memori
     * 
     */
    void clear()
    {
        (*this) = 0; //ini = nya pake operator di atas kan
    }

    //kalo medan mau ditaruh di kanan, pakai friend (kayak di luar padahal di dalam ya)
    friend Field_<T> operator* (double s, const Field_<T> &f)
    {
        Field_<T> r(f); //copy constructor
        return r *= s; //medan f berubah menjadi r*s (dengan r=f)
    }

    /**
     * @brief Me-return nilai untuk row z, column r
     * 
     * @param z Row
     * @param r Column
     * @return double 
     */
    double get(int z, int r) const 
    {
        return data[z][r];
    }

    /**
     * Operator
     * ======== 
     */

    //beberapa operator yang berguna
    T* operator[] (int i) //ini perlu dipake dua kali
    {
        return data[i];
    }

    T operator() (int z, int r)
    {
        return data[z][r];
    }

    //inisiasi awal dibutuhkan oleh clear
    void operator= (double s)
    {
        for (int z=0; z<nz; z++)
        {
            for(int r=0; r<nr; r++)
            {
                data[z][r] = s;
            }
        }
    }

    //operator pembagian (ini perlu kalo mau ngitung rho)
    void operator /= (const Field_& other)
    {
        for (int z=0; z<nz; z++)
        {
            for(int r=0; r<nr; r++)
            {
                if (other.data[z][r] != 0)
                {
                    data[z][r] /= other.data[z][r]; //lho dapat akses?
                }
                else
                {
                    data[z][r] = 0; //ini agak meragukan, pembagian dengan nol soalnya
                }
            }
        }
    }

    Field_& operator+= (const Field_& other)
    {
        for(int z=0; z<nz; z++)
        {
            for(int r=0; r<nr; r++)
            {
                data[z][r] += other.data[z][r]; 
            }
        }
        return (*this);
    }

    //perkalian dengan skalar (medan selalu ditaruh di kiri ya)
    Field_& operator*= (double s)
    {
        for(int z=0; z<nz; z++)
        {
            for(int r=0; r<nr; r++)
            {
                data[z][r] *= s; 
            }
        }
        return (*this);
    }

public:
    const int nz, nr;//jumlah grid

private:
    T **data;
};

using Field = Field_<double>;
using FieldI = Field_<int>; //ini kalo perlu aja
#endif