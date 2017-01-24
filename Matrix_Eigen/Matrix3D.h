#ifndef MATRIX3D_H_INCLUDED
#define MATRIX3D_H_INCLUDED

#include <iostream>
#include <cmath>
#include <ostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <sstream>
#include <stdlib.h>
#include <time.h>
#include <omp.h>

#include "Eigen/Dense"

using namespace std;


class Matrix3D
{
private:
    int dx,dy,dz;
    Matrix **p;

    void allocArrays()
    {
        p=new Matrix*[dz];
        for(int i=0;i<dz;++i)
        {
            p[i]=new Matrix(dx,dy);
        }
    }

public:
    Matrix3D(int sizeX, int sizeY, int sizeZ);
    Matrix3D(int sizeX, int sizeY);
    Matrix3D(int sizeX);
    Matrix3D();
    ~Matrix3D();
    Matrix3D(const Matrix3D& m);

    int sizeX () const    {return dx;}
    int sizeY () const    {return dy;}
    int sizeZ () const    {return dz;}
    double value(int i, int j, int k) const {return (*p[k])(i,j);}
    double value(int i) const {return (*p[0])(i,0);}
    Matrix** P() const {return p;}

    double& operator    ()  (int X,int Y,int Z);
    double& operator    ()  (int X);
    double& operator    ()  (Matrix3D&, double c);
    double& operator    ()  (Matrix3D&, int c);
    Matrix3D  operator    +   (const Matrix3D& m) const;
    Matrix3D  operator    -   (const Matrix3D& m) const;
    Matrix3D  operator    +   (const Matrix3D&& m) const;
    Matrix3D  operator    -   (const Matrix3D&& m) const;
    Matrix3D& operator    =   (const Matrix3D& m);
    Matrix3D& operator    =   (const Matrix3D&& m);
//    Matrix3D& operator    =   (const SpMat<double>& m);

    Matrix3D& operator    =   (const double m[]);
    Matrix3D& operator    +=  (const Matrix3D& m);
    Matrix3D& operator    -=  (const Matrix3D& m);
    Matrix3D& operator    *=  (const Matrix3D& m);
    Matrix3D& operator    /=  (const Matrix3D& m);
    Matrix3D& operator    +=  (double m);
    Matrix3D& operator    -=  (double m);
    Matrix3D& operator    *=  (double m);
    Matrix3D& operator    /=  (double m);

    Matrix3D& operator    +=  (const Matrix3D&& m);
    Matrix3D& operator    -=  (const Matrix3D&& m);
    Matrix3D& operator    *=  (const Matrix3D&& m);
    Matrix3D& operator    /=  (const Matrix3D&& m);


    friend Matrix3D operator | (const Matrix3D& m1, const Matrix3D& m2);
    friend Matrix3D operator * (const Matrix3D& m1, const Matrix3D& m2);
    friend Matrix3D operator * (double c,const Matrix3D& m2);
    friend Matrix3D operator * (const Matrix3D& m1,double c);
    friend Matrix3D operator * (int c,const Matrix3D& m2);
    friend Matrix3D operator * (const Matrix3D& m1,int c);

    friend Matrix3D operator + (double c,const Matrix3D& m);
    friend Matrix3D operator + (const Matrix3D& m, double c);
    friend Matrix3D operator + (int c,const Matrix3D& m);
    friend Matrix3D operator + (const Matrix3D& m, int c);
    friend Matrix3D operator + (const Matrix3D& m);

    friend Matrix3D operator - (const Matrix3D& m);
    friend Matrix3D operator - (const Matrix3D& m, double c);
    friend Matrix3D operator - (double c,const Matrix3D& m);
    friend Matrix3D operator - (int c,const Matrix3D& m);
    friend Matrix3D operator - (const Matrix3D& m, int c);


    friend Matrix3D operator / (const Matrix3D& m1, const Matrix3D& m2);
    friend Matrix3D operator / (const Matrix3D& m1, double c);
    friend Matrix3D operator / (double c, const Matrix3D& m1);
    friend Matrix3D operator / (const Matrix3D& m1, int c);
    friend Matrix3D operator / (int c, const Matrix3D& m1);
    friend Matrix3D operator ^ (const Matrix3D& m1, double c);
    friend Matrix3D operator ^ (const Matrix3D& m1, int c);

    /// &&
    friend Matrix3D operator | (const Matrix3D&& m1, const Matrix3D&& m2);
    friend Matrix3D operator | (const Matrix3D& m1, const Matrix3D&& m2);
    friend Matrix3D operator | (const Matrix3D&& m1, const Matrix3D& m2);
    friend Matrix3D operator * (const Matrix3D&& m1, const Matrix3D&& m2);
    friend Matrix3D operator * (const Matrix3D& m1, const Matrix3D&& m2);
    friend Matrix3D operator * (const Matrix3D&& m1, const Matrix3D& m2);
    friend Matrix3D operator * (double c,const Matrix3D&& m2);
    friend Matrix3D operator * (const Matrix3D&& m1,double c);
    friend Matrix3D operator * (int c,const Matrix3D&& m2);
    friend Matrix3D operator * (const Matrix3D&& m1,int c);

    friend Matrix3D operator + (double c,const Matrix3D&& m);
    friend Matrix3D operator + (const Matrix3D&& m, double c);
    friend Matrix3D operator + (int c,const Matrix3D&& m);
    friend Matrix3D operator + (const Matrix3D&& m, int c);
    friend Matrix3D operator + (const Matrix3D&& m);

    friend Matrix3D operator - (const Matrix3D&& m);
    friend Matrix3D operator - (const Matrix3D&& m, double c);
    friend Matrix3D operator - (double c,const Matrix3D&& m);
    friend Matrix3D operator - (int c,const Matrix3D&& m);
    friend Matrix3D operator - (const Matrix3D&& m, int c);


    friend Matrix3D operator / (const Matrix3D&& m1, const Matrix3D&& m2);
    friend Matrix3D operator / (const Matrix3D& m1, const Matrix3D&& m2);
    friend Matrix3D operator / (const Matrix3D&& m1, const Matrix3D& m2);
    friend Matrix3D operator / (const Matrix3D&& m1, double c);
    friend Matrix3D operator / (double c, const Matrix3D&& m1);
    friend Matrix3D operator / (const Matrix3D&& m1, int c);
    friend Matrix3D operator / (int c, const Matrix3D&& m1);
    friend Matrix3D operator ^ (const Matrix3D&& m1, double c);
    friend Matrix3D operator ^ (const Matrix3D&& m1, int c);
    ///

    friend Matrix3D operator > (const Matrix3D& m1, const Matrix3D& m2);
    friend Matrix3D operator > (const Matrix3D& m1, double c);
    friend Matrix3D operator > (double c, const Matrix3D& m2);

    friend Matrix3D operator >= (const Matrix3D& m1, const Matrix3D& m2);
    friend Matrix3D operator >= (const Matrix3D& m1, double c);
    friend Matrix3D operator >= (double c, const Matrix3D& m2);

    friend Matrix3D operator < (const Matrix3D& m1, const Matrix3D& m2);
    friend Matrix3D operator < (const Matrix3D& m1, double c);
    friend Matrix3D operator < (double c, const Matrix3D& m2);

    friend Matrix3D operator <= (const Matrix3D& m1, const Matrix3D& m2);
    friend Matrix3D operator <= (const Matrix3D& m1, double c);
    friend Matrix3D operator <= (double c, const Matrix3D& m2);

    friend Matrix3D operator == (const Matrix3D& m1, const Matrix3D& m2);
    friend Matrix3D operator == (const Matrix3D& m1, double c);
    friend Matrix3D operator == (double c, const Matrix3D& m2);

    friend Matrix3D operator != (const Matrix3D& m1, const Matrix3D& m2);
    friend Matrix3D operator != (const Matrix3D& m1, double c);
    friend Matrix3D operator != (double c, const Matrix3D& m2);

    friend ostream &operator << (ostream &out, const Matrix3D &m);


    Matrix3D diag() const;
    void show(int) const;
    void show(string,int) const;
    Matrix3D trans() const;
    Matrix3D det() const;
    double det2() const;
    Matrix3D index() const;
    double sum() const;
    double max() const;
    double min() const;
    void remove_nan(double) ;
    void remove_nan(int) ;
    void change(const Matrix3D &, double, double);
    void change(const Matrix3D &, double, int);
    void change(const Matrix3D &, int, double);
    void change(const Matrix3D &, int, int);
    void change(const Matrix3D&,double);
    Matrix3D row(int i);
    Matrix3D col(int j);
    Matrix&   slice(int);
    Matrix3D range(int,int,int,int);
    void swap_row(int,int);
    void swap_col(int,int);


    static Matrix3D ones(int i, int j, int k);
    static Matrix3D zeros(int i, int j, int k);
    static Matrix3D random(int i, int j, int k);
    static Matrix3D readcsv   (string filename) ;
    static void   writecsv  (string filename, const  Matrix3D &m) ;
    static Matrix3D combine(Matrix3D* M[],int);






};


///DEFULT CONSTRUCTOR
Matrix3D::Matrix3D() : dx(1),dy(1),dz(1){allocArrays();}
///CONSTRUCTOR
Matrix3D::Matrix3D(int sizeX,int sizeY,int sizeZ):dx(sizeX),dy(sizeY),dz(sizeZ)
{
    allocArrays();
}

Matrix3D::Matrix3D(int sizeX,int sizeY):dx(sizeX),dy(sizeY),dz(1)
{
    allocArrays();
}

Matrix3D::Matrix3D(int sizeX):dx(sizeX),dy(1),dz(1)
{
    allocArrays();

}
///COPY CONSTRUCTOR
Matrix3D::Matrix3D(const Matrix3D& m):dx(m.dx),dy(m.dy),dz(m.dz)
{
    allocArrays();
    ///#pragma omp parallel for
    for(int i=0;i<dz;++i) (*p[i])=(*m.p[i]);
}

///DISTRUCTOR
Matrix3D::~Matrix3D()
{
    for(int i=0;i<dz;++i)
        delete p[i];///might be a problem w/o space

    delete[] p;
}

///()
double& Matrix3D::operator()(int X,int Y,int Z)
{
    return (*p[Z])(X,Y);
}

double& Matrix3D::operator()(int X)
{
    if(dx==1) return (*p[0])(0,X);
    else if(dy==1) return (*p[0])(X,0);
    else {cout<<"Not a vector. "<<endl; return (*p[0])(0,X);}
}

///=
Matrix3D& Matrix3D::operator=(const Matrix3D& m)
{
    ///judge
    if(this==&m)
    {
        return *this;
    }
    else
    {
        if((dx!=m.dx||dy!=m.dy||dz!=m.dz))
        {
            if (dx!=1&&dy!=1&&dz!=1) cout<<"Warning(=): dimonsion not agree. New matrix will be built."<<endl;
            this->~Matrix3D();
            dx=m.dx;
            dy=m.dy;
            dz=m.dz;
            allocArrays();
        }
        for(int i=0;i<dz;++i) *p[i]=*m.p[i];
    }
    return *this;
}
Matrix3D& Matrix3D::operator=(const Matrix3D&& m)
{
    ///judge
    if(this==&m)
    {
        return *this;
    }
    else
    {
        if((dx!=m.dx||dy!=m.dy||dz!=m.dz))
        {
            if (dx!=1&&dy!=1&&dz!=1) cout<<"Warning(=): dimonsion not agree. New matrix will be built."<<endl;
            this->~Matrix3D();
            dx=m.dx;
            dy=m.dy;
            dz=m.dz;
            allocArrays();
        }
        for(int i=0;i<dz;++i) *p[i]=*m.p[i];
    }
    return *this;
}



///+
Matrix3D Matrix3D::operator+(const Matrix3D& m) const
{
    ///judge
    Matrix3D temp(*this);
        if(dx!=m.dx||dy!=m.dy||dz!=m.dz)
            cout<<"Error(+): dimonsion not agree."<<endl;
        else{

            for(int i=0;i<dz;++i) (*temp.p[i])+=(*m.p[i]);
        }

    return temp;
}
Matrix3D Matrix3D::operator+(const Matrix3D&& m) const
{
    ///judge
    Matrix3D temp(*this);
        if(dx!=m.dx||dy!=m.dy||dz!=m.dz)
            cout<<"Error(+): dimonsion not agree."<<endl;
        else{

            for(int i=0;i<dz;++i) (*temp.p[i])+=(*m.p[i]);
        }

    return temp;
}
///-
Matrix3D Matrix3D::operator-(const Matrix3D& m) const
{
    ///judge
    Matrix3D temp(*this);
        if(dx!=m.dx||dy!=m.dy||dz!=m.dz)
            cout<<"Error(-): dimonsion not agree."<<endl;
        else
        {

            for(int i=0;i<dz;++i) (*temp.p[i])-=(*m.p[i]);
        }

    return temp;
}
Matrix3D Matrix3D::operator-(const Matrix3D&& m) const
{
    ///judge
    Matrix3D temp(*this);
        if(dx!=m.dx||dy!=m.dy||dz!=m.dz)
            cout<<"Error(-): dimonsion not agree."<<endl;
        else
        {

            for(int i=0;i<dz;++i) (*temp.p[i])-=(*m.p[i]);
        }

    return temp;
}

///+=
Matrix3D& Matrix3D::operator+=(const Matrix3D& m)
{
    ///judge
        if(dx!=m.dx||dy!=m.dy||dz!=m.dz)
            cout<<"Error(+=): dimonsion does not agree."<<endl;
        else *this=(*this)+m;

    return *this;
}
Matrix3D& Matrix3D::operator+=(const Matrix3D&& m)
{
    ///judge
        if(dx!=m.dx||dy!=m.dy||dz!=m.dz)
            cout<<"Error(+=): dimonsion does not agree."<<endl;
        else *this=(*this)+m;

    return *this;
}
Matrix3D& Matrix3D::operator+=(double m)
{
    ///judge
            *this=(*this)+m;

    return *this;
}
///-=
Matrix3D& Matrix3D::operator-=(const Matrix3D& m)
{
    ///judge

        if(dx!=m.dx||dy!=m.dy||dz!=m.dz)
            cout<<"Error(-=): dimonsion not agree."<<endl;
        else *this=(*this)-m;


    return *this;
}
Matrix3D& Matrix3D::operator-=(const Matrix3D&& m)
{
    ///judge

        if(dx!=m.dx||dy!=m.dy||dz!=m.dz)
            cout<<"Error(-=): dimonsion not agree."<<endl;
        else *this=(*this)-m;


    return *this;
}
Matrix3D& Matrix3D::operator-=(double m)
{
    ///judge
            *this=(*this)-m;

    return *this;
}
/// *=
Matrix3D& Matrix3D::operator*=(const Matrix3D& m)
{
    ///judge

        if(dx!=m.dx||dy!=m.dy||dz!=m.dz)
            cout<<"Error(*=): dimonsion not agree."<<endl;
        else{

            *this=(*this)*m;

        }

    return *this;
}
Matrix3D& Matrix3D::operator*=(const Matrix3D&& m)
{
    ///judge

        if(dx!=m.dx||dy!=m.dy||dz!=m.dz)
            cout<<"Error(*=): dimonsion not agree."<<endl;
        else{

            *this=(*this)*m;

        }

    return *this;
}
Matrix3D& Matrix3D::operator*=(double m)
{
    ///judge
            *this=(*this)*m;

    return *this;
}
/// /=
Matrix3D& Matrix3D::operator/=(const Matrix3D& m)
{
    ///judge

        if(dx!=m.dx||dy!=m.dy||dz!=m.dz)
            cout<<"Error(/=): dimonsion not agree."<<endl;
        else{
            *this=(*this)/m;
        }

    return *this;
}
Matrix3D& Matrix3D::operator/=(const Matrix3D&& m)
{
    ///judge

        if(dx!=m.dx||dy!=m.dy||dz!=m.dz)
            cout<<"Error(/=): dimonsion not agree."<<endl;
        else{
            *this=(*this)/m;
        }

    return *this;
}
Matrix3D& Matrix3D::operator/=(double m)
{
    ///judge

            *this=(*this)/m;


    return *this;
}
///  |

///m1*m2
Matrix3D operator * (const Matrix3D& m1,const Matrix3D& m2)
{
    Matrix3D prod(m1);
    if(m1.dx!=m2.dx||m1.dy!=m2.dy||m1.dz!=m2.dz)
        cout<<"Error(*): dimonsion does not agree."<<endl;
    else
    {
        for(int i=0;i<m1.dz;++i) (*prod.p[i])*=(*m2.p[i]);
    }
    return prod;
}

Matrix3D operator * (double c,const Matrix3D& m2)
{
    Matrix3D prod(m2);
    for(int i=0;i<m2.dz;++i) *prod.p[i]*=c;

    return prod;
}
Matrix3D operator * (const Matrix3D& m1,double c)
{
    return c*m1;
}

Matrix3D operator * (int c1,const Matrix3D& m2)
{
    return double(c1)*m2;
}
Matrix3D operator * (const Matrix3D& m1,int c1)
{
    return double(c1)*m1;
}


///c+m1
Matrix3D operator + (double c,const Matrix3D& m)
{
    Matrix3D plus_M(m);
    plus_M+=c*Matrix3D::ones(m.dx,m.dy,m.dz);

    return plus_M;
}
Matrix3D operator + (const Matrix3D& m, double c)
{
   return (c+m);
}

Matrix3D operator + (int c1,const Matrix3D& m)
{
    return double(c1)+m;
}
Matrix3D operator + (const Matrix3D& m, int c1)
{
    double c=double(c1);
    return (c+m);
}

Matrix3D operator + (const Matrix3D& m)
{
    return m;
}
///c-m1
Matrix3D operator - (const Matrix3D& m)
{
    return -1*m;
}

Matrix3D operator - (const Matrix3D& m,  double c)
{
   return m+(-c);
}
Matrix3D operator - ( double c,const Matrix3D& m)
{
    return -(m-c);
}


Matrix3D operator - ( int c1,const Matrix3D& m)
{
    return double(c1)-m;
}
Matrix3D operator - (const Matrix3D& m,  int c1)
{
    return m-double(c1);
}


///  /
Matrix3D operator / (const Matrix3D& m1,const Matrix3D& m2)
{
    Matrix3D prod(m1);
    if(m1.dx!=m2.dx||m1.dy!=m2.dy||m1.dz!=m2.dz)
        cout<<"Error(/): dimonsion does not agree."<<endl;
    else
        for(int i=0;i<m1.dz;++i) *prod.p[i]/=(*m2.p[i]);

    return prod;
}

Matrix3D operator / (const Matrix3D& m1, double c)
{
    Matrix3D result(m1);
    for(int i=0;i<m1.dz;++i) (*result.p[i])/=c;
    return result;
}

Matrix3D operator / (double c,const Matrix3D& m1)
{
    Matrix3D prod(m1);
    for(int i=0;i<m1.dz;++i) *prod.p[i]=c/(*prod.p[i]);

    return prod;
}

Matrix3D operator / (const Matrix3D& m, int c1)
{
    return m/double(c1);
}

Matrix3D operator / (int c1,const Matrix3D& m)
{
    return double(c1)/m;
}

/// ^
Matrix3D operator ^ (const Matrix3D& m1, double c)
{
    Matrix3D prod(m1.dx,m1.dy,m1.dz);
    for(int i=0;i<m1.dz;++i) *prod.p[i]=(*m1.p[i])^c;
    return prod;
}

Matrix3D operator ^ (const Matrix3D& m1, int c1)
{
    double c=double(c1);
    return m1^c;
}
/// &&

///  |

///m1*m2
Matrix3D operator * (const Matrix3D&& m1,const Matrix3D&& m2)
{
    Matrix3D prod(m1);
    if(m1.dx!=m2.dx||m1.dy!=m2.dy||m1.dz!=m2.dz)
        cout<<"Error(*): dimonsion does not agree."<<endl;
    else
    {
        for(int i=0;i<m1.dz;++i) (*prod.p[i])*=(*m2.p[i]);
    }
    return prod;
}
Matrix3D operator * (const Matrix3D& m1,const Matrix3D&& m2)
{
    Matrix3D prod(m1);
    if(m1.dx!=m2.dx||m1.dy!=m2.dy||m1.dz!=m2.dz)
        cout<<"Error(*): dimonsion does not agree."<<endl;
    else
    {
        for(int i=0;i<m1.dz;++i) (*prod.p[i])*=(*m2.p[i]);
    }
    return prod;
}
Matrix3D operator * (const Matrix3D&& m1,const Matrix3D& m2)
{
    Matrix3D prod(m1);
    if(m1.dx!=m2.dx||m1.dy!=m2.dy||m1.dz!=m2.dz)
        cout<<"Error(*): dimonsion does not agree."<<endl;
    else
    {
        for(int i=0;i<m1.dz;++i) (*prod.p[i])*=(*m2.p[i]);
    }
    return prod;
}
Matrix3D operator * (double c,const Matrix3D&& m2)
{
    Matrix3D prod(m2);
    for(int i=0;i<m2.dz;++i) *prod.p[i]*=c;

    return prod;
}
Matrix3D operator * (const Matrix3D&& m1,double c)
{
    return c*m1;
}

Matrix3D operator * (int c1,const Matrix3D&& m2)
{
    return double(c1)*m2;
}
Matrix3D operator * (const Matrix3D&& m1,int c1)
{
    return double(c1)*m1;
}


///c+m1
Matrix3D operator + (double c,const Matrix3D&& m)
{
    Matrix3D plus_M(m);
    plus_M+=c*Matrix3D::ones(m.dx,m.dy,m.dz);

    return plus_M;
}
Matrix3D operator + (const Matrix3D&& m, double c)
{
   return (c+m);
}

Matrix3D operator + (int c1,const Matrix3D&& m)
{
    return double(c1)+m;
}
Matrix3D operator + (const Matrix3D&& m, int c1)
{
    double c=double(c1);
    return (c+m);
}

Matrix3D operator + (const Matrix3D&& m)
{
    return m;
}
///c-m1
Matrix3D operator - (const Matrix3D&& m)
{
    return -1*m;
}

Matrix3D operator - (const Matrix3D&& m,  double c)
{
   return m+(-c);
}
Matrix3D operator - ( double c,const Matrix3D&& m)
{
    return -(m-c);
}


Matrix3D operator - ( int c1,const Matrix3D&& m)
{
    return double(c1)-m;
}
Matrix3D operator - (const Matrix3D&& m,  int c1)
{
    return m-double(c1);
}


///  /
Matrix3D operator / (const Matrix3D&& m1,const Matrix3D&& m2)
{
    Matrix3D prod(m1);
    if(m1.dx!=m2.dx||m1.dy!=m2.dy||m1.dz!=m2.dz)
        cout<<"Error(/): dimonsion does not agree."<<endl;
    else
        for(int i=0;i<m1.dz;++i) *prod.p[i]/=(*m2.p[i]);

    return prod;
}
Matrix3D operator / (const Matrix3D&& m1,const Matrix3D& m2)
{
    Matrix3D prod(m1);
    if(m1.dx!=m2.dx||m1.dy!=m2.dy||m1.dz!=m2.dz)
        cout<<"Error(/): dimonsion does not agree."<<endl;
    else
        for(int i=0;i<m1.dz;++i) *prod.p[i]/=(*m2.p[i]);

    return prod;
}
Matrix3D operator / (const Matrix3D& m1,const Matrix3D&& m2)
{
    Matrix3D prod(m1);
    if(m1.dx!=m2.dx||m1.dy!=m2.dy||m1.dz!=m2.dz)
        cout<<"Error(/): dimonsion does not agree."<<endl;
    else
        for(int i=0;i<m1.dz;++i) *prod.p[i]/=(*m2.p[i]);

    return prod;
}
Matrix3D operator / (const Matrix3D&& m1, double c)
{
    Matrix3D result(m1);
    for(int i=0;i<m1.dz;++i) (*result.p[i])/=c;
    return result;
}

Matrix3D operator / (double c,const Matrix3D&& m1)
{
    Matrix3D prod(m1);
    for(int i=0;i<m1.dz;++i) *prod.p[i]=c/(*prod.p[i]);

    return prod;
}

Matrix3D operator / (const Matrix3D&& m, int c1)
{
    return m/double(c1);
}

Matrix3D operator / (int c1,const Matrix3D&& m)
{
    return double(c1)/m;
}

/// ^
Matrix3D operator ^ (const Matrix3D&& m1, double c)
{
    Matrix3D prod(m1.dx,m1.dy,m1.dz);
    for(int i=0;i<m1.dz;++i) *prod.p[i]=(*m1.p[i])^c;
    return prod;
}

Matrix3D operator ^ (const Matrix3D&& m1, int c1)
{
    double c=double(c1);
    return m1^c;
}
///
/// <<
ostream &operator << (ostream &out, const Matrix3D &m)
{
    for (int k=0;k<m.dz;k++)
        for(int j=0;j<m.dx;j++)
        {
            for (int i=0;i<m.dy;i++)
                out<<m.value(i,j,k)<<" ";
            out<<endl;
        }
    return out;
}
/// >
Matrix3D operator > (const Matrix3D& m1, const Matrix3D& m2)
{
    int Nx=m1.dx,Ny=m1.dy,Nz=m1.dz;
    Matrix3D result(Nx,Ny,Nz);
    if(m1.dx==m2.dx&&m1.dy==m2.dy&&m1.dz==m2.dz)
        for(int i=0;i<Nx;++i)
            for(int j=0;j<Ny;j++)
                for(int k=0;k<Nz;k++)
                    result(i,j,k)=(m1.value(i,j,k)>m2.value(i,j,k));
    else
        cout<<"Error(>): Dimension does not agree."<<endl;
    return result;
}
Matrix3D operator > (const Matrix3D& m1, double c)
{
    int Nx=m1.dx,Ny=m1.dy,Nz=m1.dz;
    Matrix3D result(Nx,Ny,Nz);

        for(int i=0;i<Nx;++i)
            for(int j=0;j<Ny;j++)
                for(int k=0;k<Nz;k++)
                    result(i,j,k)=(m1.value(i,j,k)>c);

    return result;
}
Matrix3D operator > (double c, const Matrix3D& m1)
{
    int Nx=m1.dx,Ny=m1.dy,Nz=m1.dz;
    Matrix3D result(Nx,Ny,Nz);

        for(int i=0;i<Nx;++i)
            for(int j=0;j<Ny;j++)
                for(int k=0;k<Nz;k++)
                    result(i,j,k)=(c>m1.value(i,j,k));

    return result;
}
/// >=

Matrix3D operator >= (const Matrix3D& m1, const Matrix3D& m2)
{
    int Nx=m1.dx,Ny=m1.dy,Nz=m1.dz;
    Matrix3D result(Nx,Ny,Nz);
    if(m1.dx==m2.dx&&m1.dy==m2.dy&&m1.dz==m2.dz)
        for(int i=0;i<Nx;++i)
            for(int j=0;j<Ny;j++)
                for(int k=0;k<Nz;k++)
                    result(i,j,k)=(m1.value(i,j,k)>=m2.value(i,j,k));
    else
        cout<<"Error(>=): Dimension does not agree."<<endl;
    return result;
}
Matrix3D operator >= (const Matrix3D& m1, double c)
{
    int Nx=m1.dx,Ny=m1.dy,Nz=m1.dz;
    Matrix3D result(Nx,Ny,Nz);

        for(int i=0;i<Nx;++i)
            for(int j=0;j<Ny;j++)
                for(int k=0;k<Nz;k++)
                    result(i,j,k)=(m1.value(i,j,k)>=c);

    return result;
}
Matrix3D operator >= (double c, const Matrix3D& m1)
{
    int Nx=m1.dx,Ny=m1.dy,Nz=m1.dz;
    Matrix3D result(Nx,Ny,Nz);

        for(int i=0;i<Nx;++i)
            for(int j=0;j<Ny;j++)
                for(int k=0;k<Nz;k++)
                    result(i,j,k)=(c>=m1.value(i,j,k));

    return result;
}

/// <
Matrix3D operator < (const Matrix3D& m1, const Matrix3D& m2)
{
    return (m2>m1);
}
Matrix3D operator < (const Matrix3D& m1, double c)
{
    return (c>m1);
}
Matrix3D operator < (double c, const Matrix3D& m2)
{
    return (m2>c);
}
/// <=
Matrix3D operator <= (const Matrix3D& m1, const Matrix3D& m2)
{
    return (m2>=m1);
}
Matrix3D operator <= (const Matrix3D& m1, double c)
{
    return (c>=m1);
}
Matrix3D operator <= (double c, const Matrix3D& m2)
{
    return (m2>=c);
}
/// ==
Matrix3D operator == (const Matrix3D& m1, const Matrix3D& m2)
{
    int Nx=m1.dx,Ny=m1.dy,Nz=m1.dz;
    Matrix3D result(Nx,Ny,Nz);
    if(m1.dx==m2.dx&&m1.dy==m2.dy&&m1.dz==m2.dz)
        for(int i=0;i<Nx;++i)
            for(int j=0;j<Ny;j++)
                for(int k=0;k<Nz;k++)
                    result(i,j,k)=(m1.value(i,j,k)==m2.value(i,j,k));
    else
        cout<<"Error(==): Dimension does not agree."<<endl;
    return result;
}
Matrix3D operator == (const Matrix3D& m1, double c)
{
    int Nx=m1.dx,Ny=m1.dy,Nz=m1.dz;
    Matrix3D result(Nx,Ny,Nz);

        for(int i=0;i<Nx;++i)
            for(int j=0;j<Ny;j++)
                for(int k=0;k<Nz;k++)
                    result(i,j,k)=(c==m1.value(i,j,k));

    return result;
}
Matrix3D operator == (double c, const Matrix3D& m2)
{
    return (m2==c);
}
/// !=
Matrix3D operator != (const Matrix3D& m1, const Matrix3D& m2)
{
    int Nx=m1.dx,Ny=m1.dy,Nz=m1.dz;
    Matrix3D result(Nx,Ny,Nz);
    if(m1.dx==m2.dx&&m1.dy==m2.dy&&m1.dz==m2.dz)
        for(int i=0;i<Nx;++i)
            for(int j=0;j<Ny;j++)
                for(int k=0;k<Nz;k++)
                    result(i,j,k)=(m1.value(i,j,k)!=m2.value(i,j,k));
    else
        cout<<"Error(==): Dimension does not agree."<<endl;
    return result;
}
Matrix3D operator != (const Matrix3D& m1, double c)
{
    int Nx=m1.dx,Ny=m1.dy,Nz=m1.dz;
    Matrix3D result(Nx,Ny,Nz);

        for(int i=0;i<Nx;++i)
            for(int j=0;j<Ny;j++)
                for(int k=0;k<Nz;k++)
                    result(i,j,k)=(c!=m1.value(i,j,k));

    return result;
}
Matrix3D operator != (double c, const Matrix3D& m2)
{
    return (m2!=c);
}
///DISPLAY
void Matrix3D::show(int n=15) const
{
    cout.precision(n);
    for(int k=0;k<dz;k++){
    for(int j=0;j<dx;j++)
    {
        for(int i=0;i<dy;i++)
            cout<<this->value(j,i,k)<<"    ";

        cout<<endl;
    }
    cout<<endl;
    }
    cout<<endl;
}
void Matrix3D::show(string name,int n) const
{

    cout<<"Matrix3D: "<<name<<endl;
    show(n);
}

///// read csv
//Matrix3D Matrix3D::readcsv(string filename)
//{
//    string line,data;
//    /// find Y,X
//    ifstream file1(filename.c_str());
//    int X=1,Y=1;
//    getline(file1,line);
//    X=count(line.begin(),line.end(),',')+1;
//    while(getline(file1,line)) Y++;
//    file1.close();
//    /// read data to matrix
//    Matrix3D dataMatrix(Y,X);
//    bool status=(*dataMatrix.p).load(filename,csv_ascii);
//
//    if(status==0) cout<<"Loading file "<<filename<<" failed."<<endl;
//    return dataMatrix;
//}
//
//void Matrix3D::writecsv (string filename, const Matrix3D &m)
//{
//    (*m.p).save(filename,csv_ascii);
//}
//
///// combine
//Matrix3D Matrix3D::combine(Matrix3D* M[],int n)
//{
//    int X=0,Y=0;
//    for(int i=0;i<n;i++){
//        Y+=(*M[i]).sizeY();
//    }
//    Matrix3D result(X,Y);
//
//    int SUM=0;
//    for(int K=0;K<n;K++)
//    {
//        for(int i=0;i<(*M[K]).sizeX();i++)
//            for(int j=0;j<(*M[K]).sizeY();j++)
//                result(SUM+i,j)=(*M[K])(i,j);
//
//        SUM+=(*M[K]).sizeX();
//    }
//
//    return result;
//
//}


/// Transpose ??
//inline Matrix3D Matrix3D::trans() const
//{
//    Matrix3D m(dy,dx);
//    (*m.p)=(*p).t();
//
//    return m;
//}

/// index
Matrix3D Matrix3D::index() const
{
    Matrix3D m(dx,dy,dz);
    for(int i=0;i<dx;i++)
        for(int j=0;j<dy;j++)
            for(int k=0;k<dz;k++)
            {
                if((*this).value(i,j,k)!=0) m(i,j,k)=1;
            }

    return m;
}

double Matrix3D::sum() const
{
    double result=0.;
    for(int i=0;i<dz;++i) result+=p[i]->sum();
    return result;
}

Matrix3D Matrix3D::ones(int i, int j, int k)
{
    Matrix3D m(i,j,k);
    for(int n=0;n<k;++n) *m.p[n]=Matrix::ones(i,j);
    return m;
}

double Matrix3D::max() const
{
    Matrix result(dz);
    for(int i=0;i<dz;++i) result(i)=p[i]->max();
    return result.max();
}

double Matrix3D::min() const
{
    Matrix result(dz);
    for(int i=0;i<dz;++i) result(i)=p[i]->min();
    return result.min();
}

void Matrix3D::remove_nan(double x)
{
    ///#pragma omp parallel for
    for(int k=0;k<dz;k++) p[k]->remove_nan(x);
}

void Matrix3D::remove_nan(int x)
{
    remove_nan(double(x));
}

Matrix3D Matrix3D::zeros(int i, int j, int k)
{
    Matrix3D m(i,j,k);
    return m;
}

Matrix3D Matrix3D::random(int I, int J, int K)
{
    Matrix3D m(I,J,K);
    for(int i=0;i<K;++i) *m.p[i]=Matrix::random(I,J);
    return m;
}

void Matrix3D::change(const Matrix3D& m, double c1, double c2)
{
    ///#pragma omp parallel for
    for(int k=0;k<dz;k++) p[k]->change(*m.P()[k],c1,c2);
}

void Matrix3D::change(const Matrix3D& m, int c1, double c2)
{
    change(m,double(c1),c2);
}

void Matrix3D::change(const Matrix3D& m, double c1, int c2)
{
    change(m,c1,double(c2));
}

void Matrix3D::change(const Matrix3D& m, int c1, int c2)
{
    change(m,double(c1),double(c2));
}
void Matrix3D::change(const Matrix3D& m,double c)
{
    change(m,1,c);
}
//Matrix3D Matrix3D::row(int I)
//{
//    Matrix3D Row(1,dy);
//    *Row.p=(*p).row(I);
//
//    return Row;
//}
//
//Matrix3D Matrix3D::col(int I)
//{
//    Matrix3D Col(dx,1);
//    *Col.p=(*p).col(I);
//
//    return Col;
//}
Matrix& Matrix3D::slice(int I)
{
    return *p[I];
}
//Matrix3D Matrix3D::range(int m, int n,int k,int q)
//{
//        Matrix3D result(n-m+1,q-k+1);
//        *result.p=(*p).submat(m,n,k,q);
//
//        return result;
//}
//
//void Matrix3D::swap_row(int i,int j)
//{
//    if (i>=dx||j>=dx) cout<<"Error(swap_row): out of range. "<<endl;
//    else{
//        (*p).swap_rows(i,j);
//    }
//}
//void Matrix3D::swap_col(int i,int j)
//{
//    if (i>=dy||j>=dy) cout<<"Error(swap_col): out of range. "<<endl;
//    else{
//        (*p).swap_cols(i,j);
//    }
//}
//Matrix3D Matrix3D::diag() const
//{
//    Matrix3D result(dx,dy);
//    if(dx!=dy)
//        cout<<"Error(diag()): Input is not a square matrix."<<endl;
//    else
//        *result.p=diagmat(*p);
//
//    return result;
//}


#endif // MATRIX3D_H_INCLUDED
