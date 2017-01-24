#ifndef MATRIX_MATH_H_INCLUDED
#define MATRIX_MATH_H_INCLUDED

#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <ostream>
#include "tk/spline.h"

Matrix log ( const Matrix &m)
{
    Matrix result(m.sizeX(),m.sizeY());
    for (int i=0;i<m.sizeX();++i)
        for (int j=0;j<m.sizeY();++j)
        {
            if (m.value(i,j)<0)
            {
                cout<<"Negative number inside the Matrix, please check. "<<endl;
                break;
            }
            else result(i,j)=std::log(m.value(i,j));
        }
    return result;
}
Matrix3D log ( const Matrix3D &m)
{
    Matrix3D result(m.sizeX(),m.sizeY(),m.sizeZ());
    for(int k=0;k<m.sizeZ();++k)
        for (int j=0;j<m.sizeY();++j)
            for (int i=0;i<m.sizeX();++i)
            {
                if (m.value(i,j,k)<0)
                {
                    cout<<"Negative number inside the Matrix, please check. "<<endl;
                    break;
                }
                else result(i,j,k)=std::log(m.value(i,j,k));
            }
    return result;
}
Matrix exp ( const Matrix &m)
{
    Matrix result(m.sizeX(),m.sizeY());
    for (int i=0;i<m.sizeX();++i)
        for (int j=0;j<m.sizeY();++j)
            result(i,j)=std::exp(m.value(i,j));

    return result;
}
Matrix3D exp ( const Matrix3D &m)
{
    int I=m.sizeX(),J=m.sizeY(),K=m.sizeZ();
    Matrix3D result(I,J,K);
    for (int i=0;i<I;++i)
        for (int j=0;j<J;++j)
            for (int k=0;k<K;++k)
            result(i,j,k)=std::exp(m.value(i,j,k));

    return result;
}
Matrix abs (const Matrix &m)
{
    Matrix result;
    result.P()=m.P().cwiseAbs();
    return result;
}
Matrix3D abs (  Matrix3D &m)
{
    Matrix3D result(m);
    for(int i=0;i<m.sizeZ();++i) result.slice(i)=abs(m.slice(i));
    return result;
}
//Matrix_Sp abs ( const Matrix_Sp &m)
//{
//    Matrix_Sp result;
//    result=arma::abs(m.P());
//    return result;
//}
Matrix sin ( const Matrix &m)
{
    Matrix result(m.sizeX(),m.sizeY());
    for (int i=0;i<m.sizeX();++i)
        for (int j=0;j<m.sizeY();++j)
            result(i,j)=std::sin(m.value(i,j));

    return result;
}


double mean(Matrix &m)
{
    double mean=m.sum()/(m.sizeX()*m.sizeY());
    return mean;
}

double mean_sq(Matrix &m)
{
    Matrix result=m^2;
    return mean(result);

}
//Matrix sample_mean ( Matrix &m)
//{
//    Matrix result(1,m.sizeY());
//    for(int i=0;i<m.sizeY();++i)
//    {
//        Matrix col=m.col(i);
//        result(i)=mean(col);
//    }
//
//    return result;
//}

//Matrix sample_mean_adj ( Matrix &m)
//{
//    Matrix result(1,m.sizeY());
//    for(int i=0;i<m.sizeY();++i)
//    {
//        Matrix col=m.col(i);
//        result(i)=mean(col)/std::abs(mean(col))*sqrt(mean_sq(col));
//    }
//
//    return result;
//}


Matrix logsig ( const Matrix &m,double beta=1.0)
{
    Matrix result=1.0/(1.0+exp(-beta*m));
    return result;
}

Matrix tanh ( const Matrix &m)
{
    Matrix result=(exp(m)-exp(-m))/(exp(m)+exp(-m));
    return result;
}

Matrix tansig(const Matrix &m)
{
    return tanh(m);

}

Matrix d_logsig(const Matrix &m)
{
    Matrix result=(logsig(m)*(1.0-logsig(m)));
    return result;
}

Matrix d_tansig(const Matrix &m)
{
    Matrix result=(1.0-(tanh(m)^2));
    return result;
}

Matrix purelin(const Matrix &m)
{
    return m;
}

Matrix normalize(const Matrix &input,double R, double r)
{
    double M=input.max();
    double m=input.min();
    return (((R-r)*input+(M*r-m*R))/(M-m));
}

Matrix de_normalize(const Matrix &input, const Matrix &origin)
{
    double R=input.max();
    double r=input.min();
    double M=origin.max();
    double m=origin.min();
    return (((M-m)*input+(R*m-r*M))/(R-r));
}

double interp1_old(const Matrix &x,const Matrix &y,const double& xq)
{
    double result=0.0;
    if (x.sizeY()!=1||y.sizeY()!=1)
        cout<<"Error(interp1): input is not col vector."<<endl;
    else if(xq==x.min()||xq==x.max())
    {
        for(int i=0;i<x.sizeX();++i)
        {
            if(x.value(i)==xq) {
                if(y.value(i)<1e-15) return 0.0;
                else return y.value(i);
            }
            break;
        }
    }
    else if(xq<x.min()||xq>x.max())
    {
        if(xq==0||xq==-0) return 0;
        else {cout<<"out of range."<<endl;return 0;}
    }
    else{
        int i;
        if(x.value(0)>x.value(1))
            {for(i=0;i<x.sizeX();++i)
                if(x.value(i)<=xq) break;}
        else
            {for(i=0;i<x.sizeX();++i)
                if(x.value(i)>=xq) break;}

        result=y.value(i-1)+(y.value(i)-y.value(i-1))*(xq-x.value(i-1))/(x.value(i)-x.value(i-1));
    }
    return result;
}

class LinearInterp
{
private:
    Matrix prime;
    Matrix X,Y;
    int N;
    bool ascend;
public:
    LinearInterp(const Matrix &x,const Matrix &y)
    {
        X=x,Y=y;
        if(X.sizeY()==1) N=X.sizeX();
        else if(X.sizeX()==1) N=X.sizeY();
        else cout<<"interp1: input is not a vector."<<endl;
        if(X(0)<X(1)) ascend=1;
        else ascend=0;

        prime=Matrix(N-1);
        for(int i=1;i<N;++i) prime(i-1)=(Y(i)-Y(i-1))/(X(i)-X(i-1));
    }
    double Interp(double xq)
    {
        if(xq<X.min()||xq>X.max())
        {
            if(xq!=0||xq!=-0) cout<<"Interp1: input out of range."<<endl;
            return 0;
        }
        else{
            int i;
            for(i=1;i<N;++i)
                if(ascend==1)
                {
                    if(X(i)>=xq) break;
                }
                else
                {
                    if(X(i)<=xq) break;
                }

            return Y(i-1)+(xq-X(i-1))*prime(i-1);

        }
    }



};
double interp1(const Matrix &x,const Matrix &y,const double& xq)
{
    LinearInterp LI(x,y);
    return LI.Interp(xq);

}


Matrix interp1(const Matrix &x,const Matrix &y,const Matrix& xq)
{
    LinearInterp LI(x,y);
    int Nx=xq.sizeX(),Ny=xq.sizeY();
    Matrix result(Nx,Ny);
    #pragma omp parallel for
    for(int j=0;j<Ny;++j)
    {
        for(int i=0;i<Nx;++i)
            result(i,j)=LI.Interp(xq.value(i,j));
    }


    return result;
}

Matrix3D interp1(const Matrix &x,const Matrix &y,const Matrix3D& xq)
{
    LinearInterp LI(x,y);
    int Nx=xq.sizeX(),Ny=xq.sizeY(),Nz=xq.sizeZ();
    Matrix3D result(Nx,Ny,Nz);

        for(int k=0;k<Nz;++k)
        {
            #pragma omp parallel for
            for(int j=0;j<Ny;++j)
            {
                for(int i=0;i<Nx;++i)
                    result(i,j,k)=LI.Interp(xq.value(i,j,k));
            }
        }


    return result;
}
Matrix3D interp3(const Matrix &x,const Matrix &y,const Matrix3D& xq)
{

    int Nx=xq.sizeX(),Ny=xq.sizeY(),Nz=xq.sizeZ();
    Matrix3D result(Nx,Ny,Nz);
    int N=x.sizeX();

    if (N==1) cout<<"Error(interp1): input is not col vector."<<endl;
    else{

        vector<double> X(N),Y(N);
        for(int i=0;i<N;++i)
        {
            X[i]=x.value(i);
            Y[i]=y.value(i);
        }
        tk::spline s;
        s.set_points(X,Y);
        for(int k=0;k<Nz;++k)
            for(int j=0;j<Ny;++j)
                for(int i=0;i<Nx;++i)
                    if(xq.value(i,j,k)>0)
                        result(i,j,k)=s(xq.value(i,j,k));
    }

    return result;
}

//Matrix interp1(Matrix x,Matrix y,const Matrix& xq)
//{
//    Matrix result(xq.sizeX(),xq.sizeY());
//    if (x.sizeY()!=1||y.sizeY()!=1)
//        cout<<"Error(interp1): input is not col vector."<<endl;
//    else{
//        for(int i=0;i<xq.sizeX();++i)
//            for(int j=0;j<xq.sizeY();++j)
//                result(i,j)=interp1_d(x,y,xq.value(i,j));
//    }
//
//    return result;
//}

Matrix pow(Matrix &m,double c)
{
    return m^c;
}

//Matrix trans_sparse(Matrix &A,int n)
//{
//    int u=A.sizeX();
//    int v=A.sizeY();
//    double C_pos=(v/n-1)/2;
//    Matrix M(u,v);
//        for (int k=0;k<=u/n-C_pos-1;++k){
//            M.P()(span((k+C_pos)*n,(k+C_pos)*n+n-1),span(0,n-1))=A.P()(span(k*n,k*n+n-1),span(v-n,v-1)).t();
//            M.P()(span(k*n,k*n+n-1),span(v-n,v-1))=A.P()(span((k+C_pos)*n,(k+C_pos)*n+n-1),span(0,n-1)).t();
//        }
//        for (int k=0;k<=u/n-1-1;++k){
//            M.P()(span((k+1)*n,(k+1)*n+n-1),span((C_pos-1)*n,C_pos*n-1))=A.P()(span(k*n,k*n+n-1),span((C_pos+1)*n,(C_pos+2)*n-1)).t();
//            M.P()(span(k*n,k*n+n-1),span((C_pos+1)*n,(C_pos+2)*n-1))=A.P()(span((k+1)*n,(k+1)*n+n-1),span((C_pos-1)*n,C_pos*n-1)).t();
//        }
//        for (int k=0;k<=u/n-1;++k){
//            M.P()(span(k*n,k*n+n-1),span((C_pos)*n,(C_pos+1)*n-1))=A.P()(span(k*n,k*n+n-1),span((C_pos)*n,(C_pos+1)*n-1)).t();
//        }
//
//
//    return M;
//
//}
void remove_res(Matrix& table)
{
    for(int i=0;i<table.sizeX();++i)
        for(int j=0;j<table.sizeY();++j)
            if(table(i,j)<1e-15) table(i,j)=0;
}



#endif // MATRIX_MATH_H_INCLUDED
