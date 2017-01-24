
#ifndef MATRIX_SP_H_INCLUDED
#define MATRIX_SP_H_INCLUDED
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



using namespace std;


typedef Eigen::SparseMatrix<double,Eigen::ColMajor> SpMat;
//typedef Eigen::SparseMatrix<double,Eigen::RowMajor> SpMat;

typedef Eigen::Triplet<double> T;

class Matrix_Sp
{
private:
    int dx,dy;
    SpMat *p;

    void allocArrays()
    {
        p=new SpMat(dx,dy);
    }

public:
    Matrix_Sp(int sizeX, int sizeY);
    Matrix_Sp(int sizeX);
    Matrix_Sp();
    ~Matrix_Sp();
    Matrix_Sp(const Matrix_Sp& m);

    int sizeX () const    {return dx;}
    int sizeY () const    {return dy;}
    double value(int i, int j) const {return (*p).coeffRef(i,j);}
    SpMat& P() const {return *p;}
    void reserve(int n){p->reserve(Eigen::VectorXi::Constant(dy,n));}

    double& operator    ()  (int X,int Y);
    double& operator    ()  (int X);
    double& operator    ()  (Matrix_Sp&, double c);
    double& operator    ()  (Matrix_Sp&, int c);
    Matrix_Sp  operator    +   (const Matrix_Sp& m);
    Matrix_Sp  operator    -   (const Matrix_Sp& m);
    Matrix_Sp& operator    =   (const Matrix_Sp& m);
    Matrix_Sp& operator    =   (const Matrix_Sp&& m);
    Matrix_Sp& operator    =   (const MatrixXd& m);
    Matrix_Sp& operator    =   (const SpMat& m);
    Matrix_Sp& operator    =   (const Matrix& m);
    Matrix_Sp& operator    =   (const double m[]);
    Matrix_Sp& operator    +=  (const Matrix_Sp& m);
    Matrix_Sp& operator    -=  (const Matrix_Sp& m);
    Matrix_Sp& operator    *=  (const Matrix_Sp& m);
    Matrix_Sp& operator    /=  (const Matrix_Sp& m);
    Matrix_Sp& operator    +=  (double m);
    Matrix_Sp& operator    -=  (double m);
    Matrix_Sp& operator    *=  (double m);
    Matrix_Sp& operator    /=  (double m);


    friend Matrix_Sp operator | (const Matrix_Sp& m1, const Matrix_Sp& m2);
    friend Matrix    operator | (const Matrix_Sp& m1, const Matrix& m2);
    friend Matrix_Sp operator * (const Matrix_Sp& m1, const Matrix_Sp& m2);
    friend Matrix_Sp operator * (double c,const Matrix_Sp& m2);
    friend Matrix_Sp operator * (const Matrix_Sp& m1,double c);
    friend Matrix_Sp operator * (int c,const Matrix_Sp& m2);
    friend Matrix_Sp operator * (const Matrix_Sp& m1,int c);

    friend Matrix_Sp operator + (double c,const Matrix_Sp& m);
    friend Matrix_Sp operator + (const Matrix_Sp& m, double c);
    friend Matrix_Sp operator + (int c,const Matrix_Sp& m);
    friend Matrix_Sp operator + (const Matrix_Sp& m, int c);
    friend Matrix_Sp operator + (const Matrix_Sp& m);

    friend Matrix_Sp operator - (const Matrix_Sp& m);
    friend Matrix_Sp operator - (const Matrix_Sp& m, double c);
    friend Matrix_Sp operator - (double c,const Matrix_Sp& m);
    friend Matrix_Sp operator - (int c,const Matrix_Sp& m);
    friend Matrix_Sp operator - (const Matrix_Sp& m, int c);


    friend Matrix_Sp operator / (const Matrix_Sp& m1, const Matrix_Sp& m2);
    friend Matrix_Sp operator / (const Matrix_Sp& m1, double c);
    friend Matrix_Sp operator / (double c, const Matrix_Sp& m1);
    friend Matrix_Sp operator / (const Matrix_Sp& m1, int c);
    friend Matrix_Sp operator / (int c, const Matrix_Sp& m1);
    friend Matrix_Sp operator ^ (const Matrix_Sp& m1, double c);
    friend Matrix_Sp operator ^ (const Matrix_Sp& m1, int c);

    friend ostream &operator << (ostream &out, const Matrix_Sp &m);




    Matrix_Sp diag();
    void show(string,int) const;
    void show(int) const;
    Matrix_Sp trans() const;
    Matrix_Sp det() const;
    double det2() const;
    Matrix_Sp index() const;
    double sum() const;
    double max() const;
    double min() const;
    void remove_nan(double) ;
    void remove_nan(int) ;
    void change(const Matrix_Sp &, double, double);
    void change(const Matrix_Sp &, double, int);
    void change(const Matrix_Sp &, int, double);
    void change(const Matrix_Sp &, int, int);
    Matrix_Sp row(int i);
    Matrix_Sp col(int j);
    Matrix_Sp range(int,int,int,int);
    void swap_row(int,int);
    void swap_col(int,int);




    static Matrix_Sp ones(int i, int j);
    static Matrix_Sp zeros(int i, int j);
    static Matrix_Sp random(int i, int j);
    static Matrix_Sp readcsv   (string filename) ;
    static void   writecsv  (string filename, const  Matrix_Sp &m) ;
    static Matrix_Sp combine(Matrix_Sp* M[],int);


    static Matrix SOLVE(Matrix_Sp &A,Matrix &B)
    {
        Matrix result(B);
        Eigen::SimplicialCholesky<SpMat> chol(A.P());  // performs a Cholesky factorization of A
        result.P() = chol.solve(B.P());

        return result;
    }
    static Matrix SOLVER_CG(Matrix_Sp &A,Matrix &B, double Tol=1e-6,int Itr=100,bool disp=0)
    {
        Matrix result(B);
        Eigen::ConjugateGradient<SpMat> cg;
        cg.setTolerance(Tol);
        cg.setMaxIterations(Itr);
        cg.compute(A.P());

        result.P() = cg.solve(B.P());
        if(disp!=0)
        {
            std::cout << "#iterations:     " << cg.iterations() << std::endl;
            std::cout << "estimated error: " << cg.error()      << std::endl;
        }

        return result;
    }
    static Matrix SOLVER_BiCGSTAB(Matrix_Sp &A,Matrix &B, double Tol=1e-6,int Itr=100,bool disp=0)
    {
        Matrix result(B);
        Eigen::BiCGSTAB<SpMat> cg;
        cg.setTolerance(Tol);
        cg.setMaxIterations(Itr);
        cg.compute(A.P());

        result.P() = cg.solve(B.P());
        if(disp!=0)
        {
            std::cout << "#iterations:     " << cg.iterations() << std::endl;
            std::cout << "estimated error: " << cg.error()      << std::endl;
        }

        return result;
    }
    static Matrix SOLVER_BiCGSTAB_ImLUT(Matrix_Sp &A,Matrix &B, double Tol=1e-6,int Itr=100,bool disp=0)
    {
        Matrix result(B);
        Eigen::BiCGSTAB<SpMat,Eigen::IncompleteLUT<double> > cg;
        cg.setTolerance(Tol);
        cg.setMaxIterations(Itr);
        cg.compute(A.P());

        result.P() = cg.solve(B.P());
        if(disp!=0)
        {
            std::cout << "#iterations:     " << cg.iterations() << std::endl;
            std::cout << "estimated error: " << cg.error()      << std::endl;
        }


        return result;
    }
    static Matrix SOLVER_DGMRES(Matrix_Sp &A,Matrix &B, double Eigenvalue=1,int restart=30,bool disp=0)
    {
        Matrix result(B);
        Eigen::DGMRES<SpMat> cg;
//        cg.setTolerance(Tol);
//        cg.setMaxIterations(Itr);
        cg.set_restart(restart); /// Set restarting value
        cg.setEigenv(Eigenvalue);
        cg.compute(A.P());

        result.P() = cg.solve(B.P());
//        if(disp!=0)
//        {
//            std::cout << "#iterations:     " << cg.iterations() << std::endl;
//            std::cout << "estimated error: " << cg.error()      << std::endl;
//        }


        return result;
    }
    static Matrix SOLVER_DGMRES_ImLUT(Matrix_Sp &A,Matrix &B, double Tol=1e-6,int Itr=100,bool disp=0,int Fillfactor=0,double Droptol=1e-6)
    {
        Matrix result(B);
//        Eigen::IncompleteLUT<double> ImLUT;

        Eigen::DGMRES<SpMat,Eigen::IncompleteLUT<double> > cg;
        cg.preconditioner().setFillfactor(Fillfactor);
        cg.preconditioner().setDroptol(Droptol);
        cg.setTolerance(Tol);
        cg.setMaxIterations(Itr);
        cg.set_restart(30); /// Set restarting value
        cg.setEigenv(1);
        cg.disp=disp;
        cg.compute(A.P());

        result.P() = cg.solve(B.P());
        if(disp!=0)
        {
            std::cout << "#iterations:     " << cg.iterations() << std::endl;
            std::cout << "estimated error: " << cg.error()      << std::endl;
        }


        return result;
    }
    static Matrix SOLVER_LU(Matrix_Sp &A,Matrix &B, double Tol=1e-6,int Itr=100)
    {
        Matrix result(B);
        Eigen::SparseLU<SpMat, Eigen::COLAMDOrdering<int> > solver;
        A.P().makeCompressed();
        solver.analyzePattern(A.P());
        solver.factorize(A.P());

        result.P() = solver.solve(B.P());


        return result;
    }
    static Matrix SOLVER_QR(Matrix_Sp &A,Matrix &B, double Tol=1e-6,int Itr=100)
    {
        Matrix result(B);
        Eigen::SparseQR<SpMat, Eigen::COLAMDOrdering<int> > solver;
        A.P().makeCompressed();
        solver.analyzePattern(A.P());
        solver.factorize(A.P());

        result.P() = solver.solve(B.P());


        return result;
    }

};


///DEFULT CONSTRUCTOR
Matrix_Sp::Matrix_Sp() : dx(1),dy(1) {allocArrays();}
///CONSTRUCTOR
Matrix_Sp::Matrix_Sp(int sizeX,int sizeY):dx(sizeX),dy(sizeY)
{
    allocArrays();
}

Matrix_Sp::Matrix_Sp(int sizeX):dx(sizeX),dy(1)
{
    allocArrays();

}
///COPY CONSTRUCTOR
Matrix_Sp::Matrix_Sp(const Matrix_Sp& m):dx(m.dx),dy(m.dy)
{
    allocArrays();
    ///#pragma omp parallel for
    (*p)=(*m.p);
}

///DISTRUCTOR
Matrix_Sp::~Matrix_Sp()
{
    delete p;
}

///()
double& Matrix_Sp::operator()(int X,int Y)
{

    return (*p).coeffRef(X,Y);
}
double& Matrix_Sp::operator()(int X)
{
    if(dx==1) return (*p).coeffRef(0,X);
    else if(dy==1) return (*p).coeffRef(X,0);
    else {cout<<"Not a vector. "<<endl; return (*p).coeffRef(0,X);}
}

///=
Matrix_Sp& Matrix_Sp::operator=(const Matrix_Sp& m)
{
    ///judge
//    if(this==&m)
//    {
//        return *this;
//    }
//    else
//    {

        if((dx!=m.dx||dy!=m.dy))
        {
            if (dx!=1&&dy!=1) cout<<"Warning(=): dimonsion not agree. New matrix will be built."<<endl;
            this->~Matrix_Sp();
            dx=m.dx;
            dy=m.dy;
            allocArrays();
        }
        (*p)=(*m.p);
//    }
    return *this;
}
Matrix_Sp& Matrix_Sp::operator=(const Matrix_Sp&& m)
{
    ///judge
//    if(this==&m)
//    {
//        return *this;
//    }
//    else
//    {
        if((dx!=m.dx||dy!=m.dy))
        {
            if (dx!=1&&dy!=1) cout<<"Warning(=): dimonsion not agree. New matrix will be built."<<endl;
            this->~Matrix_Sp();
            dx=m.dx;
            dy=m.dy;
            allocArrays();
        }
        (*p)=(*m.p);
//    }
    return *this;
}
Matrix_Sp& Matrix_Sp::operator=(const MatrixXd& m)
{
//    if(p==&m)
//    {
//        return *this;
//    }
//    else
//    {
        if((dx!=int(m.rows())||dy!=int(m.cols())))
        {
            if (dx!=1&&dy!=1) cout<<"Warning(=): dimonsion not agree. New matrix will be built."<<endl;
            this->~Matrix_Sp();
            dx=m.rows();
            dy=m.cols();
            p=new SpMat (dx,dy);
        }
        for(int i=0;i<dx;i++)
            for(int j=0;j<dy;j++)
                if(m(i,j)!=0)
                    (*p).coeffRef(i,j)=m(i,j);
//    }
    return *this;
}
Matrix_Sp& Matrix_Sp::operator=(const Matrix& m)
{
    ///judge
//    if(this==&m)
//    {
//        return *this;
//    }
//    else
//    {
        if((dx!=m.sizeX()||dy!=m.sizeY()))
        {
            if (dx!=1&&dy!=1) cout<<"Warning(=): dimonsion not agree. New matrix will be built."<<endl;
            this->~Matrix_Sp();
            dx=m.sizeX();
            dy=m.sizeY();
            p=new SpMat (dx,dy);
        }
        *this=m.P();
//    }
    return *this;
}

Matrix_Sp& Matrix_Sp::operator=(const SpMat& m)
{
//    if(p==&m)
//    {
//        return *this;
//    }
//    else
//    {
        if((dx!=int(m.rows())||dy!=int(m.cols())))
        {
            if (dx!=1&&dy!=1) cout<<"Warning(=): dimonsion not agree. New matrix will be built."<<endl;
            this->~Matrix_Sp();
            dx=m.rows();
            dy=m.cols();
            p=new SpMat (dx,dy);
        }
        *p=m;
//    }
    return *this;
}
Matrix_Sp& Matrix_Sp::operator = (const double m[])
{
        ///#pragma omp parallel for
        for(int i=0;i<dx;i++)
            for(int j=0;j<dy;j++)
                (*p).coeffRef(i,j)=m[i*dy+j];

    return *this;
}
///+
Matrix_Sp Matrix_Sp::operator+(const Matrix_Sp& m)
{
    ///judge
    Matrix_Sp temp(*this);
        if(dx!=m.dx||dy!=m.dy)
            cout<<"Error(+): dimonsion not agree."<<endl;
        else{

            (*temp.p)=(*temp.p)+(*m.p);
        }

    return temp;
}
///-
Matrix_Sp Matrix_Sp::operator-(const Matrix_Sp& m)
{
    ///judge
    Matrix_Sp temp(*this);
        if(dx!=m.dx||dy!=m.dy)
            cout<<"Error(-): dimonsion not agree."<<endl;
        else
        {

            (*temp.p)=(*temp.p)-(*m.p);
        }

    return temp;
}

///+=
Matrix_Sp& Matrix_Sp::operator+=(const Matrix_Sp& m)
{
    ///judge
        if(dx!=m.dx||dy!=m.dy)
            cout<<"Error(+=): dimonsion does not agree."<<endl;
        else{

            (*p)=(*p)+(*m.p);
        }

    return *this;
}
Matrix_Sp& Matrix_Sp::operator+=(double m)
{
    ///judge
    for(int i=0;i<dx;i++)
        for(int j=0;j<dy;j++)
            (*p).coeffRef(i,j)=(*p).coeffRef(i,j)+m;

    return *this;
}
///-=
Matrix_Sp& Matrix_Sp::operator-=(const Matrix_Sp& m)
{
    ///judge

        if(dx!=m.dx||dy!=m.dy)
            cout<<"Error(-=): dimonsion not agree."<<endl;
        else{
            (*p)=(*p)-(*m.p);
        }

    return *this;
}
Matrix_Sp& Matrix_Sp::operator-=(double m)
{
    ///judge
    for(int i=0;i<dx;i++)
        for(int j=0;j<dy;j++)
            (*p).coeffRef(i,j)=(*p).coeffRef(i,j)-m;

    return *this;
}
/// *=
Matrix_Sp& Matrix_Sp::operator*=(const Matrix_Sp& m)
{
    ///judge

        if(dx!=m.dx||dy!=m.dy)
            cout<<"Error(*=): dimonsion not agree."<<endl;
        else{
            (*p)=(*p).cwiseProduct(*m.p);
        }

    return *this;
}
Matrix_Sp& Matrix_Sp::operator*=(double m)
{
    ///judge
            (*p)=(*p)*m;

    return *this;
}
/// /=
//Matrix_Sp& Matrix_Sp::operator/=(const Matrix_Sp& m)
//{
//    ///judge
//
//        if(dx!=m.dx||dy!=m.dy)
//            cout<<"Error(/=): dimonsion not agree."<<endl;
//        else{
//        *p=p->array()/m.p->array();
//        }
//
//    return *this;
//}
Matrix_Sp& Matrix_Sp::operator/=(double m)
{
    ///judge

            (*p)=(*p)/m;


    return *this;
}
///  |
Matrix_Sp operator | (const Matrix_Sp& m1,const Matrix_Sp& m2)
{
    Matrix_Sp prod(m1.dx,m2.dy);
    if(m1.dy!=m2.dx)
        cout<<"Error(|): dimonsion does not agree."<<endl;
    else
    {
        *prod.p=(*m1.p)*(*m2.p);
    }
    return prod;
}
Matrix operator | (const Matrix_Sp& m1,const Matrix& m2)
{
    Matrix prod;
    if(m1.dy!=m2.sizeX())
        cout<<"Error(|): dimonsion does not agree."<<endl;
    else
    {
        prod=Matrix(Eigen::MatrixXd((*m1.p)*(m2.P())));
    }
    return prod;
}
///m1*m2
Matrix_Sp operator * (const Matrix_Sp& m1,const Matrix_Sp& m2)
{
    Matrix_Sp prod(m1.dx,m1.dy);
    if(m1.dx!=m2.dx||m1.dy!=m2.dy)
        cout<<"Error(*): dimonsion does not agree."<<endl;
    else
    {
        *prod.p=(*m1.p).cwiseProduct(*m2.p);
    }
    return prod;
}
Matrix_Sp operator * (double c,const Matrix_Sp& m2)
{
    Matrix_Sp prod(m2.dx,m2.dy);
    *prod.p=c*(*m2.p);

    return prod;
}
Matrix_Sp operator * (const Matrix_Sp& m1,double c)
{
    return c*m1;
}

Matrix_Sp operator * (int c1,const Matrix_Sp& m2)
{
    return double(c1)*m2;
}
Matrix_Sp operator * (const Matrix_Sp& m1,int c1)
{
    return double(c1)*m1;
}


///c+m1
Matrix_Sp operator + (double c,const Matrix_Sp& m)
{
    Matrix_Sp plus_M=m;
    for(int i=0;i<m.dx;i++)
        for(int j=0;j<m.dy;j++)
            (*plus_M.p).coeffRef(i,j)+=c;

    return plus_M;
}
Matrix_Sp operator + (const Matrix_Sp& m, double c)
{
   return (c+m);
}

Matrix_Sp operator + (int c1,const Matrix_Sp& m)
{
    return double(c1)+m;
}
Matrix_Sp operator + (const Matrix_Sp& m, int c1)
{
    double c=double(c1);
    return (c+m);
}

Matrix_Sp operator + (const Matrix_Sp& m)
{
    return m;
}
///c-m1
Matrix_Sp operator - (const Matrix_Sp& m)
{
    return -1*m;
}

Matrix_Sp operator - (const Matrix_Sp& m,  double c)
{
   return m+(-c);
}
Matrix_Sp operator - ( double c,const Matrix_Sp& m)
{
    return -(m-c);
}


Matrix_Sp operator - ( int c1,const Matrix_Sp& m)
{
    return double(c1)-m;
}
Matrix_Sp operator - (const Matrix_Sp& m,  int c1)
{
    return m-double(c1);
}


///  /
//Matrix_Sp operator / (const Matrix_Sp& m1,const Matrix_Sp& m2)
//{
//    Matrix_Sp prod=m1;
//    if(m1.dx!=m2.dx||m1.dy!=m2.dy)
//        cout<<"Error(/): dimonsion does not agree."<<endl;
//    else
//
//        *prod.p=(*prod.p)/(*m2.p);
//
//    return prod;
//}

Matrix_Sp operator / (const Matrix_Sp& m1, double c)
{
    Matrix_Sp result=m1;
    (*result.p)=(*m1.p)/c;
    return result;
}

Matrix_Sp operator / (double c,const Matrix_Sp& m1)
{
    Matrix_Sp prod(m1.dx,m1.dy);
    for(int i=0;i<m1.dx;i++)
        for(int j=0;j<m1.dy;j++)
            if((*prod.p).coeffRef(i,j)!=0)
                (*prod.p).coeffRef(i,j)=c/(*prod.p).coeffRef(i,j);

    return prod;
}

Matrix_Sp operator / (const Matrix_Sp& m, int c1)
{
    return m/double(c1);
}

Matrix_Sp operator / (int c1,const Matrix_Sp& m)
{
    return double(c1)/m;
}


///DISPLAY
void Matrix_Sp::show(int n=15) const
{

    cout.precision(n);
    for(int j=0;j<dx;j++)
    {
        for(int i=0;i<dy;i++)
            cout<<(*this).value(j,i)<<"    ";

        cout<<endl;
    }
    cout<<endl;
}
void Matrix_Sp::show(string name, int n=15) const
{
    cout<<"Matrix: "<<name<<endl;
    show(n);
}
/// Transpose
inline Matrix_Sp Matrix_Sp::trans() const
{
    Matrix_Sp m;
    *m.p=(*p).transpose();

    return m;
}

//double Matrix_Sp::max() const
//{
//    return (*p).cwiseMax();
//}
//
//double Matrix_Sp::min() const
//{
//    return (*p).cwiseMin();
//}


#endif // MATRIX_SP_H_INCLUDED
