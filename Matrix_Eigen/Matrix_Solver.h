#ifndef MATRIX_SOLVER_H_INCLUDED
#define MATRIX_SOLVER_H_INCLUDED

#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>
#include <typeinfo>

namespace SOLVER
{


Matrix SolveM(Matrix &M,Matrix &p)
{
    return (p/(M|Matrix::ones(p.sizeX(),1)));
}

Matrix SolveM_Sparse(Matrix &M,Matrix &p)
{
    return (p/M);
}

MatrixXd diag(Matrix& A)
{
    return A.P().diagonal();
}

Matrix Dot_S(Matrix& A, Matrix& X,int x,int y)
{
    Matrix result(X.sizeX());
    for(int k=0;k<X.sizeX()/x;k++){
        if(k<y)
            for(int i=0;i<x;i++)
                result(k*x+i)=(A.range(k*x,(k+1)*x-1,k*x,(k+y)*x+2)|(X.range(k*x,(k+y)*x+2,0,0)))(0,0);
        else{
        if(k>X.sizeX()-y-1)
            for(int i=0;i<x;i++)
                result(k*x+i)=(A.range(k*x,(k+1)*x-1,(k-y)*x,(k)*x+2)|(X.range((k-y)*x,(k)*x+2,0,0)))(0,0);
        else
            for(int i=0;i<x;i++)
                result(k*x+i)=(A.range(k*x,(k+1)*x-1,(k-y)*x,(k+y+1)*x-1)|(X.range((k-y)*x,(k+y)*x+2,0,0)))(0,0);
        }
    }
    return result;

}



enum DISP {ON,OFF};
Matrix BiCG_Pre_Acc(Matrix &A,Matrix &b,Matrix &X0, double tol, int max_it,DISP display=ON)
{

    string Title="BiCG_Pre_Acc";
    int flag=0;
    Matrix r,rt,g,gt,p,pt,M,Mt,X,z,zt;
    double a,beta,rho=0,rho_1;

    M=A.P().diagonal();
    X=X0;
    r=b-(A|X0);
    rt=r;

    for (int k=1;k<=max_it;k++)
    {
            z=r/M;
            zt=rt/M;


        if(k==1){
            rho=(z.trans()|rt)(0,0);
            p=z;
            pt=zt;
        }
        else{
            rho_1=rho;
            rho=(z.trans()|rt)(0,0);
            beta=rho/rho_1;
            p=z+beta*p;
            pt=zt+beta*pt;
        }
        g=(A|p);
        gt=A.trans()|(pt.P());
        a=rho/(pt.trans()|g)(0,0);
        X+=a*p;
        r-=a*g;
        rt-=a*gt;

        if(abs(r).max()<tol){
            if (display==ON) cout<<Title<<": Converged after "<<k<<" times iteration. "<<endl;
            flag=1;
            break;

         }

    }
    if(flag==0){
        cout<<Title<<": Doesn't reach tol "<<tol<<"after max iteration: "<<max_it<<". "<<endl;
    }
    if (display==ON)std::cout<<"Max Residual: "<<abs(r).max()<<endl;

    return X;
}
Matrix BiCG_Pre_Acc_Sp(Matrix_Sp &A,Matrix &b, double tol, int max_it,DISP display=ON)
{

    string Title="BiCG_Pre_Acc_Sp";
    int flag=0;
    Eigen::MatrixXd r,rt,g,gt,p,pt,M,Mt,X,z,zt,At;
    double a,beta,rho=0,rho_1,res;

//    M=Eigen::MatrixXd(A.P().diagonal());
    At=A.P().transpose();
    X=Eigen::VectorXd(b.sizeX());
    X.fill(0.);
    r=b.P()-(A.P()*X);
    rt=r;

    for (int k=1;k<=max_it;k++)
    {
            z=r,zt=rt;
            for(int i=0;i<z.rows();++i) {
                    z(i)/=A(i,i);
                    zt(i)/=A(i,i);
            }


        if(k==1){
            rho=(z.transpose()*rt)(0,0);
            p=z;
            pt=zt;
        }
        else{
            rho_1=rho;
            rho=(z.transpose()*rt)(0,0);
            beta=rho/rho_1;
            p=z+beta*p;
            pt=zt+beta*pt;
        }
        g=(A.P()*p);
        gt=At*pt;
        a=rho/(pt.transpose()*g)(0,0);
        X+=a*p;
        r-=a*g;
        rt-=a*gt;
        res=r.cwiseAbs () .maxCoeff();
        if(res<tol){
            if (display==ON) cout<<Title<<": Converged after "<<k<<" times iteration. "<<endl;
            flag=1;
            break;

         }

    }
    if(flag==0){
        cout<<Title<<": Doesn't reach tol "<<tol<<"after max iteration: "<<max_it<<". "<<endl;
    }
    if (display==ON)std::cout<<"Max Residual: "<<res<<endl;

    return X;
}

Matrix BiCG_St(Matrix &A,Matrix &b,Matrix &X0, double tol, int max_it)
{
    string Title="BiCG_St";
    int flag=0;
    double rho,rho_1,a,w=0.,beta;
    Matrix p,pt,v,s,st,t;
    Matrix M=A.diag();
    Matrix X=X0;
    Matrix r0=b-(A|X0);
    Matrix r=r0;
    Matrix rt=r0;
    for (int k=1;k<=max_it;k++)
    {

        if(k==1){
            rho=(rt.trans()|r)(0,0);
            if (rho==0) {
                cout<<"Method fails. "<<endl;
                break;
            }
            p=r;
        }
        else{
            rho_1=rho;
            rho=(rt.trans()|r)(0,0);
            if (rho==0) {
                cout<<"Method fails. "<<endl;
                break;
            }
            beta=rho/rho_1*(a/w);
            p=r+beta*(p-w*v);
        }
        pt=SolveM(M,p);
        v=A|pt;
        a=(rho/(rt.trans()|v))(0,0);
        s=r-a*v;
        st=SolveM(M,s);
        t=A|st;
        w=((t.trans()|s)/(t.trans()|t))(0,0);
        X+=a*pt+w*st;
        r=s-w*t;

        if(abs(r).max()<tol){
            cout<<Title<<": Converged after "<<k<<" times iteration. "<<endl;
            flag=1;
            break;
        }

    }
    if(flag==0){
        cout<<Title<<": Doesn't reach tol "<<tol<<"after max iteratoin: "<<max_it<<". "<<endl;
    }
    std::cout<<"Max Residual: "<<abs(r).max()<<endl;

    return X;
}

//void Solver_test(string folder)
//{
//    Matrix A=Matrix::readcsv(folder+"/A_non_sys.csv",10,10);
//    Matrix b=Matrix::readcsv(folder+"/b_non_sys.csv",10,1);
//
//    Matrix X0=Matrix::zeros(10,1);
//    A.show();
//    b.show();
//    Matrix X=BiCG_St(A,b,X0,1e-5,30);
//
//    X.show();
//
//}
}

#endif // MATRIX_SOLVER_H_INCLUDED
