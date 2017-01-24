#ifndef MATRIX_S_H_INCLUDED
#define MATRIX_S_H_INCLUDED


#include <typeinfo>
///**************   Class Body  *****************///

template<class T>
class Matrix_S
{
private:
    int dx,dy;
    T **p;

    void allocArrays()
    {
        p=new T*[dx];
        for(int i=0;i<dx;++i)
            p[i]=new T[dy];
    }

public:
    Matrix_S(int sizeX, int sizeY);
    Matrix_S(int sizeX);
    Matrix_S();
    ~Matrix_S();
    Matrix_S(const Matrix_S<T>& m);


    int sizeX () const    {return dx;}
    int sizeY () const    {return dy;}
    const T* begin(){return p[0][0];}
    const T* end(){return p[dx-1][dy-1];}

//    T& P() const {return *p;}
    T& operator    ()  (int X,int Y);
    T& operator    ()  (int X);
    T& operator    ()  (Matrix_S<T>&, T c);
    T& operator    ()  (Matrix_S<T>&, int c);
    Matrix_S& operator    =   (const Matrix_S& m);
    Matrix_S& operator    =   (T* m[]);
    Matrix_S& operator    =   (T  m[]);
    Matrix_S<T>& operator    +=  (const Matrix_S<T>& m);
    Matrix_S<T>& operator    -=  (const Matrix_S<T>& m);
    Matrix_S<T>& operator    *=  (const Matrix_S<T>& m);
    Matrix_S<T>& operator    /=  (const Matrix_S<T>& m);
    Matrix_S<T>  operator    +   (const Matrix_S<T>& m);
    Matrix_S<T>  operator    -   (const Matrix_S<T>& m);

    template<class T1> friend Matrix_S<T1> operator | ( Matrix_S<T1>& m1,  Matrix_S<T1>& m2);
     template<class T1> friend Matrix_S<T1*> operator | ( Matrix_S<T1*>& m1,  Matrix_S<T1*>& m2);

    template<class T1> friend Matrix operator | ( Matrix_S<Matrix>& m1,  Matrix& m2);
    template<class T1> friend Matrix_S<T1> operator * (const Matrix_S<T1>& m1, const Matrix_S<T1>& m2);
    template<class T1> friend Matrix_S<T1> operator * (double c,const Matrix_S<T1>& m2);
    template<class T1> friend Matrix_S<T1> operator * (const Matrix_S<T1>& m1,double c);
    template<class T1> friend Matrix_S<T1> operator * (int c,const Matrix_S<T1>& m2);
    template<class T1> friend Matrix_S<T1> operator * (const Matrix_S<T1>& m1,int c);

    template<class T1> friend Matrix_S<T1> operator + (double c,const Matrix_S<T1>& m);
    template<class T1> friend Matrix_S<T1> operator + (const Matrix_S<T1>& m, double c);
    template<class T1> friend Matrix_S<T1> operator + (int c,const Matrix_S<T1>& m);
    template<class T1> friend Matrix_S<T1> operator + (const Matrix_S<T1>& m, int c);
    template<class T1> friend Matrix_S<T1> operator + (const Matrix_S<T1>& m);
    template<class T1,class T2> friend Matrix_S<T1> operator + (const Matrix_S<T1>&, const T2& );

    template<class T1> friend Matrix_S<T1> operator - (double c,const Matrix_S<T1>& m);
    template<class T1> friend Matrix_S<T1> operator - (const Matrix_S<T1>& m, double c);
    template<class T1> friend Matrix_S<T1> operator - (int c,const Matrix_S<T1>& m);
    template<class T1> friend Matrix_S<T1> operator - (const Matrix_S<T1>& m, int c);
    template<class T1> friend Matrix_S<T1> operator - (const Matrix_S<T1>& m);
    template<class T1,class T2> friend Matrix_S<T1> operator - (const Matrix_S<T1>&, const T2& );


    template<class T1> friend Matrix_S<T1> operator / (const Matrix_S<T1>& m1, const Matrix_S<T1>& m2);
    template<class T1> friend Matrix_S<T1> operator / (const Matrix_S<T1>& m1, double c);
    template<class T1> friend Matrix_S<T1> operator / (double c, const Matrix_S<T1>& m1);
    template<class T1> friend Matrix_S<T1> operator / (const Matrix_S<T1>& m1, int c);
    template<class T1> friend Matrix_S<T1> operator / (int c, const Matrix_S<T1>& m1);
    template<class T1> friend Matrix_S<T1> operator ^ (const Matrix_S<T1>& m1, double c);
    template<class T1> friend Matrix_S<T1> operator ^ (const Matrix_S<T1>& m1, int c);

    template<class T1> friend ostream& operator << (ostream& out, const Matrix_S<T1>& m);
    template<class T1> friend istream& operator >> (istream& input, const Matrix_S<T1>& m);

    Matrix_S<T> diag();
    void show(string name="",int n=15,int c=0) const;
    Matrix_S trans() const;
    Matrix_S det() const;
    double det2() const;
    Matrix_S<T> index() const;
    T sum() const;
    double max() const;
    double min() const;
    void remove_nan(double c)
    {
        for(int i=0;i<dx;++i)
            for(int j=0;j<dy;++j)
                p[i][j].remove_nan(c);
    }
    void remove_nan(int c){remove_nan(double(c));}
    void change(const Matrix_S &, double, double);
    void change(const Matrix_S &, double, int);
    void change(const Matrix_S &, int, double);
    void change(const Matrix_S &, int, int);
    Matrix_S row(int i);
    Matrix_S col(int j);
    Matrix_S range(int m,int n);
    void swap_row(int,int);
    void swap_col(int,int);

    static Matrix_S ones(int i, int j);
    static Matrix_S zeros(int i, int j);
    static Matrix_S random(int i, int j);
    static Matrix_S readcsv   (string filename, int col, int row) ;
    static void   writecsv  (string filename, const  Matrix_S &m) ;





};
///**************   Methods defination  *****************///

///DEFULT CONSTRUCTOR
template<class T> Matrix_S<T>::Matrix_S() : dx(1),dy(1) {allocArrays();}
///CONSTRUCTOR

template<class T> Matrix_S<T>::Matrix_S(int sizeX,int sizeY):dx(sizeX),dy(sizeY)
{
    T newT;
    allocArrays();
    for(int i=0;i<dx;++i)
        for(int j=0;j<dy;++j)
            p[i][j]=newT;
}

template<class T> Matrix_S<T>::Matrix_S(int sizeX):dx(sizeX),dy(1)
{
    T newT;
    allocArrays();
    for(int i=0;i<dx;++i)
       p[i][0]=newT;
}
///COPY CONSTRUCTOR
template<class T> Matrix_S<T>::Matrix_S(const Matrix_S<T>& m):dx(m.dx),dy(m.dy)
{
    allocArrays();
    for(int i=0;i<dx;++i)
        for(int j=0;j<dy;++j)
            p[i][j]=m.p[i][j];
}

///DISTRUCTOR
template<class T> Matrix_S<T>::~Matrix_S()
{
    for(int i=0;i<dx;++i)
        delete[] p[i];///might be a problem w/o space

    delete[] p;
}

///()
template<class T> T& Matrix_S<T>::operator()(int X,int Y)
{
    return p[X][Y];
}
template<class T> T& Matrix_S<T>::operator()(int X)
{
    if(dx==1) return p[0][X];
    else if(dy==1) return p[X][0];
    else {cout<<"Not a vector. "<<endl; return p[0][X];}
}

///=
template<class T> Matrix_S<T>& Matrix_S<T>::operator=(const Matrix_S<T>& m)
{
    ///judge
    if(this==&m)
    {
        return *this;
    }
    else
    {
        if((dx!=m.dx||dy!=m.dy))
        {
            if (dx!=1&&dy!=1) cout<<"Warning(=): Dimonsion dosen't agree. New matrix will be built."<<endl;
            this->~Matrix_S();
            dx=m.dx;
            dy=m.dy;
            allocArrays();
        }

        for(int i=0;i<dx;++i)
            for(int j=0;j<dy;++j)
                p[i][j]=m.p[i][j];
    }
    return *this;
}

template<class T> Matrix_S<T>& Matrix_S<T>::operator=(T* m[])
{
        for(int i=0;i<dx;++i)
            for(int j=0;j<dy;++j)
                p[i][j]=*m[i*dy+j];

    return *this;
}

template<class T> Matrix_S<T>& Matrix_S<T>::operator=(T m[])
{
        for(int i=0;i<dx;++i)
            for(int j=0;j<dy;++j)
                p[i][j]=m[i*dy+j];

    return *this;
}

///-=
template<class T> Matrix_S<T>& Matrix_S<T>::operator-=(const Matrix_S<T>& m)
{
    ///judge

        if(dx!=m.dx||dy!=m.dy)
            cout<<"Error(-=): dimonsion not agree."<<endl;
        else{
            for(int i=0;i<dx;++i)
                for(int j=0;j<dy;++j)
                    p[i][j]-=m.p[i][j];
        }

    return *this;
}

/// *=
template<class T> Matrix_S<T>&  Matrix_S<T>::operator*=(const Matrix_S<T>& m)
{
    ///judge

        if(dx!=m.dx||dy!=m.dy)
            cout<<"Error(*=): dimonsion not agree."<<endl;
        else{
            for(int i=0;i<dx;++i)
                for(int j=0;j<dy;++j)
                    p[i][j]*=m.p[i][j];
        }

    return *this;
}
/// /=
template<class T> Matrix_S<T>& Matrix_S<T>::operator/=(const Matrix_S<T>& m)
{
    ///judge

        if(dx!=m.dx||dy!=m.dy)
            cout<<"Error(/=): dimonsion not agree."<<endl;
        else{
            for(int i=0;i<dx;++i)
                for(int j=0;j<dy;++j)
                    p[i][j]/=m.p[i][j];
        }

    return *this;
}
///+
template<class T> Matrix_S<T> Matrix_S<T>::operator+(const Matrix_S<T>& m)
{
    ///judge
    Matrix_S temp(*this);
        if(dx!=m.dx||dy!=m.dy)
            cout<<"Error(+): dimonsion not agree."<<endl;
        else{
            for(int i=0;i<dx;++i)
                for(int j=0;j<dy;++j)
                    temp.p[i][j]+=m.p[i][j];
        }

    return temp;
}
///-
template<class T> Matrix_S<T> Matrix_S<T>::operator-(const Matrix_S<T>& m)
{
    ///judge
    Matrix_S temp(*this);
        if(dx!=m.dx||dy!=m.dy)
            cout<<"Error(-): dimonsion not agree."<<endl;
        else
        {
            for(int i=0;i<dx;++i)
                for(int j=0;j<dy;++j)
                    temp.p[i][j]=temp.p[i][j]-m.p[i][j];
        }

    return temp;
}

///DISPLAY
template<class T> void Matrix_S<T>::show(string name,int n,int c) const
{
    cout.precision(n);
    cout<<"Matrix_S: "<<name<<endl;
    for(int i=0;i<dy;++i)
    {
        for(int j=0;j<dx;++j){
            if (typeid(p[j][i])==typeid(Matrix))
            {
                cout<<"Row["<<j<<"]Col["<<i<<"]: "<<endl;
                cout<<p[j][i]<<endl;
            }
            else
                cout<<p[j][i]<<"    ";

            if(c==1){int test; cin>>test;}
        }
        cout<<endl;
    }
    cout<<endl;
}


Matrix_S<Matrix> combine(Matrix* m[],int y,int x)
{
    Matrix_S<Matrix> result(y,x);
    for(int i=0;i<y;++i)
        for(int j=0;j<x;++j)
            result(i,j)=*m[i+j];

    return result;
}
template<class T> T Matrix_S<T>::sum() const
{
    T result=0.*p[0][0];
    for(int i=0;i<dx;++i)
        for(int j=0;j<dy;++j)
            result+=p[i][j];
    return result;
}

template<class T> Matrix_S<T> Matrix_S<T>::index() const
{
    Matrix_S<T> result(dx,dy);
    for(int i=0;i<dx;++i)
        for(int j=0;j<dy;++j)
            result(i,j) = p[i][j].index();
    return result;
}
  ///******************** Friends *************************///
    template<class T1>   Matrix_S<T1> operator | ( Matrix_S<T1>& m1, Matrix_S<T1>& m2)
    {
    Matrix_S<T1> m(m1.dx,m2.dy);
    if(m1.dy!=m2.dx) cout<<"Error(|): OUTTER Dimension doesn't agree. "<<endl;
    else
    {

        if(typeid(T1)==typeid(Matrix))
        {
            for(int i;i<m.dx;++i)
                for(int j=0;j<m.dy;++j)
                    m.p[i][j]=Matrix(m2(i,j).sizeX(),m2(i,j).sizeY());


//        ///#pragma omp parallel for
        for(int i=0;i<m1.dx;++i)
            for(int j=0;j<m2.dy;++j)
                for(int k=0;k<m1.dy;k++)
                {
                    if(m1.p[i][k].sizeX()==1&&m1.p[i][j].sizeY()==1&&m1.p[i][j].sizeZ==1)
                    {
                        double X=m1.p[i][j](0,0,0); /// or should use int instead?
                        if(X==0) continue;
                        if(X==1) m.p[i][j]+=m2.p[k][j];
                        else     m.p[i][j]+=(X*m2.p[k][j]);
                    }
                    else m.p[i][j]+=(m1.p[i][k]*m2.p[k][j]);
                }

        }


        return m;
    }
    }
    template<class T1>   Matrix_S<T1*> operator | ( Matrix_S<T1*>& m1, Matrix_S<T1*>& m2)
    {
        Matrix_S<T1*> m(m1.dx,m2.dy);
        if(m1.dy!=m2.dx) cout<<"Error(|*): OUTTER Dimension doesn't agree. "<<endl;
        else
        {
    //        ///#pragma omp parallel for
            for(int i=0;i<m1.dx;++i)
                for(int j=0;j<m2.dy;++j)
                    for(int k=0;k<m1.dy;k++)
                    {
                        if(typeid(T1)==typeid(Matrix3D))
                        {
                            if(m1.p[i][k]->sizeX()==1&&m1.p[i][j]->sizeY()==1&&m1.p[i][j]->sizeZ==1)
                            {
                                double X=(*m1.p[i][j])(0,0,0); /// or should use int instead?
                                if(X==0) continue;
                                if(X==1) *m.p[i][j]+=(*m2.p[k][j]);
                                else     *m.p[i][j]+=X*(*m2.p[k][j]);
                            }
                        }
                        else *m.p[i][j]+=(*m1.p[i][k])*(*m2.p[k][j]);
                    }

            }


            return m;
    }


    ///m1*m2
     template<class T1>  Matrix_S<T1> operator * (const Matrix_S<T1>& m1,const Matrix_S<T1>& m2)
    {
        Matrix_S<T1> prod(m1.dx,m2.dy);
        if(m1.dx!=m2.dx||m1.dy!=m2.dy)
            cout<<"Error(*): dimonsion does not agree."<<endl;
        else
        {
            for(int i=0;i<m1.dx;++i)
                for(int j=0;j<m2.dy;++j)
                    prod.p[i][j]=m1.p[i][j]*m2.p[i][j];
        }
        return prod;
    }
    template<class T1>   Matrix_S<T1> operator * (double c,const Matrix_S<T1>& m2)
    {
        Matrix_S<T1> prod(m2);
        for(int i=0;i<m2.dx;++i)
            for(int j=0;j<m2.dy;++j)
                prod.p[i][j]=c*m2.p[i][j];

        return prod;
    }
    template<class T1>   Matrix_S<T1> operator * (const Matrix_S<T1>& m1,double c)
    {
        return c*m1;
    }

    template<class T1>   Matrix_S<T1> operator * (int c1,const Matrix_S<T1>& m2)
    {
        double c=double(c1);
        Matrix_S<T1> prod(m2);
        for(int i=0;i<m2.dx;++i)
            for(int j=0;j<m2.dy;++j)
                prod.p[i][j]=c*m2.p[i][j];

        return prod;
    }
    template<class T1>   Matrix_S<T1> operator * (const Matrix_S<T1>& m1,int c1)
    {
        double c=double(c1);
        return c*m1;
    }


    ///c+m1
    template<class T1>   Matrix_S<T1> operator + (double c,const Matrix_S<T1>& m)
    {
        Matrix_S<T1> plus_M(m.dx,m.dy);

        for(int i=0;i<m.dx;++i)
            for(int j=0;j<m.dy;++j)
                plus_M.p[i][j]=m.p[i][j]+c;

        return plus_M;
    }
    template<class T1>   Matrix_S<T1> operator + (const Matrix_S<T1>& m, double c)
    {
       return (c+m);
    }

    template<class T1>   Matrix_S<T1> operator + (int c1,const Matrix_S<T1>& m)
    {
        double c=double(c1);
        Matrix_S<T1> plus_M(m.dx,m.dy);

        for(int i=0;i<m.dx;++i)
            for(int j=0;j<m.dy;++j)
                plus_M.p[i][j]=m.p[i][j]+c;

        return plus_M;
    }
    template<class T1>   Matrix_S<T1> operator + (const Matrix_S<T1>& m, int c1)
    {
        double c=double(c1);
        return (c+m);
    }

    template<class T1>   Matrix_S<T1> operator + (const Matrix_S<T1>& m)
    {
        return m;
    }
    template<class T1,class T2>
      Matrix_S<T1> operator + (const Matrix_S<T1>& m1,const T2& m2)
    {
        Matrix_S<T1> plus_M=m1;

        for(int i=0;i<m1.dx;++i)
            for(int j=0;j<m1.dy;++j)
                plus_M.p[i][j]+=m2;

        return plus_M;
    }
    ///c-m1
     template<class T1,class T2>
      Matrix_S<T1> operator - (const Matrix_S<T1>& m1,const T2& m2)
    {
        Matrix_S<T1> plus_M=m1;

        for(int i=0;i<m1.dx;++i)
            for(int j=0;j<m1.dy;++j)
                plus_M.p[i][j]-=m2;

        return plus_M;
    }
    template<class T1>   Matrix_S<T1> operator - ( double c,const Matrix_S<T1>& m)
    {
        Matrix_S<T1> minus_M(m.dx,m.dy);

        for(int i=0;i<m.dx;++i)
            for(int j=0;j<m.dy;++j)
                minus_M.p[i][j]=c-m.p[i][j];

        return minus_M;
    }
    template<class T1>   Matrix_S<T1> operator - (const Matrix_S<T1>& m,  double c)
    {
       return (0-(c-m));
    }

    template<class T1>   Matrix_S<T1> operator - ( int c1,const Matrix_S<T1>& m)
    {
        double c=double(c1);
        Matrix_S<T1> minus_M(m.dx,m.dy);

        for(int i=0;i<m.dx;++i)
            for(int j=0;j<m.dy;++j)
                minus_M.p[i][j]=c-m.p[i][j];

        return minus_M;
    }
    template<class T1>   Matrix_S<T1> operator - (const Matrix_S<T1>& m,  int c1)
    {
        double c=double(c1);
        return (0-(c-m));
    }

    template<class T1>   Matrix_S<T1> operator - (const Matrix_S<T1>& m)
    {
        return (0-m);
    }
    ///  /
    template<class T1>   Matrix_S<T1> operator / (const Matrix_S<T1>& m1,const Matrix_S<T1>& m2)
    {
        Matrix_S<T1> prod(m1.dx,m1.dy);
        if(m1.dx!=m2.dx||m1.dy!=m2.dy)
            cout<<"Error(/): dimonsion does not agree."<<endl;
        else
        {
            for(int i=0;i<m1.dx;++i)
                for(int j=0;j<m1.dy;++j)
                    prod.p[i][j]=m1.p[i][j]/m2.p[i][j];
        }
        return prod;
    }

    template<class T1>   Matrix_S<T1> operator / (const Matrix_S<T1>& m1, double c)
    {
        Matrix_S<T1> prod(m1.dx,m1.dy);

            for(int i=0;i<m1.dx;++i)
                for(int j=0;j<m1.dy;++j)
                    prod.p[i][j]=m1.p[i][j]/c;

        return prod;
    }

    template<class T1>   Matrix_S<T1> operator / (double c,const Matrix_S<T1>& m1)
    {
        Matrix_S<T1> prod(m1.dx,m1.dy);

            for(int i=0;i<m1.dx;++i)
                for(int j=0;j<m1.dy;++j)
                    prod.p[i][j]=c/m1.p[i][j];

        return prod;
    }

    template<class T1>   Matrix_S<T1> operator / (const Matrix_S<T1>& m1, int c1)
    {
        double c=double(c1);
        Matrix_S<T1> prod(m1.dx,m1.dy);

            for(int i=0;i<m1.dx;++i)
                for(int j=0;j<m1.dy;++j)
                    prod.p[i][j]=m1.p[i][j]/c;

        return prod;
    }

    template<class T1>   Matrix_S<T1> operator / (int c1,const Matrix_S<T1>& m1)
    {
        double c=double(c1);
        Matrix_S<T1> prod(m1.dx,m1.dy);

            for(int i=0;i<m1.dx;++i)
                for(int j=0;j<m1.dy;++j)
                    prod.p[i][j]=c/m1.p[i][j];

        return prod;
    }

    /// ^
    template<class T1>   Matrix_S<T1> operator ^ (const Matrix_S<T1>& m1, double c)
    {
        Matrix_S<T1> prod(m1.dx,m1.dy);
        if(c<0)
        {
            for(int i=0;i<m1.dx;++i)
                for(int j=0;j<m1.dy;++j)
                    prod.p[i][j]=1.0/pow(m1.p[i][j],-c);
        }
        else
        {
            for(int i=0;i<m1.dx;++i)
                for(int j=0;j<m1.dy;++j)
                    prod.p[i][j]=pow(m1.p[i][j],c);
        }
        return prod;
    }

    template<class T1>   Matrix_S<T1> operator ^ (const Matrix_S<T1>& m1, int c1)
    {
        double c=double(c1);
        return m1^c;
    }

    /// << has problem
    template<class T1>   ostream &operator << (ostream &out, const Matrix_S<T1> &m)
    {
        for (int i=0;i<m.dx;++i)
        {
            for (int j=0;j<m.dy;++j)
                if (typeid(m.p[i][j])==typeid(Matrix)||typeid(m.p[i][j])==typeid(Matrix_S<T1>))
                    out<<typeid(m.p[i][j]).name()<<"["<<m.p[i][j].sizeX()<<"x"<<m.p[i][j].sizeY()<<"]     ";
                else
                    out<<m.p[i][j]<<"   ";
            out<<endl;
        }
        return out;
    }
    ///
    template<class T1>   istream& operator >> (istream& input, const Matrix_S<T1>& m)
        {
            for(int i=0;i<m.dx;++i)
                  for(int j=0;j<m.dy;++j)
                           input>>m.p[i][j];
            return input;
         }

#endif // MATRIX_S_H_INCLUDED
