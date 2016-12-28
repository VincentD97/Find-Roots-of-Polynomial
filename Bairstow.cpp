//////////////////////////////////////////////////////////////////
//                                                              //
// To Improve :                                                 //
//   Now the program doesn't perform well when the polynomial   //
// contains several repeated roots.                             //
//                                                              //
//   e.g. (x-2)^3(x+1)^2 or (x+1)^8                             //
//                                                              //
//   Since for x in (root-∆, root+∆), f(x) is too close to 0,   //
// the program cannot precisely extract quadratic equations.    //
//                                                              //
// Reference :                                                  //
// http:// web2.uwindsor.ca/courses/engineering/ahmed/PDF%20F   //
// ILES%20FOR%20C++/Chapter%20IX%20C++.pdf                      //
//                                                              //
//////////////////////////////////////////////////////////////////

#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <complex>

using namespace std;

#define MAXDEPTH 500
#define MAXDEGREE 20
#define EPS 1e-19
#define OUTPUTEPS 1e-15

//
// Extract individual real or complex roots from list of quadratic factors
//
int roots(double *x, int n, complex<double> *root)
{
    // x = [-b +- sqrt(b^2 - 4c)] / 2
    //   = halfb +- sqrt(halfb^2 - c)    halfb = -b/2
    double halfb, c;
    int m, numroots;
    complex<double> delta;

    m = n;
    numroots = 0;

    while ( m > 1 )
    {
        halfb = -0.5 * x[m-1];
        c = x[m-2];
        delta = halfb * halfb - c;
        root[m-2] = halfb + sqrt(delta);
        root[m-1] = halfb - sqrt(delta);
        numroots += 2;
        m -= 2;
    }
    if (m == 1) {
        root[0] = -x[0];
        numroots++;
    }
    return numroots;
}


//
// Find quadratic factor using Bairstow's method
//
void find_Quad(double *a, int n, double *b, double *quad, double *err, int *iter)
{
    double *f, u, v, du = 0.0, dv = 0.0, multiplier;
    
    f = new double [n+1];
    f[0] = 1.0;
    u = 4;
    v = -3;
    *iter = 0;
    do
    {

        if (*iter > MAXDEPTH)
            break;
        
        b[n] = 1; b[n-1] = a[n-1] - u;
        for (int i = n - 2; i >= 0; i--)
            b[i] = a[i] - b[i+1] * u - b[i+2] * v;
        f[n] = 1; f[n-1] = b[n-1] - u;
        for (int i = n - 2; i >= 1; i--)
            f[i] = b[i] - u * f[i+1] - v * f[i+2];
        multiplier =  1 / ( f[1] * f[3] - f[2] * f[2] );
        du = ( b[0] * f[3] - f[2] * b[1] ) * multiplier;
        dv = ( f[1] * b[1] - b[0] * f[2] ) * multiplier;
        u += du;
        v += dv;
        
        (*iter)++;
    } while ((fabs(du)+fabs(dv)) > EPS);
    quad[1] = round(u * 100000) / 100000;
    quad[0] = round(v * 100000) / 100000;

    *err = fabs(du)+fabs(dv);
    delete [] f;
    b[n] = 0; b[n-1] = 0;
    for ( int i = n - 2; i >= 0; i-- )
        b[i] = a[i+2] - u * b[i+1] - v * b[i+2];
}

//
// extract all quadratic equations and store in x in pairs ( coeff of 2nd-order is 1 )
// last one may be a linear equation.
//
void get_All_Quads(double *a, int n, double *quad, double *x)
{
    double *quotient, *dividend, err;
    int iter, i, m;
    
    if ( a[n] != 1.0)
    {
        for ( i = n - 1; i >= 0 ;i-- )
            a[i] /= a[n];
        a[n] = 1.0;
    }
    if (n == 2)
    {
        x[1] = a[1];
        x[0] = a[0];
        return;
    }
    else if (n == 1)
    {
        x[0] = a[0];
        return;
    }
    m = n;
    quotient = new double [n+1];
    dividend = new double [n+1];
    quotient[0] = 1.0;
    for ( i = n; i >= 0 ; i-- )
    {
        dividend[i] = a[i];
        x[i] = 0.0;
    }
    do {
        find_Quad(dividend,m,quotient,quad,&err,&iter);
        x[m - 1] = quad[1];
        x[m - 2] = quad[0];
        m -= 2;
        for ( i = 0; i <= m; i++ )
            dividend[i] = quotient[i];
    } while (m > 2);

    if (m == 2)
    {
        x[1] = quotient[1];
        x[0] = quotient[0];
    }
    else x[0] = quotient[0];
    delete [] dividend;
    delete [] quotient;
}



complex<double> letstry(double *a, int n, complex<double> x)
{
    complex<double> result = 0;
    for ( int i = n; i >= 0; i-- )
    {
        result += a[i] * pow(x,i);
    }
    return result;
}

void printComplexNumber(complex<double> x, bool calledByLetstry)
{
    double re(x.real()), im(x.imag());
    if ( fabs(re) < OUTPUTEPS ) re = 0;
    if ( fabs(im) < OUTPUTEPS ) im = 0;
    
    cout.width(10);
    cout << right;
    if ( re == 0 )
        cout << "";
    else
        cout << re;
    
    if ( im > 0 )
    {
        if ( re == 0 ) cout << "       ";
        else cout << "   +   ";
        cout.width(10);
        if ( im == 1 ) cout << right << "" << " i";
        else cout << right << im << " i";
    }
    else if ( im < 0 )
    {
        if ( re == 0 ) cout << "       ";
        else cout << "   -   ";
        cout.width(10);
        if ( im == -1 ) cout << right << "" << " i";
        else cout << right << -im << " i";
    }
    if ( re == 0 && im == 0 && calledByLetstry ) cout << 0;
    cout << endl;
}

int main()
{
    int n,i,numr;
    
    cout << "Polynomial degree (1 <= n <= " << MAXDEGREE << "): ";
    cin >> n;
    if ( n < 1 || n > MAXDEGREE )
    {
        cout << "Error! Invalid degree: n = " << n << endl;
        return 1;
    }
    // get coefficients of polynomial
    double a[n+1], x[n+1], quad[2];
    complex<double> root[MAXDEGREE + 1];
    cout << "Enter coefficients, from high order to low order" << endl;
    for ( i = n; i >= 0 ; i-- )
    {
        cout << "C[" << i << "] * x^" << i << " : ";
        cin >> a[i];
        if (a[n] == 0) {
            cout << "Error! Highest coefficient cannot be 0." << endl;
            return 2;
        }
    }
    if (a[0] == 0) {
        cout << "Error! Constant term) cannot be 0." << endl;
        return 3;
    }
    
    // get roots
    get_All_Quads(a,n,quad,x);
    numr = roots(x,n,root);
    
    cout << endl << "Roots (" << numr << " found):" << endl;
    for ( i = 0; i < n; i++ )
        if ( root[i] != 0.0 )
            printComplexNumber(root[i], false);
    
    cout << "\nNow you can evaluate  f(x)  for the x you provide :       ( 0 + 0 i for quit )\n\n";
    for (;;)
    {
        double re, im;
        cout << " real(x)    =  "; cin >> re;
        cout << " imag(x)    =  "; cin >> im;
        if ( re == 0 && im == 0 )
        {
            cout << "\n\n   Thanks for using PolynomialRootFinder !\n\n";
            break;
        }
        complex<double> x(re,im);
        printComplexNumber(letstry(a, n, x), true);
        cout << endl;
    }
    
    return 0;
}
