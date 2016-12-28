#include <iostream>
#include <utility>
#include <vector>
#include <math.h>
#include <complex>

using namespace std;

/* ---------------------------------------------------- */
/*                        MACRO                         */
/* ---------------------------------------------------- */

#define PRECISION 1e-10
#define OUTPUT_PRECISION 1e-10
#define NOROOT "No Root.\n\n"
#define SINGLEROOT "The only Root is   x  =  "
#define PI 3.1415926535897932384

/* ---------------------------------------------------- */
/*                   GLOBAL VARIABLES                   */
/* ---------------------------------------------------- */

int deg;  // degree
vector<double> a(0);    // coefficients
vector<double> b(0);    // coefficients for the 1st-order derivative
vector<double> tmp_poly(0); // coefficients for a temporary polynomial
vector<complex<double>> roots(0);
vector<complex<double>> delta_roots(0);

/* ---------------------------------------------------- */
/*                 FUNCTION DECLARATION                 */
/* ---------------------------------------------------- */

inline void setup();

inline void first_deriv();
inline void calcGCD();
inline void square_free();

inline pair<double, double> bounds();
inline void initial_approx();

inline void adjust();

inline void print();

/* ---------------------------------------------------- */
/*                         MAIN                         */
/* ---------------------------------------------------- */

int main() {
    setup();
    square_free();
    
    initial_approx();
    adjust();
    print();
}

/* ---------------------------------------------------- */
/*               FUNCTION IMPLEMENTATION                */
/* ---------------------------------------------------- */

/* ---------------------- Setup ----------------------- */
inline void setup() {
    cout << "Degree ( >= 1 ) of the polynomial :    ";
    cin >> deg;
    if (deg < 1) {
        cout << "The degree must be positive.\n";
        exit(1);
    }
    cout << "Input the coefficients for  ";
    for (int i = deg; i >= 2; i--)
        cout << "a[" << i << "] x^" << i << "  +  ";
    cout << "a[1] x  +  a[0]  =  0" << endl << endl;
    a.resize(deg + 1);
    cout << "a[" << deg << "] = ";
    double highest;
    cin >> highest;
    if (!highest) {
        cout << "The highest coefficient must not equal 0.\n";
        exit(2);
    }
    a[deg] = 1;
    for (int i = deg - 1; i >= 0; i--) {
        double tmp;
        cout << "a[" << i << "] = ";
        cin >> tmp;
        a[i] = tmp / highest;
    }
    cout << endl;
    if (deg == 1) {
        double x = -a[0] / a[1];
        cout << SINGLEROOT << x << " .\n\n";
        exit(0);
    }
}

/* ------------------- Square Free -------------------- */

// 1st-order derivative
inline void first_deriv() {
    int n = (int) a.size();
    (n < 2) ? b.resize(1, 0) : b.resize(n - 1, 0);
    for (int i = 1; i < n; i++)
        b[i - 1] = i * a[i];
}

// Greatest Common Divisor of b and tmp_poly
// Finally b = 0 and tmp_poly = GCD
inline void calcGCD() {
    // assume the leading coefficients == 1
    while (b.size() != 1) {
        int b_deg = (int) b.size() - 1, tmp_deg = (int) tmp_poly.size() - 1;
        int diff = tmp_deg - b_deg;
        for (int i = b_deg; i >= 0; i--)
            tmp_poly[i + diff] -= b[i];
        int highest_deg = tmp_deg - 1;
        while (fabs(tmp_poly[highest_deg]) <= PRECISION && highest_deg > 0)
            highest_deg--;
        tmp_poly.resize(highest_deg + 1);
        double leading_coeff = tmp_poly[highest_deg];
        if (fabs(leading_coeff) > PRECISION) {
            tmp_poly[highest_deg] = 1;
            for (int i = highest_deg - 1; i >= 0; i--)
                tmp_poly[i] /= leading_coeff;
        }
        if (highest_deg < b_deg)
            b.swap(tmp_poly);
    }
    if (fabs(b[0]) > PRECISION)
        tmp_poly.resize(1, 1);
}

inline void square_free() {
    first_deriv();
    int b_deg = (int) b.size() - 1;
    double leading_coeff = b[b_deg];
    if (leading_coeff != 1) {
        b[b_deg] = 1;
        for (int i = b_deg - 1; i >= 0; i--)
            b[i] /= leading_coeff;
    }
    tmp_poly.assign(a.begin(), a.end());    // keep deg(tmp_poly) >= deg(b)
    calcGCD();
    int tmp_deg = (int) tmp_poly.size() - 1;
    if (!tmp_deg)
        return;
    // store a / tmp_poly in b
    b.resize(deg - tmp_deg + 1, 0);
    for (int i = deg - tmp_deg; i >= 0; i--) {
        double tmp = a[i + tmp_deg];
        b[i] = tmp;
        for (int j = tmp_deg - 1; j >= 0; j--)
            a[i + j] -= tmp * tmp_poly[j];
    }
    a.assign(b.begin(), b.end());
    deg = (int) a.size() - 1;
}

/* ------------------ Initial Approx ------------------ */

// Lower and Upper Bounds of absolute values of all roots
inline pair<double, double> bounds() {
    // calc max{|a[1]|, |a[2]|, ..., |a[deg-1]|}
    double max_abs = 0.0;
    for (int i = 1; i < deg; i++) {
        double tmp = fabs(a[i]);
        if (tmp > max_abs)
            max_abs = tmp;
    }
    double a0 = a[0];
    double lower_bound = a0 / (a0 + max(max_abs, a[deg]));
    double upper_bound = max(max(1 + fabs(a[deg - 1]), fabs(a0)), max_abs);
    return {lower_bound, upper_bound};
}

inline void initial_approx() {
    pair<double, double> tmp = bounds();
    double rho = (tmp.first + tmp.second) / 2.0;
    double theta = 2 * PI / deg;
    for (int i = 0; i <= deg - 2; i++)
        roots.push_back(polar(rho, theta * i));
    roots.push_back(polar(rho, -theta / 2.0));
}

/* ----------------------- Solve ---------------------- */

inline void adjust() {
    delta_roots.resize(deg, 1);
    while (1) {
        bool done = true;
        for (int i = 0; i < deg; i++) {
            first_deriv();
            complex<double> orig = {0, 0}, first_de = {0, 0}; //
            complex<double> z = roots[i];
            for (int j = deg; j >= 0; j--)
                orig += a[j] * pow(z, j);
            for (int j = deg - 1; j >= 0; j--)
                first_de += b[j] * pow(z, j);
            complex<double> frac = first_de / orig;
            for (int j = 0; j < deg; j++) {
                if (j == i)
                    continue;
                complex<double> tmp = z - roots[j];
                if (abs(tmp) < PRECISION)
                    tmp = {PRECISION, 0};
                frac -= pow(tmp, -1);
            }
            if (abs(frac) < PRECISION)
                frac = {PRECISION, 0};
            delta_roots[i] = pow(frac, -1);
            if (abs(delta_roots[i]) > PRECISION)
                done = false;
        }
        if (done)
            break;
        for (int i = 0; i < deg; i++)
            roots[i] -= delta_roots[i];
    }
}

/* ------------------- Print Result ------------------- */

inline void print() {
    int n = (int) roots.size();
    if (n == 1)
        cout << SINGLEROOT << roots[0].real() << " .\n\n";
    else {
        cout << "The " << n << " distinct Roots are\n";
        for (int i = 0; i < n; i++) {
            cout << "   x[" << i + 1 << "] = ";
            
            double re = roots[i].real(), im = roots[i].imag();
            if ( fabs(re) < OUTPUT_PRECISION )
                re = 0;
            if ( fabs(im) < OUTPUT_PRECISION )
                im = 0;
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

            cout << endl;
        }
        cout << endl;
    }
}
