#include <iostream>
#include <cmath>
#include <vector>
using namespace std;

#define PRECISION 0.000001
#define MAXTIMES 20
#define NOROOT "No Root in such range.\n\n"
#define SINGLEROOT "The only Root in such range is   x  =  "

int n;  // degree
double range_min, range_max;
vector<double> a(0);    // coefficients
vector<double> realroots(0);
vector<double> SP(0);   // Stationary Points
vector<double> com(0);  // polynomial for complex roots only

inline void setup();
void findRealRoots(vector<double>& coeff);
inline void getRemainingPoly();


int main() {
    setup();
    
    if (n == 1) {
        int x = -a[0] / a[1];
        if (x >= range_min && x <= range_max)
            cout << SINGLEROOT << x << " .\n\n";
        else
            cout << NOROOT;
        return 0;
    }
    
    findRealRoots(a);
    
    if (realroots.size() < n) {
        getRemainingPoly();
    }
    
    
    int num = (int) realroots.size();
    if (!num)
        cout << NOROOT;
    else if (num == 1)
        cout << SINGLEROOT << realroots[0] << " .\n\n";
    else {
        cout << "The " << num << " Roots in such range are\n";
        for (int i = 0; i < num; i++)
            cout << "   x[" << i + 1 << "] = " << realroots[i] << endl;
        cout << endl;
    }
    
    return 0;
}

inline void setup() {
    cout << "Degree ( >= 1 ) of the polynomial :    ";
    cin >> n;
    if (n < 1) {
        cout << "The degree must be positive.\n";
        exit(1);
    }
    cout << "Input the coefficients for  ";
    for (int i = n; i >= 2; i--)
        cout << "a[" << i << "] x^" << i << "  +  ";
    cout << "a[1] x  +  a[0]  =  0" << endl << endl;
    a.resize(n + 1);
    cout << "a[" << n << "] = ";
    cin >> a[n];
    if (!a[n]) {
        cout << "The highest coefficient must not equal 0.\n";
        exit(2);
    }
    for (int i = n - 1; i >= 0; i--) {
        cout << "a[" << i << "] = ";
        cin >> a[i];
    }
    
    cout << "\nFind roots in the range :  MIN  =  ";
    cin >> range_min;
    cout << "and MAX  =  ";
    cin >> range_max;
    if (range_max < range_min) {
        cout << "The range is invalid.\n";
        exit(3);
    }
    if (range_max == range_min) {
        double result = 0;
        for (int i = n; i >= 0 ; i--)
            result += a[i] * pow(range_min, i);
        (result) ? cout << endl << NOROOT : cout << endl << SINGLEROOT << range_min << " .\n\n";
        exit(0);
    }
    cout << endl;
}

void findRealRoots(vector<double>& coeff) {
    int degree = (int) coeff.size() - 1;
    
    /* If degree == 2 -- Start */
    if (degree == 2) {
        realroots.resize(0);
        double a = coeff[2], b = coeff[1], c = coeff[0];
        double delta = b * b - 4 * a * c;
        if (delta == 0) {
            double x = -b / (2 * a);
            if (x >= range_min && x <= range_max)
                realroots.push_back(x);
        }
        else if (delta > 0) {
            double denom = 2 * a, sqroot = sqrt(delta);   // denominator
            double x = (-b - sqroot) / denom;
            if (x >= range_min && x <= range_max)
                realroots.push_back(x);
            x = (-b + sqroot) / (2 * a);
            if (x >= range_min && x <= range_max)
                realroots.push_back(x);
        }
        return;
    }
    /* If degree == 2 -- End */
    
    /* Find the derivative polynomial -- Start */
    vector<double> coeff_tmp(0);
    for (int i = 1; i <= degree; i++)
        coeff_tmp.push_back((i) * coeff[i]);
    /* Find the derivative polynomial -- End */
    
    findRealRoots(coeff_tmp);
    
    /* Find all Stationary Points -- Start */
    SP.resize(0);
    SP.push_back(range_min);
    int SPnum = (int) realroots.size();    // Stationary Point Number
    if (SPnum) {
        if (realroots[0] > range_min)
            SP.push_back(realroots[0]);
        if (SPnum > 1)
            for (int i = 1; i < SPnum; i++)
                SP.push_back(realroots[i]);
    }
    if (SP.back() < range_max)
        SP.push_back(range_max);
    SPnum = (int) SP.size();
    realroots.resize(0);
    /* Find all Stationary Points -- End */
    
    /* Generate roots -- Start */
    /* Find f(Stationary Point Numbers) */
    vector<double> f(0);
    for (int i = 0; i < SPnum; i++) {
        double result = 0;
        double x = SP[i];
        for (int i = degree; i >= 0 ; i--)
            result += coeff[i] * pow(x, i);
        f.push_back(result);
    }
    /* Find at most 1 root between every 2 successive Stationary Points ( [a, b) ) */
    double left = range_min, right = SP[1];
    for (int i = 0; i < SPnum - 1; i++) {
        double f_left = f[i], f_right = f[i + 1];
        if (f_left * f_right >= 0) {
            if (!f_left)
                realroots.push_back(left);
            if (i < SPnum - 2) {
                left = right;
                right = SP[i + 2];
            }
            continue;
        }
        double left_tmp = left, right_tmp = right;
        bool done = false;
        while (right_tmp - left_tmp > PRECISION) {
            double x = (left_tmp + right_tmp) / 2;
            double f_mid = 0;
            for (int i = degree; i >= 0 ; i--)
                f_mid += coeff[i] * pow(x, i);
            if (fabs(f_mid) < PRECISION) {
                done = true;
                realroots.push_back(x);
                break;
            }
            f_mid * f_left > 0 ? left_tmp = x : right_tmp = x;
        }
        if (!done)
            realroots.push_back(left_tmp);
        if (i < SPnum - 2) {
            left = right;
            right = SP[i + 2];
        }
    }
    if (f[SPnum - 1] < PRECISION && realroots.back() < range_max)
        realroots.push_back(range_max);
    /* Generate roots -- End */
}

inline void getRemainingPoly() {
    vector<double> divisor(1,1);
    int RRnum = (int) realroots.size();   // Real Roots Number
    for (int i = 0; i < RRnum; i++) {
        divisor.push_back(1);
        double root = realroots[i];
        for (int j = i; j > 0; j--)
            divisor[j] = divisor[j - 1] - divisor[j] * root;
        divisor[0] *= -root;
    }
    int CRnum = n - RRnum;  // Complex Roots Number
    com.resize(CRnum + 1, 0);
    double first_coeff = a[n];
    for (int i = CRnum; i >= 0; i--) {
        com[i] = a[i + RRnum] / first_coeff;
        double k = a[i + RRnum];
        for (int j = RRnum - 1; j >= 0; j--)
            a[i + j] -= k * divisor[j];
    }
}




/*

1
-4
6
-26
65
-42
 
 
 */
