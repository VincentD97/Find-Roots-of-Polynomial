# Find Roots of Polynomial

A program to find all complex roots for a real-coefficient polynomial equation.

## Attempt I. Bisection Method

Advantage :

- Easy implementation.

Disadvantage :

- Only work for real roots.

## Attempt II. Bairstow's Method

Find all quadratic factors of the polynomial and solve each of them. See [here](http://web2.uwindsor.ca/courses/engineering/ahmed/PDF%20FILES%20FOR%20C++/Chapter%20IX%20C++.pdf).

Advantage :

- Work for complex roots.

Disadvantage :

- Behave badly at repeated roots (multiple roots or close roots). e.g. (x-2)^3(x+1)^2 or (x+1)^8. Since for x in (root-∆, root+∆), f(x) is too close to 0, the program cannot precisely extract quadratic equations.

## Final Work. Aberth Method

See [here](https://en.wikipedia.org/wiki/Aberth_method).

We first initialize the approximations of zeros within bounds given by [Rouché theorem](https://en.wikipedia.org/wiki/Properties_of_polynomial_roots).

Then we adjust their values until the offsets approach 0. Suppose the current approximations of the zeros are

<img src="http://latex.codecogs.com/gif.latex?z_1,\dots,z_n\in\mathbb{C}"/>

We calculate a set of better approximations by

<img src="http://latex.codecogs.com/gif.latex?z_k'=z_k-\frac1{\frac{p'(z_k)}{p(z_k)}-\sum_{j=0;\,j!=k}^n\frac1{z_k-z_j}}"/>

[derived](https://en.wikipedia.org/wiki/Aberth_method) from the univariate Newton method iteration on an optimized function.


## Test Case

Test Case 1
```
 deg = 4
 a[4] = 2 
 a[3] = 15
 a[2] = 11
 a[1] = 13
 a[0] = -105
```
 
Test Case 2
```
 deg = 10
 a[10] = 5
 a[9] = -50
 a[8] = 175
 a[7] = -180
 a[6] = -325
 a[5] = 830
 a[4] = -135
 a[3] = -920
 a[2] = 520
 a[1] = 320
 a[0] = -240
```

Test Case 3
```
 deg = 5
 a[5] = 1
 a[4] = -4
 a[3] = 6
 a[2] = -26
 a[1] = 65
 a[0] = -42
```
