Suppose the current approximations of the zeros are

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
