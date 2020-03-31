Multi-dimensional Polynomials
==============================

Peform operations with multivariate polynomials such as replacement and multiplication in rust.
The supported type for the coefficients are f64` and `f128`. 

Here is a list of the function implemented for the structure `MPolynomial`:
 - `add`: add a monomial to the polynomial
 - `update`: update the coefficient of a monomial
 - `replace`: replace one variable with a linear polynomial (e.g. `x1-> c+2x1+3x2`)
 - `mult`: multiply the current polynomial by another and overwrite
