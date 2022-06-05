Multi-dimensional Polynomials
==============================

Peform operations with multivariate polynomials such as replacement and multiplication in **rust**.
With this library one can import polynomial functions into rust in an easy way and perform simple operations end evaluations.


The supported type for the coefficients are:
 - `f64` (real and complex)
 - `f128` (real and complex)
 - `num::BigRational`

There is a parser that imports expressions from string into `num::BigRational`:

```rust
extern crate mpolynomial;

use mpolynomial::MPolynomial;
use mpolynomial::parser::parse_expression;
use num::{BigInt, BigRational};

fn main() {
    let var_names = vec![
        String::from("x1"),
        String::from("x2"),
        String::from("x3"),
        String::from("x4"),
    ];

    // Geometric Series
    let mpoly_str = "(1-x1^11)^2";
    let factor_str = "(1-x1)^2";
    let res_str = "(1+x1+x1^2+x1^3+x1^4+x1^5+x1^6+x1^7+x1^8+x1^9+x1^10)^2";

    let mut mpoly = parse_expression(mpoly_str, &var_names);
    let factor = parse_expression(factor_str, &var_names);
    let res = parse_expression(res_str, &var_names);

    mpoly.exact_division(&factor).unwrap();
    assert_eq!(res.to_str(&var_names), mpoly.to_str(&var_names));

    // Multivariate
    let mpoly_str = "(11 x1^32 + 12 x2 + 99 x3 + 21 x1*x4)^10";
    let factor_str = "(11 x1**32 + 12 x2 + 99 x3 + 21 x1*x4)^8";

    let factor = parse_expression(factor_str, &var_names);
    let mut mpoly = parse_expression(mpoly_str, &var_names);

    mpoly.exact_division(&factor).unwrap();
    //println!("{}", mpoly);
    let mut res_str = String::new();
    res_str += "+9801*x3^2+2376*x2^1*x3^1+144*x2^2+4158*x1^1*x3^1*x4^1";
    res_str += "+504*x1^1*x2^1*x4^1+441*x1^2*x4^2+2178*x1^32*x3^1";
    res_str += "+264*x1^32*x2^1+462*x1^33*x4^1+121*x1^64";
    assert_eq!(res_str, format!("{}", mpoly));
}
```


Here is a list of the function implemented for the structure `MPolynomial`:
 - `add`: add a monomial to the polynomial
 - `update`: update the coefficient of a monomial
 - `replace`: replace one variable with a linear polynomial (e.g. `x1-> c+2x1+3x2`)
 - `mult`: multiply the current polynomial by another and overwrite
 - `exact_division`: divide by one of the factors of the polynomial (otherwise fails)
 - `euclidean_division`: returns both the quotient and remainder of the division
 - `multivariate_gcd`: compute the gcd of two polynomial
 - ...
 - ...
