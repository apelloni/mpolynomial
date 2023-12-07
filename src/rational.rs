use crate::gcd::multivariate_gcd;
use crate::{Field, MPolynomial, NumberLike};
use num::BigRational;
use regex::Regex;
use std::fmt;
use std::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign};

/// MPolyRat needs absolute precision in order to
/// simplify the expression for this reason is implemented
/// only for BigRationlals
#[derive(Debug, Clone)]
pub struct MPolyRat {
    pub num: MPolynomial<BigRational>,
    pub den: MPolynomial<BigRational>,
    pub n_var: usize,
}

impl MPolyRat {
    pub fn new(n_var: usize) -> MPolyRat {
        MPolyRat {
            num: MPolynomial::new(n_var),
            den: MPolynomial::one(n_var),
            n_var,
        }
    }
    // TODO: implement from string
    pub fn from_mpolynomials(
        num: &MPolynomial<BigRational>,
        den: &MPolynomial<BigRational>,
    ) -> MPolyRat {
        assert!(
            num.n_var == den.n_var,
            "Numerator and Denominator must have the same number of variables"
        );
        MPolyRat {
            num: num.clone(),
            den: den.clone(),
            n_var: num.n_var,
        }
    }

    pub fn drop_zeros(&mut self) {
        self.num.drop_zeros();
        self.den.drop_zeros();
    }

    pub fn clear(&mut self) {
        self.num.clear();
        self.den = MPolynomial::one(self.n_var);
    }

    pub fn is_constant(&mut self) -> bool {
        self.reduce();
        self.num.is_constant() && self.den.is_constant()
    }

    pub fn reduce(&mut self) {
        let factor = multivariate_gcd(&self.num, &self.den);
        //println!("num = {};",self.num);
        //println!("den = {};",self.den);
        //println!("factor = {};",factor);
        self.num.exact_division(&factor).unwrap();
        self.den.exact_division(&factor).unwrap();
    }

    pub fn scale(&mut self, factor: &BigRational) {
        self.num.scale(factor);
    }

    pub fn pown(&mut self, n: usize) {
        self.num.pown(n);
        self.den.pown(n);
    }

    /// Print the polynomail to string using the given names for the variables
    /// The variable names must be of the form <char:1 or more><digits:0 or more>
    pub fn to_str(&self, var_names: &[String]) -> String {
        format!(
            "({})/({})",
            self.num.to_str(var_names),
            self.den.to_str(var_names)
        )
    }
}

impl fmt::Display for MPolyRat {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "({})/({})", self.num, self.den,)
    }
}

impl<'a> Mul<&'a MPolyRat> for &'a MPolyRat {
    type Output = MPolyRat;
    fn mul(self, other: Self) -> MPolyRat {
        let mut out = self.clone();
        out *= other;
        out
    }
}
impl<'a> Add<&'a MPolyRat> for &'a MPolyRat {
    type Output = MPolyRat;
    fn add(self, other: Self) -> MPolyRat {
        let mut out = self.clone();
        out += other;
        out
    }
}
impl<'a> Sub<&'a MPolyRat> for &'a MPolyRat {
    type Output = MPolyRat;
    fn sub(self, other: Self) -> MPolyRat {
        let mut out = self.clone();
        out -= other;
        out
    }
}

impl<'a> MulAssign<&'a MPolyRat> for MPolyRat {
    fn mul_assign(&mut self, other: &'a Self) {
        self.num.mult(&other.num);
        self.den.mult(&other.den);
    }
}

impl<'a> MulAssign<BigRational> for MPolyRat {
    fn mul_assign(&mut self, other: BigRational) {
        self.num.scale(&other);
    }
}

impl<'a> DivAssign<&'a MPolyRat> for MPolyRat {
    fn div_assign(&mut self, other: &'a Self) {
        self.num *= &other.den;
        self.den *= &other.num;
    }
}

impl<'a> DivAssign<BigRational> for MPolyRat {
    fn div_assign(&mut self, other: BigRational) {
        self.num.scale(&(num::one::<BigRational>() / other.clone()));
    }
}

impl<'a> AddAssign<&'a MPolyRat> for MPolyRat {
    fn add_assign(&mut self, other: &'a Self) {
        self.num = &(&self.num * &other.den) + &(&other.num * &self.den);
        self.den *= &other.den;
    }
}

impl<'a> SubAssign<&'a MPolyRat> for MPolyRat {
    fn sub_assign(&mut self, other: &Self) {
        self.num = &(&self.num * &other.den) - &(&other.num * &self.den);
        self.den *= &other.den;
    }
}

impl PartialEq<MPolyRat> for MPolyRat {
    fn eq(&self, other: &MPolyRat) -> bool {
        self.num == other.num && self.den == other.den
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parser::parse_polynomial;

    #[test]
    fn check_univariate() {
        let variable = &[String::from("x")];
        let poly_a = parse_polynomial("(1-x)(1+x)^2", variable);
        let poly_b = parse_polynomial("(1-x)", variable);
        let red_a = parse_polynomial("(1+x)^2", variable);
        let red_b = parse_polynomial("1", variable);
        let mut mpolyrat = MPolyRat::from_mpolynomials(&poly_a, &poly_b);

        println!("num = {}", poly_a.to_str(variable));
        println!("den = {}", poly_b.to_str(variable));

        let start = std::time::Instant::now();
        mpolyrat.reduce();
        println!("\tReduce in {:?}", start.elapsed());
        println!(
            "a/b = ({}) / ({})",
            mpolyrat.num.to_str(variable),
            mpolyrat.den.to_str(variable)
        );
        assert_eq!(red_a.to_str(variable), mpolyrat.num.to_str(variable));
        assert_eq!(red_b.to_str(variable), mpolyrat.den.to_str(variable));

        println!("========");
        //let mpoly_str = "(11 x1^32 + 12 x2 + 99 x3 + 21 x1*x4)^10";
        //let factor_str = "(11 x1**32 + 12 x2 + 99 x3 + 21 x1*x4)^8";

        let poly_a = parse_polynomial("(1-x)(1+x)^2", variable);
        let poly_b = parse_polynomial("(1-x)(1+3x)", variable);
        let red_a = parse_polynomial("(1+x)^2", variable);
        let red_b = parse_polynomial("(1+3x)", variable);
        //panic!();
        let mut mpolyrat = MPolyRat::from_mpolynomials(&poly_a, &poly_b);

        println!("num = {}", poly_a.to_str(variable));
        println!("den = {}", poly_b.to_str(variable));

        let start = std::time::Instant::now();
        mpolyrat.reduce();
        println!("\tReduce in {:?}", start.elapsed());
        println!("a/b = {}", mpolyrat.to_str(variable));
        assert_eq!(red_a.to_str(variable), mpolyrat.num.to_str(variable));
        assert_eq!(red_b.to_str(variable), mpolyrat.den.to_str(variable));
    }

    #[test]
    fn check_multivariate() {
        let var_names = vec![
            String::from("x1"),
            String::from("x2"),
            String::from("x3"),
            String::from("x4"),
        ];
        // Multivariate
        let mpoly_a_str = "(x3-x2^2)^2";
        let mpoly_b_str = "(x3-x2^2)^1(1-x2)";
        let mpoly_a = parse_polynomial(mpoly_a_str, &var_names);
        let mpoly_b = parse_polynomial(mpoly_b_str, &var_names);
        let red_a = parse_polynomial("(x3-x2^2)", &var_names);
        let red_b = parse_polynomial("(1-x2)", &var_names);
        let mut mpolyrat = MPolyRat::from_mpolynomials(&mpoly_a, &mpoly_b);

        println!("num = {}", mpoly_a.to_str(&var_names));
        println!("den = {}", mpoly_b.to_str(&var_names));

        let start = std::time::Instant::now();
        mpolyrat.reduce();
        println!("\tReduce in {:?}", start.elapsed());
        println!("a/b = {}", mpolyrat.to_str(&var_names));
        assert_eq!(red_a.to_str(&var_names), mpolyrat.num.to_str(&var_names));
        assert_eq!(red_b.to_str(&var_names), mpolyrat.den.to_str(&var_names));

        // Multivariate
        let mpoly_a_str = "(11 x1^32 + 12 x2 + 99 x3 + 21 x1*x4)^3";
        let mpoly_b_str = "(11 x1^32 + 12 x2 + 99 x3 + 21 x1*x4)^2*(1-x1)";
        let mpoly_a = parse_polynomial(mpoly_a_str, &var_names);
        let mpoly_b = parse_polynomial(mpoly_b_str, &var_names);
        let red_a = parse_polynomial("(11 x1^32 + 12 x2 + 99 x3 + 21 x1*x4)", &var_names);
        let red_b = parse_polynomial("(1-x1)", &var_names);
        let mut mpolyrat = MPolyRat::from_mpolynomials(&mpoly_a, &mpoly_b);

        println!("num = {}", mpoly_a.to_str(&var_names));
        println!("den = {}", mpoly_b.to_str(&var_names));

        let start = std::time::Instant::now();
        mpolyrat.reduce();
        println!("\tReduce in {:?}", start.elapsed());
        println!("a/b = {}", mpolyrat.to_str(&var_names));

        assert_eq!(red_a.to_str(&var_names), mpolyrat.num.to_str(&var_names));
        assert_eq!(red_b.to_str(&var_names), mpolyrat.den.to_str(&var_names));
    }

    #[test]
    fn operation_mult() {
        let var_names = vec![
            String::from("x1"),
            String::from("x2"),
            String::from("x3"),
            String::from("x4"),
        ];
        // Multivariate
        let mpoly_a_str = "(x3-x2^2)^2";
        let mpoly_b_str = "(x3-x2^2)^1(1-x2)";
        let mpoly_a = parse_polynomial(mpoly_a_str, &var_names);
        let mpoly_b = parse_polynomial(mpoly_b_str, &var_names);
        let red_a = parse_polynomial("(x3-x2^2)", &var_names);
        let red_b = parse_polynomial("(1-x2)", &var_names);
        let mut mpolyrat = MPolyRat::from_mpolynomials(&mpoly_a, &mpoly_b);

        println!("num = {}", mpoly_a.to_str(&var_names));
        println!("den = {}", mpoly_b.to_str(&var_names));

        let start = std::time::Instant::now();
        mpolyrat.reduce();
        println!("\tReduce in {:?}", start.elapsed());
        println!("a/b = {}", mpolyrat.to_str(&var_names));
        assert_eq!(red_a.to_str(&var_names), mpolyrat.num.to_str(&var_names));
        assert_eq!(red_b.to_str(&var_names), mpolyrat.den.to_str(&var_names));

        // Multivariate
        let mpoly_a_str = "(11 x1^32 + 12 x2 + 99 x3 + 21 x1*x4)^3";
        let mpoly_b_str = "(11 x1^32 + 12 x2 + 99 x3 + 21 x1*x4)^2*(1-x1)";
        let mpoly_a = parse_polynomial(mpoly_a_str, &var_names);
        let mpoly_b = parse_polynomial(mpoly_b_str, &var_names);
        let red_a = parse_polynomial("(11 x1^32 + 12 x2 + 99 x3 + 21 x1*x4)", &var_names);
        let red_b = parse_polynomial("(1-x1)", &var_names);
        let mut mpolyrat = MPolyRat::from_mpolynomials(&mpoly_a, &mpoly_b);

        println!("num = {}", mpoly_a.to_str(&var_names));
        println!("den = {}", mpoly_b.to_str(&var_names));

        let start = std::time::Instant::now();
        mpolyrat.reduce();
        println!("\tReduce in {:?}", start.elapsed());
        println!("a/b = {}", mpolyrat.to_str(&var_names));

        assert_eq!(red_a.to_str(&var_names), mpolyrat.num.to_str(&var_names));
        assert_eq!(red_b.to_str(&var_names), mpolyrat.den.to_str(&var_names));
    }
}
