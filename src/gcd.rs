use crate::{Field, MPolynomial, NumberLike};
use num::BigRational;

/// Euclidean division
/// Same as long division but it reports the quotient and remainder
pub fn euclidean_division<T: Field + NumberLike>(
    poly_a: &MPolynomial<T>,
    poly_b: &MPolynomial<T>,
) -> (MPolynomial<T>, MPolynomial<T>) {
    assert_eq!(poly_a.n_var, poly_b.n_var);
    let mut r = poly_a.clone();
    let mut q = poly_b.clone();
    if r.is_zero() {
        return (poly_b.clone(), MPolynomial::zero(1));
    }
    if q.is_zero() {
        return (poly_a.clone(), MPolynomial::zero(1));
    }

    let mut s = MPolynomial::zero(poly_a.n_var);
    let (c, deg) = poly_b.leading_coeff().expect("polynomial is zero!");

    //println!("deg : {:?}", deg);
    //println!("c = {}", &c);
    //q.clear();
    q = MPolynomial::zero(poly_a.n_var);
    //    while r.powers.last().unwrap() >= &deg {
    while r
        .powers
        .last()
        .unwrap()
        .iter()
        .zip(deg.iter())
        .all(|(p_to, p_from)| p_to >= p_from)
    {
        match r.leading_coeff() {
            Some((leading_c, leading_deg)) => {
                // Update rescaling monomial
                s.coeffs[0] = leading_c.clone() / c.clone();
                for (pn, p) in leading_deg.iter().zip(deg.iter()).enumerate() {
                    s.powers[0][pn] = p.0 - p.1;
                }

                // Append to the quotient and subtract from remainder
                q += &s;
                r -= &(&s * poly_b);
                r.drop_zeros();
            }
            None => {
                return (q.clone(), r.clone());
            }
        }
    }
    //println!("q={};", q);
    //println!("r={};", r);
    return (q.clone(), r.clone());
}

/// Univariate GCD
/// Input: a and b ≠ 0 two polynomials in the variable x;
/// Output: q, the quotient, and r, the remainder;
/// Begin
///     r_0:=a
///     r_1:=b
///     r_2:=0
///     for ( i := 1 ; r1 ≠ 0; i := i + 1) do
///         r2 := rem( r0 , r1 )
///         r0=r1
///         r1=r2
///     end do
///     return r1
///  end
pub fn univariate_gcd<T: Field + NumberLike>(
    poly_a: &MPolynomial<T>,
    poly_b: &MPolynomial<T>,
) -> MPolynomial<T> {
    assert!(
        poly_a.n_var == 1,
        "univariate_gcd only for univariate polynomials"
    );
    assert!(
        poly_b.n_var == 1,
        "univariate_gcd only for univariate polynomials"
    );
    //if poly_a.powers.last().unwrap()[0] < poly_b.powers.last().unwrap()[0] {
    if poly_a.powers.last().unwrap() < poly_b.powers.last().unwrap() {
        return univariate_gcd(poly_b, poly_a);
    }

    let mut r0 = poly_a.clone();
    let mut r1 = poly_b.clone();
    let mut r2: MPolynomial<T>;
    while !r1.is_zero() {
        r2 = euclidean_division(&r0, &r1).1;
        r0 = r1.clone();
        r1 = r2.clone();
    }
    // Before return the GCD normalized it w.r.t. the second polynomial
    r0.drop_zeros();
    r0.scale(&(poly_b.coeffs[0].clone() / r0.coeffs[0].clone()))
        .clone()
}

pub fn multivariate_gcd(
    poly_a: &MPolynomial<BigRational>,
    poly_b: &MPolynomial<BigRational>,
) -> MPolynomial<BigRational> {
    assert!(
        poly_a.n_var == poly_b.n_var,
        "GCD required two polynomial with the same number of variables"
    );
    //if poly_a.powers.last().unwrap()[0] < poly_b.powers.last().unwrap()[0] {
    if poly_a.powers.last().unwrap() < poly_b.powers.last().unwrap() {
        return multivariate_gcd(poly_b, poly_a);
    }

    let mut r0 = poly_a.clone();
    let mut r1 = poly_b.clone();
    let mut r2: MPolynomial<BigRational>;
    let mut q: MPolynomial<BigRational>;
    while !r1.is_zero() {
        (q, r2) = euclidean_division(&r0, &r1);
        // If q is zero it means that the reduction of r0 by r1 is complete
        // -> r0 and r1 make up our Grobner basis and the divisor of the
        //  principle part (w.r.t. the leading variable in the lexicographic order)
        //  of the one with the lowest degree corresponds to the commond divisor
        if q.is_zero() {
            //println!("r0: {}", r0);
            //println!("r1: {}", r1);
            if r0.leading_coeff().unwrap().1 > r1.leading_coeff().unwrap().1 {
                r0 = r1.clone();
            }
            // Obtain the highest rank that can appear in the principle part
            r0.ranks_update();
            let (pp_var, pp_deg) = r0
                .max_rank
                .iter()
                .enumerate()
                .find(|(_, &r)| r > 0)
                .unwrap();
            // Extract the coefficient of the principle part (up to a monomial)
            r1 = MPolynomial::one(r0.n_var);
            r1.powers[0][pp_var] = *pp_deg;
            (q, _) = euclidean_division(&r0, &r1);
            // Any overall monomial will be part of the principle part and not of
            // its coefficient
            if !q.is_monomial() {
                q.monomial_rescale();
            }
            // Get the principle part and store it in r0
            (r0, _) = euclidean_division(&r0, &q);
            break;
        }

        r0 = r1.clone();
        r1 = r2.clone();
    }
    // Before return the GCD normalized it w.r.t. the second polynomial
    r0.drop_zeros();
    r0.scale(&(poly_b.coeffs[0].clone() / r0.coeffs[0].clone()))
        .clone()
}

#[cfg(test)]
mod tests {
    //use super::*;
    use crate::parser::parse_expression;

    #[test]
    fn check_euclidean_division() {
        let variable = &[String::from("x")];
        let mut poly_a = parse_expression("(1-x)(1+x)^2", variable);
        let poly_b = parse_expression("(1-x)", variable);
        let res = parse_expression("(1-x)", variable);

        println!("a = {}", poly_a.to_str(variable));
        println!("b = {}", poly_b.to_str(variable));

        let start = std::time::Instant::now();
        let q = super::euclidean_division(&poly_a, &poly_b);
        println!("\tEuclidean Division in {:?}", start.elapsed());

        let start = std::time::Instant::now();
        poly_a.exact_division(&q.0).unwrap();
        println!("\tLong Division  in {:?}", start.elapsed());

        println!("q = {}", q.0.to_str(variable));
        println!("r = {}", q.1.to_str(variable));
        println!("a/b = {}", poly_a.to_str(variable));
        assert_eq!(res.to_str(variable), poly_a.to_str(variable));

        // Geometric Series
        println!("========");
        let poly_str = "(1-x^11)^2";
        let factor_str = "(1-x)^2";
        let res_str = "(1+x+x^2+x^3+x^4+x^5+x^6+x^7+x^8+x^9+x^10)^2";

        let mut poly = parse_expression(poly_str, variable);
        let factor = parse_expression(factor_str, variable);
        let res = parse_expression(res_str, variable);

        println!("a = {}", poly.to_str(variable));
        println!("b = {}", factor.to_str(variable));

        let start = std::time::Instant::now();
        let q = super::euclidean_division(&poly, &factor);
        println!("\tEuclidean Division in {:?}", start.elapsed());

        let start = std::time::Instant::now();
        poly.exact_division(&factor).unwrap();
        println!("\tLong Division  in {:?}", start.elapsed());
        println!("q = {}", q.0.to_str(variable));
        println!("r = {}", q.1.to_str(variable));
        println!("a/b = {}", poly_a.to_str(variable));
        assert_eq!(res.to_str(variable), poly.to_str(variable));
        assert_eq!(res.to_str(variable), q.0.to_str(variable));
    }

    #[test]
    fn check_univariate_gcd() {
        let variable = &[String::from("x")];
        let mut poly_a = parse_expression("(1-x)(1+x)^2", variable);
        let mut poly_b = parse_expression("(1-x)", variable);
        let red_a = parse_expression("(1+x)^2", variable);
        let red_b = parse_expression("1", variable);

        println!("a = {}", poly_a.to_str(variable));
        println!("b = {}", poly_b.to_str(variable));

        let start = std::time::Instant::now();
        let gcd = super::univariate_gcd(&poly_a, &poly_b);
        println!("\tGCD in {:?}", start.elapsed());

        let start = std::time::Instant::now();
        poly_a.exact_division(&gcd).unwrap();
        poly_b.exact_division(&gcd).unwrap();
        println!("\tLD  in {:?}", start.elapsed());

        println!("gcd = {}", gcd.to_str(variable));
        println!(
            "a/b = ({}) / ({})",
            poly_a.to_str(variable),
            poly_b.to_str(variable)
        );
        assert_eq!(red_a.to_str(variable), poly_a.to_str(variable));
        assert_eq!(red_b.to_str(variable), poly_b.to_str(variable));

        println!("========");
        //let mpoly_str = "(11 x1^32 + 12 x2 + 99 x3 + 21 x1*x4)^10";
        //let factor_str = "(11 x1**32 + 12 x2 + 99 x3 + 21 x1*x4)^8";

        let mut poly_a = parse_expression("(1-x)(1+x)^2", variable);
        let mut poly_b = parse_expression("(1-x)(1+3x)", variable);
        let red_a = parse_expression("(1+x)^2", variable);
        let red_b = parse_expression("(1+3x)", variable);

        println!("a = {}", poly_a.to_str(variable));
        println!("b = {}", poly_b.to_str(variable));

        let start = std::time::Instant::now();
        let gcd = super::univariate_gcd(&poly_a, &poly_b);
        println!("\tGCD in {:?}", start.elapsed());

        let start = std::time::Instant::now();
        poly_a.exact_division(&gcd).unwrap();
        poly_b.exact_division(&gcd).unwrap();
        println!("\tLD  in {:?}", start.elapsed());

        println!("gcd = {}", gcd.to_str(variable));
        println!(
            "a/b = ({}) / ({})",
            poly_a.to_str(variable),
            poly_b.to_str(variable)
        );
        assert_eq!(red_a.to_str(variable), poly_a.to_str(variable));
        assert_eq!(red_b.to_str(variable), poly_b.to_str(variable));
        //panic!();
    }

    #[test]
    fn multivariate_gcd() {
        let var_names = vec![
            String::from("x1"),
            String::from("x2"),
            String::from("x3"),
            String::from("x4"),
        ];
        // Multivariate
        let mpoly_a_str = "(x3-x2^2)^2";
        let mpoly_b_str = "(x3-x2^2)^1(1-x2)";
        let mut mpoly_a = parse_expression(mpoly_a_str, &var_names);
        let mut mpoly_b = parse_expression(mpoly_b_str, &var_names);
        println!("b = {}", mpoly_b.to_str(&var_names));

        println!("a = {}", mpoly_a_str);
        println!("b = {}", mpoly_b_str);
        let start = std::time::Instant::now();
        let gcd = super::multivariate_gcd(&mpoly_a, &mpoly_b);
        println!("\tGCD in {:?}", start.elapsed());

        let start = std::time::Instant::now();
        mpoly_a.exact_division(&gcd).unwrap();
        mpoly_b.exact_division(&gcd).unwrap();
        println!("\tLD  in {:?}", start.elapsed());

        println!("gcd = {}", gcd.to_str(&var_names));
        println!(
            "a/b = ({}) / ({})",
            mpoly_a.to_str(&var_names),
            mpoly_b.to_str(&var_names)
        );

        // Multivariate
        let mpoly_a_str = "(11 x1^32 + 12 x2 + 99 x3 + 21 x1*x4)^3";
        let mpoly_b_str = "(11 x1^32 + 12 x2 + 99 x3 + 21 x1*x4)^2*(1-x1)";
        let mut mpoly_a = parse_expression(mpoly_a_str, &var_names);
        let mut mpoly_b = parse_expression(mpoly_b_str, &var_names);
        println!("b = {}", mpoly_b.to_str(&var_names));

        println!("a = {}", mpoly_a_str);
        println!("b = {}", mpoly_b_str);
        let start = std::time::Instant::now();
        let gcd = super::multivariate_gcd(&mpoly_a, &mpoly_b);
        println!("\tGCD in {:?}", start.elapsed());

        let start = std::time::Instant::now();
        mpoly_a.exact_division(&gcd).unwrap();
        mpoly_b.exact_division(&gcd).unwrap();
        println!("\tLD  in {:?}", start.elapsed());

        println!("gcd = {}", gcd.to_str(&var_names));
        println!(
            "a/b = ({}) / ({})",
            mpoly_a.to_str(&var_names),
            mpoly_b.to_str(&var_names)
        );
    }
}
