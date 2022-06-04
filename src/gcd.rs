use crate::{Field, MPolynomial};

/// Euclidean division
/// Input: a and b ≠ 0 two polynomials in the variable x;
/// Output: q, the quotient, and r, the remainder;
/// Begin
///     q := 0
///     r := a
///     d := deg(b)
///     c := lc(b)
///     while deg(r) ≥ d do
///         s := lc(r)/c x^(deg(r)−d)
///         q := q + s
///         r := r − sb
///     end do
///     return (q, r)
/// end
pub fn euclidean_division<T: Field>(
    poly_a: &MPolynomial<T>,
    poly_b: &MPolynomial<T>,
) -> (MPolynomial<T>, MPolynomial<T>) {
    assert!(
        poly_a.n_var == 1,
        "univariate_gcd only for univariate polynomials"
    );
    assert!(
        poly_b.n_var == 1,
        "univariate_gcd only for univariate polynomials"
    );
    if poly_a.powers.last().unwrap()[0] < poly_b.powers.last().unwrap()[0] {
        return euclidean_division(poly_b, poly_a);
    }

    let mut r = poly_a.clone();
    let mut q = poly_b.clone();
    if r.is_zero() {
        return (poly_b.clone(), MPolynomial::zero(1));
    }
    if q.is_zero() {
        return (poly_a.clone(), MPolynomial::zero(1));
    }

    let mut s = MPolynomial::zero(1);
    //let mut leading_coeff = q.coeffs.last().expect("polynomial is zero!");
    let c = q.leading_coeff().expect("polynomial is zero!").0.clone();
    let deg = poly_b.powers.last().unwrap()[0];

    q.clear();
    while r.powers.last().unwrap()[0] >= deg {
        //println!("r: {:?}",r.coeffs.iter().map(|x| format!("{}",x)).collect::<Vec<String>>());
        //println!("c = {}",&c);
        s.coeffs[0] = r.leading_coeff().unwrap().0.clone() / c.clone();
        s.powers[0][0] = r.powers.last().unwrap()[0] - deg;

        q += &s;
        r -= &(&s * poly_b);
        r.drop_zeros();
    }

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
pub fn univariate_gcd<T: Field>(
    poly_a: &MPolynomial<T>,
    poly_b: &MPolynomial<T>,
) -> MPolynomial<T> {
    let mut r0 = poly_a.clone();
    let mut r1 = poly_b.clone();
    let mut r2 : MPolynomial<T>;
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
    }
}
