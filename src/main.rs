extern crate mpolynomial;
use color_eyre::eyre::Result;
use mpolynomial::parser::parse_polynomial;
use mpolynomial::MPolynomial;

#[allow(unused_variables, unused_mut, unused_assignments)]
fn main() {
    color_eyre::install().unwrap();

    let mut mpoly1 = MPolynomial::new(2);
    mpoly1.add_coeff(&[0, 1], 1.0);
    mpoly1.add_coeff(&[1, 0], 1.0);
    let mut mpoly2 = mpoly1.clone();
    println!("mpoly1: {}", mpoly1);

    // -=
    mpoly1 -= &mpoly2;
    println!("mpoly1-mpoly2: {}", mpoly1);
    println!("clean -> : {}", mpoly1.drop_zeros());

    // +=
    mpoly1 = mpoly2.clone();
    mpoly1 += &mpoly2;
    println!("mpoly1+mpoly2: {}", mpoly1);

    // Rescale
    mpoly1 = mpoly2.clone();
    println!("mpoly1 * 10: {}", mpoly1.scale(&10.0));

    // Multiply
    mpoly1 = mpoly2.clone();
    println!("mpoly1: {}", mpoly1);
    println!("mpoly2: {}", mpoly2);
    mpoly1.mult(&mpoly2);
    mpoly1.mult(&mpoly2);
    println!("mpoly1*mpoly2*mpoly2: {}", mpoly1);

    // Powers
    println!("\nTEST POWERS");
    let start = std::time::Instant::now();
    println!("Linear Power: (11 x1 * 12 x2 + 99 x3 + 21x4)^30");
    let mpoly = MPolynomial::linear_pown(&[11.0, 12.0, 99.0, 21.0], &[1, 2, 3, 4], 4, 30);
    println!("\t#{} in {:?}", mpoly.coeffs.len(), start.elapsed());
    //for (pows, c) in mpoly.powers.iter().zip(mpoly.coeffs.iter()) {
    //    print!("{:+}", c);
    //    for (n, p) in pows.iter().enumerate() {
    //        if *p > 0 {
    //            print!("*x{}^{}", n+1, p);
    //        }
    //    }
    //}
    //println!("");
    // test the time
    let n_var = 3;
    let n = 5;
    // Create two copies of (3x1 +5x2+ x3)
    let mut mpoly1 = mpolynomial::MPolynomial::new(n_var);
    mpoly1.add_coeff(&[1, 0, 0], 3.0);
    mpoly1.add_coeff(&[0, 1, 0], 5.0);
    mpoly1.add_coeff(&[0, 0, 1], 1.0);
    mpoly1.add_coeff(&[1, 1, 0], 3.0);
    mpoly1.add_coeff(&[0, 1, 1], 5.0);
    mpoly1.add_coeff(&[1, 0, 1], 1.0);
    mpoly1.add_coeff(&[1, 0, 1], 3.0);
    mpoly1.add_coeff(&[1, 1, 0], 5.0);
    mpoly1.add_coeff(&[0, 1, 1], 1.0);
    mpoly1.add_coeff(&[1, 1, 1], 3.0);
    //    mpoly1.add_coeff(&vec![0, 2, 0], 5.0);
    //    mpoly1.add_coeff(&vec![0, 0, 2], 1.0);
    //    mpoly1.add_coeff(&vec![2, 0, 0], 3.0);
    //    mpoly1.add_coeff(&vec![1, 2, 0], 5.0);
    //    mpoly1.add_coeff(&vec![1, 0, 2], 1.0);
    //    mpoly1.add_coeff(&vec![2, 1, 0], 3.0);
    //    mpoly1.add_coeff(&vec![0, 2, 1], 5.0);
    //    mpoly1.add_coeff(&vec![0, 1, 2], 1.0);
    //    mpoly1.add_coeff(&vec![2, 0, 1], 3.0);
    //    mpoly1.add_coeff(&vec![1, 2, 1], 5.0);
    //    mpoly1.add_coeff(&vec![1, 1, 2], 1.0);

    let mut mpoly2 = mpoly1.clone();
    let mut mpoly3 = mpoly1.clone();
    let mut mpoly4 = mpoly1.clone();
    println!("General: {} pow {}", mpoly1, n);

    // multiply
    let start = std::time::Instant::now();
    for _ in 1..n {
        mpoly1.mult(&mpoly2);
    }
    println!("\t[mult] #{} in {:?}", mpoly1.coeffs.len(), start.elapsed());
    let start = std::time::Instant::now();
    mpoly2.pown(n);
    println!("\t[pown] #{} in {:?}", mpoly2.coeffs.len(), start.elapsed());
    assert_eq!(mpoly1, mpoly2);
    let start = std::time::Instant::now();
    mpoly4.pown3(n);
    println!(
        "\t[pown3] #{} in {:?}",
        mpoly4.coeffs.len(),
        start.elapsed()
    );
    assert_eq!(mpoly2, mpoly4);
    let start = std::time::Instant::now();
    mpoly3.pown2(n);
    println!(
        "\t[pown2] #{} in {:?}",
        mpoly3.coeffs.len(),
        start.elapsed()
    );
    assert_eq!(mpoly2, mpoly3);

    // power of linear expression
    println!("Linear:");
    let mut mpoly1 = mpolynomial::MPolynomial::new(n_var);
    mpoly1.add_coeff(&[1, 0, 0], 3.0);
    mpoly1.add_coeff(&[0, 1, 0], 5.0);
    mpoly1.add_coeff(&[0, 0, 1], 1.0);
    let mut mpoly2 = mpoly1.clone();

    let start = std::time::Instant::now();
    for _ in 1..n {
        mpoly1.mult(&mpoly2);
    }
    println!("\t[mult]        in {:?}", start.elapsed());
    let start = std::time::Instant::now();
    mpoly2.pown2(n);
    println!("\t[pown2]       in {:?}", start.elapsed());
    assert_eq!(mpoly1, mpoly2);
    let start = std::time::Instant::now();
    mpoly2 = mpolynomial::MPolynomial::linear_pown(&[1.0, 3.0, 5.0], &[3, 1, 2], n_var, n);
    println!("\t[linear_pown] in {:?}", start.elapsed());
    println!("#{}", mpoly1.coeffs.len());
    assert_eq!(mpoly1, mpoly2);

    let mut mpoly: MPolynomial<f64> = MPolynomial::new(3);
    mpoly.coeffs.clear();
    mpoly.powers.clear();
    mpoly.add_coeff(&[1, 0, 1], 1.0);
    mpoly.add_coeff(&[2, 1, 0], 3.0);
    mpoly.add_coeff(&[1, 3, 0], 5.0);
    println!("{}", mpoly);
    mpoly.replace(1, &[7.0], &[0]);
    println!("{}", mpoly);

    // Scalar
    let n_var = 3;
    let n = 5;
    // Create two copies of (3x1 +5x2+ x3)
    let mut mpoly1 = mpolynomial::MPolynomial::new(n_var);
    mpoly1.add_coeff(&[0, 0, 1], 3.0);
    let mut mpoly2 = mpoly1.clone();
    let mut mpoly3 = mpoly1.clone();
    let mut mpoly4 = mpoly1.clone();
    println!("Scalar: {} pow {} = {}", mpoly1, n, 3.0_f64.powi(n as i32));

    // multiply
    let start = std::time::Instant::now();
    for _ in 1..n {
        mpoly1.mult(&mpoly2);
    }
    println!(
        "\t[mult]  #{} in {:?}",
        mpoly1.coeffs.len(),
        start.elapsed()
    );
    let start = std::time::Instant::now();
    mpoly2.pown(n);
    println!(
        "\t[pown]  #{} in {:?}",
        mpoly2.coeffs.len(),
        start.elapsed()
    );
    assert_eq!(mpoly1, mpoly2);
    let start = std::time::Instant::now();
    mpoly4.pown3(n);
    println!(
        "\t[pown3] #{} in {:?}",
        mpoly4.coeffs.len(),
        start.elapsed()
    );
    assert_eq!(mpoly2, mpoly4);
    let start = std::time::Instant::now();
    mpoly3.pown2(n);
    println!(
        "\t[pown2] #{} in {:?}",
        mpoly3.coeffs.len(),
        start.elapsed()
    );
    assert_eq!(mpoly2, mpoly3);
    let start = std::time::Instant::now();
    mpoly4 = mpolynomial::MPolynomial::linear_pown(&[3.0], &[3], n_var, n);
    println!("\t[linear_pown] {:?}", start.elapsed());
    assert_eq!(mpoly3, mpoly4);
    let start = std::time::Instant::now();
    let mpoly5 = 3.0_f64.powi(n as i32);
    println!("\t[powi]        {:?}", start.elapsed());
    let start = std::time::Instant::now();
    let mpoly5 = MPolynomial::powi(&3.0, n as i32);
    println!("\t[powi]        {:?}", start.elapsed());

    // Use BigRational
    let start = std::time::Instant::now();
    let c1 = num::BigInt::from(-10 as i64);
    let c2 = num::BigInt::from(-231 as i64);
    let c = num::BigRational::from(c1) / num::BigRational::from(c2);
    let mut mpoly: MPolynomial<num::BigRational> = MPolynomial::new(3);
    println!("{}", c);
    mpoly.add_coeff(&[0, 0, 0], c.clone());
    let mut poly = mpoly.clone();
    let start = std::time::Instant::now();
    poly.pown(100);
    println!("\t[pown]        {:?}", start.elapsed());
    poly = mpoly.clone();
    let start = std::time::Instant::now();
    mpoly.pown2(100);
    println!("\t[pow2]        {:?}", start.elapsed());
    poly = mpoly.clone();
    let start = std::time::Instant::now();
    MPolynomial::powi(&c, 100);
    println!("\t[powi]        {:?}", start.elapsed());
    //println!("{}", mpoly);

    //
    let mut mpoly: MPolynomial<num::BigRational> = MPolynomial::new(3);
    let c1 = num::BigRational::from(num::BigInt::from(123 as i64))
        / num::BigRational::from(num::BigInt::from(-145 as i64));
    let c2 = num::BigRational::from(num::BigInt::from(58 as i64))
        / num::BigRational::from(num::BigInt::from(85 as i64));
    let v1 = num::BigRational::from(num::BigInt::from(1 as i64))
        / num::BigRational::from(num::BigInt::from(2 as i64));

    mpoly.add_coeff(&[1, 1, 0], c1);
    mpoly.add_coeff(&[1, 2, 0], c2);
    let var_names = vec![String::from("z"), String::from("w2"), String::from("x")];
    println!("{}", mpoly);
    println!("{}", mpoly.to_str(&var_names));

    mpoly.replace(2, &[v1], &[0]);
    println!("{}", mpoly);
    println!("{:?}", 1_i32.checked_div(2));

    // Try factorizing
    let var_names = vec![
        String::from("x1"),
        String::from("x2"),
        String::from("x3"),
        String::from("x4"),
    ];

    let mpoly_str = "(11 x1 + 12 x2 + 99 x3 + 21x4)^30";
    let factor_str = "(11 x1 + 12 x2 + 99 x3 + 21x4)^28";
    println!("Check factorization of : {}", mpoly_str);
    println!("  factor: {}", factor_str);
    let c1 = num::BigRational::from(num::BigInt::from(11 as i64));
    let c2 = num::BigRational::from(num::BigInt::from(12 as i64));
    let c3 = num::BigRational::from(num::BigInt::from(99 as i64));
    let c4 = num::BigRational::from(num::BigInt::from(21 as i64));
    let start = std::time::Instant::now();
    let mut mpoly = MPolynomial::linear_pown(
        &[c1.clone(), c2.clone(), c3.clone(), c4.clone()],
        &[1, 2, 3, 4],
        4,
        30,
    );
    let factor = MPolynomial::linear_pown(
        &[c1.clone(), c2.clone(), c3.clone(), c4.clone()],
        &[1, 2, 3, 4],
        4,
        28,
    );
    //println!("{}", mpoly);
    println!(
        "\t[linear_pown] {} in {:?}",
        mpoly.coeffs.len(),
        start.elapsed()
    );
    //    factor.add_coeff(&[2,2,2,2],1.0);
    let start = std::time::Instant::now();
    let mut factor_b = parse_polynomial(factor_str, &var_names);
    let mut mpoly_b = parse_polynomial(mpoly_str, &var_names);
    println!(
        "\t[parser] {} in {:?}",
        mpoly_b.coeffs.len(),
        start.elapsed()
    );

    //println!("{}",mpoly);
    //println!("{}",factor);
    let start = std::time::Instant::now();
    mpoly.exact_division(&factor).unwrap();
    //println!("{}", mpoly);
    assert_eq!("+441*x4^2+4158*x3*x4+9801*x3^2+504*x2*x4+2376*x2*x3+144*x2^2+462*x1*x4+2178*x1*x3+264*x1*x2+121*x1^2", mpoly.to_str(&var_names));
    println!("\t[exact division] in {:?}", start.elapsed());

    factor_b -= &factor;
    factor_b.drop_zeros();
    //println!("{}", mpoly_b);
    //println!("{}", factor_b);
    //println!("{}", factor_b.to_str(&var_names));
}
