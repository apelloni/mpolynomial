use mpolynomial::MPolynomial;

fn main() {
    let mut mpoly: MPolynomial<f64> = MPolynomial::new(3);
    mpoly.coeffs.clear();
    mpoly.powers.clear();
    println!("{:?}", mpoly);
    mpoly.add(&vec![0, 0, 0], 7.0);
    mpoly.add(&vec![0, 1, 0], 1.0);
    mpoly.add(&vec![1, 1, 0], 2.0);
    mpoly.add(&vec![0,1], 100.0);
    mpoly.add(&vec![3, 1, 0], 3.0);
    mpoly.add(&vec![1, 2, 0], 4.0);
    mpoly.add(&vec![1, 1, 3], 5.0);
    mpoly.add(&vec![0, 0, 1], 6.0);
    mpoly.add(&vec![0, 0], 100.0);
    mpoly.add(&vec![0, 0, 0], 7.0);
    mpoly.add(&vec![3,1], 100.0);
    println!("{:?}", mpoly);
    println!("{}", mpoly);
    let start = std::time::Instant::now();
    //let mpoly = MPolynomial::linear_pown(&[1.0, 1.0], &[1, 2], 2, 2);
    //println!("{}", mpoly);
    let mpoly = MPolynomial::linear_pown(&[11.0, 12.0, 99.0, 21.0,], &[1, 2, 3, 4], 4, 30);
    println!("#{} in {:?}", mpoly.coeffs.len(), start.elapsed());
    //for (pows, c) in mpoly.powers.iter().zip(mpoly.coeffs.iter()) {
    //    print!("{:+}", c);
    //    for (n, p) in pows.iter().enumerate() {
    //        if *p > 0 {
    //            print!("*x{}^{}", n+1, p);
    //        }
    //    }
    //}
    //println!("");
    let mut mpoly: MPolynomial<f64> = MPolynomial::new(3);
    mpoly.coeffs.clear();
    mpoly.powers.clear();
    mpoly.add(&vec![1,0,1], 1.0);
    mpoly.add(&vec![2,1,0], 3.0);
    mpoly.add(&vec![1,3,0], 5.0);
    println!("{}", mpoly);
    mpoly.replace(1, &[7.0,11.0,-13.0,17.0], &[0,1,2,3]);
    println!("{}", mpoly);

    println!("Hello, world!");
}
