use num::Integer;
use std::fmt;
use crate::{MPolynomial, FloatLike};

/* Utility functions */
pub fn binomial<T: Integer + Copy + std::fmt::Debug>(n: T, k: T) -> T {
    if k > n {
        return T::zero();
    }
    if k > n - k {
        return binomial(n, n - k);
    }
    let mut res = T::one();
    let mut i = T::one();
    loop {
        if i > k {
            break;
        }
        res = res * (n + T::one() - i) / i;
        i = i + T::one();
    }
    res
}

/// Calculate the multinomial coefficient.
pub fn multinomial<T:Integer + Copy + std::fmt::Debug>(k: &[T]) -> T {
    let mut res = T::one();
    let mut p = T::zero();
    for &i in k {
        p = p + i;
        res = res * binomial(p, i);
    }
    res
}

pub fn next_combination_with_replacement(state: &mut [usize], max_entry: usize) -> bool {
    for i in (0..state.len()).rev() {
        if state[i] < max_entry {
            state[i] += 1;
            for j in i + 1..state.len() {
                state[j] = state[i]
            }
            return true;
        }
    }
    false
}

impl fmt::Display for MPolynomial<f64> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for (pows, c) in self.powers.iter().zip(self.coeffs.iter()) {
            write!(f, "{:+}", c)?;
            for (n, p) in pows.iter().enumerate() {
                if *p > 0 {
                    write!(f, "*x{}^{}", n + 1, p)?;
                }
            }
        }
        write!(f, "")
    }
}

impl<T: FloatLike> PartialEq<MPolynomial<T>> for MPolynomial<T> {
    fn eq(&self, other: &MPolynomial<T>) -> bool {
        if self.coeffs.len() != other.coeffs.len() {
            return false;
        }
        for ((p1, c1), (p2, c2)) in self
            .powers
            .iter()
            .zip(self.coeffs.iter())
            .zip(other.powers.iter().zip(other.coeffs.iter()))
        {
            //println!(
            //    "{:?} {:?}  vs {:?} {:?} :: {:?}",
            //    p1,
            //    c1,
            //    p2,
            //    c2,
            //    ((*c1 - *c2) / *c2).abs() * T::from_f64(1e16).unwrap()
            //);
            if p1 != p2 || ((*c1 - *c2) / *c2).abs() > T::from_f64(1e-14).unwrap() {
                //c1 != c2 {
                return false;
            }
        }
        true
    }
}
