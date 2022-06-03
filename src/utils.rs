extern crate regex;

use crate::{Field, FloatLike, MPolynomial, NumberLike, RealNumberLike};
use num::Integer;
use regex::Regex;
use std::fmt;

/* Utility functions */
pub fn binomial<T: Integer + Clone + std::fmt::Debug>(n: &T, k: &T) -> T {
    if k > n {
        return T::zero();
    }
    if *k > n.clone() - k.clone() {
        return binomial(n, &(n.clone() - k.clone()));
    }
    let mut res = T::one();
    let mut i = T::one();
    loop {
        if i > *k {
            break;
        }

        res = res * (n.clone() + T::one() - i.clone()) / i.clone();
        i = i.clone() + T::one();
    }
    res
}


/// Calculate the multinomial coefficient.
pub fn multinomial<T: Integer + Clone + std::fmt::Debug>(k: &[T]) -> T {
    let mut res = T::one();
    let mut p = T::zero();
    for i in k {
        p = p.clone() + i.clone();
        res = res * binomial(&p, i);
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

impl<T: NumberLike> fmt::Display for MPolynomial<T> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mut is_zero = true;
        for (pows, c) in self.powers.iter().zip(self.coeffs.iter()) {
            if T::zero() != *c {
                is_zero = false;
                write!(f, "{:+}", c)?;
                for (n, p) in pows.iter().enumerate() {
                    if *p > 0 {
                        write!(f, "*x{}^{}", n + 1, p)?;
                    }
                }
            }
        }
        if is_zero {
            write!(f, "+0")
        } else {
            write!(f, "")
        }
    }
}

impl<T: Field> MPolynomial<T> {
    /// Print the polynomail to string using the given names for the variables
    /// The variable names must be of the form <char:1 or more><digits:0 or more>
    pub fn to_str(&self, var_names: &[String]) -> String {
        let rg = Regex::new(r"^[[:alpha:]]+\d*$").unwrap();
        let rg2 = Regex::new(r"([+-])1\*").unwrap();
        assert_eq!(var_names.len(), self.n_var);
        assert!(var_names.iter().all(|var| rg.is_match(var)));

        let mut spoly = String::new();
        for (pows, c) in self.powers.iter().zip(self.coeffs.iter()) {
            if T::zero() != *c {
                std::fmt::write(&mut spoly, format_args!("{:+}", c)).unwrap();
                for (n, p) in pows.iter().enumerate() {
                    if *p > 0 {
                        std::fmt::write(&mut spoly, format_args!("*{}^{}", var_names[n], p))
                            .unwrap();
                    }
                }
            }
        }
        if spoly.len() == 0 {
            spoly = String::from("+0");
        }
        format!("{}", rg2.replace_all(spoly.as_str(), "$1"))
    }
}

impl<T: RealNumberLike> fmt::Display for MPolynomial<num::Complex<T>> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for (i, (pows, c)) in self.powers.iter().zip(self.coeffs.iter()).enumerate() {
            if i > 0 {
                write!(f, "+({:})", c)?;
            } else {
                write!(f, "({:})", c)?;
            }
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
