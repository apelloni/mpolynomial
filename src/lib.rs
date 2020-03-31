use f128;
use std::fmt;

use num::Integer;
use num_traits::{Float, FloatConst, FromPrimitive, Num, One, Signed, ToPrimitive, Zero};
//use utils::Signum;

pub trait FloatLike:
From<f64>
+ Num
+ FromPrimitive
+ Float
+ Signed
+ FloatConst
+ std::fmt::LowerExp
+ std::fmt::Debug
//+ num_traits::float::FloatCore
+ 'static
//+ Signum
{
}

impl FloatLike for f64 {}
impl FloatLike for f128::f128 {}

pub const MAX_VARIABLE: usize = 7;

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
pub fn multinomial<T: Integer + Copy + std::fmt::Debug>(k: &[T]) -> T {
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

#[derive(Default, Debug, Clone)]
pub struct MulinomaialCache<T: FloatLike> {
    powers: Vec<Vec<usize>>,
    coeffs: Vec<T>,
}
impl<T: FloatLike> MulinomaialCache<T> {
    pub fn new() -> MulinomaialCache<T> {
        MulinomaialCache {
            powers: vec![],
            coeffs: vec![],
        }
    }
}

#[derive(Default, Debug, Clone)]
pub struct MPolynomial<T: FloatLike> {
    //pub powers: Vec<(usize,Vec<usize>)>,
    pub powers: Vec<Vec<usize>>,
    pub coeffs: Vec<T>,
    n_var: usize,
    max_rank: usize,
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

impl<T: FloatLike> MPolynomial<T> {
    pub fn new(n_var: usize) -> MPolynomial<T> {
        assert!(n_var <= MAX_VARIABLE, "Increase the value of MAX_VARIABLE");
        MPolynomial {
            powers: vec![],
            coeffs: vec![],
            //powers: vec![vec![0; n_var]],
            //coeffs: vec![T::zero()],
            n_var: n_var,
            max_rank: 0,
        }
    }

    /// Add a new coefficient to the polynomial keeping the list sorted
    pub fn add(&mut self, pows: &Vec<usize>, coeff: T) -> bool {
        //match self.powers.binary_search(&(pows.len(),pows.clone())) {
        match self.powers.binary_search(pows) {
            Ok(pos) => {
                self.coeffs[pos] = self.coeffs[pos] + coeff;
                false
            }
            Err(pos) => {
                self.powers.insert(pos, pows.clone());
                self.coeffs.insert(pos, coeff);
                let rank = pows.iter().sum();
                if rank > self.max_rank {
                    self.max_rank = rank;
                }
                true
            }
        }
    }

    /// Add a new coefficient to the polynomial keeping the list sorted
    pub fn update(&mut self, pows: &Vec<usize>, coeff: T) -> usize {
        //match self.powers.binary_search(&(pows.len(),pows.clone())) {
        match self.powers.binary_search(pows) {
            Ok(pos) => {
                self.coeffs[pos] = coeff;
                pos
            }
            Err(pos) => {
                self.powers.insert(pos, pows.clone());
                self.coeffs.insert(pos, coeff);
                let rank = pows.iter().sum();
                if rank > self.max_rank {
                    self.max_rank = rank;
                }
                pos
            }
        }
    }
    
    // Take the n-th power of a multivariate linear function
    // (c + a1 * x1 + a2 * x2 + .. + ai * xi )^n
    // ids: are the ids coerresponding to the variables and constant
    // TODO: allow to give a container to avoid always creating a new mltnom
    pub fn linear_pown(coeffs: &[T], ids: &[usize], n_var: usize, n: usize) -> MPolynomial<T> {
        let mut mpoly: MPolynomial<T> = MPolynomial::new(n_var);
        // Extract information
        let mut coeff: T;
        let mut pows = vec![0; n_var];
        let mut coeff_pows = [0; MAX_VARIABLE+1];
        let mut vars = vec![0; n];
        loop {
            // Reset values
            for n in 0..n_var {
                pows[n] = 0;
                coeff_pows[n] = 0;
            }
            coeff_pows[n_var] = 0;
            coeff = T::one();

            // Extract info
            for v in vars.iter() {
                coeff = coeff * coeffs[*v];
                // TODO: merge these two 
                if *v > 0{
                    pows[ids[*v]-1] += 1;
                }
                coeff_pows[ids[*v]]+=1;
            }
            // Multiply by multinomial
            mpoly.add(
                &pows,
                coeff * T::from_usize(multinomial::<usize>(&coeff_pows[..=n_var])).unwrap(),
            );
            if !next_combination_with_replacement(&mut vars, ids.len() - 1) {
                break;
            }
        }
        mpoly
    }

    // Replace one of the variables by a multinomial which is at most linear in all
    // the variables
    pub fn replace(&mut self, var_id: usize, coeffs: &[T], ids: &[usize]) {
        // using the fact the the coefficient are sorted we can simply iterate over them
        // and update the values
        assert_ne!(
            var_id, 0,
            "var_id starts from 1 with 0 reserved for the constant"
        );
        // Check if the replace variable depends on itself
        let reduce = ids.iter().all(|&x| x != var_id);
        println!("REDUCE: {}  for {:?} with {}", reduce, ids , var_id);

        let n_var = if reduce { self.n_var - 1 } else { self.n_var };

        let mut new_pows: Vec<usize>;
        let mut new_coeff: T;
        let mut var_pow: usize;
        let mut coeff: T;
        let mut mpoly: MPolynomial<T>;
        let mut pos = 0;
        loop {
            println!("---");
            if pos >= self.powers.len() {
                break;
            }
            // store coefficient
            coeff = self.coeffs[pos];
            var_pow = self.powers[pos][var_id - 1];
            // Take the corresponding power of the transformation
            mpoly = MPolynomial::linear_pown(coeffs, ids, self.n_var, var_pow);
            println!("{:?}", mpoly);
            for (rep_pows, rep_coeff) in mpoly.powers.iter().zip(mpoly.coeffs.iter()) {
                println!("pos: {}", pos);
                new_pows = rep_pows
                    .iter()
                    .zip(self.powers[pos].iter())
                    .enumerate()
                    .map(|(n, (x, y))| if n == var_id - 1 { *x } else { x + y })
                    .collect();
                new_coeff = *rep_coeff * coeff;
                print!(
                    "\t{:?} -> {:?} :: {:?}",
                    self.powers[pos], new_pows, new_coeff
                );
                if new_pows[var_id - 1] == var_pow {
                    self.update(&new_pows, new_coeff);
                    print!(" UPDATE!");
                } else {
                    if self.add(&new_pows, new_coeff) {
                        print!(" NEW!");
                        pos += 1;
                    }
                }
                println!("");
            }
            if reduce && var_pow > 0{
                self.powers.remove(pos);
                self.coeffs.remove(pos);
            }else{
            pos += 1;
            }
        }
    }
}
