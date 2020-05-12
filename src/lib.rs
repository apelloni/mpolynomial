extern crate f128;

use num_traits::{Float, FloatConst, FromPrimitive, Num, One, Signed, ToPrimitive, Zero};
//use utils::Signum;
//pub mod test;

// TODO: use it
pub const MAX_VARIABLE: usize = 7;

pub mod utils;
use utils::{binomial, multinomial, next_combination_with_replacement};

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

#[derive(Default, Debug, Clone)]
pub struct MPolynomialCache<
    T: From<f64>
        + Num
        + FromPrimitive
        + Float
        + Signed
        + FloatConst
        + std::fmt::LowerExp
        + std::fmt::Debug
        + 'static,
>
//<T: FloatLike>
{
    pub powers: Vec<Vec<usize>>,
    pub coeffs: Vec<T>,
    size: usize,
}
impl<
        T: From<f64>
            + Num
            + FromPrimitive
            + Float
            + Signed
            + FloatConst
            + std::fmt::LowerExp
            + std::fmt::Debug
            + 'static,
    > MPolynomialCache<T>
{
    pub fn new() -> MPolynomialCache<T> {
        MPolynomialCache {
            powers: vec![],
            coeffs: vec![],
            size: 0,
        }
    }
}

#[derive(Default, Debug, Clone)]
pub struct MPolynomial<
    T: From<f64>
        + Num
        + FromPrimitive
        + Float
        + Signed
        + FloatConst
        + std::fmt::LowerExp
        + std::fmt::Debug
        + 'static,
>
//<T: FloatLike>
{
    //pub powers: Vec<(usize,Vec<usize>)>,
    pub powers: Vec<Vec<usize>>,
    pub coeffs: Vec<T>,
    n_var: usize,
    max_rank: usize,
    pub cache: MPolynomialCache<T>,
}

impl<
        T: From<f64>
            + Num
            + FromPrimitive
            + Float
            + Signed
            + FloatConst
            + std::fmt::LowerExp
            + std::fmt::Debug
            + 'static,
    > MPolynomial<T>
{
    pub fn new(n_var: usize) -> MPolynomial<T> {
        assert!(n_var <= MAX_VARIABLE, "Increase the value of MAX_VARIABLE");
        MPolynomial {
            powers: vec![],
            coeffs: vec![],
            n_var: n_var,
            max_rank: 0,
            cache: MPolynomialCache::new(),
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
        let mut coeff_pows = [0; MAX_VARIABLE + 1];
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
                if ids[*v] > 0 {
                    pows[ids[*v] - 1] += 1;
                }
                coeff_pows[ids[*v]] += 1;
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

    /// Replace one of the variables by a multinomial which is at most linear in all
    /// the variables
    pub fn replace(&mut self, var_id: usize, coeffs: &[T], ids: &[usize]) {
        // using the fact the the coefficient are sorted we can simply iterate over them
        // and update the values
        assert_ne!(
            var_id, 0,
            "var_id starts from 1 with 0 reserved for the constant"
        );
        // Check if the replace variable depends on itself
        let reduce = ids.iter().all(|&x| x != var_id);
        //println!("REDUCE: {}  for {:?} with {}", reduce, ids , var_id);

        // Containers
        let mut new_pows: Vec<usize>;
        let mut new_coeff: T;
        let mut var_pow: usize;
        let mut coeff: T;
        let mut mpoly: MPolynomial<T>;
        let mut pos = 0;
        loop {
            //println!("---");
            if pos >= self.powers.len() {
                break;
            }
            // store coefficient
            coeff = self.coeffs[pos];
            var_pow = self.powers[pos][var_id - 1];
            // Take the corresponding power of the transformation
            mpoly = MPolynomial::linear_pown(coeffs, ids, self.n_var, var_pow);
            for (rep_pows, rep_coeff) in mpoly.powers.iter().zip(mpoly.coeffs.iter()) {
                //println!("pos: {}", pos);
                new_pows = rep_pows
                    .iter()
                    .zip(self.powers[pos].iter())
                    .enumerate()
                    .map(|(n, (x, y))| if n == var_id - 1 { *x } else { x + y })
                    .collect();
                new_coeff = *rep_coeff * coeff;
                //print!(
                //    "\t{:?} -> {:?} :: {:?}",
                //    self.powers[pos], new_pows, new_coeff
                //);
                // If we are replacing the original coefficient we need to replace the
                // coefficient instead then add the contribution
                if new_pows[var_id - 1] == var_pow {
                    self.update(&new_pows, new_coeff);
                //print!(" UPDATE!");
                } else {
                    // Whenever a new element is added we need to take into account
                    // the shift in position
                    if self.add(&new_pows, new_coeff) {
                        //print!(" NEW!");
                        pos += 1;
                    }
                }
                //println!("");
            }
            // When we have a reduction we need to get rid of all the
            // coefficients that now don't exist anymore
            if reduce && var_pow > 0 {
                self.powers.remove(pos);
                self.coeffs.remove(pos);
            } else {
                pos += 1;
            }
        }
    }

    fn to_cache(&mut self) {
        // If cache is too small then resize
        if self.cache.coeffs.len() < self.coeffs.len() {
            self.cache
                .powers
                .resize(self.powers.len(), vec![0; self.n_var]);
            self.cache.coeffs.resize(self.coeffs.len(), T::zero());
        }
        // Store actual used size
        self.cache.size = self.coeffs.len();

        // Store coefficients
        for ((p1, c1), (p2, c2)) in self.powers.iter().zip(self.coeffs.iter()).zip(
            self.cache
                .powers
                .iter_mut()
                .zip(self.cache.coeffs.iter_mut()),
        ) {
            *c2 = *c1;
            for (p1_i, p2_i) in p1.iter().zip(p2.iter_mut()) {
                *p2_i = *p1_i;
            }
        }
    }

    pub fn mult(&mut self, other: &MPolynomial<T>) {
        // In oder to multiply our polynomial by another we need to store
        // its coefficients
        self.to_cache();
        //self.coeffs.clear();
        //self.powers.clear();
        //Shift all the entries by the first entry in other
        for p1 in self.powers.iter_mut() {
            for (p1_i, p2_i) in p1.iter_mut().zip(other.powers[0].iter()) {
                *p1_i += p2_i;
            }
        }

        let mut pows = vec![0; self.n_var];
        let mut first = true;
        for (p2, &c2) in other.powers.iter().zip(other.coeffs.iter()) {
            for n in 0..self.cache.size {
                for (i, p2_i) in p2.iter().enumerate() {
                    pows[i] = p2_i + self.cache.powers[n][i];
                }
                if first {
                    self.update(&pows, c2 * self.cache.coeffs[n]);
                } else {
                    self.add(&pows, c2 * self.cache.coeffs[n]);
                }
            }
            first = false;
        }
    }
    pub fn square(&mut self) {
        // In oder to multiply our polynomial by another we need to store
        // its coefficients
        self.to_cache();
        //self.coeffs.clear();
        //self.powers.clear();
        //Shift all the entries by the first entry in other
        for n in 0..self.cache.size {
            for (p1_i, p2_i) in self.powers[n].iter_mut().zip(self.cache.powers[0].iter()) {
                *p1_i += p2_i;
            }
            self.coeffs[n] = self.coeffs[n] * self.cache.coeffs[0];
        }
        //
        let mut pows = vec![0; self.n_var];
        for n1 in 0..self.cache.size {
            for n2 in 1..self.cache.size {
                for i in 0..self.n_var {
                    pows[i] = self.cache.powers[n1][i] + self.cache.powers[n2][i];
                }
                self.add(&pows, self.cache.coeffs[n1] * self.cache.coeffs[n2]);
            }
        }
    }

    pub fn cube(&mut self) {
        // In oder to multiply our polynomial by another we need to store
        // its coefficients
        self.to_cache();
        self.coeffs.clear();
        self.powers.clear();
        //Shift all the entries by the first entry in other
        let mut pows = vec![0; self.n_var];
        for n1 in 0..self.cache.size {
            for n2 in 0..self.cache.size {
                for n3 in 0..self.cache.size {
                    for i in 0..self.n_var {
                        pows[i] = self.cache.powers[n1][i]
                            + self.cache.powers[n2][i]
                            + self.cache.powers[n3][i];
                    }
                    self.add(
                        &pows,
                        self.cache.coeffs[n1] * self.cache.coeffs[n2] * self.cache.coeffs[n3],
                    );
                }
            }
        }
    }
    pub fn pown2(&mut self, n: usize) {
        // In oder to multiply our polynomial by another we need to store
        // its coefficients
        self.to_cache();
        self.coeffs.clear();
        self.powers.clear();
        //Shift all the entries by the first entry in other
        let mut pows = vec![0; self.n_var];
        let mut coeff = T::one();
        let mut v = vec![0; n];
        let mut pick = 0;
        let mut k_i: usize;
        let mut k = vec![0; n];
        loop {
            // reset values
            for i in 0..n {
                k[i] = 0;
            }
            pick = v[0];
            k_i = 0;

            // read out infos
            for (which, &n) in v.iter().enumerate() {
                if pick != n {
                    pick = n;
                    k_i += 1;
                }
                k[k_i] += 1;
                if which == 0 {
                    for i in 0..self.n_var {
                        pows[i] = self.cache.powers[n][i]
                    }
                    coeff = self.cache.coeffs[n]
                } else {
                    for i in 0..self.n_var {
                        pows[i] += self.cache.powers[n][i]
                    }
                    coeff = coeff * self.cache.coeffs[n]
                }
            }
            //println!("> {:?} | {:?}", &k[..=k_i], pows);
            self.add(
                &pows,
                coeff * T::from_usize(multinomial(&k[..=k_i])).unwrap(),
            );
            if !next_combination_with_replacement(&mut v, self.cache.size - 1) {
                break;
            }
        }
    }

    pub fn pown(&mut self, n: usize) {
        match n {
            0 => {
                self.coeffs.resize(1, T::one());
                self.powers.resize(1, vec![0; self.n_var]);
                self.max_rank = 0;
            }
            1 => return,
            _ => {
                if n % 2 == 0 {
                    self.pown(n / 2);
                    self.square();
                } else {
                    let old = self.clone();
                    self.pown(n / 2);
                    self.square();
                    self.mult(&old);
                }
            }
        }
    }
}
