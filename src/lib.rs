extern crate f128;
extern crate itertools;
extern crate pest;
#[macro_use]
extern crate pest_derive;

use num_traits::{Float, FloatConst, FromPrimitive, Num, Signed, ToPrimitive, Zero};

use itertools::Itertools;

use num::traits::{NumAssign, NumOps, NumRef};
use num::NumCast;
use num::{BigInt, BigRational};
use std::fmt::{Debug, Display, LowerExp};
use std::iter::Sum;
use std::ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign};
//use utils::Signum;
//pub mod test;

// TODO: use it
pub const MAX_VARIABLE: usize = 7;

pub mod gcd;
pub mod parser;
pub mod rational;
pub mod utils;
use utils::{multinomial, next_combination_with_replacement};

pub trait FloatLike:
    From<f64>
    + Field
    + Num
    + FromPrimitive
    + ToPrimitive
    + Float
    + Signed
    + FloatConst
    + LowerExp
    + Debug
    + 'static
{
}
pub trait NumberLike
where
    Self: Field,
    //Self: Num,
    //Self: Zero,
    Self: PartialEq,
    Self: PartialOrd,
    Self: Signed,
    Self: Debug,
{
}
pub trait Field
where
    Self: Num,
    Self: Mul<Self, Output = Self>,
    Self: MulAssign<Self>,
    Self: AddAssign<Self>,
    Self: SubAssign<Self>,
    Self: Add<Self, Output = Self>,
    Self: Sub<Self, Output = Self>,
    Self: Neg<Output = Self>,
    Self: Sum<Self>,
    Self: Clone,
    Self: Debug,
    Self: Display,
    Self: ToPrimitive,
    Self: FromPrimitive,
{
}

pub trait RealNumberLike
where
    Self: Field,
    //Self: Num,
    Self: NumCast,
    Self: Float,
    Self: NumAssign,
    Self: NumOps,
    Self: NumRef,
{
}

impl Field for f64 {}
impl Field for f128::f128 {}
impl Field for BigRational {}
impl Field for BigInt {}

impl FloatLike for f64 {}
impl FloatLike for f128::f128 {}

impl NumberLike for f64 {}
impl NumberLike for f128::f128 {}
impl NumberLike for BigRational {}
impl NumberLike for BigInt {}

impl RealNumberLike for f64 {}
impl RealNumberLike for f128::f128 {}
impl<T: RealNumberLike> Field for num::Complex<T> {}

#[derive(Default, Debug, Clone)]
pub struct MPolynomialCache<T: Field> {
    pub powers: Vec<Vec<u16>>,
    pub coeffs: Vec<T>,
    pub size: usize,
    pub coeff_pows: Vec<u16>,
}
impl<T: Field> MPolynomialCache<T> {
    pub fn new(n_var: usize) -> MPolynomialCache<T> {
        MPolynomialCache {
            powers: vec![],
            coeffs: vec![],
            size: 0,
            coeff_pows: vec![0; n_var],
        }
    }
}

#[derive(Default, Debug, Clone)]
pub struct MPolynomial<T: Field> {
    //pub powers: Vec<(usize,Vec<usize>)>,
    pub powers: Vec<Vec<u16>>,
    pub coeffs: Vec<T>,
    pub n_var: usize,
    pub max_rank: Vec<u16>,
    pub cache: MPolynomialCache<T>,
}

impl<T: Field> MPolynomial<T> {
    /// Create an empty polynomial with n variables
    pub fn new(n_var: usize) -> MPolynomial<T> {
        assert!(n_var <= MAX_VARIABLE, "Increase the value of MAX_VARIABLE");
        MPolynomial {
            powers: vec![],
            coeffs: vec![],
            //powers: vec![vec![0; n_var]],
            //coeffs: vec![T::zero()],
            n_var: n_var,
            max_rank: vec![0; n_var],
            cache: MPolynomialCache::new(n_var),
        }
    }

    /// Create a zero polynomial wiht n variables
    pub fn zero(n_var: usize) -> MPolynomial<T> {
        assert!(n_var <= MAX_VARIABLE, "Increase the value of MAX_VARIABLE");
        MPolynomial {
            powers: vec![vec![0; n_var]],
            coeffs: vec![T::zero()],
            n_var: n_var,
            max_rank: vec![0; n_var],
            cache: MPolynomialCache::new(n_var),
        }
    }

    /// Create a monomial
    pub fn monomial(pows: &[u16], coeff: T) -> MPolynomial<T> {
        let mut monomial = MPolynomial::new(pows.len());
        monomial.add_coeff(pows, coeff);
        monomial
    }

    /// Create a polynomial equal to a constant 1 with n variables
    pub fn one(n_var: usize) -> MPolynomial<T> {
        assert!(n_var <= MAX_VARIABLE, "Increase the value of MAX_VARIABLE");
        MPolynomial {
            powers: vec![vec![0; n_var]],
            coeffs: vec![T::one()],
            n_var: n_var,
            max_rank: vec![0; n_var],
            cache: MPolynomialCache::new(n_var),
        }
    }

    /// Check if the polynomial is constant
    pub fn is_constant(&mut self) -> bool {
        self.drop_zeros();
        return self.coeffs.len() == 1 && self.powers[0].iter().all(|&x| x.is_zero());
    }

    /// Check if the polynomial is zero
    pub fn is_zero(&mut self) -> bool {
        return self.is_constant() && self.coeffs[0].is_zero();
    }

    /// Check if the polynomial is a monomial
    pub fn is_monomial(&mut self) -> bool {
        self.drop_zeros();
        return self.coeffs.len() == 1;
    }

    /// Return the highest non-zero coefficients and its powers according
    /// to the internal sorting
    pub fn leading_coeff(&self) -> Option<(&T, &Vec<u16>)> {
        self.coeffs
            .iter()
            .zip(self.powers.iter())
            .rev()
            .find(|(c, _)| *c != &T::zero())
    }

    /// Rescale all powers by the highest common monomial
    /// and return the powers corresponding powers
    /// Example:
    ///     x1*x2+x1 -> x2+1 (return [1,0])
    pub fn monomial_rescale(&mut self) -> Vec<u16> {
        let mut factor_powers = vec![u16::MAX; self.n_var];
        if !self.is_zero() {
            // Calling self::is_zero() also drops all zero
            // coefficients from the polynomial
            for pows in self.powers.iter() {
                for vn in 0..self.n_var {
                    if factor_powers[vn] > pows[vn] {
                        factor_powers[vn] = pows[vn];
                    }
                }
                if factor_powers.iter().all(|p| *p == 0) {
                    return factor_powers;
                }
            }
            for pows in self.powers.iter_mut() {
                for vn in 0..self.n_var {
                    pows[vn] -= factor_powers[vn];
                }
            }
        }
        factor_powers
    }

    /// Return the highest non-zero coefficients and its powers according
    /// to the internal sorting
    pub fn ranks_update(&mut self) {
        for vn in 0..self.n_var {
            self.max_rank[vn] = 0;
        }

        if !self.is_zero() {
            // Calling self::is_zero() also drops all zero
            // coefficients from the polynomial
            for pows in self.powers.iter() {
                for vn in 0..self.n_var {
                    if self.max_rank[vn] < pows[vn] {
                        self.max_rank[vn] = pows[vn];
                    }
                }
            }
        }
    }

    /// Add a new coefficient to the polynomial keeping the list sorted
    pub fn add_coeff(&mut self, pows: &[u16], coeff: T) -> bool {
        //match self.powers.binary_search(&(pows.len(),pows.clone())) {
        for (pow, c_pow) in pows.iter().zip(self.cache.coeff_pows.iter_mut()) {
            *c_pow = *pow;
        }
        match self.powers.binary_search(&self.cache.coeff_pows) {
            Ok(pos) => {
                self.coeffs[pos] += coeff;
                false
            }
            Err(pos) => {
                self.powers.insert(pos, self.cache.coeff_pows.clone());
                self.coeffs.insert(pos, coeff);
                for (pow, max_pow) in pows.iter().zip_eq(self.max_rank.iter_mut()) {
                    if pow > max_pow {
                        *max_pow = *pow;
                    }
                }
                true
            }
        }
    }

    /// Add a new coefficient to the polynomial keeping the list sorted
    pub fn update(&mut self, pows: &Vec<u16>, coeff: T) -> usize {
        //match self.powers.binary_search(&(pows.len(),pows.clone())) {
        match self.powers.binary_search(pows) {
            Ok(pos) => {
                self.coeffs[pos] = coeff;
                pos
            }
            Err(pos) => {
                self.powers.insert(pos, pows.clone());
                self.coeffs.insert(pos, coeff);
                for (pow, max_pow) in pows.iter().zip_eq(self.max_rank.iter_mut()) {
                    if pow > max_pow {
                        *max_pow = *pow;
                    }
                }
                pos
            }
        }
    }

    /// Clear the polynomial as if it was newly creted keeping the
    /// number of variables constant
    pub fn clear(&mut self) {
        self.powers.clear();
        self.coeffs.clear();
        for max_pow in self.max_rank.iter_mut() {
            *max_pow = 0;
        }
    }

    /// Clear the polynomial cache as if it was newly creted keeping the
    /// number of variables constant
    pub fn clear_cache(&mut self) {
        self.cache.coeffs.clear();
        self.cache.powers.clear();
        self.cache.size = 0;
    }

    /// Remove all the entries that are exactely zero and return a constant equal
    /// zero if all the coefficients are vanishing
    pub fn drop_zeros(&mut self) -> &Self {
        for i in (1..self.powers.len()).rev() {
            if self.coeffs[i] == T::zero() {
                self.coeffs.remove(i);
                self.powers.remove(i);
            }
        }

        if self.coeffs[0] == T::zero() {
            if self.coeffs.len() == 1 {
                for p in self.powers[0].iter_mut() {
                    *p = 0;
                }
            } else {
                self.coeffs.remove(0);
                self.powers.remove(0);
            }
        }
        self
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
        if n == 0 {
            mpoly.add_coeff(&pows, T::one());
            mpoly
        } else {
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
                    // TODO: merge these two
                    if ids[*v] > 0 {
                        pows[ids[*v] - 1] += 1;
                    }
                    coeff_pows[*v] += 1;
                }
                for (i, p) in coeff_pows[..=vars[n - 1]].iter().enumerate() {
                    coeff *= MPolynomial::powi(&coeffs[i], *p as i32);
                }
                // Multiply by multinomial
                mpoly.add_coeff(
                    &pows,
                    coeff * T::from_usize(multinomial::<usize>(&coeff_pows[..=n_var])).unwrap(),
                );
                if !next_combination_with_replacement(&mut vars, ids.len() - 1) {
                    break;
                }
            }
            mpoly
        }
    }

    /// Replace one of the variables by a multinomial which is at most linear in all
    /// the variables
    pub fn replace(&mut self, var_id: usize, coeffs: &[T], ids: &[usize]) -> &Self {
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
        let mut new_pows: Vec<u16>;
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
            coeff = self.coeffs[pos].clone();
            var_pow = self.powers[pos][var_id - 1] as usize;
            // Take the corresponding power of the transformation
            //mpoly = MPolynomial::linear_pown(coeffs, ids, self.n_var, var_pow);
            mpoly = MPolynomial::linear_pown(coeffs, ids, self.n_var, 1);
            mpoly.pown(var_pow);

            for (rep_pows, rep_coeff) in mpoly.powers.iter().zip(mpoly.coeffs.iter()) {
                //println!("pos: {}", pos);
                new_pows = rep_pows
                    .iter()
                    .zip(self.powers[pos].iter())
                    .enumerate()
                    .map(|(n, (x, y))| if n == var_id - 1 { *x } else { x + y })
                    .collect();
                new_coeff = rep_coeff.clone() * coeff.clone();
                //print!(
                //    "\t{:?} -> {:?} :: {:?}",
                //    self.powers[pos], new_pows, new_coeff
                //);
                // If we are replacing the original coefficient we need to replace the
                // coefficient instead then add the contribution
                if new_pows[var_id - 1] as usize == var_pow {
                    self.update(&new_pows, new_coeff);
                //print!(" UPDATE!");
                } else {
                    // Whenever a new element is added we need to take into account
                    // the shift in position
                    if self.add_coeff(&new_pows, new_coeff) {
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
        self
    }

    /// Store the current information of the polynomial into its cache
    /// This is used during multiplications in order to not lose any info
    pub fn to_cache(&mut self) {
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
            *c2 = c1.clone();
            for (p1_i, p2_i) in p1.iter().zip(p2.iter_mut()) {
                *p2_i = *p1_i;
            }
        }
    }

    /// Rescale all the coefficeints by a scalar
    pub fn scale(&mut self, scalar: &T) -> &Self {
        for c in self.coeffs.iter_mut() {
            *c *= scalar.clone()
        }
        self
    }

    /// Multiply the current polynomial with another and overwrite the original content
    pub fn mult(&mut self, other: &MPolynomial<T>) -> &Self {
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

        // TODO: use slice
        let mut pows = vec![0; self.n_var];
        let mut first = true;
        for (p2, c2) in other.powers.iter().zip(other.coeffs.iter()) {
            for n in 0..self.cache.size {
                for (i, p2_i) in p2.iter().enumerate() {
                    pows[i] = p2_i + self.cache.powers[n][i];
                }
                if first {
                    self.update(&pows, c2.clone() * self.cache.coeffs[n].clone());
                } else {
                    self.add_coeff(&pows, c2.clone() * self.cache.coeffs[n].clone());
                }
            }
            first = false;
        }
        self
    }

    pub fn square(&mut self) -> &Self {
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
            self.coeffs[n] = self.coeffs[n].clone() * self.cache.coeffs[0].clone();
        }
        //
        let mut pows = vec![0; self.n_var];
        for n1 in 0..self.cache.size {
            for n2 in 1..self.cache.size {
                for i in 0..self.n_var {
                    pows[i] = self.cache.powers[n1][i] + self.cache.powers[n2][i];
                }
                self.add_coeff(
                    &pows,
                    self.cache.coeffs[n1].clone() * self.cache.coeffs[n2].clone(),
                );
            }
        }
        self
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
                    self.add_coeff(
                        &pows,
                        self.cache.coeffs[n1].clone()
                            * self.cache.coeffs[n2].clone()
                            * self.cache.coeffs[n3].clone(),
                    );
                }
            }
        }
    }

    /// Try to factorize the given factor from the polynomial
    /// It does it by building the factorized polynomial by matching
    /// the highest coefficients in both
    /// e.g.
    ///     self:   10*x + 2*x*y
    ///     factor: x
    ///     result: 10*x + 2*x*y
    ///     STEP 1:
    ///          monomial : 2*y
    ///          result -= monomial * factor
    ///          result: 10*x
    ///     STEP 2:
    ///         monomial: 10
    ///         result -= monomial * factor
    ///         result: 0
    ///
    /// If the result==0 then the factorization is succesfull
    ///
    pub fn exact_division(&mut self, factor: &MPolynomial<T>) -> Result<(), &str> {
        self.clear_cache();
        let mut result = self.clone();

        let mut monomial = MPolynomial::<T>::new(self.n_var);
        let mut factorized = MPolynomial::<T>::new(self.n_var);

        // Extract from the factor polynomial the highest non-zero coefficient and powers
        let factor_last_non_zero_pos = factor.coeffs.len()
            - factor
                .coeffs
                .iter()
                .rev()
                .position(|c| !c.is_zero())
                .expect("Factor is Zero!")
            - 1;
        let highest_factor_monomial_powers = factor.powers[factor_last_non_zero_pos].clone();
        let highest_factor_monomial_coeff = factor.coeffs[factor_last_non_zero_pos].clone();

        //dbg!(&highest_factor_monomial_powers);

        //let mut step = 1;
        while !result.is_zero() {
            //dbg!(&step);
            //step += 1;

            // Check if we can factorize
            if !result.powers[result.coeffs.len() - 1]
                .iter()
                .zip(highest_factor_monomial_powers.iter())
                .all(|(p1, p2)| p1 >= p2)
            {
                return Err("Cannot Factorize!");
            }

            // Construct the monomial used for this step subtraction
            monomial.clear();
            monomial.powers.push(
                result.powers[result.coeffs.len() - 1]
                    .iter()
                    .zip(highest_factor_monomial_powers.iter())
                    .map(|(p1, &p2)| p1.checked_sub(p2).expect("Cannot Factorize!"))
                    .collect(),
            );
            monomial.coeffs.push(
                result.coeffs[result.coeffs.len() - 1].clone()
                    / highest_factor_monomial_coeff.clone(),
            );
            // Add monomial to factorized solution
            factorized += &monomial;

            // Subtract monomial * factor from result
            monomial.mult(factor);
            result -= &monomial;
            result.drop_zeros();
        }
        *self = factorized.clone();
        Ok(())
    }

    #[inline]
    pub fn powi(x: &T, n: i32) -> T {
        if n == 0 {
            T::one()
        } else {
            let r = MPolynomial::powi(x, n / 2);
            if n % 2 == 0 {
                r.clone() * r.clone()
            } else {
                r.clone() * r.clone() * x.clone()
            }
        }
    }

    /// Makes use of the binomial expression which is compute over u64 integers.
    /// This means that for large powers (>64) can potentially generate overflow
    pub fn pown2(&mut self, n: usize) -> &Self {
        if n == 0 {
            self.coeffs.resize(1, T::one());
            self.coeffs[0] = T::one();
            self.powers.resize(1, vec![0; self.n_var]);
            for p in self.powers[0].iter_mut() {
                *p = 0;
            }
            for max_pow in self.max_rank.iter_mut() {
                *max_pow = 0;
            }
            return self;
        }

        if self.coeffs.len() == 1 {
            self.coeffs[0] = MPolynomial::powi(&self.coeffs[0], n as i32);
            for p in self.powers[0].iter_mut() {
                *p *= n as u16;
            }
            return self;
        }

        // In oder to multiply our polynomial by another we need to store
        // its coefficients
        self.to_cache();
        self.coeffs.clear();
        self.powers.clear();
        //Shift all the entries by the first entry in other
        let mut pows = vec![0; self.n_var];
        let mut coeff = T::one();
        let mut v = vec![0; n];
        let mut pick: usize;
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
                    coeff = self.cache.coeffs[n].clone()
                } else {
                    for i in 0..self.n_var {
                        pows[i] += self.cache.powers[n][i]
                    }
                    coeff = coeff.clone() * self.cache.coeffs[n].clone()
                }
            }
            //println!("> {:?} | {:?}", &k[..=k_i], pows);
            self.add_coeff(
                &pows,
                coeff.clone() * T::from_usize(multinomial(&k[..=k_i])).unwrap(),
            );
            if !next_combination_with_replacement(&mut v, self.cache.size - 1) {
                break;
            }
        }
        self
    }

    pub fn pown3(&mut self, n: usize) -> &Self {
        if n == 0 {
            self.coeffs.resize(1, T::one());
            self.coeffs[0] = T::one();
            self.powers.resize(1, vec![0; self.n_var]);
            for p in self.powers[0].iter_mut() {
                *p = 0;
            }
            for max_pow in self.max_rank.iter_mut() {
                *max_pow = 0;
            }
        } else {
            let tmp = self.clone();
            for _ in 1..n {
                self.mult(&tmp);
            }
        }
        self
    }

    pub fn pown(&mut self, n: usize) -> &Self {
        match n {
            0 => {
                self.coeffs.resize(1, T::one());
                self.coeffs[0] = T::one();
                self.powers.resize(1, vec![0; self.n_var]);
                for p in self.powers[0].iter_mut() {
                    *p = 0;
                }
                for max_pow in self.max_rank.iter_mut() {
                    *max_pow = 0;
                }
            }
            1 => return self,
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
        self
    }
}

impl<'a, T: Field> Mul<&'a MPolynomial<T>> for &'a MPolynomial<T> {
    type Output = MPolynomial<T>;
    fn mul(self, other: Self) -> MPolynomial<T> {
        let mut out = self.clone();
        out.mult(other);
        out
    }
}

impl<'a, T: Field> MulAssign<&'a MPolynomial<T>> for MPolynomial<T> {
    fn mul_assign(&mut self, other: &'a Self) {
        self.mult(other);
    }
}

impl<'a, T: Field> MulAssign<T> for MPolynomial<T> {
    fn mul_assign(&mut self, other: T) {
        self.scale(&other);
    }
}

impl<'a, T: Field> Add<&'a MPolynomial<T>> for &'a MPolynomial<T> {
    type Output = MPolynomial<T>;
    fn add(self, other: Self) -> MPolynomial<T> {
        let mut out = self.clone();
        out += &other;
        out
    }
}

impl<'a, T: Field> AddAssign<&'a MPolynomial<T>> for MPolynomial<T> {
    fn add_assign(&mut self, other: &'a Self) {
        for (pow, c) in other.powers.iter().zip(other.coeffs.iter()) {
            self.add_coeff(pow, c.clone());
        }
    }
}

impl<'a, T: Field> SubAssign<&'a MPolynomial<T>> for MPolynomial<T> {
    fn sub_assign(&mut self, other: &Self) {
        for (pow, c) in other.powers.iter().zip(other.coeffs.iter()) {
            self.add_coeff(pow, -c.clone());
        }
    }
}

impl<'a, T: Field> Sub<&'a MPolynomial<T>> for &'a MPolynomial<T> {
    type Output = MPolynomial<T>;
    fn sub(self, other: Self) -> MPolynomial<T> {
        let mut out = self.clone();
        out -= &other;
        out
    }
}
