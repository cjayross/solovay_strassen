//! # Quick Start
//!
//! Contained within this module are two functions:
//!   * `solovay_strassen`
//!   * `is_prime`
//!
//! The function `solovay_strassen` performs a single iteration of the Solovay-Strassen
//! primality test.
//!
//! On the other hand, `is_prime` is a routine that performs the Solovay-Strassen
//! primality test a given number of times in parallel, exiting as soon as the iterator
//! encounters a witness for the compositeness of the tested integer.

extern crate num_bigint as bigint;
extern crate num_traits as traits;
extern crate num_integer as integer;
extern crate rand;
extern crate rayon;

use bigint::*;
use rayon::prelude::*;
use std::iter::repeat_with;
use std::mem;
use traits::{One, Zero};
use integer::Integer;

macro_rules! bigint {
    ($e:expr) => {
        ($e).to_bigint().unwrap()
    };
}

macro_rules! biguint {
    ($e:expr) => {
        ($e).to_biguint().unwrap()
    };
}

fn legendre_symbol(mut a: BigUint, mut n: BigUint) -> BigInt {
    let mut res: BigInt = One::one();
    let (three, five) = (biguint!(3), biguint!(5));

    while a != Zero::zero() {
        while &a % 2u8 == Zero::zero() {
            a /= 2u8;
            if &n % 8u8 == three || &n % 8u8 == five {
                res = -res;
            }
        }
        mem::swap(&mut a, &mut n);
        if &a % 4u8 == three && &n % 4u8 == three {
            res = -res;
        }
        a = &a % &n;
    }

    if n == One::one() {
        res
    } else {
        Zero::zero()
    }
}

fn __solovay_strassen(a: &BigUint, n: &BigUint) -> bool {
    let x: BigInt = legendre_symbol(a.clone(), n.clone());
    x == Zero::zero() || x.mod_floor(&bigint!(n)) != bigint!(a.modpow(&((n - 1u8) / 2u8), n))
}

/// Test whether an integer `a` is a witness for the compositeness of `n`.
///
/// # Examples
///
/// ```
/// extern crate rand;
/// use rand::distributions::{Distribution, Uniform};
/// use solovay_strassen::solovay_strassen;
///
/// let n: u64 = 27;
/// let dist = Uniform::new(2, n);
/// let mut rng = rand::thread_rng();
/// // A random integer in [2..n]
/// let a: u64 = dist.sample(&mut rng);
/// assert!(solovay_strassen(&a, &n));
/// ```
pub fn solovay_strassen<T: ToBigUint>(a: &T, n: &T) -> bool {
    let (ref a, ref n) = (biguint!(a), biguint!(n));
    __solovay_strassen(a, n)
}

/// Test whether an integer `n` is likely prime.
///
/// # Examples
///
/// ```
/// use solovay_strassen::is_prime;
///
/// // Mersenne Prime (2^31 - 1)
/// let n: u64 = 0x7FFF_FFFF;
/// // Try the solovay-strassen test 100 times in parallel
/// // Fails with a probability of at most `2_f64.pow(-100_f64)`
/// assert!(is_prime(&n, 100));
/// ```
pub fn is_prime<T: ToBigUint>(n: &T, k: usize) -> bool {
    let n = biguint!(n);
    if n < biguint!(3) {
        return true;
    }

    let mut rng = rand::thread_rng();
    let samples: Vec<BigUint> = repeat_with(|| rng.gen_biguint(n.bits()))
        .filter(|m| m < &n)
        .take(k)
        .collect();

    !samples
        .par_iter()
        .find_any(|&a| __solovay_strassen(a, &n))
        .is_some()
}

#[cfg(test)]
mod tests {
    const K: usize = 1;

    use super::*;
    use std::io;

    #[test]
    fn test_prime() -> io::Result<()> {
        let prime: u64 = 0x7FFF_FFFF;
        assert!(is_prime(&prime, K));
        Ok(())
    }

    #[test]
    fn test_composite() -> io::Result<()> {
        let composite: u64 = 0x7FFF_FFFE;
        assert!(!is_prime(&composite, K));
        Ok(())
    }
}
