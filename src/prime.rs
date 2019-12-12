extern crate rand;
extern crate rayon;

use rand::distributions::{Distribution, Uniform};
use rayon::prelude::*;
use std::{mem, u64};

fn mod_pow(mut x: u64, mut exp: u64, n: u64) -> u64 {
    let mut res = 1;

    if n == 1 {
        return 0;
    }
    assert!(n - 1 < u64::MAX / (n - 1));
    x %= n;

    while exp != 0 {
        if exp & 1 != 0 {
            res = res * x % n;
        }
        exp >>= 1;
        x = x * x % n;
    }

    res
}

fn legendre_symbol(mut a: u64, mut n: u64) -> i64 {
    let mut res: i64 = 1;

    while a != 0 {
        while a & 1 == 0 {
            a /= 2;
            res = match n % 8 {
                3 | 5 => -res,
                _ => res,
            }
        }
        mem::swap(&mut a, &mut n);
        if a % 4 == 3 && n % 4 == 3 {
            res = -res;
        }
        a = a % n;
    }

    match n {
        1 => res,
        _ => 0,
    }
}

/// Test whether an integer `a` is a witness for the compositeness of `n`.
///
/// # Examples
///
/// ```
/// extern crate rand;
/// use rand::distributions::{Distribution, Uniform};
/// use solovay_strassen::prime;
///
/// let n: u64 = 27;
/// let dist = Uniform::new(2, n);
/// let mut rng = rand::thread_rng();
/// // A random integer in [2..n]
/// let a: u64 = dist.sample(&mut rng);
/// let res: bool = prime::solovay_strassen(a, n);
/// assert_eq!(res, true);
/// ```
pub fn solovay_strassen(a: u64, n: u64) -> bool {
    let x: i64 = legendre_symbol(a, n);

    if x == 0 || x.rem_euclid(n as i64) as u64 != mod_pow(a, (n - 1) / 2, n) {
        true
    } else {
        false
    }
}

/// Test whether an integer `n` is likely prime.
///
/// # Examples
///
/// ```
/// use solovay_strassen::prime;
///
/// // Mersenne Prime (2^31 - 1)
/// let n: u64 = 0x7FFF_FFFF;
/// // Try the solovay-strassen test 100 times in parallel
/// let res: bool = prime::is_prime(n, 100);
/// // Fails with a probability of at most 2_f64.pow(-100_f64)
/// assert_eq!(res, true);
/// ```
pub fn is_prime(n: u64, k: usize) -> bool {
    let dist = Uniform::new(2, n);
    let rng = rand::thread_rng();
    let samples: Vec<u64> = dist.sample_iter(rng).take(k).collect();

    !samples
        .par_iter()
        .find_any(|&&a| solovay_strassen(a, n))
        .is_some()
}

#[cfg(test)]
mod tests {
    const K: usize = 100;

    use super::*;
    use std::io;

    #[test]
    fn test_prime() -> io::Result<()> {
        let prime: u64 = 0x7FFF_FFFF;
        assert_eq!(is_prime(prime, K), true);
        Ok(())
    }

    #[test]
    fn test_composite() -> io::Result<()> {
        let composite: u64 = 0x7FFF_FFFE;
        assert_eq!(is_prime(composite, K), false);
        Ok(())
    }
}
