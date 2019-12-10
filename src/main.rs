extern crate rand;
extern crate rayon;
extern crate clap;

use std::{u64, mem};
use clap::{Arg, App};
use rand::distributions::{Distribution, Uniform};
use rayon::prelude::*;

const K: usize = 100;

fn mod_pow(mut x: u64, mut exp: u64, n: u64) -> u64 {
  let mut res = 1;

  if n == 1 { return 0; }
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

fn solovay_strassen(a: u64, n: u64) -> bool {
  let x: i64 = legendre_symbol(a, n);

  if x == 0 || x.rem_euclid(n as i64) as u64 != mod_pow(a, (n - 1) / 2, n) {
    true
  }
  else {
    false
  }
}

fn is_prime(n: u64, k: usize) -> bool {
  let dist = Uniform::new(2, n);
  let rng = rand::thread_rng();
  let samples: Vec<u64> = dist.sample_iter(rng).take(k).collect();

  !samples.par_iter().find_any(|&&a| solovay_strassen(a, n)).is_some()
}

fn main() {
  let matches = App::new("Solovay-Strassen Test")
    .version("0.1.0")
    .author("cjayross <calvinjayross@gmail.com>")
    .about("Simple Solovay-Strassen algorithm for testing the primality of integers.")
    .arg(Arg::with_name("integer")
      .short("n")
      .long("integer")
      .takes_value(true)
      .help("Integer to test"))
    .get_matches();

  let n;
  let n_str = matches.value_of("integer");

  match n_str {
    None => panic!("Missing integer."),
    Some(s) => {
      match s.parse::<u64>() {
        Ok(_n) => n = _n,
        Err(_) => panic!("Invalid number, `{}`", s),
      }
    },
  }

  let result: bool = is_prime(n, K);

  println!("{} {}", n, match result { false => "is composite", _ => "is likely prime" });
}

#[cfg(test)]
mod tests {
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
