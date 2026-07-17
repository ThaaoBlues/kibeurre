use crate::{core::encrypt, math_utils::Vector};

pub mod math_utils;
pub mod ntt;
pub mod core;
pub mod parameters;
pub mod tui;
pub mod format_utils;
pub mod benchmark;



fn main() {
    //let _ = tui::mainloop();

    
    let mut results : Vec<(i32, f64)> = Vec::new();

    let (pub_key, priv_key) = core::generate_key_pair();

    let r = core::generate_noise_polyvector(parameters::ETA_1);
    
    let c = core::encrypt(pub_key.A, pub_key.t, r, Vector::new(&[1;256], parameters::m));
    // println!("u : {:?}", c.u);
    // println!("v : {:?}", c.v);

    let msg = core::decrypt(&c, priv_key.s);

    // println!("Decrypted message: {:?}", msg);

    // let mut result   = benchmark::benchmark_encryption(1000);
    // benchmark::results_to_csv(result,"encryption_benchmark.csv");

    // result = benchmark::benchmark_decryption(1000);
    // benchmark::results_to_csv(result,"decryption_benchmark.csv");



}
