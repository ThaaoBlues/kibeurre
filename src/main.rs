pub mod math_utils;
pub mod ntt;
pub mod core;
pub mod parameters;
pub mod tui;
pub mod format_utils;
pub mod benchmark;



fn main() {
    let _ = tui::mainloop();

    
    // let mut result   = benchmark::benchmark_encryption(1000);
    // benchmark::results_to_csv(result,"encryption_benchmark.csv");

    // result = benchmark::benchmark_decryption(1000);
    // benchmark::results_to_csv(result,"decryption_benchmark.csv");



}
