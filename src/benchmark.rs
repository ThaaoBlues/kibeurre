use crate::core::{decrypt_string, encrypt_string, generate_key_pair, generate_noise_polyvector};
use crate::parameters::{eta_1};


pub fn benchmark_encryption(max_input_len : i32) -> Vec<(i32, f64)> {
    let mut results = Vec::new();
    let (pub_key, _) = generate_key_pair();

    for input_len in 1..=max_input_len {
        let input_str = "A".repeat(input_len as usize);
        let start_time = std::time::Instant::now();
        let r = generate_noise_polyvector(eta_1);
        let _ = encrypt_string(&input_str, &pub_key,r);
        let duration = start_time.elapsed();
        results.push((input_len, duration.as_secs_f64()));
    }
    results
}


pub fn benchmark_decryption(max_input_len : i32)->Vec<(i32,f64)>{


    let mut results = Vec::new();

    let (pub_key, priv_key) = generate_key_pair();

    let mut encrypted_strings = Vec::new();
    for input_len in 1..=max_input_len {
        let input_str = "A".repeat(input_len as usize);
        let r = generate_noise_polyvector(eta_1);
        encrypted_strings.push(encrypt_string(&input_str, &pub_key,r));
    }



    for input_len  in 0..max_input_len{
        let start_time = std::time::Instant::now();
        
        decrypt_string(&mut encrypted_strings[input_len as usize], &priv_key);

        let duration = start_time.elapsed();

        results.push((input_len,duration.as_secs_f64()));
    }

    return results;
    
}



pub fn results_to_csv(results : Vec<(i32,f64)>,filename : &str){
    let mut wtr = csv::Writer::from_path(filename).unwrap();
    wtr.write_record(&["input_length", "time_taken"]).unwrap();
    for (input_length, time_taken) in results {
        wtr.write_record(&[input_length.to_string(), time_taken.to_string()]).unwrap();
    }
    wtr.flush().unwrap();
}