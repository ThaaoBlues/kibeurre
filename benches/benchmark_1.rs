use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use criterion_cycles_per_byte::CyclesPerByte;
use kibeurre::core::{decrypt_string, encrypt_string, generate_key_pair, generate_noise_polyvector};
use kibeurre::parameters::ETA_1;


const MAX_INPUT_LEN : usize = 100;

fn benchmark_encryption(c: &mut Criterion<CyclesPerByte>) {
    let mut group = c.benchmark_group("Encryption");


    let (pub_key, _) = generate_key_pair();
    
    for input_len in 1..=MAX_INPUT_LEN {
        let input_str = "A".repeat(input_len);
        
        group.bench_with_input(
            BenchmarkId::from_parameter(input_len), 

            &input_str, 
            |b, input| {
                // Measures both noise generation and encryption together

                b.iter(|| {
                    let r = generate_noise_polyvector(ETA_1);
                    let _ = encrypt_string(input, &pub_key, r);
                });
            }
        );
    }
    group.finish();
}

fn benchmark_decryption(c: &mut Criterion<CyclesPerByte>) {
    let mut group = c.benchmark_group("Decryption");

    let (pub_key, priv_key) = generate_key_pair();
    
    for input_len in 1..= MAX_INPUT_LEN {
        let input_str = "A".repeat(input_len);
        
        group.bench_with_input(
            BenchmarkId::from_parameter(input_len),
            &input_str,
            |b, input| {

                b.iter_with_setup(
                    || {
                        let r = generate_noise_polyvector(ETA_1);
                        encrypt_string(input, &pub_key, r)
                    },
                    |mut encrypted_string| {
                        let _ = decrypt_string(&mut encrypted_string, &priv_key);
                    }
                );
            }
        );
    }
    group.finish();
}


criterion_group!(
    name=benches;
    config=Criterion::default().with_measurement(CyclesPerByte);
    targets=benchmark_encryption, benchmark_decryption
);
criterion_main!(benches);