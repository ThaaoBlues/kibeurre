use crate::core::*;
use crate::ntt;
use crate::parameters::*;
use crate::format_utils::*;
use crate::math_utils::*;
use nist_pqc_seeded_rng::{NistPqcAes256CtrRng, Seed, SeedableRng};

#[cfg(test)]
mod tests {

use super::*;

    #[test]
    fn test_generate_seed_vector() {
        let seed_vector = generate_seed_vector();
        println!("Seed vector: {:?}", seed_vector);
    }

    #[test]
    #[allow(non_snake_case)]
    fn test_generate_A_from_seed() {
        let seed_vector = generate_seed_vector();
        let A = generate_A_from_seed(&seed_vector);
        println!("A: {:?}", A);
    
    }

    #[test]
    fn test_compute_t() {

        let mut nonce = 0;
        let sigma = Vec::new();
        let seed_vector = generate_seed_vector();
        #[allow(non_snake_case)]
        let A = generate_A_from_seed(&seed_vector);
        let s : PolyVector<k>;
        (s, nonce) = generate_noise_polyvector(ETA_1, nonce, &sigma);
        let e : PolyVector<k>;
        (e, _) = generate_noise_polyvector(ETA_2, nonce, &sigma);
        let t = compute_t(A, s, e);
        println!("t: {:?}", t);
    }

    #[test]
    fn test_compress_decompress() {
        let noise_vector = generate_noise_vector(ETA_1,0,&(0..32).collect()).0;
        let compressed = compress(noise_vector, D_U);
        let decompressed = decompress(compressed, D_U);
        println!("Original: {:?}", noise_vector);
        println!("Compressed: {:?}", compressed);
        println!("Decompressed: {:?}", decompressed);

    }



    #[test]
    fn test_encrypt_decrypt() {
        #[allow(non_snake_case)]
        let mut rng = NistPqcAes256CtrRng::from_seed([3u8; 48].into());
        let d = generate_d_and_z(&mut rng).0;
        let (rho,sigma) = generate_rho_and_sigma(d);
        let A = generate_A_from_seed(&rho); // generated in ntt form

        let mut nonce = 0;
        let mut s : PolyVector<k>;
        (s, nonce) = generate_noise_polyvector(ETA_1, nonce, &sigma);
        s.c = s.c.map(ntt::ntt);
        let e : PolyVector<k>;
        (e, nonce) = generate_noise_polyvector(ETA_1, nonce, &sigma);
        
        
        let t = compute_t(A, s, e);

        let hashed_pk = hash_public_key(&PublicKey { A, t, rho });




        let msg = Vector::new(&[1; 256],Q); // using seed vector as a message for testing


        let (kb,r) = derive_coins(msg.encode(1), hashed_pk);
        
        
        let encrypted_message = encrypt(A, t, msg,r);

        // at this stage, u and v are compressed 


        //println!("Encrypted message: u = {:?}, v = {:?}", encrypted_message.u, encrypted_message.v);
        let decrypted_msg = decrypt(&encrypted_message, s);

        assert_eq!(msg.c, decrypted_msg.c,"decryption failed");
    }




    use crate::format_utils::{KyberTestCase,parse_kyber_test_vectors};

    #[test]
    fn kat_test_hash_pk(){
        let test_cases : Vec<KyberTestCase> = parse_kyber_test_vectors();
        for tc in test_cases {
            println!("Running test case: {}", tc.count);

            let A  = generate_A_from_seed(&tc.rho);
            let pub_key = PublicKey { A, t: tc.t, rho: tc.rho };
            let h_pk = hash_public_key(&pub_key);
            assert_eq!(tc.pk_hash, h_pk,"hashes are not equal");
        }
    }

    #[test]
    fn kat_test_z_generation(){
        let test_cases : Vec<KyberTestCase> = parse_kyber_test_vectors();
        for tc in test_cases {
            println!("Running test case: {}", tc.count);

            let mut rng = NistPqcAes256CtrRng::from_seed(tc.rng_seed[..48].try_into().unwrap());
            let z = generate_d_and_z(&mut rng).1;
            assert_eq!(tc.z, z,"z values are not equal");
        }
    }



    #[test]
    fn kat_test_t_and_s_generation(){
        let test_cases : Vec<KyberTestCase> = parse_kyber_test_vectors();
        for tc in test_cases {
            println!("Running test case: {}", tc.count);

            let A  = generate_A_from_seed(&tc.rho);
            let mut nonce = 0;
            let mut s : PolyVector<k>;
            let mut rng = NistPqcAes256CtrRng::from_seed(tc.rng_seed[..48].try_into().unwrap());
            let d = generate_d_and_z(&mut rng).0;
            let (_,sigma) = generate_rho_and_sigma(d);


            (s,nonce) = generate_noise_polyvector(ETA_1,nonce,&sigma);
            
            let mut sk_intt = empty_polyvector();
            sk_intt.c = tc.sk.c.map(ntt::intt);

            assert_eq!(sk_intt, s,"s normal values are not equal");
            s.c = s.c.map(ntt::ntt);
            assert_eq!(tc.sk, s,"s ntt values are not equal");

            let e : PolyVector<k>;
            (e, nonce) = generate_noise_polyvector(ETA_1, nonce, &sigma);
            let t: PolyVector<k> = compute_t(A, s, e);

            assert_eq!(tc.t, t,"t values are not equal");
        }
    }




    #[test]
    fn kat_test_decapsulation(){

        let test_cases : Vec<KyberTestCase> = parse_kyber_test_vectors();
        for tc in test_cases {
            //println!("Running test case: {:?}", tc);

            let A  = generate_A_from_seed(&tc.rho);

            let mut encrypted_message = EncryptedMessage { u: tc.ct.0, v: tc.ct.1 };
            let pub_key = PublicKey { A, t: tc.t, rho: tc.rho };
            let mut rng = NistPqcAes256CtrRng::from_seed(tc.rng_seed[..48].try_into().unwrap());
            let z = generate_d_and_z(&mut rng).1;


            //good up until here
            // hash_public_key is correct
            let K = decapsulate( encrypted_message, 
                PrivateKey { s: tc.sk,
                    pk: pub_key.clone(),
                    hashed_pk: hash_public_key(&pub_key),
                    z
                }
            );
            assert_eq!(tc.ss, K,"decryption failed");
        }
    }
    
    #[test]
    fn run_enc_dec_test_on_known_vectors(){

        let test_cases = parse_kyber_test_vectors();
        for tc in test_cases {
            //println!("Running test case: {:?}", tc);
            println!("============= TESTING ENCAPSULATION ==============");

            let A  = generate_A_from_seed(&tc.rho);
            let pk = PublicKey { A, t: tc.t, rho: tc.rho };
            let mut rng = NistPqcAes256CtrRng::from_seed(tc.rng_seed[..48].try_into().unwrap());
            let (msg, msg_bytes) = generate_message_from_seed(&mut rng);
            let (encrypted_message, K) = encapsulate(&pk, msg);

            println!("{:x?}",tc.ct.0.encode(D_U as usize));

            assert_eq!(encrypted_message.u, tc.ct.0, "encryption failed for u");
            assert_eq!(encrypted_message.v, tc.ct.1, "encryption failed for v");

            println!("============= TESTING DECAPSULATION ==============");
            let ref_enc_msg = EncryptedMessage { u: tc.ct.0, v: tc.ct.1 };
            let K = decapsulate(ref_enc_msg, 
                PrivateKey { s: tc.sk,
                    pk: pk.clone(),
                    hashed_pk: hash_public_key(&pk),
                    z: generate_d_and_z(&mut rng).1
                }
            );
            assert_eq!(tc.ss, K,"decryption failed");
        }


    }
}


