use crate::math_utils::{Vector,empty_vector,PolyVector,empty_polyvector};
use crate::parameters::{n,k};
use std::fs::File;
use std::io::{BufRead, BufReader};

/*
Converts a string into a vector of Vector<n> where each Vector<n> represents a chunk of the string.
The string is split into bytes chunks of size n.
So each vector actually contains n bytes of the string.
*/
pub fn string_to_vectors(input: &str) -> Vec<Vector<n>> {
    let bytes = input.as_bytes();
    let mut vectors = Vec::new();


    for chunk in bytes.chunks(n/8) {
        let mut vector = empty_vector(); // filled with zeros, will be used as padding
        for (i, &byte) in chunk.iter().enumerate() {

            for j in 0..8 {
                let bit = (byte >> j) & 1;
                vector.set(i * 8 + j, bit as i32);
            }
        }
        vectors.push(vector);
    }

    
    vectors
}




pub fn vectors_to_string(vectors: Vec<Vector<n>>) -> String {
    let mut bytes = Vec::new();
    for vector in vectors {

        for i in 0..(n/8) {
            let mut byte = 0;
        
            
            for j in 0..8 {
                
                let bit = vector.get(i * 8 + j) as u8;
                if bit != 0 {
                    byte |= bit << j;
                }
            }

            // NULL byte => string end
            // do not append it as String constructor will
            if byte == 0 {
                break;
            }else{
                bytes.push(byte);
            }
            
        }
    }

    String::from_utf8(bytes).unwrap_or_else(|_| String::from("Invalid UTF-8"))
}







/*
Functions used to parse KAT test vectors from the NIST test files.
*/

// coin to derivate all the randomness needed for the key generation and encryption
// s/o fujisaki okamoto transformation
pub fn parse_AES_CTR_DRBG_seed(string: &str) -> Vec<u8> {
    hex::decode(string).expect("Invalid hex string")
}

fn parse_seed(seed_str: &str) -> Vec<u8> {

    // 32 bytes
    let seed_bytes = hex::decode(seed_str).expect("Invalid hex string");
    // let mut seed_vector = empty_vector();



    // for (i, &byte) in seed_bytes.iter().enumerate() {
    //     for j in 0..8 {
    //         let bit = (byte >> j) & 1;
    //         seed_vector.set(i * 8 + j, bit as i32);
    //     }
    // }


    seed_bytes
}


fn parse_compressed_polyvector(s : &str)->PolyVector<k>{

    let bytes = hex::decode(s).expect("Invalid hex string");
    let mut vector = empty_polyvector();

    let mut tmp = empty_vector();

    let mut poly_index = 0; // index of polynomial being parsed
    let mut rel_coef_index = 0; // coefficient index in polynomial we are currently parsing
    for byte_index in (0..bytes.len()).step_by(5) {

        // 5*8=40 => easy to parse 4 coefficients as they are 10 bits each
        // LITTLE ENDIAN converted to BIG ENDIAN
        let coef_1: i32  = (bytes[byte_index] as i32) | (((bytes[byte_index+1] & 0x03) as i32) << 8);
        let coef_2: i32 = ((bytes[byte_index+1] & 0xFC) >> 2) as i32 | (((bytes[byte_index+2] & 0x0F) as i32) << 6);
        let coef_3: i32 = ((bytes[byte_index+2] & 0xF0) >> 4) as i32 | (((bytes[byte_index+3] & 0x3F) as i32) << 4);
        let coef_4: i32 = ((bytes[byte_index+3] & 0xC0) >> 6) as i32 | ((bytes[byte_index+4]  as i32) << 2);

        tmp.c[rel_coef_index] = coef_1;
        tmp.c[rel_coef_index + 1] = coef_2;
        tmp.c[rel_coef_index + 2] = coef_3;
        tmp.c[rel_coef_index + 3] = coef_4;

        // 256 is not odd
        rel_coef_index += 4;

        if rel_coef_index == 256 {
            vector.c[poly_index] = tmp;
            poly_index += 1;
            rel_coef_index = 0;
            tmp = empty_vector();
        }


    }
    vector
}
fn parse_polyvector(s : &str)->PolyVector<k>{

    // t is represented as hex string of 3*256*12 bits long

    let bytes = hex::decode(s).expect("Invalid hex string");
    let mut vector = empty_polyvector();

    let mut tmp = empty_vector();

    let mut poly_index = 0; // index of polynomial being parsed
    let mut rel_coef_index = 0; // coefficient index in polynomial we are currently parsing


    for byte_index in (0..bytes.len()).step_by(3) {

        // 3*8 = 24 => easy to parse 2 coefficients on 3 bytes as they are 12 bits each
        let coef_1: i32 = (bytes[byte_index] as i32) | (((bytes[byte_index + 1] & 0x0F) as i32) << 8);
        let coef_2: i32 = ((bytes[byte_index + 1] >> 4) as i32) | ((bytes[byte_index + 2] as i32) << 4);
        
        tmp.c[rel_coef_index] = coef_1;

        // 256 is not odd, so no need to check polynome boundaries at first coefficient

        tmp.c[rel_coef_index+1] = coef_2;

        rel_coef_index += 2;

        // 256*12 = 3072 bits => 384 bytes, so we can check if we are at the end of a polynomial
        if rel_coef_index == 256 {
            vector.c[poly_index] = tmp;
            poly_index += 1;
            rel_coef_index = 0;
            tmp = empty_vector();
        }


    }
    vector
}

pub fn parse_polyvector_bytes(bytes : &[u8])->PolyVector<k>{

    let mut vector = empty_polyvector();

    let mut tmp = empty_vector();

    let mut poly_index = 0; // index of polynomial being parsed
    let mut rel_coef_index = 0; // coefficient index in polynomial we are currently parsing


    for byte_index in (0..bytes.len()).step_by(3) {

        // 3*8 = 24 => easy to parse 2 coefficients on 3 bytes as they are 12 bits each
        let coef_1: i32 = (bytes[byte_index] as i32) | (((bytes[byte_index + 1] & 0x0F) as i32) << 8);
        let coef_2: i32 = ((bytes[byte_index + 1] >> 4) as i32) | ((bytes[byte_index + 2] as i32) << 4);
        
        tmp.c[rel_coef_index] = coef_1;

        // 256 is not odd, so no need to check polynome boundaries at first coefficient

        tmp.c[rel_coef_index+1] = coef_2;

        rel_coef_index += 2;

        // 256*12 = 3072 bits => 384 bytes, so we can check if we are at the end of a polynomial
        if rel_coef_index == 256 {
            vector.c[poly_index] = tmp;
            poly_index += 1;
            rel_coef_index = 0;
            tmp = empty_vector();
        }


    }
    vector
}

fn parse_compressed_vector(string : &str)->Vector<256>{
    // 4 bit per coefficient

    let mut rel_coef_index: usize = 0; // coefficient index in polynomial we are currently parsing

    let bytes: Vec<u8> = hex::decode(string).expect("Invalid hex string");
    let mut v: Vector<256> = empty_vector();

    for byte in bytes {

        
        // 8 bits per byte, 4 bits per coefficient => 2 coefficients per byte


        let coef_1 = ((byte & 0x0F) as i32);
        let coef_2 = ((byte >> 4) as i32);


        v.c[rel_coef_index] = coef_1;
        v.c[rel_coef_index + 1] = coef_2;

        rel_coef_index += 2;
    }

    v
}





pub fn parse_public_key(string : &str)->(Vec<u8>,PolyVector<k>){
    // parse polymatrix seed A then parse t or the other way around i don't remember
    // t is first then seed A is second, both are represented as hex strings

    // a char represents 4 bits even if coded on 8 bits
    let t = parse_polyvector(string.trim().chars().take(k*256*12/4).collect::<String>().as_str());
    
    let seed = parse_seed(string.trim().chars().skip(k*256*12/4).collect::<String>().as_str());
    (seed,t)
}

pub fn parse_cyphertext(string : &str)->(PolyVector<k>, Vector<n>){
    // parse u then v, both are represented as hex strings
    let u = parse_compressed_polyvector(string.trim().chars().take(k*256*10/4).collect::<String>().as_str());
    let v = parse_compressed_vector(string.trim().chars().skip(k*256*10/4).collect::<String>().as_str());
    (u,v)
}


pub fn parse_shared_secret(string : &str)->Vec<u8>{

    let bytes: Vec<u8> = hex::decode(string).expect("Invalid hex string");
    bytes
}





pub fn parse_secret_key(string : &str)->PolyVector<k>{

    parse_polyvector(string.trim().chars().take(k*256*12/4).collect::<String>().as_str())

}

pub fn parse_pk_hash(string : &str)->Vec<u8>{
    // s | t | rho | H(pk) | z
    let hash_str = string.trim().chars().skip((k*256*12/4)*2+32*2).take(32*2).collect::<String>();
    let bytes: Vec<u8> = hex::decode(hash_str).expect("Invalid hex string");
    bytes
}

pub fn parse_z(string : &str)->Vec<u8>{
    // s | t | s | H(pk) | z
    let z_str = string.trim().chars().skip((k*256*12/4)*2+32*2+32*2).take(32*2).collect::<String>();
    let bytes: Vec<u8> = hex::decode(z_str).expect("Invalid hex string");
    bytes
}

#[derive(Debug, Clone)]
pub struct KyberTestCase {
    pub count: usize,

    pub rho: Vec<u8>,
    pub t: PolyVector<k>,
    pub sk: PolyVector<k>,
    pub pk_hash: Vec<u8>,
    pub z : Vec<u8>,
    pub rng_seed: Vec<u8>,
    // u,v 
    pub ct: (PolyVector<k>, Vector<n>),
    // derived key from m
    pub ss: Vec<u8>,
}

pub fn parse_kyber_test_vectors()->Vec<KyberTestCase>{
    let file = File::open("src/PQCkemKAT_2400.rsp").expect("Failed to open test vector file");
    let reader = BufReader::new(file);

    let mut test_cases = Vec::new();
    let mut current_case : KyberTestCase;

    current_case = KyberTestCase {
        count: 0,
        rho: Vec::new(),
        rng_seed: Vec::new(),
        t: empty_polyvector(),
        sk: empty_polyvector(),
        pk_hash: Vec::new(),
        z : Vec::new(),
        ct: (empty_polyvector(), empty_vector()),
        ss: Vec::new(),
    };
    for line in reader.lines() {
        let line = line.expect("Failed to read line");
        let trimmed = line.trim();

        // skip comments and empty lines
        if trimmed.starts_with('#') || trimmed.is_empty() {
            continue;
        }


        // split line into key and value
        if let Some((key, value)) = trimmed.split_once('=') {
            let key = key.trim();
            let value = value.trim().to_string();
            
            match key {
                "count" => {

                    current_case.count = value.parse().unwrap_or(0);
                }
                "seed" => current_case.rng_seed = parse_AES_CTR_DRBG_seed(&value),
                "pk" => {
                    let (seed, t) = parse_public_key(&value);
                    current_case.rho = seed;
                    current_case.t = t;

                },
                "sk" => {
                    current_case.sk = parse_secret_key(&value);
                    current_case.pk_hash = parse_pk_hash(&value);
                    current_case.z = parse_z(&value);
                },
                "ct" => current_case.ct = parse_cyphertext(&value),
                "ss" => {
                    current_case.ss = parse_shared_secret(&value);
                    test_cases.push(current_case.clone());
                },
                _ => {} // Ignore unknown keys
            }
        }
    }

    test_cases

}