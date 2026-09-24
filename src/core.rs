
use crate::math_utils::{PolyMatrix,PolyVector,Vector,empty_polymatrix,empty_polyvector,empty_vector};

use crate::ntt;
use rand::{Rng, random};
use shake::{ExtendableOutput, Update, XofReader,Shake256};
use crate::parameters::{D_U, D_V, ETA_1, ETA_2, k, N, Q};
use crate::format_utils::{parse_polyvector_bytes, string_to_vectors, vectors_to_string};
use nist_pqc_seeded_rng::{NistPqcAes256CtrRng, Seed, SeedableRng};
use sha3::{Digest, Sha3_256,Sha3_512};


fn cbd_sample_from_bits(bytes: &[u8], eta: usize, coeff_index: usize) -> i32 {
    // bytes: bit buffer from PRF output, LSB-first, 2*eta bits per coefficient
    let mut a = 0i32;
    let mut b = 0i32;
    let base = coeff_index * 2 * eta;
    for i in 0..eta {
        a += extract_bit(bytes, base + i) as i32;
    }
    for i in 0..eta {
        b += extract_bit(bytes, base + eta + i) as i32;
    }
    a - b
}


fn extract_bit(bytes: &[u8], bit_index: usize) -> u8 {
    (bytes[bit_index / 8] >> (bit_index % 8)) & 1
}



pub fn generate_seed_vector() -> Vec<u8>{

    let r  :[u8; 32] = [0; 32];
    let r = r.map(|_| rand::random::<u8>());

    r.to_vec()

}

#[allow(non_snake_case)]
pub fn generate_A_from_seed(seed : &Vec<u8>) -> PolyMatrix<k,k>{

    /*
    
    Generates the NTT version of the matrix A from a seed vector.
    /!\ Regarding the specification, this is A_hat (hashing coefs order j,i)
    (used in KeyGen and compute_t functions)
    NOT A_hat^T (hashing coefs order i,j)
    (Used in encrypt function)

    */

    // A is used for debugging and will be removed,
    // A_ntt is used for the actual encryption/decryption process
    
    //let mut A : PolyMatrix<k,k> = empty_polymatrix();
    let mut A_ntt : PolyMatrix<k,k> = empty_polymatrix();
    let mut generated_polynomial : Vector<N> = empty_vector();


    for i in 0..k {
        for j in 0..k {
            let mut hasher = shake::Shake128::default();

            // 2. Feed data into the hasher
            // log2(q) = 12 bits par coef de polynome ? 
            hasher.update(seed);

            // beware of the order of the indices (yes I made the mistake)
            hasher.update(&[j as u8, i as u8]);
            // log2(q) round = 12 bits => u16
            // one polynomial = 12*256 bits = 12*64 bytes => 2*u8*256

            let mut reader = hasher.finalize_xof();
      
            let mut buf: [u8;3] = [0;3]; // 24 bits =>  2 12 bits coefficients

            // Kyber assume reading LSB-FIRST
            let mut l : usize = 0;
            while l < N {

                reader.read(&mut buf);

                let c1 : i32  = buf[0] as i32 | ((buf[1] & 0x0F) as i32) << 8;
                let c2 : i32 = (((buf[1] & 0xF0) as i32) >> 4) | ((buf[2] as i32) << 4);
                

                // rejection sampling
                if c1 < Q {
                   generated_polynomial.c[l] = c1;
                   l += 1;
                }

                if c2 < Q && l < N {
                    generated_polynomial.c[l] = c2;
                    l += 1;
                }                
            }

            // no need to ntt the polynomial as the distribution would be the same in the NTT domain,
            // we can consider its output as NTT domain values
            A_ntt.set_coef(i, j, generated_polynomial);
        }

        
    }


    A_ntt

}



pub fn generate_noise_vector(eta : i16,nonce : u8,sigma_or_r : &Vec<u8>) -> (Vector<N>, u8){

    /*
    noise vectors are considered "small" vectors,
    meaning that their coefficients are symetrical modulus ( mods q )
    and additionally this modulus is constrained between -eta and eta
    */

    // wrapping mul as it should never be in a case where any overflow could happen anyway.
    // eta*64 bytes = 2*eta*256 bits
    let random_bytes = PRF(sigma_or_r,nonce,eta.wrapping_mul(64) as usize);

    let mut random_vector : Vector<N> = empty_vector();
    for i in 0..N {
        let mut sample = cbd_sample_from_bits(&random_bytes,eta as usize,i);

        // in kyber, the noise is then directly reduced modulo Q to map to [0;Q]
        sample = sample.rem_euclid(Q);
        random_vector.set(i, sample);
    }
    (random_vector, nonce + 1)
}


pub fn generate_noise_polyvector(eta : i16,mut nonce : u8,sigma_or_r : &Vec<u8>) -> (PolyVector<k>, u8){

    let mut random_polyvector : PolyVector<k> = empty_polyvector();
    let mut coef : Vector<N>;
    for i in 0..k {
        (coef, nonce) = generate_noise_vector(eta,nonce,sigma_or_r);
        random_polyvector.set(i, coef);
    }

    (random_polyvector, nonce)
}


#[allow(non_snake_case)]
pub fn compute_t(mut A : PolyMatrix<k,k>, s : PolyVector<k>, mut e : PolyVector<k>) -> PolyVector<k>{

    // t = As+e
    
    let mut t : PolyVector<k> = A.ntt_mult_vec(s);

    // as e is not in NTT domain
    e.c = e.c.map(ntt::ntt);

    t.add(e);

    t
}





pub fn compress(mut u : Vector<N>, d : i32) -> Vector<N>{
    
    for i in 0..u.c.len() {
        // Rust uses truncature for integer division, so we add Q/2 to round to the next integer instead of flooring
        //println!("compression coefficient : {:?}",((u.c[i]*(1 << d)+Q/2)/Q) % (1 << d));
        u.c[i] = ((u.c[i]*(1 << d)+Q/2)/Q) % (1 << d);

        // positive symetrical modulus
        if u.c[i] < 0 {
            u.c[i] += 1 << d;
        }
        // here, u is null ???? why ???
    }

    u

}



fn compress_polyvector(mut u : PolyVector<k>, d : i32) -> PolyVector<k>{

    for i in 0..k {
        // wtf, compress returns 0
        u.c[i] = compress(u.c[i],d);
    }
    
    u

}
pub fn decompress(mut u : Vector<N>, d : i32) -> Vector<N>{

    for i in 0..u.c.len() {
        // Rust uses truncature for integer division, so we add Q/2 to round to the next integer instead of flooring
        u.c[i] = (u.c[i] * Q + (1 << (d - 1))) / (1 << d);    }

    u

}

fn decompress_polyvector(mut u : PolyVector<k>, d : i32) -> PolyVector<k>{

    for i in 0..k {
        u.c[i] = decompress(u.c[i],d);
    }

    u
}

fn round(v : Vector<N>) -> Vector<N>{

    /*
    Map all values on the south emisphere of the circle mods Q to 1
    And all the north emisphere to 0
     */


    let mut r = empty_vector();
    for i in 0..N {

        let val = v.c[i];

        if val > Q/4 && val < 3*Q/4  {
            r.set(i,1);

        }else{
            r.set(i,0);
        }

    }

    r
}



#[derive(Debug, Clone)]
#[derive(PartialEq)]
pub struct EncryptedMessage{
    // u and v are returned in their NTT form from the encrypt function
    pub u : PolyVector<k>,
    pub v : Vector<N>
}


#[allow(non_snake_case)]
pub fn encrypt(A : PolyMatrix<k,k>,t : PolyVector<k>, msg : Vector<N>,r : Vec<u8>) -> EncryptedMessage{


    // v = t^T.r + e_2 + round(q/2)*m

    let mut  nonce = 0;
    // generation order is important because of the nonce
    let mut r_ntt : PolyVector<k>;
    (r_ntt,nonce) = generate_noise_polyvector(ETA_1,nonce,&r);
    r_ntt.c = r_ntt.c.map(ntt::ntt); 
    let mut e1 : PolyVector<k>; 
    (e1,nonce) = generate_noise_polyvector(ETA_2,nonce,&r);
    let mut e2 : Vector<N>; 
    (e2,nonce) = generate_noise_vector(ETA_2,nonce,&r);


    //println!("r_ntt = {:?}",r_ntt);
    let mut v : Vector<N> = t.ntt_dot(r_ntt);                                
    //println!("t after ntt dot with rntt : \n {:?}",v);

    v = ntt::intt(v);
    v.add(e2);

    let mut tmp = msg;
    tmp.scalar_mult((Q+1)/2); // kyber use round to next int division, thus the +1
    v.add(tmp);


    // u := NTT−1(AT◦r) + e1 
    let mut u: PolyVector<k> = A.transpose().ntt_mult_vec(r_ntt);     

    //println!("u before ntt : {:?}",u);

    u.c = u.c.map(ntt::intt);
    u.add(e1);

    // u and v returned in the normal domain
    EncryptedMessage { u:compress_polyvector(u,D_U), v:compress(v,D_V) }     // CHECK COMPRESSION ???

}

pub fn decrypt(EncryptedMessage { u, v }: & EncryptedMessage, s : PolyVector<k>) -> Vector<N>{
    // assume u and v are in the NORMAL-domain as they should be compressed from a normal domain vector 
    // s is in the NTT domain
    
    // m = round(v - s^T.u)
    let s_copy = s; // do not modify s, as it could be used for multiple decryption operations

    let mut u_decompressed = decompress_polyvector(*u,D_U);
    let mut v_decompressed = decompress(*v,D_V);

    v_decompressed = ntt::ntt(v_decompressed);

    u_decompressed.c = u_decompressed.c.map(ntt::ntt);

    v_decompressed.sub(    s_copy.ntt_dot(u_decompressed));
    //
    v_decompressed = ntt::intt(v_decompressed);
    round(v_decompressed)
}


pub fn hash_public_key(pk: &PublicKey) -> Vec<u8> {
    let mut bytes = pk.t.encode(12); // 12 bits/coeff packing, same layout as parse_polyvector expects
    bytes.extend_from_slice(&pk.rho);
    Sha3_256::digest(&bytes).to_vec()
}

pub fn KDF(K_b : Vec<u8>, hashed_c : Vec<u8>) -> Vec<u8> {
    // a non standard key derivation function is used in kyber ? why ?
    // it looks a bit like the high entropy one with HMAC
    let mut hasher : Shake256 = Shake256::default();
    hasher.update(&K_b);
    hasher.update(&hashed_c);

    let mut reader = hasher.finalize_xof();
    let mut hash = vec![0u8; 32];
    reader.read(&mut hash);
    hash.to_vec()
}


// H(c) = SHA3-256(u||v) where u is compressed to 10 bits/coeff and v is compressed to 4 bits/coeff
fn hash_c(EncryptedMessage { u, v }: EncryptedMessage) -> Vec<u8> {

    // assume u and v are already compressed
    let mut bytes = u.encode(D_U as usize); // 10 bits/coeff packing, 
    bytes.extend_from_slice(&v.encode(D_V as usize)); // 4 bits/coeff packing
    Sha3_256::digest(&bytes).to_vec()
}

pub fn decapsulate(c : EncryptedMessage, PrivateKey { s, pk, hashed_pk, z }: PrivateKey) -> Vec<u8>{

    // m = round(v - s^T.u)

    let msg = decrypt(&c, s); 
    // G(m || H(pk))
    let mut G = sha3::Sha3_512::new();
    sha3::Digest::update(&mut G, msg.encode(1));
    sha3::Digest::update(&mut G, &hashed_pk);
    let buf = G.finalize().to_vec();
    let K_b = buf[..32].to_vec();

    let sigma = buf[32..].to_vec();


    let c_p: EncryptedMessage = encrypt(pk.A, pk.t, msg,sigma);

    let hashed_c = hash_c(c.clone());

    if c_p == c {
        KDF(K_b, hashed_c)
    } else {
        KDF(z, hashed_c)
    }
}

pub fn encapsulate(pk: &PublicKey, msg : Vector<N>) -> (EncryptedMessage, Vec<u8>){

    let hashed_pk = hash_public_key(pk);
    let msg_bytes: Vec<u8> = msg.encode(1);
    let (K_b,r) = derive_coins(msg_bytes, hashed_pk);
    let c = encrypt(pk.A, pk.t, msg,r);
    let hashed_c = hash_c(c.clone());

    let K = KDF(K_b, hashed_c);

    (c,K)
}



/*
Uses encrypt function to directly encrypt a string into a vector of EncryptedMessage,
where each EncryptedMessage represents a chunk of 256 bits of the string.
*/
#[allow(non_snake_case)]
pub fn encrypt_string(input: &str,PublicKey { A, t, rho }: &PublicKey) -> Vec<EncryptedMessage> {
    // let input_chunks = string_to_vectors(input);
    // let mut encrypted_chunks = Vec::new();


    // for chunk in input_chunks {

    //     let encrypted_chunk = encrypt(*A, *t, chunk,);
    //     encrypted_chunks.push(encrypted_chunk);
    
    // }

    // encrypted_chunks

    Vec::new()
}


/*
Uses decrypt function to directly decrypt a vector of EncryptedMessage into a string,
where each EncryptedMessage represents a chunk of 256 bits of the string.
*/
pub fn decrypt_string(encrypted_chunks: &mut Vec<EncryptedMessage>, PrivateKey { s, pk: _, hashed_pk: _, z: _ }: &PrivateKey) -> String {
    let mut decrypted_chunks = Vec::new();

    for encrypted_chunk in encrypted_chunks {
        let decrypted_chunk = decrypt(encrypted_chunk, *s);
        decrypted_chunks.push(decrypted_chunk);
    }

    vectors_to_string(decrypted_chunks)
}


#[derive(Clone)]
#[allow(non_snake_case)]
pub struct PublicKey{
    pub A : PolyMatrix<k,k>,
    pub t : PolyVector<k>,
    #[allow(unused)]
    pub rho : Vec<u8>
}

impl PublicKey {

    pub fn new(rho : Vec<u8>,t : PolyVector<k>)->PublicKey{
        PublicKey { A: generate_A_from_seed(&rho),t, rho }
    }
}

#[derive(Clone)]
pub struct PrivateKey{
    pub s : PolyVector<k>,
    pub pk : PublicKey,
    pub hashed_pk : Vec<u8>,
    pub z : Vec<u8>
}



pub fn generate_key_pair(sigma : &Vec<u8>) -> (PublicKey, PrivateKey) {

    
    let seed_vector = generate_seed_vector();
    #[allow(non_snake_case)]
    let A = generate_A_from_seed(&seed_vector);
    let mut nonce = 0;
    let mut s : PolyVector<k>; 
    (s,nonce) = generate_noise_polyvector(ETA_1,nonce,sigma);
    s.c = s.c.map(ntt::ntt);
    let e : PolyVector<k>;
    (e,nonce) = generate_noise_polyvector(ETA_1,nonce,sigma);
    let t = compute_t(A, s, e);

    let public_key = PublicKey::new(seed_vector, t);
    let private_key = PrivateKey{s, pk: public_key.clone(), hashed_pk: Vec::new(), z: Vec::new()};

    (public_key, private_key)
}



pub fn generate_d_and_z(rng : &mut NistPqcAes256CtrRng)->(Vec<u8>,Vec<u8>){
    let mut d : [u8; 32] = [0; 32];
    let mut z : [u8; 32] = [0; 32];

    rng.fill_bytes(d.as_mut_slice());
    rng.fill_bytes(z.as_mut_slice());

    (d.to_vec(),z.to_vec())

}

pub fn generate_message_from_seed(rng : &mut NistPqcAes256CtrRng)->(Vector<N>, Vec<u8>){

    let mut b_msg : [u8; 32] = [0; 32];
    let mut msg_vector : Vector<N> = empty_vector();
    rng.fill_bytes(b_msg.as_mut_slice());

    for(i,byte) in b_msg.iter_mut().enumerate() {
        for j in 0..8{
            msg_vector.c[i*8+j] = ((*byte >> j) & 1) as i32;
        }
    }
    (msg_vector,b_msg.to_vec())

}


pub fn generate_rho_and_sigma(d : Vec<u8>)->(Vec<u8>,Vec<u8>){

    let mut hasher : Sha3_512 = Sha3_512::new();
    Digest::update(&mut hasher, d.as_slice());
    let hash = hasher.finalize();
    
    (hash[..32].to_vec(),hash[32..].to_vec())
}

pub fn derive_coins(_m : Vec<u8>, hashed_pk : Vec<u8>) -> (Vec<u8>, Vec<u8>) {
    let mut G : Sha3_512 = Sha3_512::new();
    Digest::update(&mut G, _m);
    Digest::update(&mut G, hashed_pk.as_slice());
    let hash = G.finalize();
    let tmp = hash.to_vec();

    (tmp[..32].to_vec(), tmp[32..].to_vec())
}

pub fn PRF(sigma_or_r : &Vec<u8>, nonce : u8, output_len : usize) -> Vec<u8> {
    let mut hasher : Shake256 = Shake256::default();
    hasher.update(sigma_or_r.as_slice());
    hasher.update(&[nonce]);

    let mut reader = hasher.finalize_xof();
    let mut hash = vec![0u8; output_len];
    reader.read(&mut hash);
    hash.to_vec()
}


