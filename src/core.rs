use std::ops::Shl;

use crate::math_utils::{PolyMatrix,PolyVector,Vector,empty_polymatrix,empty_polyvector,empty_vector};

use crate::ntt;
use rand;
use shake::{ExtendableOutput, Update, XofReader};
use crate::parameters::{n, k, q,eta_1, eta_2, d_u, d_v,m};
use crate::format_utils::{string_to_vectors,vectors_to_string};



fn cbd_sample(eta : i16) -> i16 {

    // get a probability from a uniform distribution between 0 and 1 
    // and return a sample from the centered binomial distribution

    /*
        Rounds are made by sampling 16 bits at a time from a uniform distribution 
        and then computing the sum of the differences between the first and second halves of the bits. 
        This is done until the desired number of samples (b) is reached. 
        The final sum is returned as the sample from the centered binomial distribution.
     */



    let mut sum : i16 = 0;

    let p: i32 = rand::random::<i32>();



    for i in 0..eta {
        sum += ((p & (1 << i)) >> i) as i16 - ((p & (1 << i+16)) >> i+16) as i16;

    }


    return  sum; // [-eta, eta]
    
}




pub fn generate_seed_vector() -> Vector<n>{

    let mut r  :[i32; n] = [0; n];
    for i in 0..n{
        // random bits
        r[i] = (rand::random::<u8>() % 2) as i32;
    }

    return Vector::new(&r,m);

}

#[allow(non_snake_case)]
pub fn generate_A_from_seed(seed : Vector<n>) -> PolyMatrix<k,k>{

    /*
    
    Generates the NTT version of the matrix A from a seed vector.
     */

    // A is used for debugging and will be removed,
    // A_ntt is used for the actual encryption/decryption process
    
    //let mut A : PolyMatrix<k,k> = empty_polymatrix();
    let mut A_ntt : PolyMatrix<k,k> = empty_polymatrix();
    let mut generated_polynomial : Vector<n> = empty_vector();


    for i in 0..k {
        for j in 0..k {

            for l in 0..n {
                let mut hasher = shake::Shake128::default();

                // 2. Feed data into the hasher
                // log2(m) = 12 bits par coef de polynome ? 
                hasher.update(&seed.get_bytes());
                hasher.update(&[i as u8, j as u8]);
                
                let mut reader = hasher.finalize_xof();

                // log2(m) round = 12 bits => u16
                let mut buf: [u8; 2] = [0; 2];
                reader.read(&mut buf);

                generated_polynomial.set(l,((buf[0] as i32).shl(8) as i32).wrapping_add(buf[1] as i32) % m as i32);
                
            }
            
            //A.set_coef(i,j,generated_polynomial);
            A_ntt.set_coef(i, j, ntt::ntt(generated_polynomial));
        }

        
    }


    return A_ntt;

}


fn generate_noise_vector(eta : i16) -> Vector<n>{

    /*
    noise vectors are considered "small" vectors,
    meaning that their coefficients are symetrical modulus ( mods q )
    and additionally this modulus is constrained between -eta and eta
    */

    let mut random_vector : Vector<n> = empty_vector();
    for i in 0..n {
        let sample = cbd_sample(eta);
        random_vector.set(i, sample as i32);
    }
    return random_vector;
}


pub fn generate_noise_polyvector(eta : i16) -> PolyVector<k>{

    let mut random_polyvector : PolyVector<k> = empty_polyvector();
    for i in 0..k {
        random_polyvector.set(i, generate_noise_vector(eta));
    }
    return random_polyvector;
}


#[allow(non_snake_case)]
pub fn compute_t(mut A : PolyMatrix<k,k>, s : PolyVector<k>, e : PolyVector<k>) -> PolyVector<k>{

    // t = As+e
    
    let mut t : PolyVector<k> = A.ntt_mult_vec(s);
    t.add(e);

    return t;
}




fn compress(mut u : Vector<n>, d : i32) -> Vector<n>{

    
    // TODO : find a more efficient storage than i32, as we discard log2(q)-d bits of information
    for i in 0..u.c.len() {
        u.set(i, ((u.c[i] as f32 * ((1 << d) as f32 / (q as f32)) ).round() as i32) % (1 << d));
    }


    return u;

}



fn compress_polyvector(mut u : PolyVector<k>, d : i32) -> PolyVector<k>{

    for i in 0..k {
        u.set(i, compress(u.get(i),d));
    }

    return u;

}

fn decompress(mut u : Vector<n>, d : i32) -> Vector<n>{

    for i in 0..u.c.len() {
        u.set(i, ((u.c[i] as f32 * ((q as f32)/ (1 << d) as f32) ).round() as i32) % (1 << d));
    }


    return u;

}

fn decompress_polyvector(mut u : PolyVector<k>, d : i32) -> PolyVector<k>{

    for i in 0..k {
        u.set(i, decompress(u.get(i),d));
    }

    return u;

}


fn round(v : Vector<n>) -> Vector<n>{


    /*
    Map all values on the south emisphere of the circle mods q to 1
    And all the north emisphere to 0
     */


    let mut r = empty_vector();
    for i in 0..n {

        let val = v.get(i);

        if val.abs() < q/4 {
            r.set(i,0);

        }else{
            r.set(i,1);
        }

    }



    return r;

}

fn round_polyvector(v : PolyVector<k>) -> PolyVector<k>{

    let mut r = empty_polyvector();
    for i in 0..k {
        r.set(i, round(v.get(i)));
    }

    return r;

}


#[derive(Debug)]
pub struct EncryptedMessage{
    u : PolyVector<k>,
    v : Vector<n>
}
/*
returns [u,v]
*/
#[allow(non_snake_case)]
pub fn encrypt(A : PolyMatrix<k,k>,t : PolyVector<k>, r : PolyVector<k>, msg : Vector<n>) -> EncryptedMessage{


    // v = t^T.r + e_2 + round(q/2)*m


    let e_2 : Vector<n>= generate_noise_vector(eta_2);

    
    let mut r_ntt : PolyVector<k> = empty_polyvector();

    for i in 0..k {
        r_ntt.set(i, ntt::ntt(r.get(i)));
    }
    

    let mut v : Vector<n> = t.ntt_dot(r_ntt);
    v.add(e_2);

    let mut tmp = msg.clone();
    tmp.scalar_mult(q/2);
    v.add(tmp);


    // u := NTT−1(AT◦r) + e1 
    let mut u: PolyVector<k> = A.transpose().ntt_mult_vec(r);

    for i in 0..k {
        u.set(i, ntt::intt(u.get(i)));
    }

    let e1 = generate_noise_polyvector(eta_1);

    u.add(e1);

    return EncryptedMessage { u, v };

}

pub fn decrypt(EncryptedMessage { u, v }: &mut EncryptedMessage, s : PolyVector<k>) -> Vector<n>{

    // m = round(v - s^T.u)

    // do not modify s, as it could be used for multiple decryption operations
    let s_copy = s.clone();
    v.sub(s_copy.ntt_dot(*u));
    return round(*v);

}



/*
Uses encrypt function to directly encrypt a string into a vector of EncryptedMessage,
where each EncryptedMessage represents a chunk of 256 bits of the string.
*/
#[allow(non_snake_case)]
pub fn encrypt_string(input: &str,PublicKey { A, t, seed : _ }: &PublicKey,r : PolyVector<k>) -> Vec<EncryptedMessage> {
    let input_chunks = string_to_vectors(input);
    let mut encrypted_chunks = Vec::new();

    for chunk in input_chunks {
        let encrypted_chunk = encrypt(*A, *t, r, chunk);
        encrypted_chunks.push(encrypted_chunk);
    }

    return encrypted_chunks;
}


/*
Uses decrypt function to directly decrypt a vector of EncryptedMessage into a string,
where each EncryptedMessage represents a chunk of 256 bits of the string.
*/
pub fn decrypt_string(encrypted_chunks: &mut Vec<EncryptedMessage>, PrivateKey { s }: &PrivateKey) -> String {
    let mut decrypted_chunks = Vec::new();

    for encrypted_chunk in encrypted_chunks {
        let decrypted_chunk = decrypt(encrypted_chunk, *s);
        decrypted_chunks.push(decrypted_chunk);
    }

    return vectors_to_string(decrypted_chunks);
}


#[derive(Clone)]
pub struct PublicKey{
    #[allow(non_snake_case)]
    A : PolyMatrix<k,k>,
    t : PolyVector<k>,
    seed : Vector<n>
}

impl PublicKey {

    pub fn new(seed : Vector<n>,t : PolyVector<k>)->PublicKey{
        return PublicKey { A: generate_A_from_seed(seed), t: t, seed }
    }
}

#[derive(Clone)]
pub struct PrivateKey{
    s : PolyVector<k>
}

impl PrivateKey {

    pub fn new(s : PolyVector<k>)->PrivateKey{
        return PrivateKey { s }
    }
}


pub fn generate_key_pair() -> (PublicKey, PrivateKey) {
    let seed_vector = generate_seed_vector();
    let A = generate_A_from_seed(seed_vector.clone());
    let s = generate_noise_polyvector(eta_1);
    let e = generate_noise_polyvector(eta_2);
    let t = compute_t(A, s, e);

    let public_key = PublicKey::new(seed_vector, t);
    let private_key = PrivateKey::new(s);

    return (public_key, private_key);
}




#[cfg(test)]
mod tests {
    use super::*;


    #[test]
    fn test_cbd_sample() {
        let sample = cbd_sample(2);
        assert!(sample >= -2 && sample <= 2);
    
    }
    
    #[test]
    fn test_generate_seed_vector() {
        let seed_vector = generate_seed_vector();
        println!("Seed vector: {:?}", seed_vector);
    }

    #[test]
    #[allow(non_snake_case)]
    fn test_generate_A_from_seed() {
        let seed_vector = generate_seed_vector();
        let A = generate_A_from_seed(seed_vector);
        println!("A: {:?}", A);
    
    }


    #[test]
    fn test_generate_noise_vector() {
        let noise_vector = generate_noise_vector(eta_1);
        println!("Noise vector : {:?}",noise_vector);
    }

    #[test]
    fn test_generate_noise_polyvector() {
        let noise_polyvector = generate_noise_polyvector(eta_1);
        println!("Noise polyvector: {:?}", noise_polyvector);
    }

    #[test]
    fn test_encrypt_decrypt() {
        let seed_vector = generate_seed_vector();
        #[allow(non_snake_case)]
        let A = generate_A_from_seed(seed_vector);
        let s = generate_noise_polyvector(eta_1);
        let e = generate_noise_polyvector(eta_2);
        let t = compute_t(A, s.clone(), e);
        let r = generate_noise_polyvector(eta_1);
        let msg = generate_seed_vector(); // using seed vector as a message for testing
        let mut encrypted_message = encrypt(A, t, r, msg.clone());
        println!("Encrypted message: u = {:?}, v = {:?}", encrypted_message.u, encrypted_message.v);
        let decrypted_msg = decrypt(&mut encrypted_message, s);



        println!("\n\n=======\nDecrypted message: {:?}\n\n=========Original message: {:?}========", decrypted_msg, msg);
        assert_eq!(msg.c, decrypted_msg.c,"decryption failed");
    }
}