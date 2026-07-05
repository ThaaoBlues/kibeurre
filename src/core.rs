use std::ops::Shl;

use crate::math_utils::{PolyMatrix,PolyVector,Vector,Matrix,empty_polymatrix,empty_polyvector,empty_vector};

use crate::ntt;
use rand;
use shake;
use shake::ExtendableOutput;
use shake::Update;
use shake::XofReader;
use crate::parameters::{n, k, q,eta_1, eta_2, d_u, d_v,m};



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
        r[i] = (rand::random::<i32>() % 2).abs();
    }

    return Vector::new(&r,m);

}

pub fn generate_A_from_seed(seed : Vector<n>) -> PolyMatrix<k,k>{

    /*
    
    Generates the NTT version of the matrix A from a seed vector.
     */

    // A is used for debugging and will be removed,
    // A_ntt is used for the actual encryption/decryption process
    
    let mut A : PolyMatrix<k,k> = empty_polymatrix();
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
            
            A.set_coef(i,j,generated_polynomial);
            A_ntt.set_coef(i, j, ntt::ntt(generated_polynomial));
        }

        
    }


    return A_ntt;

}


fn generate_small_vector(j : i32, eta : i32) -> Vector<n>{
    return Vector::new(&[0; n],m);
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


fn decode_pk(){}

fn encode_pk(){}


#[derive(Debug)]
pub struct EncryptedMessage{
    u : PolyVector<k>,
    v : Vector<n>
}
/*
returns [u,v]
*/
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

pub fn decrypt(EncryptedMessage { u, v }: &mut EncryptedMessage, s : &mut PolyVector<k>) -> Vector<n>{

    // m = round(v - s^T.u)
    v.sub(s.ntt_dot(*u));
    return round(*v);

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
        let A = generate_A_from_seed(seed_vector);
        let s = generate_noise_polyvector(eta_1);
        let e = generate_noise_polyvector(eta_2);
        let t = compute_t(A, s.clone(), e);
        let r = generate_noise_polyvector(eta_1);
        let msg = generate_seed_vector(); // using seed vector as a message for testing
        let mut encrypted_message = encrypt(A, t, r, msg.clone());
        println!("Encrypted message: u = {:?}, v = {:?}", encrypted_message.u, encrypted_message.v);
        let decrypted_msg = decrypt(&mut encrypted_message, &mut s.clone());



        println!("\n\n=======\nDecrypted message: {:?}\n\n=========Original message: {:?}========", decrypted_msg, msg);
        assert_eq!(msg.c, decrypted_msg.c,"decryption failed");
    }
}