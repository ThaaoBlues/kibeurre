use crate::math_utils::Vector;
use crate::math_utils::Matrix;
use crate::ntt;
use rand;
use shake;
use shake::ExtendableOutput;
use shake::Update;
use shake::XofReader;


// ML-KEM-768
const q : i32 = 3329;
const n : usize = 256;
const k : usize = 3;
const eta_1 : i32 = 2;
const eta_2 : i32 = 2;
const d_u : i32 = 10;
const d_v : i32 = 4;

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



fn generate_seed_vector() -> Vector<n>{

    let mut r  :[i32; n] = [0; n];
    for i in 0..n{
        r[i] = rand::random::<i32>() % 2;
    }

    return Vector::new(&r,3329);

}

fn generate_A_from_seed(seed : Vector<n>) -> Matrix<k,n>{

    /*
    
    Generates the NTT version of the matrix A from a seed vector.
     */

    
    let mut A : Matrix<k,n> = Matrix::new([[0; n]; k],3329);
    let mut A_ntt : Matrix<k,n> = Matrix::new([[0; n]; k],3329);


    for i in 0..k {
        for j in 0..n {
            let mut hasher = shake::Shake128::default();

            // 2. Feed data into the hasher
            // log2(3329) = 12 bits par coef de polynome ? 
            hasher.update(&seed.get_bytes());
            hasher.update(&[i as u8, j as u8]);
            
            let mut reader = hasher.finalize_xof();

            // log2(3329) round = 12 bits => u16
            let mut buf: [u8; 2] = [0; 2];
            reader.read(&mut buf);
            A.set_coef(i,j,((buf[0] as i32) << 8 + buf[1] as i32) % 3329)

        }

        A_ntt.set_row(i,ntt::ntt(A.get_row(i)));

        
    }


    return A_ntt;

}


fn generate_small_vector(j : i32, eta : i32) -> Vector<n>{
    return Vector::new(&[0; n],3329);
}


fn generate_noise_vector(eta : i16) -> Vector<n>{

    /*
    noise vectors are considered "small" vectors,
    meaning that their coefficients are symetrical modulus ( mods q )
    and additionally this modulus is constrained between -eta and eta
    */

    let mut random_vector : Vector<n> = Vector::new(&[0; n],3329);
    for i in 0..n {
        let sample = cbd_sample(eta);
        random_vector.set(i, sample as i32);
    }
    return random_vector;
}

fn compute_t(mut A : Matrix<k,n>, s : Vector<n>, e : Vector<k>) -> Vector<k>{

    // t = As+e
    
    let mut t : Vector<k> = A.mult_vec(s);
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

fn decompress(mut u : Vector<n>, d : i32) -> Vector<n>{

    for i in 0..u.c.len() {
        u.set(i, ((u.c[i] as f32 * ((q as f32)/ (1 << d) as f32) ).round() as i32) % (1 << d));
    }


    return u;

}


fn round(v : Vector<n>) -> Vector<n>{


    /*
    Map all values on the south emisphere of the circle mods q to 1
    And all the north emisphere to 0
     */


    let mut r = Vector::new(&[0; n],3329);
    for i in 0..v.c.len() {

        let mut val = v.get(i);

        if val.abs() < q/4 {
            r.set(i,0);

        }else{
            r.set(i,1);
        }

    }



    return r;

}


fn decode_pk(){}

fn encode_pk(){}


fn encrypt(v : Vector<n>,t : Vector<n>, r : Vector<n>, m : Vector<n>) -> Vector<n>{
    // v = t^T.r + e_2 + round(q/2)*m
    let r_ntt = ntt::ntt(r);

    // u := NTT−1(AT◦r) + e1 
    let u = ntt:intt()
    return t.

}

fn decrypt(){

    // m = s^T - u

}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_normal_sample() {
        println!("Normal sample: {}", normal_sample());
    }

}