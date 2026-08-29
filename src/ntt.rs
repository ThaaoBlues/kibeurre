use crate::math_utils::MontgomeryForm;
use crate::math_utils::Vector;
use crate::math_utils::square_and_mult;

// zeta powers arranged in bit-reversed order
// /!\ in the kyber reference C code, the array contains the values directly in Montgomery domain
// Here, we store the values in normal representation, and convert to Montgomery form when needed
const ZETA_TABLE : [i32;128] = [1, 1729, 2580, 3289, 2642, 630, 1897, 848, 1062, 1919, 193, 797, 2786, 3260, 569, 1746, 296, 2447, 1339, 1476, 3046, 56, 2240, 1333, 1426, 2094, 535, 2882, 2393, 2879, 1974, 821, 289, 331, 3253, 1756, 1197, 2304, 2277, 2055, 650, 1977, 2513, 632, 2865, 33, 1320, 1915, 2319, 1435, 807, 452, 1438, 2868, 1534, 2402, 2647, 2617, 1481, 648, 2474, 3110, 1227, 910, 17, 2761, 583, 2649, 1637, 723, 2288, 1100, 1409, 2662, 3281, 233, 756, 2156, 3015, 3050, 1703, 1651, 2789, 1789, 1847, 952, 1461, 2687, 939, 2308, 2437, 2388, 733, 2337, 268, 641, 1584, 2298, 2037, 3220, 375, 2549, 2090, 1645, 1063, 319, 2773, 757, 2099, 561, 2466, 2594, 2804, 1092, 403, 1026, 1143, 2150, 2775, 886, 1722, 1212, 1874, 1029, 2110, 2935, 885, 2154];

const ZETA_INV_TABLE : [i32;128] = [1, 1600, 40, 749, 2481, 1432, 2699, 687, 1583, 2760, 69, 543, 2532, 3136, 1410, 2267, 2508, 1355, 450, 936, 447, 2794, 1235, 1903, 1996, 1089, 3273, 283, 1853, 1990, 882, 3033, 2419, 2102, 219, 855, 2681, 1848, 712, 682, 927, 1795, 461, 1891, 2877, 2522, 1894,1010, 1414, 2009, 3296, 464, 2697, 816, 1352, 2679, 1274, 1052, 1025, 2132, 1573, 76, 2998, 3040, 1175, 2444, 394, 1219, 2300, 1455, 2117,1607, 2443, 554, 1179, 2186, 2303, 2926, 2237, 525, 735, 863, 2768, 1230, 2572, 556, 3010, 2266, 1684, 1239, 780, 2954, 109, 1292, 1031, 1745, 2688, 3061, 992, 2596, 941, 892, 1021, 2390, 642, 1868, 2377, 1482, 1540, 540, 1678, 1626, 279, 314, 1173, 2573, 3096, 48, 667, 1920, 2229, 1041, 2606, 1692, 680, 2746, 568, 3312];
const ZETA_0 : i32= 17;
const ZETA_INV_0  : i32 = 1175;


// Kyber style NTT, stops at degree 1 to use nth unity root
// returns 128 polynomials => vector still size 256 
pub fn ntt(p : Vector<256>) -> Vector<256>{
    
        

    // built Montgomery rpz of the vector

    let mut mtg_p: [MontgomeryForm; 256] = [MontgomeryForm::new(0);256];

   for (i, mtg) in mtg_p.iter_mut().enumerate(){
        mtg.set_n(p.c[i]);
        mtg.set_a((p.c[i]<<12) % 3329);
    }


    let mut k = 0;

    for n_sur_2 in [128,64,32,16,8,4,2]{


        
        // get chunks of length n/2 + n/2
        for chunk in (0..256).step_by(2*n_sur_2){



            k+=1;
            // perform symetric calculations on j and j+n/2
            // make use (probably unwise) of Montgomery reduction
            for i in chunk..(chunk+n_sur_2){

                let b : i32 =  mtg_p[i+n_sur_2].mult(ZETA_TABLE[k]).get_a();

                let a : i32 = mtg_p[i].get_a();


                let mut a_plus_b = a + b;
                if(a+b) > 3329 {
                    a_plus_b = a+b - 3329;
                }

                let mut a_moins_b = a - b;
                if(a-b) < 0 {
                    a_moins_b = a-b + 3329;
                }
                
                mtg_p[i].set_a(a_plus_b);

                mtg_p[i+n_sur_2].set_a(a_moins_b);

            }
            
        }
    }


    let mut reduced_p = [0;256];

    for i in 0..256{
        reduced_p[i] = mtg_p[i].reduction();
    }


    Vector::new(&reduced_p, 3329)

    // recursive version
    /*let pe = p.get_even_indexes();
    let po = p.get_odd_indexes();
    let ye = ntt(pe);
    let mut yo = ntt(po);

    let y = [0;POL];

    let n_sur_2 = DEGREE/2;
    for j in 0..(n_sur_2){
        
        let mut yo_mult = yo.scalar_mult(psi**2*j+1);
        let opp_yo_mult = yo_mult.scalar_mult(-1);
        y[j] = ye.add(yo_mult);
        y[j+n_sur_2] = ye.add(opp_yo_mult);
    }*/
}


// see end of page 14
pub fn intt(p : Vector<256>)->Vector<256>{

    

    let mut mtg_p = [MontgomeryForm::new(0);256];

    for (i,mtg) in mtg_p.iter_mut().enumerate(){
        mtg.set_n(p.c[i]);
        mtg.set_a((p.c[i]<<12) % 3329);
    }


    //let mut k = 127;
    for n_sur_2 in [2,4,8,16,32,64,128]{


        // as we do in reverse, we must start k from its value at the end of each previous ntt chunk

        // for example, for n=2 we started with k=64 during ntt 
        // the following recreates it

        // we cannot start from 0 and increment as we did during ntt
        // as in reverse, k starting value is not linear
        let mut k: usize = 256/(2*n_sur_2) - 1;

        // get 2 chunks of length n/2
        for chunk in (0..256).step_by(2*n_sur_2){

            
            // we still increment k as we are using inverse zeta table
            k+=1;
            
            
            // perform symetric calculations on j and j+n/2
            // make use (probably unwise) of Montgomery reduction
            for i in chunk..(chunk+n_sur_2){

                let b : i32 =  mtg_p[i+n_sur_2].get_a();

                let a : i32 = mtg_p[i].get_a();

                let mut a_plus_b = a + b;
                if(a+b) > 3329 {
                    a_plus_b = a+b - 3329;
                }

                let mut a_moins_b = a - b;
                if(a-b) < 0 {
                    a_moins_b = a-b + 3329;
                }
                
                mtg_p[i].set_a(a_plus_b);


                mtg_p[i+n_sur_2].set_a(a_moins_b);
                mtg_p[i+n_sur_2] = mtg_p[i+n_sur_2].mult(ZETA_INV_TABLE[k])

            }

            
        }

    }


    let mut reduced_p = [0;256];

    for i in 0..256{
        reduced_p[i] = mtg_p[i].reduction();
    }
    let mut ret = Vector::new(&reduced_p, 3329);

    // 128^-1 % 3329 = 3303
    ret.scalar_mult(3303);

    ret
}


/*
Uses NTT to quickly compute the product of two degree 256 polynomials 
ASSUME BOTH A AND B ARE ALREADY IN NTT FORM
*/
pub fn poly_mult(ntt_a : Vector<256>, ntt_b : Vector<256>)->Vector<256>{


    let mut ntt_c: Vector<256> = Vector::new(&[0;256], 3329);
    // step by 2 as standardised kyber NTT returns degree 1 polynomials
    // instead of an array of constants

    let mut i    = 0;
    for zeta in ZETA_TABLE.iter(){
        let a0 = ntt_a.c[i];
        let a1 = ntt_a.c[i+1];
        let b0 = ntt_b.c[i];
        let b1 = ntt_b.c[i+1];

        ntt_c.c[i] = a0.wrapping_mul(b0).wrapping_add(zeta.wrapping_mul(a1.wrapping_mul(b1)));

        // x^2 % (x^2+1) = -1
        // -1 = zeta^n/2 
        ntt_c.c[i+1] = a0.wrapping_mul(b1).wrapping_add(a1.wrapping_mul(b0));


        i += 2;
    }

    ntt_c
}



pub fn compute_zeta_table() -> [i32;128]{

    let mut z : [i32;128] = [0;128];
    for i in 0..128u8{

        let i_r = i.reverse_bits() >> 1; // /!\ we only need 7 bits to compute powers of 17, as 17^128 = 1 mod 3329

        let res = square_and_mult(ZETA_0 as u32, i_r as u32, 3329) as i32;
        //println!("17^{i_r}%3329 = {res}");
        z[ i as usize] = res; // bit reversed order
    

    }

    z
}


pub fn compute_inv_zeta_table() -> [i32;128]{

    let mut z : [i32;128] = [0;128];
    for i in 0..128u8{

        // /!\ we only need 7 bits to compute powers of 17, as 17^128 = 1 mod 3329
        let i_r = i.reverse_bits() >> 1;

        let res = square_and_mult(ZETA_INV_0 as u32, i_r as u32, 3329) as i32;
        println!("17^-{i_r}%3329 = {res}");
        z[ i as usize] = res; // bit reversed order
    

    }

    z
}



#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn print_zeta_tables() {
        let zeta_table = compute_zeta_table();
        let zeta_inv_table = compute_inv_zeta_table();  
        println!("ZETA_TABLE: {:?}", zeta_table);
        println!("ZETA_INV_TABLE: {:?}", zeta_inv_table);
    }


    #[test]
    fn test_ntt_constant_one() {
        // NTT of f(x) = 1 should result in [1, 0, 1, 0, 1, 0...] 
        // because we split into degree-1 polynomials (a + bx)
        let mut coefficients = [0i32; 256];
        coefficients[0] = 1;
        let p = Vector::new(&coefficients, 3329);

        let ntt_p = ntt(p);

        //println!("{:?}",ntt_p);

        // In Kyber NTT, f(x)=1 maps to constant 1 in every small block
        for i in (0..256).step_by(2) {
            assert_eq!(ntt_p.c[i], 1, "Even index {} should be 1", i);
            assert_eq!(ntt_p.c[i+1], 0, "Odd index {} should be 0", i+1);
        }
    }

        #[test]
    fn test_ntt_constant() {
        // NTT of f(x) = 1 should result in [2, 0, 2, 0, 2, 0...] 
        // because we split into degree-1 polynomials (a + bx)
        let mut coefficients = [0i32; 256];
        coefficients[0] = 1325;
        let p = Vector::new(&coefficients, 3329);

        let ntt_p = ntt(p);

        //println!("{:?}",ntt_p);

        // In Kyber NTT, f(x)=1 maps to constant 1 in every small block
        for i in (0..256).step_by(2) {
            assert_eq!(ntt_p.c[i], 1325, "Even index {} should be 1325", i);
            assert_eq!(ntt_p.c[i+1], 0, "Odd index {} should be 0", i+1);
        }
    }


    
    #[test]
    fn test_square_and_mult(){
        let mut tmp: u32;
        
        for i in 0..31{
            tmp = square_and_mult(17, i, 3329);
            let real_result = (ZETA_0 as u128).pow(i) % 3329;
            //println!("17^{i}={tmp} \treal result = {real_result}");
            assert_eq!(tmp,real_result as u32);
        }
    }

    #[test]
    fn test_ntt_intt_random_polynomials() {
        // Deterministic PRNG so failures are reproducible.
        // We don't need a cryptographic RNG for a correctness test.
        let mut seed: u64 = 0x1234_5678_9ABC_DEF0;

        fn next_random(seed: &mut u64) -> i32 {
            // xorshift64*
            *seed ^= *seed << 13;
            *seed ^= *seed >> 7;
            *seed ^= *seed << 17;

            (*seed % 3329) as i32
        }

        const NB_TESTS: usize = 100;

        for test_number in 0..NB_TESTS {
            let mut coefficients = [0i32; 256];

            for coefficient in coefficients.iter_mut() {
                *coefficient = next_random(&mut seed);
            }

            let p = Vector::new(&coefficients, 3329);

            let ntt_p = ntt(p);
            let back_p = intt(ntt_p);

            for i in 0..256 {
                assert_eq!(
                    back_p.c[i],
                    coefficients[i],
                    "NTT/INTT roundtrip failed: test {}, coefficient {}, expected {}, got {}",
                    test_number,
                    i,
                    coefficients[i],
                    back_p.c[i]
                );
            }
        }
}

}
