use crate::math_utils::MontgomeryForm;
use crate::math_utils::Vector;
use crate::math_utils::square_and_mult;

// zeta powers arranged in bit-reversed order
    const ZETA_TABLE : [i32;128] = [1, 3328, 1729, 1600, 2580, 749, 3289, 40, 2642, 687, 630, 2699, 1897, 1432, 848, 2481, 1062, 2267, 1919, 1410, 193, 3136, 797, 2532, 2786, 543, 3260, 69, 569, 2760, 1746, 1583, 296, 3033, 2447, 882, 1339, 1990, 1476, 1853, 3046, 283, 56, 3273, 2240, 1089, 1333, 1996, 1426, 1903, 2094, 1235, 535, 2794, 2882, 447, 2393, 936, 2879, 450, 1974, 1355, 821, 2508, 289, 3040, 331, 2998, 3253, 76, 1756, 1573, 1197, 2132, 2304, 1025, 2277, 1052, 2055, 1274, 650, 2679, 1977, 1352, 2513, 816, 632, 2697, 2865, 464, 33, 3296, 1320, 2009, 1915, 1414, 2319, 1010, 1435, 1894, 807, 2522, 452, 2877, 1438, 1891, 2868, 461, 1534, 1795, 2402, 927, 2647, 682, 2617, 712, 1481, 1848, 648, 2681, 2474, 855, 3110, 219, 1227, 2102, 910, 2419];
    /*[
    2285, 2571, 2970, 1812, 1493, 1422, 287, 202, 3158, 622, 1577, 182, 962,
    2127, 1855, 1468, 573, 2004, 264, 383, 2500, 1458, 1727, 3199, 2648, 1017,
    732, 608, 1787, 411, 3124, 1758, 1223, 652, 2777, 1015, 2036, 1491, 3047,
    1785, 516, 3321, 3009, 2663, 1711, 2167, 126, 1469, 2476, 3239, 3058, 830,
    107, 1908, 3082, 2378, 2931, 961, 1821, 2604, 448, 2264, 677, 2054, 2226,
    430, 555, 843, 2078, 871, 1550, 105, 422, 587, 177, 3094, 3038, 2869, 1574,
    1653, 3083, 778, 1159, 3182, 2552, 1483, 2727, 1119, 1739, 644, 2457, 349,
    418, 329, 3173, 3254, 817, 1097, 603, 610, 1322, 2044, 1864, 384, 2114, 3193,
    1218, 1994, 2455, 220, 2142, 1670, 2144, 1799, 2051, 794, 1819, 2475, 2459,
    478, 3221, 3021, 996, 991, 958, 1869, 1522, 1628
    ];*/

    const ZETA_INV_TABLE : [i32;128] = [1, 3328, 1600, 1729, 40, 3289, 749, 2580, 2481, 848, 1432, 1897, 2699, 630, 687, 2642, 1583, 1746, 2760, 569, 69, 3260, 543, 2786, 2532, 797, 3136, 193, 1410, 1919, 2267, 1062, 2508, 821, 1355, 1974, 450, 2879, 936, 2393, 447, 2882, 2794, 535, 1235, 2094, 1903, 1426, 1996, 1333, 1089, 2240, 3273, 56, 283, 3046, 1853, 1476, 1990, 1339, 882, 2447, 3033, 296, 2419, 910, 2102, 1227, 219, 3110, 855, 2474, 2681, 648, 1848, 1481, 712, 2617, 682, 2647, 927, 2402, 1795, 1534, 461, 2868, 1891, 1438, 2877, 452, 2522, 807, 1894, 1435, 1010, 2319, 1414, 1915, 2009, 1320, 3296, 33, 464, 2865, 2697, 632, 816, 2513, 1352, 1977, 2679, 650, 1274, 2055, 1052, 2277, 1025, 2304, 2132, 1197, 1573, 1756, 76, 3253, 2998, 331, 3040, 289];
    const ZETA_0 : i32= 17;
    const ZETA_INV_0  : i32 = 1175;

// Kyber style NTT, stops at degree 1 to use nth unity root
// returns 128 polynomials => vector still size 256 
pub fn ntt(p : Vector<256>) -> Vector<256>{
    
        



    // built Montgomery rpz of the vector

    let mut mtg_p = [MontgomeryForm::new(0);256];

   for i in 0..256{
        mtg_p[i].set_n(p.c[i]);
        mtg_p[i].set_a((p.c[i]<<12) % 3329);
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


    return Vector::new(&reduced_p, 3329)

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

    for i in 0..256{
        mtg_p[i].set_n(p.c[i]);
        mtg_p[i].set_a((p.c[i]<<12) % 3329);
    }



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

    return ret;
}


/*
Uses NTT to quickly compute the product of two degree 256 polynomials 
*/
pub fn poly_mult(a : Vector<256>, b : Vector<256>)->Vector<256>{


    

    let ntt_a = ntt(a);
    let ntt_b = ntt(b);



    let mut ntt_c: Vector<256> = Vector::new(&[0;256], 3329);
    // step by 2 as standardised kyber NTT returns degree 1 polynomials
    // instead of an array of constants

    let mut zeta_index : usize = 0;
    for i in (0..256).step_by(2){
        let a0 = ntt_a.c[i];
        let a1 = ntt_a.c[i+1];
        let b0 = ntt_b.c[i];
        let b1 = ntt_b.c[i+1];

        ntt_c.c[i] = a0.wrapping_mul(b0).wrapping_add(ZETA_TABLE[zeta_index as usize].wrapping_mul(a1.wrapping_mul(b1)));

        // x^2 % (x^2+1) = -1
        // -1 = zeta^n/2 
        ntt_c.c[i+1] = a0.wrapping_mul(b1).wrapping_add(a1.wrapping_mul(b0));
        zeta_index+=1;
    }

    return intt(ntt_c);
}



pub fn compute_zeta_table() -> [i32;128]{

    let mut z : [i32;128] = [0;128];
    for i in 0..128u8{

        let i_r = i.reverse_bits();

        let res = square_and_mult(ZETA_0 as u32, i_r as u32, 3329) as i32;
        //println!("17^{i_r}%3329 = {res}");
        z[ i as usize] = res; // bit reversed order
    

    }

    println!("{:?}",z);

    return z;
}


pub fn compute_inv_zeta_table() -> [i32;128]{

    let mut z : [i32;128] = [0;128];
    for i in 0..128u8{

        let i_r = i.reverse_bits();

        let res = square_and_mult(ZETA_INV_0 as u32, i_r as u32, 3329) as i32;
        println!("17^-{i_r}%3329 = {res}");
        z[ i as usize] = res; // bit reversed order
    

    }

    println!("{:?}",z);

    return z;
}



#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ntt_intt_roundtrip() {
        // Create a simple polynomial: f(x) = 1 + 2x + 3x^2 + ...
        let mut coefficients = [0; 256];
        for i in 0..256 {
            coefficients[i] = i as i32;
        }
        let p = Vector::new(&coefficients, 3329);

        // Perform NTT then iNTT
        let ntt_p = ntt(p.clone());
        let back_p = intt(ntt_p);
        println!("BACK_P = {:?}",back_p);

        // Check if we got the original back
        for i in 0..256 {
            assert_eq!(
                back_p.c[i], 
                coefficients[i], 
                "Roundtrip failed at index {}. Expected {}, got {}", 
                i, coefficients[i], back_p.c[i]
            );
        }
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
        let mut tmp = 0;
        
        for i in 0..31{
            tmp = square_and_mult(17, i, 3329);
            let real_result = (ZETA_0 as u128).pow(i) % 3329;
            //println!("17^{i}={tmp} \treal result = {real_result}");
            assert_eq!(tmp,real_result as u32);
        }
    }


}
