use crate::math_utils::MontgomeryForm;
use crate::math_utils::Vector;


// Kyber style NTT, stops at degree 1 to use nth unity root
// returns 128 polynomials => vector still size 256 
fn ntt(p : Vector<256>) -> Vector<256>{
    
        
    const zeta_0 : i32 = 17;

// zeta powers arranged in bite-reversed order
    const ZETA_TABLE : [i32;128] = [ 
0, 64, 32, 96, 16, 80, 48, 112, 8, 72, 40, 104, 24, 88, 56, 120,
  4, 68, 36, 100, 20, 84, 52, 116, 12, 76, 44, 108, 28, 92, 60, 124,
  2, 66, 34, 98, 18, 82, 50, 114, 10, 74, 42, 106, 26, 90, 58, 122,
  6, 70, 38, 102, 22, 86, 54, 118, 14, 78, 46, 110, 30, 94, 62, 126,
  1, 65, 33, 97, 17, 81, 49, 113, 9, 73, 41, 105, 25, 89, 57, 121,
  5, 69, 37, 101, 21, 85, 53, 117, 13, 77, 45, 109, 29, 93, 61, 125,
  3, 67, 35, 99, 19, 83, 51, 115, 11, 75, 43, 107, 27, 91, 59, 123,
  7, 71, 39, 103, 23, 87, 55, 119, 15, 79, 47, 111, 31, 95, 63, 127]; 
    

    // built Montgomery rpz of the vector

    let mut mtg_p = [MontgomeryForm::new(0);256];

    for i in 0..256{
        mtg_p[i].set_n(p.c[i]);
    }


    let mut k = 0;

    for n_sur_2 in [128,64,32,16,8,4,2]{


        
        // get chunks of length n/2 + n/2
        for chunk in (0..256).step_by(2*n_sur_2){



            k+=1;
            // perform symetric calculations on j and j+n/2
            // make use (probably unwise) of Montgomery reduction
            for i in chunk..(chunk+n_sur_2){
                let b : i32 =  mtg_p[i+n_sur_2].mult(ZETA_TABLE[k]).reduction();

                let a : i32 = mtg_p[i].get_a();
                
                mtg_p[i].set_a((a + b) % 3329);

                mtg_p[i+n_sur_2].set_a((a - b + 3329) % 3329);

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


fn intt(){

}