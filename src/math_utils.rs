use crate::parameters::{n, k,m};
use crate::ntt;


const POLYNOMIAL_SIZE : usize = n;
pub const EMPTY_VECTOR : [i32;n] = [0;n];
pub fn empty_vector() -> Vector<n>{
    Vector::new(&EMPTY_VECTOR,m)
}

pub fn empty_polyvector() -> PolyVector<k>{
    PolyVector::new(&[Vector::new(&EMPTY_VECTOR, m); k],m)
}

pub fn empty_polymatrix() -> PolyMatrix<k,k>{
    PolyMatrix::new([PolyVector::new(&[Vector::new(&EMPTY_VECTOR, m); k],m); k],m)
}

#[derive(Copy)]
#[derive(Clone)]
#[derive(Debug)]
pub struct Vector<const VECTOR_SIZE : usize> {
    pub c : [i32;VECTOR_SIZE],
    m : i32
}

impl<const VECTOR_SIZE : usize> Vector<{VECTOR_SIZE}> {

    pub fn new( val : &[i32;VECTOR_SIZE],_m : i32 ) -> Vector<VECTOR_SIZE>{
        Vector { c : *val, m : _m}
    }

    pub fn set(&mut self,index : usize,v : i32){
        self.c[index] = v;
    }

    pub fn get(&self, index : usize)->i32{
        self.c[index]
    }


    pub fn add(&mut self, v2 : Vector<VECTOR_SIZE>){
        

        //let mut tmp : [i32;VECTOR_SIZE] = [0;VECTOR_SIZE];

        for i in 0..VECTOR_SIZE {
            self.c[i] = (self.c[i] + v2.c[i]) % self.m;
        }

        //Vector::new( &tmp);

    }


    pub fn sub(&mut self, v2 : Vector<VECTOR_SIZE>){
        

        //let mut tmp : [i32;VECTOR_SIZE] = [0;VECTOR_SIZE];

        for i in 0..VECTOR_SIZE {
            self.c[i] = (self.c[i] + v2.c[i]) % self.m;
        }

        //Vector::new( &tmp);

    }




    pub fn dot(&mut self, v2 : Vector<VECTOR_SIZE>) -> i32{
        let mut ret : i32 = 0;


        for i in 0..VECTOR_SIZE {
            let tmp = MontgomeryForm::new(self.c[i]);
            ret += tmp.mult(v2.c[i]).to_standard();
        }

        //println!("{:?} DOT {:?} = {}",self, v2,ret);
        ret % self.m
    }



    pub fn scalar_mult(&mut self,v : i32){

        for i in 0..VECTOR_SIZE {
            let tmp = MontgomeryForm::new(v);
            self.c[i] = tmp.mult(self.c[i]).to_standard();
        } 

    }



    pub fn get_bytes(&self)->Vec<u8>{
        let mut ret : Vec<u8> = Vec::new();

        for i in 0..VECTOR_SIZE {
           ret.extend_from_slice(&self.c[i].to_le_bytes());
        }

        ret
    }

}


/*

Vectors of size 256 polynomials
*/


#[derive(Copy)]
#[derive(Clone)]
#[derive(Debug)]
pub struct PolyVector<const VECTOR_SIZE : usize> {
    pub c : [Vector<POLYNOMIAL_SIZE>;VECTOR_SIZE],
    m : i32
}

impl<const VECTOR_SIZE : usize> PolyVector<{VECTOR_SIZE}> {

    pub fn new( val : &[Vector<POLYNOMIAL_SIZE>;VECTOR_SIZE],_m : i32 ) -> PolyVector<VECTOR_SIZE>{
        PolyVector { c : *val, m : _m}
    }

    pub fn set(&mut self,index : usize,v : Vector<POLYNOMIAL_SIZE>){
        self.c[index] = v;
    }
    pub fn get(&self, index : usize)->Vector<POLYNOMIAL_SIZE>{
        self.c[index]
    }


    pub fn add(&mut self, v2 : PolyVector<VECTOR_SIZE>){
        

        //let mut tmp : [i32;VECTOR_SIZE] = [0;VECTOR_SIZE];

        for i in 0..VECTOR_SIZE {
            self.c[i].add(v2.c[i]);
        }

        //Vector::new( &tmp);

    }


    pub fn ntt_dot(&self, v2 : PolyVector<VECTOR_SIZE>) -> Vector<POLYNOMIAL_SIZE>{
        /* assume both vectors are in NTT form, and compute the dot product in NTT form
         assume the first vector is line and the second one is column
         so the result is a single polynomial

        */

        

        let mut ret : Vector<POLYNOMIAL_SIZE> = Vector::new(&[0;POLYNOMIAL_SIZE],self.m);

        for poly_index in 0..VECTOR_SIZE {
            let tmp : Vector<POLYNOMIAL_SIZE> = ntt::poly_mult(self.c[poly_index],v2.c[poly_index]);
            //println!("tmp : {:?}",tmp);
            for poly_coef_index in 0..POLYNOMIAL_SIZE {
                
                ret.c[poly_index] += tmp.c[poly_coef_index];
            }
            
        }

        for poly_coef_index in 0..POLYNOMIAL_SIZE {    
            ret.c[poly_coef_index] %= n as i32;
        }

        //println!("ret : {:?}",ret);
        ret
    
    }



    pub fn scalar_mult(&mut self,v : i32){
        // apply multiplication on each polynomial by n

        for i in 0..VECTOR_SIZE {
           self.c[i].scalar_mult(v);
        } 

    }

    


}


/*
Matrices of polynomials

*/

#[derive(Copy)]
#[derive(Clone)]
#[derive(Debug)]
pub struct PolyMatrix<const RN : usize, const CN : usize>{
    row_vecs : [PolyVector<CN>;RN],
    col_vecs : [PolyVector<RN>;CN],
    m : i32
}


impl <const RN : usize, const CN : usize> PolyMatrix<{RN},{CN}>{

    pub fn new(val : [PolyVector<CN>;RN], _m : i32)->PolyMatrix<RN, CN>{


        // the array of arrays we get as parameter
        // is a rows representation of the matrix


        // build two representations to ease matrix multiplications
        // one where we store matrix by its rows
        // and the second by its columns

        // rows
        let r_v : [PolyVector<CN>;RN] = val;

        // columns
        let mut c_v : [PolyVector<RN>;CN] = [PolyVector::new(&[Vector::new(&[0;POLYNOMIAL_SIZE],_m);RN],_m);CN];


        for (i, col) in c_v.iter_mut().enumerate(){
            for (j,new_row) in val.iter().enumerate()  {
                col.c[j] = new_row.c[i];
            }
        }
        
        PolyMatrix { col_vecs : c_v, row_vecs : r_v, m : _m}

    }



    pub fn ntt_mult_vec(&mut self, v : PolyVector<CN>)->PolyVector<RN>{
        /*
        Assume matrix and vector are in NTT form, and compute the product in NTT form
        
        */

        // A = Av
        let mut ret : PolyVector<RN> = PolyVector::new(&[Vector::new(&[0;POLYNOMIAL_SIZE],self.m);RN],self.m);
        for i in 0..RN{

            
            // NTT PolyVector multiplication returns a Vector
            // as we sum every result, like in classic vector multiplication
            // (it is in fact, a classic vector)
            ret.c[i] = self.row_vecs[i].ntt_dot(v);
        }

        ret
    }


    #[allow(non_snake_case)]
    pub fn matmul(&mut self, B : PolyMatrix<CN,RN>)->PolyMatrix<CN,RN>{
        
        let mut ret : PolyMatrix<CN,RN> = PolyMatrix::new([PolyVector::new(&[Vector::new(&[0;POLYNOMIAL_SIZE],self.m);RN],self.m);CN],self.m);

        for i in 0..RN{
            
            for j in 0..CN {
                ret.row_vecs[i].c[j] = self.row_vecs[i].ntt_dot(B.col_vecs[j]);
                ret.col_vecs[j].c[i] = ret.row_vecs[i].c[j];
            }
        }


        ret


    }

    pub fn transpose(& self) -> PolyMatrix<CN,RN>{
        let mut ret : PolyMatrix<CN,RN> = PolyMatrix::new([PolyVector::new(&[Vector::new(&[0;POLYNOMIAL_SIZE],self.m);RN],self.m);CN],self.m);
        ret.col_vecs = self.row_vecs;
        ret.row_vecs = self.col_vecs;

        ret
    }

    pub fn set_row(&mut self, row_index : usize, v : PolyVector<CN>){
        self.row_vecs[row_index] = v;

        // update the columns based representation
        for col_index in 0..RN {
            self.col_vecs[col_index].set(row_index,self.row_vecs[row_index].get(col_index));
        }
        
        
    }

    pub fn get_row(&self, row_index : usize) -> PolyVector<CN>{
        self.row_vecs[row_index]
    }

    pub fn set_col(&mut self, col_index : usize, v : PolyVector<RN>){
        self.col_vecs[col_index] = v;

        // update the columns based representation
        for row_index in 0..RN {
            self.row_vecs[col_index].set(row_index,self.col_vecs[col_index].get(row_index));
        }    
    }

    pub fn set_coef(&mut self, col_index : usize, row_index : usize, v : Vector<POLYNOMIAL_SIZE>){
        self.row_vecs[row_index].set(col_index,v);

        self.col_vecs[col_index].set(row_index,v);
    }

}



// matrix of integers
#[derive(Copy)]
#[derive(Clone)]
#[derive(Debug)]
pub struct Matrix<const RN : usize, const CN : usize>{
    row_vecs : [Vector<CN>;RN],
    col_vecs : [Vector<RN>;CN],
    m : i32
}

impl <const RN : usize, const CN : usize> Matrix<{RN},{CN}>{

    pub fn new(val : [[i32;CN];RN], _m : i32)->Matrix<RN, CN>{


        // the array of arrays we get as parameter
        // is a rows representation of the matrix


        // build two representations to ease matrix multiplications
        // one where we store matrix by its rows
        // and the second by its columns

        // rows
        let mut r_v : [Vector<CN>;RN] = [Vector::new(&[0;CN],_m);RN];

        // columns
        let mut c_v : [Vector<RN>;CN] = [Vector::new(&[0;RN],_m);CN];

        for i in 0..CN{ 
            r_v[i] = Vector::new(&val[i],_m);
        }



        for i in 0..CN{
            for (j,new_row) in val.iter().enumerate() {
                c_v[i].c[j] = new_row[i];
            }
        }
        
        Matrix { col_vecs : c_v, row_vecs : r_v, m : _m}

    }


    pub fn mult_vec(&mut self, v : Vector<CN>)->Vector<RN>{
        // A = Av

        let mut ret : Vector<RN> = Vector::new(&[0;RN],self.m);

        for i in 0..RN{
            ret.c[i] = self.row_vecs[i].dot(v);
        }

        ret
    }

    #[allow(non_snake_case)]
    pub fn matmul(&mut self, B : Matrix<CN,RN>)->Matrix<CN,RN>{
        
        let mut ret : Matrix<CN,RN> = Matrix::new([[0;RN];CN],self.m);

        for i in 0..RN{
            
            for j in 0..CN {
                ret.row_vecs[i].c[j] = self.row_vecs[i].dot(B.col_vecs[j]);
                ret.col_vecs[j].c[i] = ret.row_vecs[i].c[j];
            }
        }


        ret


    }



    pub fn set_row(&mut self, row_index : usize, v : Vector<CN>){
        self.row_vecs[row_index] = v;

        // update the columns based representation
        for col_index in 0..RN {
            self.col_vecs[col_index].set(row_index,self.row_vecs[row_index].get(col_index));
        }
        
        
    }

    pub fn get_row(&self, row_index : usize) -> Vector<CN>{
        self.row_vecs[row_index]
    }

    pub fn set_col(&mut self, col_index : usize, v : Vector<RN>){
        self.col_vecs[col_index] = v;

        // update the columns based representation
        for row_index in 0..RN {
            self.row_vecs[col_index].set(row_index,self.col_vecs[col_index].get(row_index));
        }    
    }

    pub fn set_coef(&mut self, col_index : usize, row_index : usize, v : i32){
        self.row_vecs[row_index].set(col_index,v);

        self.col_vecs[col_index].set(row_index,v);
    }

}




/*
Montgomery form enables us to speed up contiguous modular multoplications

*/

#[derive(Copy)]
#[derive(Clone)]
#[derive(Debug)]
pub struct MontgomeryForm{

    a : i32,
    n : i32,
    q : i32,
    r_pow : u32,
    r : i32,
    r_1 : i32,
    k : i32
}


impl MontgomeryForm{

    // we assume q is prime as it will be in this usage ( q = 3329)
    // thus, we can hardcode values of r, r⁻¹ and k that will always work
    // n < q so i16 is more than sufficient here
    pub fn new(v: i32) -> MontgomeryForm{

        
        MontgomeryForm { a: (v<<12) % 3329, n:v, q :3329, r_pow:12, r : 4096,r_1 : 2704, k: 3327 }

    }

    pub fn set_a(&mut self, a : i32){
        self.a = a;
    }

    pub fn set_n(&mut self, v : i32){
        self.n = v;
    }

    pub fn get_a(&self)-> i32{
        self.a
    }

    pub fn mult(& self, b : i32) -> MontgomeryForm{

        
        // montgomery form of b
        let b_b : i32 = (b<<self.r_pow) % self.q;
        
        //let c_b : i32 = (i32::from(self.a) * i32::from(b_b)) % i32::from(self.q);

        let mut ab = MontgomeryForm { 
            a:self.a * b_b, 
            n: 0, 
            q: self.q, 
            r_pow:self.r_pow,
            r: self.r, 
            r_1: self.r_1, 
            k: self.k
        };
        
        let c_b = ab.reduction();

        MontgomeryForm { a:c_b, n: 0, q: self.q, r_pow:self.r_pow,r: self.r, r_1: self.r_1, k: self.k}
        
    }

    // goes back to [0;q] from montgomery form
    pub fn reduction(&mut self) -> i32{
        let s : i32  = i32::try_from((i64::from(self.a) * i64::from(self.k)) %  i64::from(self.r)).unwrap();


        let t : i32= self.a + s*self.q;


        let u : i32 = t >> self.r_pow; // u = t/r

        //println!("{0} = {1} ?",(u % self.q), i32::try_from((i64::from(self.a)*i64::from(self.r_1)) % i64::from(self.q)).unwrap() );

        
        if u < self.q {
            self.n = u;
            u
        }else{
            self.n = u-self.q;
            self.n
        }

        
    }



    pub fn to_standard(&mut self)->i32{
        (self.a * self.r_1) % self.q
    }
}


pub fn square_and_mult(x : u32, e : u32, _m : u32) -> u32{
    if e == 0 {
        return 1;
    }else if e == 1 {
        return x;
    }
    let mut ret : u32 = 1;
    let mut iter : u32= x ;
    let mut exp = e;
    // we decompose the exponent in binary
    // so we can multiply x by the iterator everytime have to
    // and thus, multiply by x^k, k being a power of 2
    // this boils down to summing powers of 2 until we recreate the exponent we wanted
    while exp > 1 {

        iter = iter *iter % _m;
        exp >>= 1;
        // current_pow = current_pow + k
        if exp % 2 == 1 {
            ret = ret * iter % _m;
        }
    }
    // final pow+1 if power was odd 
    if e % 2 == 1 {
        ret = ret*x % _m;
    }

    ret
}



pub fn mod_inv(zeta : i32,_n : i32) -> i32{
    let mut r : i32 = 0;
    let mut r_1 : i32;
    let mut u : i32 = 0;
    let mut v : i32 = 1;
    let mut u_1 : i32 = 1;
    let mut v_1 : i32 = 0;
    let mut b = _n;
    let mut a = zeta;

    // first round

    // au + bv = pgcd(a,b) = 1
    // b = n ; a = zeta; u = zeta⁻¹
    let mut tmp:i32;
    while r != 1 {
        
        r_1 = r;
        r = a % b;

        tmp = u;
        u = u_1 - (r_1/r)*u;
        u_1 = tmp;

        tmp = v;
        v = v_1 - (r_1/r)*v;
        v_1 = tmp;
        
        a = b;
        b = r;
        
    }
    println!("u={u}, v={v}");
    v
}






#[cfg(test)]
mod tests {
    use super::*;
    use rand::Rng;

#[test]
fn test_vector_add_with_modulo() {
    // (2000 + 1500) % 3329 = 3500 % 3329 = 171
    // (3328 + 1) % 3329 = 0
    let mut v1 = Vector::new(&[2000, 3328], 3329);
    let v2 = Vector::new(&[1500, 1], 3329);
    v1.add(v2);
    
    assert_eq!(v1.c, [171, 0], "L'addition modulaire a échoué");
}

#[test]
fn test_vector_dot_product_with_modulo() {
    // (1000 * 4 % 3329) + (2 * 500 % 3329)
    // (4000 % 3329) + (1000 % 3329) = 671 + 1000 = 1671
    let mut v1 = Vector::new(&[1000, 2], 3329);
    let v2 = Vector::new(&[4, 500], 3329);
    let result = v1.dot(v2);
    
    assert_eq!(result, 1671, "Le produit scalaire est faux !");
}

#[test]
fn test_vector_scalar_mult_modulo() {
    // (2000 * 2) % 3329 = 4000 % 3329 = 671
    // (3000 * 2) % 3329 = 6000 % 3329 = 2671
    let mut v = Vector::new(&[2000, 3000], 3329);
    v.scalar_mult(2);

    assert_eq!(v.c, [671, 2671], "La multiplication par un scalaire n'applique pas encore le modulo");
}

#[test]
fn test_matrix_matmul_modulo() {
    // Matrice 2x2 avec modulo 3329
    // Row 0, Col 0: (2000*2 + 1000*1) % 3329 = 5000 % 3329 = 1671
    // Row 0, Col 1: (2000*1 + 1000*2) % 3329 = 4000 % 3329 = 671
    let mut m1 = Matrix::new([[2000, 1000], [0, 1]], 3329);
    let m2 = Matrix::new([[2, 1], [1, 2]], 3329); 
    let res = m1.matmul(m2);
    
    assert_eq!(res.row_vecs[0].c[0], 1671, "Matmul R0C0 incorrect");
    assert_eq!(res.row_vecs[0].c[1], 671,  "Matmul R0C1 incorrect");
    assert_eq!(res.row_vecs[1].c[0], 1,    "Matmul R1C0 incorrect");
    assert_eq!(res.row_vecs[1].c[1], 2,    "Matmul R1C1 incorrect");
}


#[test]

fn test_montgomery_form(){

    let mut rng = rand::rng();
    for _ in 0..1000{

        let a : i32 = (rng.next_u32() % (1<<14) ).try_into().unwrap();
        let _m : i32 = (rng.next_u32() % (1<<14) ).try_into().unwrap();

        let mtg = MontgomeryForm::new(a);
        let res = mtg.mult(_m).to_standard();
        let true_res = (a*_m)%3329;
        assert_eq!(res,true_res, "reduction failed : {res} != {true_res}");
    }



}
}

