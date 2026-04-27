#[derive(Copy)]
#[derive(Clone)]
#[derive(Debug)]
pub struct Vector<const VECTOR_SIZE : usize> {
    pub c : [i32;VECTOR_SIZE],
    m : i32
}

impl<const VECTOR_SIZE : usize> Vector<{VECTOR_SIZE}> {

    pub fn new( val : &[i32;VECTOR_SIZE],m : i32 ) -> Vector<VECTOR_SIZE>{
        Vector { c : *val, m : m}
    }

    fn size(&self)->usize{
        return self.c.len();
    }


    fn add(&mut self, v2 : Vector<VECTOR_SIZE>){
        

        //let mut tmp : [i32;VECTOR_SIZE] = [0;VECTOR_SIZE];

        for i in 0..VECTOR_SIZE {
            self.c[i] = (self.c[i] + v2.c[i]) % self.m;
        }

        //Vector::new( &tmp);

    }



    fn dot(&mut self, v2 : Vector<VECTOR_SIZE>) -> i32{
        let mut ret : i32 = 0;


        for i in 0..VECTOR_SIZE {
            let tmp = MontgomeryForm::new(self.c[i]);
            ret += tmp.mult(v2.c[i]).to_standard();
        }

        //println!("{:?} DOT {:?} = {}",self, v2,ret);
        ret % self.m
    }



    fn scalar_mult(&mut self,n : i32){

        for i in 0..VECTOR_SIZE {
            let tmp = MontgomeryForm::new(n);
            self.c[i] = tmp.mult(self.c[i]).to_standard();
        } 

    }


    fn get_odd_indexes<const HALF_SIZE: usize>(&self) -> Vector<HALF_SIZE>{
        
        assert_eq!(HALF_SIZE, VECTOR_SIZE/2);

        let mut v = Vector::new(&[0;HALF_SIZE], self.m);


        for i in 0..HALF_SIZE{ 

            v.c[i] = self.c[2*i+1];
            
        }

        return v;
    }


    fn get_even_indexes<const HALF_SIZE: usize>(&self) -> Vector<HALF_SIZE>{
        
        // rust is an useless piece of shit in its current state
        // you can't even perform calculations on genericity parameter
        // and put it as return
        assert_eq!(HALF_SIZE, VECTOR_SIZE/2);

        let mut v = Vector::new(&[0;HALF_SIZE], self.m);


        for i in 0..HALF_SIZE{ 

            v.c[i] = self.c[2*i];
            
        }

        return v;
    }

}



#[derive(Copy)]
#[derive(Clone)]
#[derive(Debug)]
pub struct Matrix<const RN : usize, const CN : usize>{
    row_vecs : [Vector<CN>;RN],
    col_vecs : [Vector<RN>;CN],
    m : i32
}


impl <const RN : usize, const CN : usize> Matrix<{RN},{CN}>{

    fn new(val : [[i32;CN];RN], m : i32)->Matrix<RN, CN>{


        // the array of arrays we get as parameter
        // is a rows representation of the matrix


        // build two representations to ease matrix multiplications
        // one where we store matrix by its rows
        // and the second by its columns

        // rows
        let mut r_v : [Vector<CN>;RN] = [Vector::new(&[0;CN],m);RN];

        // columns
        let mut c_v : [Vector<RN>;CN] = [Vector::new(&[0;RN],m);CN];

        for i in 0..CN{ 
            r_v[i] = Vector::new(&val[i],m);
        }



        for i in 0..CN{
            for j in 0..RN {
                c_v[i].c[j] = val[j][i];
            }
        }
        
        Matrix { col_vecs : c_v, row_vecs : r_v, m : m}

    }

    fn matmul(&mut self, B : Matrix<CN,RN>)->Matrix<CN,RN>{
        
        let mut ret : Matrix<CN,RN> = Matrix::new([[0;RN];CN],self.m);

        for i in 0..RN{
            
            for j in 0..CN {
                ret.row_vecs[i].c[j] = self.row_vecs[i].dot(B.col_vecs[j]);
                ret.col_vecs[j].c[i] = ret.row_vecs[i].c[j];
            }
        }


        ret


    }

    fn transpose(& self) -> Matrix<CN,RN>{
        let mut ret : Matrix<CN,RN> = Matrix::new([[0;RN];CN],self.m);
        ret.col_vecs = self.row_vecs;
        ret.row_vecs = self.col_vecs;

        ret
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
    pub fn new(n: i32) -> MontgomeryForm{

        
        MontgomeryForm { a: (n<<12) % 3329, n, q :3329, r_pow:12, r : 4096,r_1 : 2704, k: 3327 }

    }

    pub fn set_a(&mut self, a : i32){
        self.a = a;
    }

    pub fn set_n(&mut self, n : i32){
        self.n = n;
    }

    pub fn get_a(&self)-> i32{
        return self.a;
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

        return MontgomeryForm { a:c_b, n: 0, q: self.q, r_pow:self.r_pow,r: self.r, r_1: self.r_1, k: self.k};
        
    }

    // goes back to [0;q] from montgomery form
    pub fn reduction(&mut self) -> i32{
        let s : i32  = i32::try_from((i64::from(self.a) * i64::from(self.k)) %  i64::from(self.r)).unwrap();


        let t : i32= dbg!(self.a + s*self.q);


        let u : i32 = t >> self.r_pow; // u = t/r

        println!("{0} = {1} ?",(u % self.q), i32::try_from((i64::from(self.a)*i64::from(self.r_1)) % i64::from(self.q)).unwrap() );

        
        if(u < self.q){
            self.n = u;
            return dbg!(u);
        }else{
            self.n = u-self.q;
            return self.n;
        }

        
    }



    pub fn to_standard(&mut self)->i32{
        (self.a * self.r_1) % self.q
    }
}







    #[cfg(test)]
    mod tests {
        use super::*;
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
        let a = MontgomeryForm::new(3948);

        assert_eq!(a.mult(124).to_standard(),189);

    }
}