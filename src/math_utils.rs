#[derive(Copy)]
#[derive(Clone)]
#[derive(Debug)]
pub struct Vector<const VECTOR_SIZE : usize> {
    c : [i32;VECTOR_SIZE],
    m : i32
}

impl<const VECTOR_SIZE : usize> Vector<{VECTOR_SIZE}> {

    pub fn new( val : &[i32;VECTOR_SIZE],m : i32 ) -> Vector<VECTOR_SIZE>{
        Vector { c : *val, m : m}
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
            ret += self.c[i] * v2.c[i];
        }

        //println!("{:?} DOT {:?} = {}",self, v2,ret);

        ret % self.m
    }


    fn scalar_mult(&mut self,n : i32){

        for i in 0..VECTOR_SIZE {
           self.c[i] = (self.c[i]*n) % self.m;
        }

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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_vector_add_with_modulo() {
        // Test du dépassement de modulo : (8 + 5) % 10 = 3
        let mut v1 = Vector::new(&[8, 15], 10);
        let v2 = Vector::new(&[5, 2], 10);
        v1.add(v2);
        
        assert_eq!(v1.c, [3, 7], "L'addition modulaire a échoué");
    }

    #[test]
    fn test_vector_dot_product_with_modulo() {
        // (3*4 % 10) + (5*5 % 10) = 2 + 5 = 7
        let mut v1 = Vector::new(&[3, 5], 10);
        let v2 = Vector::new(&[4, 5], 10);
        let result = v1.dot(v2);
        
        assert_eq!(result, 7, "Le produit scalaire est faux !");
    }

    #[test]
    fn test_vector_scalar_mult_modulo() {
        let mut v = Vector::new(&[6, 2], 10);
        v.scalarMult(2);

        assert_eq!(v.c, [2, 4], "La multiplication par un scalaire n'applique pas encore le modulo");
    }

    #[test]
    fn test_matrix_matmul_modulo() {
        // Matrice 2x2 avec modulo 5
        // [2, 3] * [1, 4] -> (2*1 % 5) + (3*4 % 5) = 2 + 2 = 4
        let mut m1 = Matrix::new([[2, 3], [0, 1]], 5);
        let m2 = Matrix::new([[1, 4], [0, 1]], 5); 
        let res = m1.matmul(m2);
        
        assert_eq!(res.row_vecs[0].c[0], 2, "La multiplication matricielle est fausse !");
        assert_eq!(res.row_vecs[0].c[1], 1, "La multiplication matricielle est fausse !");
        assert_eq!(res.row_vecs[1].c[0], 0, "La multiplication matricielle est fausse !");
        assert_eq!(res.row_vecs[1].c[1], 1, "La multiplication matricielle est fausse !");

    }
}