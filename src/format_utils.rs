use crate::math_utils::{Vector,empty_vector};
use crate::parameters::{n, m};


/*
Converts a string into a vector of Vector<n> where each Vector<n> represents a chunk of the string.
The string is split into bytes chunks of size n.
So each vector actually contains n bytes of the string.
*/
pub fn string_to_vectors(input: &str) -> Vec<Vector<n>> {
    let bytes = input.as_bytes();
    let mut vectors = Vec::new();


    for chunk in bytes.chunks(n) {
        let mut vector = empty_vector(); // filled with zeros, will be used as padding
        for (i, &byte) in chunk.iter().enumerate() {

            for j in 0..7 {
                let bit = (byte >> j) & 1;
                vector.set(i * 8 + j, bit as i32);
            }
        }
        vectors.push(vector);
    }
    return vectors;
}



pub fn vectors_to_string(vectors: Vec<Vector<n>>) -> String {
    let mut bytes = Vec::new();
    for vector in vectors {

        for i in 0..32 {
            let mut byte = 0;
        
            
            for j in 0..7 {
                
                let bit = vector.get(i * 8 + j) as u8;
                if bit != 0 {
                    byte |= bit << j;
                }
            }

            // NULL byte => string end
            // do not append it as String constructor will
            if byte == 0 {
                break;
            }else{
                bytes.push(byte);
            }
            
        }
    }
    return String::from_utf8(bytes).unwrap_or_else(|_| String::from("Invalid UTF-8"));
}