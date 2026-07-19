use kibeurre::core::{encrypt, decrypt, generate_key_pair, generate_noise_polyvector};
use kibeurre::parameters::{ETA_1, m};
use kibeurre::math_utils::Vector;

fn enc(){
    let (pub_key, priv_key) = generate_key_pair();

    let r = generate_noise_polyvector(ETA_1);

    let c = encrypt(pub_key.A, pub_key.t, r, Vector::new(&[1;256], m));
}


fn enc_dec(){
    let (pub_key, priv_key) = generate_key_pair();

    let r = generate_noise_polyvector(ETA_1);

    let c = encrypt(pub_key.A, pub_key.t, r, Vector::new(&[1;256], m));
    println!("u : {:?}", c.u);
    println!("v : {:?}", c.v);

    let msg = decrypt(&c, priv_key.s);
    println!("Decrypted message: {:?}", msg);
}

fn main(){
    enc();
    //enc_dec();
}

