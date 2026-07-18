
fn main(){
    let (pub_key, priv_key) = core::generate_key_pair();

    let r = core::generate_noise_polyvector(parameters::ETA_1);

    let c = core::encrypt(pub_key.A, pub_key.t, r, Vector::new(&[1;256], parameters::m));
    println!("u : {:?}", c.u);
    println!("v : {:?}", c.v);

    let msg = core::decrypt(&c, priv_key.s);
    println!("Decrypted message: {:?}", msg);
}

