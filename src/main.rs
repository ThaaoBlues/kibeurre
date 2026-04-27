pub mod math_utils;
pub mod ntt;
fn main() {
    let v = math_utils::Vector::new(&[0;5],10);

    println!("{:?}",v);

}
