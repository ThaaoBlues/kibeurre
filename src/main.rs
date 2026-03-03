pub mod math_utils;

fn main() {
    let v = math_utils::Vector::new(&[0;5],10);

    println!("{:?}",v);

}
