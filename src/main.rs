pub mod math_utils;
pub mod ntt;
fn main() {
    ntt::compute_zeta_table();
    ntt::compute_inv_zeta_table();
}
