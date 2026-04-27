# KIBEURRE

Kibeurre is a naive kyber implementation written in Rust. The goal of this projet is to make me work on lattices/LWE and Rust ofc :)



## math_utils.rs

Contains all basic vector and matrices operations in Z/qZ :
- vectors addition
- vector scalar multiplication
- matrix multiplication (not optimised)
- Montgomery form

### Multiplications in Z/3329Z 
- I use Montgomery method to speedup modulus computation when multiplying two numbers (sadly not by much)