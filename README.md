# KIBEURRE

Kibeurre is a naive kyber implementation written in Rust. The goal of this projet is to make me work on lattices/LWE and Rust ofc :)



## math_utils.rs

Contains all basic vector and matrices operations in Z/qZ :
- vectors addition
- vector scalar multiplication
- matrix multiplication (not optimised)
- Montgomery form
- 


# Kyber core principles

- [Cryptograhy 101 Kyber Course](https://www.youtube.com/watch?v=9NKm84vKALc&list=PLA1qgQLL41SSUOHlq8ADraKKzv47v2yrF)

- [Kyber repository](https://github.com/pq-crystals/kyber)



# NTT

I first learned NTT throught theses ressources : 
- [Satriawan, Ardianto, et al. « A Complete Beginner Guide to the Number Theoretic Transform (NTT) ». nᵒ 2024/585, 2024, Cryptology ePrint Archive. Cryptology ePrint Archive (eprint.iacr.org), https://eprint.iacr.org/2024/585.
](2024-585.pdf)

- [Cryptograhy 101 Kyber Course](https://youtu.be/ey1ND_xPITw)

- [Reducible FFT video](https://youtu.be/h7apO7q16V0)

Then implemented a naive recursive version and finally derived my algorithm from the "in place" implemetation from kyber round 3 NIST submission.
In the reference implementation pre-computed arrays of zeta powers and inverse powers are in Montgomery form. 
I prefered to let them be in the standard one and compute the Montgomery form on the fly, as i don't have another goal than self-learning and remembering easily is. Which speed optimisation can sometimes obstruct.



### Multiplications in Z/3329Z 
- I use Montgomery method to speedup modulus computation when multiplying two numbers (sadly not by much)