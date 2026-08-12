def PRNG(seed=42):                         // Linear Congruential Generator
    a = 1664525.0
    c = 1013904223.0
    m = 2**32
    seed = (a * seed + c) % m
    return float(seed/m)
