import numpy as np

def print_step(name, poly):
    print(f"\n{name}:")
    print(f"First few coeffs: {poly.coeffs[:4]}")
    print(f"Scale: {poly.scale}")
    magnitude = np.max(np.abs(poly.coeffs))
    print(f"Max magnitude: {magnitude}")
    print(f"log2(max): {np.log2(magnitude) if magnitude > 0 else '-inf'}")

class RingPolynomial:
    def __init__(self, coeffs, N, modulus):
        self.N = N
        self.modulus = modulus
        # Safely handle negative numbers in modulo
        self.coeffs = np.remainder(np.array(coeffs, dtype=np.int64)[:N], modulus)
        self.scale = 1.0
    
    def __add__(self, other):
        if self.modulus != other.modulus:
            raise ValueError(f"Modulus mismatch: {self.modulus} vs {other.modulus}")
        result = np.remainder(self.coeffs + other.coeffs, self.modulus)
        poly = RingPolynomial(result, self.N, self.modulus)
        poly.scale = self.scale
        return poly
    
    def __mul__(self, other):
        if self.modulus != other.modulus:
            raise ValueError(f"Modulus mismatch: {self.modulus} vs {other.modulus}")
            
        # Use remainder instead of mod for proper handling of negatives
        conv = np.zeros(2 * self.N, dtype=np.int64)
        for i in range(len(self.coeffs)):
            for j in range(len(other.coeffs)):
                prod = np.remainder(self.coeffs[i] * other.coeffs[j], self.modulus)
                conv[i + j] = np.remainder(conv[i + j] + prod, self.modulus)
        
        result = np.zeros(self.N, dtype=np.int64)
        for i in range(len(conv)):
            if i < self.N:
                result[i] = np.remainder(result[i] + conv[i], self.modulus)
            else:
                result[i - self.N] = np.remainder(result[i - self.N] - conv[i], self.modulus)
        
        poly = RingPolynomial(result, self.N, self.modulus)
        poly.scale = self.scale * other.scale
        return poly

def encode_complex(values, scale, N, modulus):
    slots = N // 2
    evals = np.zeros(N, dtype=np.complex128)
    
    # Map values to evaluation points
    for i, v in enumerate(values):
        if i >= slots:
            break
        evals[i] = v
        evals[N-1-i] = np.conj(v)
    
    print("\nEvaluations:", evals)
    
    # Get polynomial coefficients
    coeffs = np.fft.ifft(evals) * N
    print("\nUnscaled coefficients:", coeffs)
    
    # Scale and round to integers
    scaled_coeffs = np.round(np.real(coeffs) * scale).astype(np.int64)
    print("\nScaled coefficients:", scaled_coeffs)
    
    poly = RingPolynomial(scaled_coeffs, N, modulus)
    poly.scale = float(scale)
    
    print_step("After encoding", poly)
    return poly

def generate_keys(N, q):
    # Generate sparse ternary secret key
    h = N // 4  # Hamming weight
    s = np.zeros(N, dtype=np.int64)
    positions = np.random.choice(N, h, replace=False)
    s[positions] = np.random.choice([-1, 1], h)
    s_poly = RingPolynomial(s, N, q)
    print_step("Secret key", s_poly)
    
    # Sample a and e
    a = np.random.randint(-q//4, q//4, size=N, dtype=np.int64)
    e = np.random.normal(0, 0.001, size=N).astype(np.int64)  # 0.01 -> 0.001
    
    a_poly = RingPolynomial(a, N, q)
    e_poly = RingPolynomial(e, N, q)
    
    # -a*s + e
    neg_one = RingPolynomial([-1] + [0]*(N-1), N, q)
    temp = a_poly * s_poly
    print_step("a*s", temp)
    temp = temp * neg_one
    print_step("-a*s", temp)
    b_poly = temp + e_poly
    print_step("Public key b=-a*s+e", b_poly)
    
    return s_poly, (a_poly, b_poly)

def encrypt(m_poly, public_key, q):
    a_poly, b_poly = public_key
    
    # Small error terms
    e1 = np.random.normal(0, 0.01, size=m_poly.N).astype(np.int64)
    e2 = np.random.normal(0, 0.01, size=m_poly.N).astype(np.int64)
    e1_poly = RingPolynomial(e1, m_poly.N, q)
    e2_poly = RingPolynomial(e2, m_poly.N, q)
    e1_poly.scale = m_poly.scale
    e2_poly.scale = m_poly.scale
    
    print_step("Message to encrypt", m_poly)
    
    c0 = b_poly * m_poly + e1_poly
    print_step("c0", c0)
    
    c1 = a_poly * m_poly + e2_poly
    print_step("c1", c1)
    
    return (c0, c1)

def decrypt(cipher, s_poly, q):
    c0, c1 = cipher
    print_step("c0 in decrypt", c0)
    print_step("c1 in decrypt", c1)
    
    temp = c1 * s_poly
    print_step("c1*s", temp)
    
    result = c0 + temp
    print_step("c0 + c1*s", result)
    
    return result

def decode_complex(poly, modulus):
    # Center coefficients around zero
    coeffs = poly.coeffs.astype(np.float64)
    coeffs = np.where(coeffs > modulus//2, coeffs - modulus, coeffs)
    print("\nCentered coefficients:", coeffs)
    
    # Unscale
    complex_coeffs = coeffs / poly.scale
    print("\nUnscaled coefficients:", complex_coeffs)
    
    # FFT to get values back
    evals = np.fft.fft(complex_coeffs) / poly.N
    print("\nEvaluations:", evals)
    
    # Take first half of slots
    result = []
    slots = len(coeffs) // 2
    for i in range(slots):
        result.append(evals[i])
    return result

if __name__ == "__main__":
    # Parameters
    N = 8
    # Use smaller prime modulus initially
    q = 2**20 - 3  # Mersenne-like prime
    scale = 2**10  # Smaller scale to match modulus
    
    print("\nParameters:")
    print(f"N = {N}")
    print(f"log2(q) ≈ {np.log2(q):.2f}")
    print(f"log2(scale) = {np.log2(scale):.2f}")
    
    values = [1 + 2j, 3 - 1j]
    print("\nOriginal values:", values)
    
    s_poly, public_key = generate_keys(N, q)
    print("\nKey generation successful")
    
    m_poly = encode_complex(values, scale, N, q)
    print("\nEncoding successful")
    
    cipher = encrypt(m_poly, public_key, q)
    print("\nEncryption successful")
    
    decrypted_poly = decrypt(cipher, s_poly, q)
    print("\nDecryption successful") 
    
    decrypted_values = decode_complex(decrypted_poly, q)
    print("\nDecoded values:", decrypted_values)
    
    errors = [abs(complex(1,2) - decrypted_values[0]),
              abs(complex(3,-1) - decrypted_values[1])]
    print("\nAbsolute errors:", errors)
