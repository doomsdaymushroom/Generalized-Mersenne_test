#include <iostream>
#include <cstdint>
#include <cmath>
#include <stdexcept>

// Use safer fixed-width integer types
using int32 = int32_t;
using int64 = int64_t;
using uint32 = uint32_t;
using uint64 = uint64_t;

// Prime decomposition struct: stores parameters for generalized Mersenne prime decomposition
struct PrimeDecomposition {
    int exponent_p;    // 2^p term exponent
    int coefficient_k; // linear coefficient k
    int shift_q;       // shift parameter q
    uint32 modulus_R;  // modulus base R=2^p
    bool is_valid;     // decomposition validity flag

    explicit PrimeDecomposition(int p = -1, int k = -1, int q = -1,
        uint32 R = 0, bool valid = false)
        : exponent_p(p), coefficient_k(k), shift_q(q),
        modulus_R(R), is_valid(valid) {}
};

// Helper function declarations
static bool IsPowerOfTwo(uint32 num) noexcept;
uint32 CalculateMontgomeryInverse(uint32 q, uint32 R);
uint64 CalculateBarrettParameter(uint32 q, uint32 R);

/* Generalized Mersenne prime decomposition
 * Parameters: x - prime to decompose
 * Returns: decomposition parameters struct
 * Algorithm: Decompose prime into the form 2^p - k*2^q + 1
 */
PrimeDecomposition DecomposePrime(uint32 x) {
    constexpr uint32 MIN_PRIME = 2;
    if (x == MIN_PRIME) { // Handle special case for the smallest prime
        return PrimeDecomposition(1, 1, 0, 2, true);
    }

    uint32 temp = x - 1;
    if (IsPowerOfTwo(temp)) { // Special case of 2^m +1 form
        int m = log2(temp);
        return PrimeDecomposition(m, 0, 1, 1U << (m + 1), true);
    }

    // Standard decomposition process
    int q = 0;
    uint32 s = temp;
    while ((s & 1) == 0) { // Extract power of two factors
        ++q;
        s >>= 1;
    }

    // Find smallest power of two greater than s
    uint32 t = 1U << static_cast<uint32>(ceil(log2(s)));
    int p = q + static_cast<int>(log2(t));

    return PrimeDecomposition(
        p,
        static_cast<int>(t - s),
        q,
        1U << p,
        true
    );
}

/* Generalized Mersenne modulus reduction algorithm
 * Parameters: a,b - input operands, Q - modulus
 * Returns: (a*b) mod Q
 * Features: Optimize modulus operation using special number form
 */
uint32 GeneralizedMersenneReduce(uint32 a, uint32 b, uint32 Q) {
    const PrimeDecomposition params = DecomposePrime(Q);
    if (!params.is_valid) {
        throw std::invalid_argument("Invalid prime decomposition");
    }

    uint64 product = static_cast<uint64>(a) * b;
    uint64 residual = product;

    // Optimized reduction process
    while (residual > 2 * Q) {                                  // Hardware-friendly design completed via paper analysis
        uint64 step1 = 0; 
        uint64 Cres;
        const int shift1 = params.exponent_p;
        const int shift2 = 2 * params.exponent_p - params.shift_q;

        if (Q > (1 << params.exponent_p)){
            Cres = (residual >> shift1) - (residual >> shift2);
        }else{
            Cres = (residual >> shift1) + params.coefficient_k * (residual >> shift2); 
        }/* This part of the hardware has been optimized and replaced with truncation in the hardware, without the need for multiplication and branching conditions 
          * (as this step only requires obtaining an approximate result)
          */
         
        step1 += Cres;

        const uint64 step2 = ((step1 * (Q >> params.shift_q)) << params.shift_q);
        residual -= step2 + step1;
    }

    return static_cast<uint32>((residual >= Q) ? (residual - Q) : residual);
}

/* Montgomery modulus reduction core function
 * Parameters: res_m - pre-multiplied result, q - modulus, inv - precomputed inverse, R - Montgomery base
 * Returns: reduction result
 * Complexity: O(1)
 */
uint32 MontgomeryReduce(uint64 res_m, uint32 q, uint32 inv, uint32 R_shift) {
    const uint32 R = 1U << R_shift;
    const uint32 m = (static_cast<uint32>(res_m) * inv) % R;
    const uint64 y = (res_m + static_cast<uint64>(m) * q) >> R_shift;
    return (y >= q) ? y - q : static_cast<uint32>(y);
}

// Calculate Montgomery inverse
uint32 CalculateMontgomeryInverse(uint32 q, uint32 R) {
    for (uint32 i = 1; i < R; ++i) {
        if ((static_cast<uint64>(i) * q) % R == R - 1) {
            return i;
        }
    }
    throw std::domain_error("Montgomery inverse not found");
}

/* Montgomery modular multiplication algorithm
 * Parameters: a,b - operands, q - modulus
 * Returns: (a*b) mod q
 */
uint32 MontgomeryMultiply(uint32 a, uint32 b, uint32 q) {
    const PrimeDecomposition params = DecomposePrime(q);
    if (!params.is_valid || params.modulus_R == 0) {
        throw std::invalid_argument("Invalid modulus for Montgomery");
    }

    const uint32 R = params.modulus_R;
    const uint32 R_shift = params.exponent_p;
    const uint32 inv = CalculateMontgomeryInverse(q, R);

    // Convert to Montgomery form
    const uint32 a_m = (static_cast<uint64>(a) * R) % q;
    const uint32 b_m = (static_cast<uint64>(b) * R) % q;

    // Core computation
    const uint64 product = static_cast<uint64>(a_m) * b_m;
    uint32 res = MontgomeryReduce(product, q, inv, R_shift);

    // Second reduction and convert back to normal domain
    res = MontgomeryReduce(res, q, inv, R_shift);
    return res % q;
}

/* Barrett modulus reduction algorithm
 * Parameters: a,b - operands, q - modulus, param_p - precomputed parameter
 * Returns: (a*b) mod q
 */
uint32 BarrettReduce(uint32 a, uint32 b, uint32 q, uint64 param_p) {
    const uint64 product = static_cast<uint64>(a) * b;
    const uint64 quotient = (product * param_p) >> (2 * static_cast<uint32>(log2(q) + 1));
    uint64 result = product - quotient * q;

    // Final adjustment
    while (result >= q) result -= q;
    return static_cast<uint32>(result);
}

// Helper function implementation
bool IsPowerOfTwo(uint32 num) noexcept {
    return num && !(num & (num - 1));
}

uint64 CalculateBarrettParameter(uint32 q, uint32 R) {
    return (R * R) / q;
}

// Validation function
void RunVerification(uint32 x, uint32 y, uint32 Q) {
    try {
        const uint64 golden = (static_cast<uint64>(x) * y) % Q;
        std::cout << "Golden reference: " << golden << "\n";

        const uint32 mersenne = GeneralizedMersenneReduce(x, y, Q);
        std::cout << "Generalized Mersenne: " << mersenne
            << (mersenne == golden ? " √ " : " × ") << "\n";

        const uint32 montgomery = MontgomeryMultiply(x, y, Q);
        std::cout << "Montgomery: " << montgomery
            << (montgomery == golden ? " √ " : " × ") << "\n";

        const PrimeDecomposition params = DecomposePrime(Q);
        const uint64 barrett_param = CalculateBarrettParameter(Q, params.modulus_R);
        const uint32 barrett = BarrettReduce(x, y, Q, barrett_param);
        std::cout << "Barrett: " << barrett
            << (barrett == golden ? " √ " : " × ") << "\n\n";
    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
    }
}

int main() {
    // Test cases
    constexpr uint32 TEST_Q = 1073479681; // Typical security primes  Kyber:3329/7681 NewHope:12289 NTRU:65537 Dilithum:8380417 qTESLA v2.0:8404993 HPS:1073479681
    constexpr uint32 TEST_X = 412223;
    constexpr uint32 TEST_Y = 412132;

    std::cout << "=== Modular Reduction Validation ===\n";
    RunVerification(TEST_X, TEST_Y, TEST_Q);

    // Boundary testing
    std::cout << "=== Boundary Case Testing ===\n";
    RunVerification(TEST_Q - 1, TEST_Q - 1, TEST_Q); // Maximum input test
    RunVerification(0, 12345, TEST_Q);          // Zero input test

    return 0;
}
