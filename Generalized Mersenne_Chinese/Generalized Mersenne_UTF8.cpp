#include <iostream>
#include <cstdint>
#include <cmath>
#include <stdexcept>

// 使用更安全的固定宽度整数类型
using int32 = int32_t;
using int64 = int64_t;
using uint32 = uint32_t;
using uint64 = uint64_t;

// 质数分解结构体：存储广义梅森素数分解参数
struct PrimeDecomposition {
    int exponent_p;    // 2^p 项指数
    int coefficient_k; // 线性系数k
    int shift_q;       // 位移参数q
    uint32 modulus_R;  // 模数基数R=2^p
    bool is_valid;     // 分解有效性标志

    explicit PrimeDecomposition(int p = -1, int k = -1, int q = -1,
        uint32 R = 0, bool valid = false)
        : exponent_p(p), coefficient_k(k), shift_q(q),
        modulus_R(R), is_valid(valid) {}
};

// 辅助函数声明
static bool IsPowerOfTwo(uint32 num) noexcept;
uint32 CalculateMontgomeryInverse(uint32 q, uint32 R);
uint64 CalculateBarrettParameter(uint32 q, uint32 R);

/* 广义梅森素数分解
 * 参数: x - 需要分解的质数
 * 返回值: 分解参数结构体
 * 算法思想: 将质数分解为 2^p - k*2^q + 1 的形式
 */
PrimeDecomposition DecomposePrime(uint32 x) {
    constexpr uint32 MIN_PRIME = 2;
    if (x == MIN_PRIME) { // 处理最小质数特殊情况
        return PrimeDecomposition(1, 1, 0, 2, true);
    }

    uint32 temp = x - 1;
    if (IsPowerOfTwo(temp)) { // 形式为2^m +1的特殊情况
        int m = log2(temp);
        return PrimeDecomposition(m, 0, 1, 1U << (m + 1), true);
    }

    // 常规分解流程
    int q = 0;
    uint32 s = temp;
    while ((s & 1) == 0) { // 提取2的幂因子
        ++q;
        s >>= 1;
    }

    // 寻找大于s的最小2的幂
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

/* 广义梅森模约简算法
 * 参数: a,b - 输入操作数, Q - 模数
 * 返回值: (a*b) mod Q
 * 算法特点: 利用特殊数形式优化模运算
 */
uint32 GeneralizedMersenneReduce(uint32 a, uint32 b, uint32 Q) {
    const PrimeDecomposition params = DecomposePrime(Q);
    if (!params.is_valid) {
        throw std::invalid_argument("Invalid prime decomposition");
    }

    uint64 product = static_cast<uint64>(a) * b;
    uint64 residual = product;

    // 优化约简过程
    while (residual > 2 * Q) {                                  //已通过论文完成硬件友好化设计
        uint64 step1 = 0; 
        uint64 Cres;
        const int shift1 = params.exponent_p;
        const int shift2 = 2 * params.exponent_p - params.shift_q;

        if (Q > (1 << params.exponent_p)){
            Cres = (residual >> shift1) - (residual >> shift2);
        }else{
            Cres = (residual >> shift1) + params.coefficient_k * (residual >> shift2); //硬件中由于采用截位方案，无需该步骤的判断，且该步骤已硬件友好为并行方案。
        }
        step1 += Cres;

        const uint64 step2 = ((step1 * (Q >> params.shift_q)) << params.shift_q);
        residual -= step2 + step1;
    }

    return static_cast<uint32>((residual >= Q) ? (residual - Q) : residual);
}

/* 蒙哥马利模约简核心函数
 * 参数: res_m - 预乘结果, q - 模数, inv - 预计算逆元, R - 蒙哥马利基数
 * 返回值: 约简结果
 * 算法复杂度: O(1)
 */
uint32 MontgomeryReduce(uint64 res_m, uint32 q, uint32 inv, uint32 R_shift) {
    const uint32 R = 1U << R_shift;
    const uint32 m = (static_cast<uint32>(res_m) * inv) % R;
    const uint64 y = (res_m + static_cast<uint64>(m) * q) >> R_shift;
    return (y >= q) ? y - q : static_cast<uint32>(y);
}

// 计算蒙哥马利逆元
uint32 CalculateMontgomeryInverse(uint32 q, uint32 R) {
    for (uint32 i = 1; i < R; ++i) {
        if ((static_cast<uint64>(i) * q) % R == R - 1) {
            return i;
        }
    }
    throw std::domain_error("Montgomery inverse not found");
}

/* 蒙哥马利模乘算法
 * 参数: a,b - 操作数, q - 模数
 * 返回值: (a*b) mod q
 */
uint32 MontgomeryMultiply(uint32 a, uint32 b, uint32 q) {
    const PrimeDecomposition params = DecomposePrime(q);
    if (!params.is_valid || params.modulus_R == 0) {
        throw std::invalid_argument("Invalid modulus for Montgomery");
    }

    const uint32 R = params.modulus_R;
    const uint32 R_shift = params.exponent_p;
    const uint32 inv = CalculateMontgomeryInverse(q, R);

    // 转换为蒙哥马利形式
    const uint32 a_m = (static_cast<uint64>(a) * R) % q;
    const uint32 b_m = (static_cast<uint64>(b) * R) % q;

    // 核心计算
    const uint64 product = static_cast<uint64>(a_m) * b_m;
    uint32 res = MontgomeryReduce(product, q, inv, R_shift);

    // 二次约简并转换回普通域
    res = MontgomeryReduce(res, q, inv, R_shift);
    return res % q;
}

/* Barrett模约简算法
 * 参数: a,b - 操作数, q - 模数, param_p - 预计算参数
 * 返回值: (a*b) mod q
 */
uint32 BarrettReduce(uint32 a, uint32 b, uint32 q, uint64 param_p) {
    const uint64 product = static_cast<uint64>(a) * b;
    const uint64 quotient = (product * param_p) >> (2 * static_cast<uint32>(log2(q) + 1));
    uint64 result = product - quotient * q;

    // 最终调整
    while (result >= q) result -= q;
    return static_cast<uint32>(result);
}

// 辅助函数实现
bool IsPowerOfTwo(uint32 num) noexcept {
    return num && !(num & (num - 1));
}

uint64 CalculateBarrettParameter(uint32 q, uint32 R) {
    return (R * R) / q;
}

// 验证函数
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
    // 测试用例
    constexpr uint32 TEST_Q = 1073479681; // 典型安全素数  Kyber:3329/7681 NewHope:12289 NTRU:65537 Dilithum:8380417 qTESLA v2.0:8404993 HPS:1073479681
    constexpr uint32 TEST_X = 412223;     // TEST_X < TEST_Q
    constexpr uint32 TEST_Y = 412132;     // TEST_Y < TEST_Q

    std::cout << "=== Modular Reduction Validation ===\n";
    RunVerification(TEST_X, TEST_Y, TEST_Q);

    // 边界测试
    std::cout << "=== Boundary Case Testing ===\n";
    RunVerification(TEST_Q - 1, TEST_Q - 1, TEST_Q); // 最大输入测试
    RunVerification(0, 12345, TEST_Q);          // 零输入测试

    return 0;
}