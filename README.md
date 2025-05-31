## 代码介绍 | Code Overview  
本算法针对广义梅森素数（形如 `2^p - k*2^q + 1`）提出一种**通用模约简算法**，支持所有素数的模运算。相较于传统蒙哥马利（Montgomery）和巴雷特（Barrett）算法，本方案采用更小规模的乘法器，在**领域专用架构（DSA）**中可实现更优的面积效率与能效表现。详细实验数据请参考论文：  
[**Fast Modular Reduction Algorithm and Reconfigurable Domain-Specific Architecture Design Based on Generalized Mersenne Primes**] 

**核心优势**  
- ✅ 支持广义梅森素数的高效模运算  
- ✅ 硬件友好设计：减少乘法器尺寸，优化面积与能效  
- ✅ 提供可重构架构设计，适配不同参数场景  

---  
This algorithm proposes a **generic modular reduction method** for Generalized Mersenne Primes (form: `2^p - k*2^q + 1`), supporting all prime numbers. Compared to traditional Montgomery and Barrett algorithms, our approach utilizes smaller multipliers, achieving superior area efficiency and energy performance in **Domain-Specific Architectures (DSA)**. For detailed benchmarks, refer to the paper:  
[**Fast Modular Reduction Algorithm and Reconfigurable Domain-Specific Architecture Design Based on Generalized Mersenne Primes**]  

**Key Advantages**  
- ✅ Efficient modular operations for Generalized Mersenne Primes  
- ✅ Hardware-friendly design: Reduced multiplier size, optimized area and energy  
- ✅ Reconfigurable architecture for diverse parameter scenarios  

---

## 文件结构 | File Structure  
### 目录说明 | Directory Description  
1. **`Generalized Mersenne_Chinese`**  
   - **中文注释代码版本**，包含以下内容：  
     - **C 模型** (`Generalized Mersenne_UTF8.cpp`)  
       - 实现广义梅森算法、蒙哥马利算法、巴雷特算法  
       - 通过修改主函数中的 `TEST_Q`, `TEST_X`, `TEST_Y` 可验证不同模数与输入数据的正确性  
       - ⚠️ 注意：需防范 32 位整型溢出问题  
     - **Matlab 遍历测试工程** (`kyber_test_UTF8.m`)  
       - 使用素数 `3329`（Kyber 算法参数）进行广义梅森算法遍历测试  
       - 测试时长：约 25 秒（AMD 5600X CPU）  

---  
2. **`Generalized Mersenne_English`**  
   - **English-annotated code version**, including:  
     - **C Model** (`Generalized Mersenne.cpp`)  
       - Implements Generalized Mersenne, Montgomery, and Barrett algorithms  
       - Modify `TEST_Q`, `TEST_X`, `TEST_Y` in the main function to validate correctness under different moduli and inputs  
       - ⚠️ Note: Prevent 32-bit integer overflow  
     - **Matlab Test Suite** (`kyber_test.m`)  
       - Exhaustive test with prime `3329` (Kyber algorithm parameter)  
       - Execution time: ~25 seconds (AMD 5600X CPU)  

---  
3. **`time_comparison.cpp`**  
   - Time comparison of a single operation for outputting Generalized Mersenne, Montgomery, and Barrett algorithms
   
---

## 作者 | Author  
- 邮箱：<xinyu_wang@smail.nju.edu.cn>  

---  
- Email: <xinyu_wang@smail.nju.edu.cn>  

---

**License**  
本代码允许商业使用、修改与私有化部署，但需事先联系版权方（xinyu_wang@smail.nju.edu.cn）并获得书面授权。未经许可不得移除原始声明或用于衍生品分发。
This code permits commercial use, modification, and proprietary deployment, provided that prior written authorization is obtained from the copyright holder (xinyu_wang@smail.nju.edu.cn). Removal of original notices or unlicensed distribution of derivatives is prohibited.

---

**更新日期 | Last Updated**  
2025年5月 | MAY 2025 # Generalized-Mersenne_test
