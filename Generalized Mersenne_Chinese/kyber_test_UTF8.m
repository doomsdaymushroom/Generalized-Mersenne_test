% Kyber q = 2^12 - 3*2^8 + 1 (reference value)
clear
clc
close all

% 算法参数初始化
m = 6;          % 模数参数1
n = 8;          % 模数参数2
k = 3;          % 系数参数
shift_para = floor(log(k)/log(2));  % 预计算移位参数（当前未使用）
Q = 2^(2*m) - k*2^n + 1;           % 广义梅森素数模数

% 性能统计初始化
r2 = 0;         % 保留变量（当前未使用）
tic
loop_max = 0;   % 记录最大循环次数
value_min = 0;  % 记录最小余数值
data_index = 0; % 数据索引（当前未使用）

% 主计算循环
for a = 1:3328
    for b = 1:3328
        x = a * b;
        r1 = x;
        loop_cnt = 0;
        
        % 快速约简阶段
        while (r1 > 3 * Q)
            % 移位操作组
            t1 = bitshift(r1, -(2*m));        % 第一级移位：2^2m
            t2 = bitshift(t1, -(2*m - n));    % 第二级移位：2^(2m-n)
            
            % 系数合成
            t = t1 + 3 * t2;                 % 带系数的移位组合
            
            % 快速乘法实现
            tmp1 = t * bitshift(Q, -n);       % 低位乘法项
            t_mul_q = bitshift(tmp1, n) + t;  % 重组乘法结果
            
            % 余数更新
            r1 = r1 - t_mul_q;
            loop_cnt = loop_cnt + 1;
        end
        
        % 最终模约简
        res = r1;
        while res > Q
            res = res - Q;
        end
        
        % 统计信息更新
        if loop_cnt > loop_max
            loop_max = loop_cnt;
        end
        if res < value_min
            value_min = res;
        end
        
        % 结果验证
        demo = mod(x, Q);
        if demo ~= res
            error("校验错误 at a=%d b=%d", a, b);
        end
    end
end
toc

% Results display
fprintf('\n所有测试用例验证成功! All test cases verified successfully!\n');
fprintf('最大循环次数: %d (Maximum loop count)\n', loop_max);