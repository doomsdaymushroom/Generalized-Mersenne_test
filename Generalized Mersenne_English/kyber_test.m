% Kyber q = 2^12 - 3*2^8 + 1 (reference value)
clear
clc
close all

% Algorithm parameters initialization
m = 6;          % Modulus parameter 1
n = 8;          % Modulus parameter 2
k = 3;          % Coefficient parameter
shift_para = floor(log(k)/log(2));  % Precomputed shift parameter (currently unused)
Q = 2^(2*m) - k*2^n + 1;           % Generalized Mersenne-like prime modulus

% Performance tracking initialization
r2 = 0;         % Reserved variable (currently unused)
tic
loop_max = 0;   % Track maximum loop count
value_min = 0;  % Track minimum remainder value
data_index = 0; % Data index (currently unused)

% Main computation loop
for a = 1:3328
    for b = 1:3328
        x = a * b;
        r1 = x;
        loop_cnt = 0;
        
        % Fast reduction phase
        while (r1 > 2 * Q)
            % Shift operations group
            t1 = bitshift(r1, -(2*m));       % First-stage shift: 2^2m bits
            t2 = bitshift(t1, -(2*m - n));   % Second-stage shift: 2^(2m-n) bits
            
            % Coefficient synthesis
            t = t1 + t2;                % Shift combination with coefficients
            
            % Fast multiplication implementation
            tmp1 = t * bitshift(Q, -n);      % Lower bits multiplication term
            t_mul_q = bitshift(tmp1, n) + t; % Reconstruct full multiplication result
            
            % Remainder update
            r1 = r1 - t_mul_q;
            loop_cnt = loop_cnt + 1;
        end
        
        % Final modulus reduction
        res = r1;
        while res > Q
            res = res - Q;
        end
        
        % Statistics update
        if loop_cnt > loop_max
            loop_max = loop_cnt;
        end
        if res < value_min
            value_min = res;
        end
        
        % Result validation
        demo = mod(x, Q);
        if demo ~= res
            error("Verification failed at a=%d b=%d", a, b);
        end
    end
end
toc

% Results display
fprintf('\n所有测试用例验证成功! All test cases verified successfully!\n');
fprintf('最大循环次数: %d (Maximum loop count)\n', loop_max);