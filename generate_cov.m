function C = generate_cov(n, y_log, k)
%GENERATE_COV 通过高斯分布生成模型协方差矩阵
%   k 强相关层/点间距（暂定为1σ）

l = [y_log(1); y_log(2:end)-y_log(1:end-1)] ./ 2;
C = diag(ones(1, n));
for i = 1:n
    for j = i-1:2:i+1
        if j > 0 && j <= n
            if i > j
                C(i, j) = exp(-((l(i)+l(j)+sum(l(j+1:i-1)*2))^2)/(2*k^2));
            else
                C(i, j) = exp(-((l(i)+l(j)+sum(l(i+1:j-1)*2))^2)/(2*k^2));
            end
        end
    end
end