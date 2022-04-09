function Cm = model_cov(n, z_log, z_smooth_log)
%MODEL_COV 通过高斯分布生成模型协方差矩阵
%   z_smooth 模型平滑参数：强相关层间距（暂定为1σ）

l = [z_log(1); z_log(2:end)-z_log(1:end-1)] ./ 2; % 假设最后一层与上面完全不相关
Cm = diag(ones(1, n));
for i = 1:n
    for j = i-1:2:i+1
        if j > 0 && j <= n
            if i > j
                Cm(i, j) = exp(-(((l(i)+l(j)+sum(l(j+1:i-1)*2))^2)/(2*z_smooth_log^2)));
            else
                Cm(i, j) = exp(-(((l(i)+l(j)+sum(l(i+1:j-1)*2))^2)/(2*z_smooth_log^2)));
            end
        end
    end
end