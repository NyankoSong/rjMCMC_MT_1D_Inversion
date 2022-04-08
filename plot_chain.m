% 绘制指定深度的电阻率迭代变化曲线

vec = zeros(N, 2);

for i = 1:N
    zk = model_cell{i, 1} - z_mesh(133);
    zk2 = model_cell{i, 1} - z_mesh(157);
    z_ind = find(zk >= 0, 1);
    z_ind2 = find(zk2 >= 0, 1);
    vec(i, 1) = model_cell{i, 2}(z_ind);
    vec(i, 2) = model_cell{i, 2}(z_ind2);
end

figure(1)
semilogy(vec(:, 1));

figure(2)
semilogy(vec(:, 2));

figure(3)
semilogy([model_cell{:, 5}])
