function ovar = mesh_func(ivar, var_mesh)
%MESH_FUNC 将数据约化到网格上
% 

len = length(ivar);
ovar = zeros(len, 1);
for i = 1:len
    [~, ind] = min(abs(var_mesh-ivar(i)));
    ovar(i) = var_mesh(ind);
end

end

