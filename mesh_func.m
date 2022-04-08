function ovar = mesh_func(ivar, var_mesh)
%MESH_FUNC ������Լ����������
% 

len = length(ivar);
ovar = zeros(len, 1);
for i = 1:len
    [~, ind] = min(abs(var_mesh-ivar(i)));
    ovar(i) = var_mesh(ind);
end

end

