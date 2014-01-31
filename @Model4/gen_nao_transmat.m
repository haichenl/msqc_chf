function res = gen_nao_transmat(obj)
% use nao class to generate transformation matrix 
obj.nao = nao(obj);
res = obj.nao.c_nmb;
end