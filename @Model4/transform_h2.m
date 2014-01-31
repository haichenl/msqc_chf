function res = transform_h2(obj, H2, transmat)
n_basisuse = size(transmat,2);
lined_transmat = reshape(transmat,obj.nbasis.*n_basisuse,1);
transmath2 = reshape(permute(reshape(lined_transmat*lined_transmat',obj.nbasis,n_basisuse,obj.nbasis,n_basisuse),[1 3 2 4]),obj.nbasis.^2,n_basisuse.^2);

res = reshape(transmath2'*reshape(H2,obj.nbasis.^2,obj.nbasis.^2)*transmath2,n_basisuse,n_basisuse,n_basisuse,n_basisuse);

end