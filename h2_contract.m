function unique_h2_elements_vector = h2_contract(h2_4d_array, nbasis)
% contract h2 elements into unique_h2_elements_vector
% vector format: [ {cn2plusn diag terms} {cn2plusn+cn2plusn.*(cn2plusn-1)./2 offdiagterms} ]
if (nargin < 2)
    nbasis = size(h2_4d_array, 1);
end
nsq = nbasis.^2;
num_cn2plusn = nbasis + nbasis.*(nbasis-1)./2;
ordervec = zeros(1,nsq);
for i=1:nsq
    p = ceil(i./nbasis);
    q = mod(i,nbasis);
    if(q==0)
        q = nbasis; % to make mod(i,nbasis) yield 1~nbasis rather than 0~(nbasis-1)
    end
    if(p==q)
        ordervec(i) = p;
    elseif(p<q)
        ordervec(i) = nbasis + nbasis.*(nbasis-1)./2 + (2.*nbasis-p).*(p-1)./2 + (q-p);
    else
        ordervec(i) = nbasis + (2.*nbasis-q).*(q-1)./2 + (p-q);
    end
end
nsq_mat= zeros(nsq,nsq);
nsq_mat(ordervec, ordervec) = reshape(h2_4d_array,nsq,nsq);
cn2plusn_mat = nsq_mat(1:num_cn2plusn,1:num_cn2plusn);
unique_h2_elements_vector(1:num_cn2plusn) = diag(cn2plusn_mat);
unique_h2_elements_vector(num_cn2plusn+1:num_cn2plusn+num_cn2plusn.*(num_cn2plusn-1)./2) = cn2plusn_mat(tril(true(num_cn2plusn),-1));
end