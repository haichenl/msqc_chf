function h2_4d_array = h2_distribute(unique_h2_elements_vector, nbasis)
% distribute h2 elements from unique_h2_elements_vector
% vector format: [ {cn2plusn diag terms} {cn2plusn+cn2plusn.*(cn2plusn-1)./2 offdiagterms} ]

% generate a cn2+n by cn2+n matrix
% diagonal terms
nsq = nbasis.^2;
nsq_mat = zeros(nsq,nsq);
num_cn2plusn = nbasis + nbasis.*(nbasis-1)./2;
diag_h2vec = unique_h2_elements_vector(1:num_cn2plusn);
offdiag_h2vec = unique_h2_elements_vector(num_cn2plusn+1:num_cn2plusn+num_cn2plusn.*(num_cn2plusn-1)./2);
nsq_mat_diag = diag(diag_h2vec);
% off-diagonal terms
nsq_mat_offdiag = zeros(num_cn2plusn, num_cn2plusn);
% generate triangular matrix then the final cn2+n by cn2+n matrix
nsq_mat_offdiag(tril(true(num_cn2plusn),-1)) = offdiag_h2vec ;
nsq_mat(1:num_cn2plusn,1:num_cn2plusn) = nsq_mat_offdiag + nsq_mat_offdiag.' + nsq_mat_diag;

% duplicate part of the cn2+n by cn2+n matrix to yield a n^2 by n^2 matrix
nsq_mat(num_cn2plusn+1:nsq, num_cn2plusn+1:nsq) = nsq_mat(nbasis+1:num_cn2plusn, nbasis+1:num_cn2plusn);
nsq_mat(1:num_cn2plusn, num_cn2plusn+1:nsq) = nsq_mat(1:num_cn2plusn, nbasis+1:num_cn2plusn);
nsq_mat(num_cn2plusn+1:nsq, 1:num_cn2plusn) = nsq_mat(nbasis+1:num_cn2plusn, 1:num_cn2plusn);

% reshape the n^2 by n^2 matrix into a n by n by n by n quad-matrix

nums = reshape(1:nsq,nbasis,nbasis);
pickvec = nums(tril(true(nbasis)));

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
h2_4d_array = reshape(nsq_mat(ordervec,ordervec),nbasis,nbasis,nbasis,nbasis); % final result
end