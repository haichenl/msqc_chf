function p2hf = dens2p(obj)
% Input:

% output:
%     density2p:  (nbasis,nbasis,nbasis,nbasis) 2-particle density matrix
%                  valid only for restricted hartree fock theory

% For Restricted Hartree Fock theory, the 2-particle density is derivable 
% from the 1-particle density, as follows:

p1hf = obj.dens;
nb = obj.nbasis;
lined_p1hf = reshape(p1hf,nb.^2,1);
p2hf = lined_p1hf*lined_p1hf';
% p2hf = zeros(nb, nb, nb, nb);
% for a=1:nb
%    for b=1:nb
%       for c=1:nb
%          for d=1:nb
% %             p2hf(a,b,c,d) = 2.*p1hf(a,b).*p1hf(c,d) - p1hf(a,d).*p1hf(c,b);
%             p2hf(a,b,c,d) =  - p1hf(a,d).*p1hf(c,b);
% %             p2hf(a,b,c,d) = 2.*p1hf(a,b).*p1hf(c,d);
%          end
%       end
%    end
% end
end


