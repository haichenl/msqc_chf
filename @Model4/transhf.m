function [orb,Eorb,Ehf,Eelec] = transhf(obj,transmat,eps,maxIter,minIter)
% Solve Hartree Fock equations
% Input:
%   obj:  Holds hamiltonian information (H1,H2,S,nelec,Hnuc)
%          [see Fragment class for definitions of these properties]
%   eps:   convergence criteria
%          [defaults to 1e-8]
% Output:
%   orb:   (nbasis,nbasis)
%            orbital coefficient matrix (atom basis #, mol orb #)
%   Eorb:  (nbasis,1)  molecular orbital energies
%   Ehf:    total Hartree-Fock energy

if (nargin < 3)
    eps = 1.0e-10;
end
if (nargin < 4)
    maxIter = 5000;
end
if (nargin < 5)
    minIter = 5;
end

H1 = obj.H1();
H2 = obj.H2();
S = obj.S();

n_basisuse = size(transmat,2);

H1trans = transmat'*H1*transmat;
H2trans = obj.transform_h2(H2,transmat);

% obj.h1trans = H1trans;
% obj.h2trans = H2trans;

% we can modify H1trans and H2trans here 

H2jktrans = obj.genh2jk(H2trans);
Strans = transmat'*S*transmat;
Enuc = obj.Hnuc;
Nelec = obj.frag.nelec;

% Step 3 -- Calculate transformation matrix (eq. 3.167).
X = real(inv(sqrtm(Strans)));

% Step 4 -- Guess at density matrix -- all zeros right now.
if (isempty(obj.density_transhf))
    Pn = zeros(n_basisuse);
else
    Pn = obj.density_transhf;
end
% Pn = zeros(n_basisuse);

% how to get a better transformed density?

% Pn = transmat'*obj.densitySave*transmat; % failed

iter = 0;
finished = false;

% Begin iteration through.
while (~finished)  % Step 11 -- Test convergence
    
    P = Pn;
    
    % Step 5 -- Build 2-electron components of Fock matrix.
    G = reshape(H2jktrans*reshape(P,n_basisuse.^2,1),n_basisuse,n_basisuse);
    
    % Step 6 -- Obtain F (fock matrix).
    F = H1trans + G;
    
    % Step 7 -- Calculate the transformed F matrix.
    Ft = X'*F*X;
    
    % Step 8 -- Find e and the transformed expansion coefficient matrices.
    [Ct1,e1] = eig(Ft);
    Ct1 = real(Ct1);
    e1 = real(e1);
    e2 = diag(e1);
    [e, i1] = sort(e2);
    Ct = Ct1(:,i1);
    
    % Step 9 -- Transform Ct back to C.
    C = X*Ct;
    
    % Step 10 -- Calculate the new density matrix.
    filled = 1:(Nelec/2);
    Pn = 2* C(:,filled)*( C(:,filled)');
    iter = iter + 1;

    changeInDensity = max(max(abs(P - Pn)));
    if (iter > maxIter)
        finished = true;
    elseif (iter > minIter)
        if (changeInDensity < eps)
            finished = true;
        end
    end
end
% End of iteration of steps 5-11.

P = Pn;  % For convenience.
obj.density_transhf = P;

% Step 12: Output.

%Total energy
%3.184: E0 = 1/2 Sum(i,j) {P(j,i)[H1(i,j) + F(i,j)]}
Ee = sum(sum(P.*(H1trans+F)));
Ehf = Ee/2 + Enuc;

% electronic energy
Eelec = Ee/2;

% Orbital energies.
Eorb = e;

% Molecular orbital components.
orb = C;

% disp(iter);

if (iter+1 > maxIter)
    disp('You are living on the edge.. hartree fock didn''t converge');
end
%{
Adapted from "Modern quantum chemistry", by Attila Szab? Neil S. Ostlund
Numbered equations also adapted from here.
1. Specify a molecule
2. Calculate S(i,j), H^core (H1), and (i j|k l)(H2)
    -These first two steps are done by Gaussian
3. Diagonalize overlap matrix S and obtain X from 3.167
    3.167: X = S^(-1/2)
4. Guess the density matrix P (first guess is zeros here)
5. Calculate matrix G of 3.154 from P and H2
    G(i,j) = Sum(k, l){P(k,l)[(i j|l k)-1/2(i k|l j)]}
6. Add G to core-Hamiltonian  to get Fock matrix
    3.154: F(i,j) = H1(i,j) + G(i,j)
7. Calculate transformed Fock matrix F' = X'(t)FX
8. Diagonalize F' to obtain C' and epsilon
9. Calculate C = XC'
10. Form new density matrix P from C w/ 3.145
    3.145: P(i,j) = 2 Sum(1-Nelec/2){C(i,a) C*(j,a)}
11. Has P converged to within eps?
    No? -> Step 5 w/ new P from 10.
    Yes? -> Step 12
12. Use resultant solution, represented by C,P,F to calculate outputs
%}
