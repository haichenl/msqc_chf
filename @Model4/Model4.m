classdef Model4 < handle
    %{
   Interpolate between narrow and diffuse STO-3G matrix elements.
    %}
    properties
        % Input to the class
        frag;  % Fragment with regular STO-3G basis set
        
        % Most recent predictions of the model
        Ehf     % Hartree Fock energy
        newEhf
        Eorb    % (nbasis,1)      molecular orbital energies
        orb     % (nbasis,nbasis) molecular orbital coefficients
        Eelec   % electronic energy
        
        singlezeta
        nao
        
        dens    % density matrix 
        denst   % transformed density matrix 
        
        % Useful properties initialized in the constructor
        natom   % number of atoms in fragment
        nelec   % number of electrons in the fragment
        Z       % (1,natom) atomic numbers of the molecules
        aType   % (1,natom) atom type (initialized to Z)
        rcart   % (3,natom) cartesian coordinates of the atoms
        nbasis  % number of atomic (and molecular) basis functions
        
        %end
        %properties (Transient)
        densitySave   % most recent density matrix
        density_transhf
        % used to start HF iterations
    end
    methods (Static)
        h2 = H2slater(F0, G1, F2)
    end
    methods
        
        function res = Model4(frag_)
            if (nargin ~= 0)
                res.frag = frag_;
                res.natom = frag_.natom;
                res.nelec = frag_.nelec;
                res.Z     = frag_.Z;
                res.aType = res.Z;
                res.rcart = frag_.rcart;
                res.nbasis = frag_.nbasis;
                res.singlezeta = res.initialize_singlezeta();
            end
        end
        
        function res = initialize_singlezeta(obj)
            if( strcmpi(obj.frag.config.basisSet(1:3), 'sto') )
                res = 1;
            elseif ( strcmpi(obj.frag.config.basisSet, '6-31g') || strcmpi(obj.frag.config.basisSet, '3-21g') )
                res = 0;
            else
                res = 0;
            end
        end
        
        function res = H1ini(obj)
            res = obj.frag.KE + sum(obj.frag.H1en,3);
        end
        
        function res = H1(obj)
            res = obj.frag.KE + sum(obj.frag.H1en,3);
        end
        
        function res = KE(obj)  % modifying matrix
            res   = obj.frag.KE;
        end
        
        function res = Hnuc(obj)
            res = obj.frag.Hnuc;
        end
        
        function res = H2(obj)
            res = obj.frag.H2;
        end
        
        function res = S(obj)
            res = obj.frag.S;
        end
        
        function res = pars(obj)
            res = obj.frag.config.zmat.pars;
        end
        
        function res = transform_h2(~, H2, transmat)
            % transmath2 = zeros(n_basisori^2,n_basisuse^2);
            % for i=1:n_basisuse
            %     for j=1:n_basisuse
            %         transmath2(:,(i-1).*n_basisuse+j) = reshape(transmat(:,j)*transmat(:,i)',n_basisori.^2,1);
            %     end
            % end
            n_basisori = size(transmat,1);
            n_basisuse = size(transmat,2);
            lined_transmat = reshape(transmat,n_basisori.*n_basisuse,1);
            transmath2 = reshape(permute(reshape(lined_transmat*lined_transmat',n_basisori,n_basisuse,n_basisori,n_basisuse),[1 3 2 4]),n_basisori.^2,n_basisuse.^2);
            res = reshape(transmath2'*reshape(H2,n_basisori.^2,n_basisori.^2)*transmath2,n_basisuse,n_basisuse,n_basisuse,n_basisuse);
        end
        
        function res = genh2jk(~, H2)
            nbasisuse = size(H2,1);
            res = reshape(H2,nbasisuse.^2,nbasisuse.^2) - reshape(permute(H2,[1 4 2 3]),nbasisuse.^2,nbasisuse.^2)./2;
        end
        
    end % methods
    
end %
