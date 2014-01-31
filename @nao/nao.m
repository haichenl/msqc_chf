classdef nao < handle
    properties
        nbasis
        densd
        s
        almblocks
        unsymmblocks
        n_nmb
        
        c_nmb
        w_nmb
    end
    methods
        function res = nao(model4)
            if(isempty(model4.dens))
                error('please solve hf first');
            end
            res.nbasis = model4.nbasis;
            res.densd = model4.dens;
            res.s = model4.S;
            res.almblocks = nao.gen_almblocks(model4);
            res.unsymmblocks = nao.gen_unsymmblocks(model4);
            res.n_nmb = nao.calc_n_nmb(model4);
            [res.c_nmb, res.w_nmb] = res.donao;
        end
        
        function [c_nmb, w_nmb] = donao(obj)
            
            densp = obj.s*obj.densd*obj.s;
            n_almblocks = size(obj.almblocks,2);
            pnaot = zeros(obj.nbasis);
            pnao = zeros(obj.nbasis); % transforms AO to pre-NAO
            n_red = zeros(obj.nbasis); % transforms AO to NAO
            w_pnao = zeros(obj.nbasis); % weighted occupation of pre-NAOs
            w_nao = zeros(obj.nbasis); % weighted occupation of NAOs
            
            % symmetrically average over dens and s
            symm_avg_densp = densp;
            symm_avg_s = obj.s;
%             for i=1:length(obj.unsymmblocks)
%                 avgdens = nao.blockavg(densp,obj.almblocks,obj.unsymmblocks{i});
%                 avgs = nao.blockavg(obj.s,obj.almblocks,obj.unsymmblocks{i});
%                 for j=obj.unsymmblocks{i}
%                     symm_avg_densp(obj.almblocks{j},obj.almblocks{j}) = avgdens;
%                     symm_avg_s(obj.almblocks{j},obj.almblocks{j}) = avgs;
%                 end
%             end
            
            % diagonalize each (Alm) diagonal block to form pnao
            for i=1:n_almblocks
                x = inv(sqrtm(symm_avg_s(obj.almblocks{i}, obj.almblocks{i})));
                [pnaot(obj.almblocks{i}, obj.almblocks{i}), w_pnao(obj.almblocks{i}, obj.almblocks{i})] = eig(x*symm_avg_densp(obj.almblocks{i},obj.almblocks{i})*x);
                pnao(obj.almblocks{i}, obj.almblocks{i}) = x*pnaot(obj.almblocks{i}, obj.almblocks{i});
            end
            
            % separate nmb with nrb while keeping their original order
            [~, sep_ordervec] = sort(diag(w_pnao), 'descend');
            nmb_indices = sort(sep_ordervec(1:obj.n_nmb));
            nrb_indices = sort(sep_ordervec(obj.n_nmb+1:obj.nbasis));
            
            % symmetrically orthogonalize nmb
            orpnao = pnao;
            w_pnao_nmb = w_pnao(nmb_indices, nmb_indices);
            s_pnao_nmb = pnao(:, nmb_indices)'*symm_avg_s*pnao(:, nmb_indices);
            orpnao(:, nmb_indices) = pnao(:, nmb_indices)*(w_pnao_nmb/sqrtm(w_pnao_nmb*s_pnao_nmb*w_pnao_nmb));
            
            % schmidt orthogonalize each nrb to all nmb's
            os = eye(obj.nbasis) - orpnao(:,nmb_indices)*orpnao(:,nmb_indices)'*symm_avg_s;
            orpnao(:, nrb_indices) = os*orpnao(:, nrb_indices);
            
            % diagonalize each (Anlm) diagonal block only in nrb space
            % not needed for 2-zeta basis sets
            w_orpnao = diag(diag(os'*w_pnao*os));
            w_orpnao_nrb = w_orpnao(nrb_indices, nrb_indices);
            
            % symmetrically orthogonalize nrb
            s_orpnao_nrb = orpnao(:, nrb_indices)'*symm_avg_s*orpnao(:, nrb_indices);
            orpnao(:, nrb_indices) = orpnao(:, nrb_indices)*(w_orpnao_nrb/sqrtm(w_orpnao_nrb*s_orpnao_nrb*w_orpnao_nrb));
            
            % re-diagonalization
            new_symm_avg_densp = orpnao'*symm_avg_densp*orpnao;
            for i=1:n_almblocks
                [n_red(obj.almblocks{i}, obj.almblocks{i}), w_nao(obj.almblocks{i}, obj.almblocks{i})] = eig(new_symm_avg_densp(obj.almblocks{i},obj.almblocks{i}));
            end
            
            ao_to_nao_tr = orpnao*n_red;
            c_nmb = ao_to_nao_tr(:, nmb_indices);
            w_nmb = diag(ao_to_nao_tr(:, nmb_indices)'*symm_avg_densp*ao_to_nao_tr(:, nmb_indices));
            
        end
        
    end
    
    methods (Static)
        
        function res = gen_almblocks(model4)
            % generate a almblocks cell array
            if( strcmpi(model4.frag.config.basisSet(1:3), 'sto') ) % minimal basis
                curr_ll_bas_start = 1;
                for atom_z=model4.Z
                    if(atom_z>1) % heavy atom
                        res{curr_ll_bas_start} = curr_ll_bas_start; % 1s
                        res{curr_ll_bas_start+1} = curr_ll_bas_start+1; % 2s
                        res{curr_ll_bas_start+2} = curr_ll_bas_start+2; % 2p1
                        res{curr_ll_bas_start+3} = curr_ll_bas_start+3; % 2p2
                        res{curr_ll_bas_start+4} = curr_ll_bas_start+4; % 2p3
                        curr_ll_bas_start = curr_ll_bas_start + 5; % ll jump to next atom
                    else % hydrogen
                        res{curr_ll_bas_start} = curr_ll_bas_start; % 1s
                        curr_ll_bas_start = curr_ll_bas_start + 1; % ll jump to next atom
                    end
                end
            elseif( strcmpi(model4.frag.config.basisSet, '6-31g') || strcmpi(model4.frag.config.basisSet, '3-21g') ) % split valence
                curr_ll_bas_start = 1;
                curr_hl_bas_start = 1;
                for atom_z=model4.Z
                    if(atom_z>1) % heavy atom
                        res{curr_ll_bas_start} = [curr_hl_bas_start curr_hl_bas_start+1 curr_hl_bas_start+5]; % 1s & 2s
                        res{curr_ll_bas_start+1} = [curr_hl_bas_start+2 curr_hl_bas_start+6]; % 2p1
                        res{curr_ll_bas_start+2} = [curr_hl_bas_start+3 curr_hl_bas_start+7]; % 2p2
                        res{curr_ll_bas_start+3} = [curr_hl_bas_start+4 curr_hl_bas_start+8]; % 2p3
                        curr_ll_bas_start = curr_ll_bas_start + 4; % ll jump to next atom
                        curr_hl_bas_start = curr_hl_bas_start + 9; % hl jump to next atom
                    else % hydrogen
                        res{curr_ll_bas_start} = [curr_hl_bas_start curr_hl_bas_start+1]; % 1s
                        curr_ll_bas_start = curr_ll_bas_start + 1; % ll jump to next atom
                        curr_hl_bas_start = curr_hl_bas_start + 2; % hl jump to next atom
                    end
                end
            end
        end
        
        function res = gen_unsymmblocks(model4)
            % generate an unsymmblocks cell array
            res = {};
            if( strcmpi(model4.frag.config.basisSet, '6-31g') || strcmpi(model4.frag.config.basisSet, '3-21g') ) % split valence
                curr_ll_bas_start = 1;
                for atom_z=model4.Z
                    if(atom_z>1) % heavy atom
                        res{end+1} = [curr_ll_bas_start+1 curr_ll_bas_start+2 curr_ll_bas_start+3];
                        curr_ll_bas_start = curr_ll_bas_start + 4; % ll jump to next atom
                    else % hydrogen
                        curr_ll_bas_start = curr_ll_bas_start + 1; % ll jump to next atom
                    end
                end
            else
                error('basis set not supported');
            end
        end
        
        function res = calc_n_nmb(model4)
            res = 0;
            if( strcmpi(model4.frag.config.basisSet, '6-31g') || strcmpi(model4.frag.config.basisSet, '3-21g') ) % split valence
                for atom_z=model4.Z
                    if(atom_z>1) % 2nd row heavy atom
                        res = res + 5;
                    else % hydrogen
                        res = res + 1;
                    end
                end
            else
                error('basis set not supported');
            end
        end
        
        function res = blockavg(mat,almblocks,avgindices)
            % avgindices must be an array composed of consecutive numbers
            n = length(avgindices);
            res = [];
            for i=avgindices
                if(isempty(res))
                    res = mat(almblocks{i}, almblocks{i});
                else
                    res = res + mat(almblocks{i}, almblocks{i});
                end
            end
            res = res./n;
        end
        
    end
end



