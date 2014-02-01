classdef int < handle
    properties
        basis_set
        cgf % cell array of cgf structures 
        % cgf.gfnum
        % cgf.d(i)
        % cgf.alpha(i)
        % cgf.pos
        % cgf.l
        % cgf.m
    end
    
    methods
        
        function res = gf_s(alpha1, alpha2, distsq)
            res = (sqrt(pi./(alpha1+alpha2)).^3).*exp(-alpha1.*alpha2.*distsq./(alpha1+alpha2));
        end
        
        function res = cgf_s(cgf1, cgf2, distdiff_x, distdiff_y, distdiff_z)
            % distdiff values are for computing ke (and other derivatives)
            if(cgf1.l==cgf2.l && cgf1.m==cgf2.m)
                distvec = [distdiff_x, distdiff_y, distdiff_z];
                distsq = distvec*distvec';
                to_be_summed_mat = zeros(cgf1.gfnum,cgf2.gfnum);
                for i=1:cgf1.gfnum
                    for j=1:cgf2.gfnum
                        to_be_summed_mat(i,j) = cgf1.d(i).*cgf2.d(j).*gf_s(cgf1.alpha(i), cgf2.alpha(j), distsq);
                    end
                end
                res = sum(sum(to_be_summed_mat));
            else
                res = 0;
                return;
            end
        end
        
        function res = cgf_ke(cgf1, cgf2)
            
        end
        
        function res = cgf_int(alpha1, alpha2, distsq, int_type)
            if(int_type==0) % 0: s
                res = cgf_s(alpha1, alpha2, distsq);
            else
                res = 0;
            end
        end
    end
end