classdef FitSingle < handle
    properties
        ll  % low level Model3
        hl  % high level Model3
    end
    methods
        
        function res = FitSingle(fileName, read_index)
            if (nargin < 2)
                read_index = 1;
            end
            load(fileName,'LL','HL');
            res.ll = Model4(LL{read_index,1});
            res.hl = Model4(HL{read_index,1});

            % MATLAB derping.
            warning('off', 'MATLAB:mir_warning_maybe_uninitialized_temporary');
        end
        
        function solvehf(obj)
%             tic
            obj.ll.solvehf();
            obj.hl.solvehf();
%             toc
        end
        
%         function solvehft(obj)
%             tic
%             obj.ll.solvehft();
%             obj.hl.solvehft();
%             toc
%         end
        
    end
end