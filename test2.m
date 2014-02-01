filename = '631ppgss.txt';

bsfid = fopen(filename, 'r');

tline = fgets(bsfid);

while ischar(tline)
    if(length(tline)>3)
        if(strcmp(tline(1:4), '****')) % entry of basis set region
            break;
        end
    end
    tline = fgets(bsfid);
end

tline = fgets(bsfid);
element_ind = 0;
bs_cgf_ind = 0;
while ischar(tline)
    z = name2z(strtrim(tline(1:2)));
    if(z~=0) % element exist
        element_ind = element_ind + 1;
        bs_cgf_ind = 0;
    end
    tline = fgets(bsfid); % this line has angular momentum type and gfnum
    if(ischar(tline))
        while(~strcmp(tline(1:4), '****'))
            element{element_ind}.z = z;
            bs_cgf_ind = bs_cgf_ind + 1;
            ltypecell = regexpi(tline, '[A-z]', 'match');
            l_components_vec = 0;
            for l_ind=1:length(ltypecell)
                l_components_vec(l_ind) = ltype2l(ltypecell{l_ind});
            end
            gfnumcell = strtrim(regexpi(tline, '[1-9][\d]*[^\.]', 'match'));
            element{element_ind}.bs_cgf{bs_cgf_ind}.gfnum = str2double(gfnumcell{1});
            for gf_ind=1:element{element_ind}.bs_cgf{bs_cgf_ind}.gfnum
                tline = fgets(bsfid);
                alpha_and_d_cell = regexpi(tline, '-?([1-9]\d*\.\d*|0\.\d*[1-9]\d*|0?\.0+|0)', 'match');
                element{element_ind}.bs_cgf{bs_cgf_ind}.alpha(gf_ind) = str2double(alpha_and_d_cell{1});
                element{element_ind}.bs_cgf{bs_cgf_ind}.d(gf_ind) = str2double(alpha_and_d_cell{2});
                element{element_ind}.bs_cgf{bs_cgf_ind}.l_compo = l_components_vec;
            end
            tline = fgets(bsfid);
        end
    else
        break;
    end

    
    tline = fgets(bsfid);
end
%             tline = fgets(bsfid); % this line has element name
%             z = name2z(strtrim(tline(1:2)));
% %             disp(z)
%             if(z~=0) % element exist
%                 tline = fgets(bsfid); % this line has angular momentum type and gfnum
%                 while(~strcmpi(tline(1), ' ')) % entry to a new cgf 
%                     bs_cgf_ind = bs_cgf_ind + 1;
%                     ltypecell = regexpi(tline, '[A-z]+', 'match');
%                     
%                     gfnumcell = strtrim(regexpi(tline, '[1-9][\d]*[^\.]', 'match'));
%                     bs_cgf{bs_cgf_ind}.gfnum = str2double(gfnumcell{1});
%                     for gf_ind=1:bs_cgf{bs_cgf_ind}.gfnum
%                         tline = fgets(bsfid);
%                         alpha_and_d_cell = regexpi(tline, '-?([1-9]\d*\.\d*|0\.\d*[1-9]\d*|0?\.0+|0)', 'match');
%                         bs_cgf{bs_cgf_ind}.alpha(gf_ind) = str2double(alpha_and_d_cell{1});
%                         bs_cgf{bs_cgf_ind}.d(gf_ind) = str2double(alpha_and_d_cell{2});
%                     end
%                     
%                 end
%                 
% %                 disp(ltype)
% %                 disp(gfnum)
% %                 no_m_cgf.l = type2l(ltype);
%             end
%         end
%         
%     end
    
    
%     tline = fgets(bsfid);
% end
fclose(bsfid);