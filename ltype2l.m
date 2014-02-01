function l = ltype2l(ltype)
l_list = {'S','P','D','F','G'};
for i=1:length(l_list)
    if(strcmpi(ltype, l_list{i}))
        l = i-1;
        return;
    end
end
end