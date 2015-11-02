function str = calcmathexpression(str)

if isempty(regexp(str,'\(|\)|+|/|*|-','once'))
    return
end

pos = regexp(str,'=');
if isempty(pos)
    return
end
pos = [pos,length(str)+1];
for n_equal=length(pos)-1:-1:1
    str_tmp = strtrim(str(pos(n_equal)+1:pos(n_equal+1)-1));
    if isempty(regexp(str_tmp,'\(|\)|+|/|*|-','once'))
        continue
    end
    while isempty(str2num(str_tmp))
        str_tmp = str_tmp(1:end-1);
        if isempty(str_tmp)
            break
        end
    end
    if ~isempty(str_tmp)
        str = strrep(str,str_tmp,[num2str(str2num(str_tmp),'%20.15g'),' ']);
    end
end

end

