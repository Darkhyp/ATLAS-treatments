function val = getvar(varname,variables)
val = NaN;
if ~isempty(variables)
    for n_var=1:size(variables,1)
        if strcmpi(variables{n_var,1},varname)
            val = variables{n_var,3};
            return
        end
    end
end
end
