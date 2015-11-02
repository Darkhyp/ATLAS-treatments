function variables = addvar(varname,varstr,variables,externvars)
% chek if variable is external
if ~isempty(externvars)
    for n_var=1:size(externvars,1)
        if strcmpi(externvars{n_var,1},varname)
            % skip if variable is in the list of external variables
            return
        end
    end
end
% add or replace previous value
if isempty(variables)
    [var,type] = vartype(varstr);
    variables{1,4} = type;
    variables{1,3} = var;
    variables{1,2} = varstr;
    variables{1,1} = varname;
    return
else
    n = size(variables,1)+1;
    for n_var=1:size(variables,1)
        if strcmpi(variables{n_var,1},varname)
            n = n_var;
            break
        end
    end
    var_tmp = multipledef(varstr,variables);
    if ~isempty(regexpi(var_tmp,'\$','once'))
        fprintf('\nUnknown variable: %s',var_tmp)
        return
    %     exit(1)
    end

    
    [var,type] = vartype(var_tmp);
    variables{n,4} = type;
    variables{n,3} = var;
    variables{n,2} = var_tmp;
    variables{n,1} = varname;
    variables = sortvars(variables);
end
end

function [var,type] = vartype(var)
    if isempty(str2num(var))
        var  = regexprep(var,'"','');
        type = 'char';
    else
        var  = str2num(var);
        type = 'num';
    end
end
