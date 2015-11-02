function variables2 =  sortvars(variables)
if ~isempty(variables)
    for n_var=size(variables,1):-1:1
%         var_name_tmp = variables{n_var,1};
        var_size(n_var,1) = length(variables{n_var,1});
%         var_name{n_var,1} = variables{n_var,1};
    end
%     [~,index] = sortrows(table(var_name,var_size),{'var_size','var_name'},{'descend','descend'});
    [~,index] = sortrows(var_size,-1);

    for n_var=size(variables,1):-1:1
        variables2{n_var,4} = variables{index(n_var),4};
        variables2{n_var,3} = variables{index(n_var),3};
        variables2{n_var,2} = variables{index(n_var),2};
        variables2{n_var,1} = variables{index(n_var),1};
    end
    
end
end
