function varval = multipledef(varval,variables,key)
% key = 0 - put corrected string
% key = 1 - put original string
if ~isempty(regexpi(varval,'\$','once'))
    if ~exist('key','var')
        key = 0;
    end
    for n_var=1:size(variables,1)
        if isempty(regexpi(varval,'\$','once'))
            break
        end
        pos = regexpi(varval,variables{n_var,1});
        if isempty(pos)
            continue
        end
        if strcmpi(varval(pos-1),'''')
            if key==1
                varval = regexprep(varval,['\$''',variables{n_var,1},''''],variables{n_var,2});
            else
                varval = regexprep(varval,['\$''',variables{n_var,1},''''],num2str(variables{n_var,3},'%20.15g'));
            end
        else
            if key==1
                varval = regexprep(varval,['\$',variables{n_var,1}],variables{n_var,2});
            else
                varval = regexprep(varval,['\$',variables{n_var,1}],num2str(variables{n_var,3},'%20.15g'));
            end
        end
    end
end


