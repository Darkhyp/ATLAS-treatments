function cond = analysecondition(condstr,variables)
condstr = strtrim(condstr);
% find expression in bracats
% bra = regexp(condstr,'(');
% cat = regexp(condstr,')');

switch true
    case ~isempty(regexpi(condstr,'\|', 'once'))
        % check for logical OR 
        pos = regexpi(condstr,'\|', 'once');
        cond = analysecondition(condstr(1:pos-1),variables) || analysecondition(condstr(pos+1:end),variables);
    case ~isempty(regexpi(condstr,'\&', 'once'))
        % check for logical AND 
        pos = regexpi(condstr,'\&', 'once');
        cond = analysecondition(condstr(1:pos-1),variables) && analysecondition(condstr(pos+1:end),variables);
    case ~isempty(regexpi(condstr,'\^=', 'once')) || ~isempty(regexpi(condstr,'!=', 'once')) || ~isempty(regexpi(condstr,'<>', 'once'))
        % check for logical = or comparison two strings 
        cond_tmp = strsplit(condstr,{'^=','!=','<>'},'CollapseDelimiters',true);
        if length(cond_tmp)==2
            if isempty(str2num(cond_tmp{1})) && isempty(str2num(cond_tmp{2}))
                fprintf('\nError ''^=''');
            else
                cond = analysecondition(cond_tmp{1},variables) ~= analysecondition(cond_tmp{2},variables);
%                 cond = abs(analysecondition(cond_tmp{1},variables) - analysecondition(cond_tmp{2},variables))>1e-12;
            end
        end
    case ~isempty(regexpi(condstr,'=', 'once'))
        % check for logical = or comparison two strings 
        cond_tmp = strsplit(condstr,'=','CollapseDelimiters',true);
        if length(cond_tmp)==2
            if isempty(str2num(cond_tmp{1})) && isempty(str2num(cond_tmp{2}))
                cond = strcmpi(analysecondition(cond_tmp{1},variables),analysecondition(cond_tmp{2},variables));
            else
                cond = analysecondition(cond_tmp{1},variables) == analysecondition(cond_tmp{2},variables);
            end
        end
    case ~isempty(regexpi(condstr,'<', 'once'))
        cond_tmp = strsplit(condstr,'<','CollapseDelimiters',true);
        if length(cond_tmp)==2
            if isempty(str2num(cond_tmp{1})) && isempty(str2num(cond_tmp{2}))
                fprintf('\nError ''<''');
            else
                cond = analysecondition(cond_tmp{1},variables) < analysecondition(cond_tmp{2},variables);
            end
        end
    case ~isempty(regexpi(condstr,'>', 'once'))
        cond_tmp = strsplit(condstr,'>','CollapseDelimiters',true);
        if length(cond_tmp)==2
            if isempty(str2num(cond_tmp{1})) && isempty(str2num(cond_tmp{2}))
                fprintf('\nError ''>''');
            else
                cond = analysecondition(cond_tmp{1},variables) > analysecondition(cond_tmp{2},variables);
            end
        end
    case ~isempty(regexpi(condstr,'<=', 'once'))
        cond_tmp = strsplit(condstr,'<=','CollapseDelimiters',true);
        if length(cond_tmp)==2
            if isempty(str2num(cond_tmp{1})) && isempty(str2num(cond_tmp{2}))
                fprintf('\nError ''<=''');
            else
                cond = analysecondition(cond_tmp{1},variables) <= analysecondition(cond_tmp{2},variables);
            end
        end
    case ~isempty(regexpi(condstr,'>=', 'once'))
        cond_tmp = strsplit(condstr,'>=','CollapseDelimiters',true);
        if length(cond_tmp)==2
            if isempty(str2num(cond_tmp{1})) && isempty(str2num(cond_tmp{2}))
                fprintf('\nError ''>=''');
            else
                cond = analysecondition(cond_tmp{1},variables) >= analysecondition(cond_tmp{2},variables);
            end
        end
    otherwise
        % this variable
        if isempty(str2num(condstr))
            cond = regexprep(condstr,'"','');
        else
            cond = str2num(condstr);
            cond = fix(cond*1e12)/1e12;
        end
end
end



%{
function out = varfind(str,variables)
    out = 0; 
    if ~strcmpi(str(1),'$')
        return
    end
    str = regexprep(str,'\$|''','');
    for n_var=1:size(variables,1)
        if strcmpi(str,variables{n_var,1})
            out = n_var; 
            return
        end
    end
end
%}
