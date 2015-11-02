function out = SCperformance(infilename,width,isShow)
if ~exist('isShow','var')
    isShow=false;
end

out.IV          = NaN;
out.V_OC        = NaN;
out.J_SC        = NaN;
out.CE          = NaN;
out.FF          = NaN;
out.V_maxpower  = NaN;
out.J_maxpower  = NaN;

% fin = fopen(infilename);
% if fin<=0
%     fprintf('\nFile %s is not found',infilename)
%     return
% end
% load file to memory and split on strings
out.filename = infilename;
out.width = width;
% check file existing
if ~exist(infilename,'file')
    fprintf('\nFile %s is not found',infilename)
    return
end
% check file size
tmp = dir(infilename);
if tmp.bytes==0
    fprintf('\nFile %s is empty',infilename)
    return
end
m = memmapfile(infilename);
instrings = strsplit(char(m.Data.'),'\n','CollapseDelimiters',true).';

s1='d ';

% N_data              = length(regexp(s,['\n',s1]));
% out_tmp(N_data,1)   = 0;
out_tmp	= [];
index = 0;

% while ~feof(fin)
%     getstr = fgetl(fin);
for str_num=1:length(instrings)
    getstr = strtrim(instrings{str_num});

% only for test
% fprintf('\nRead %i string',str_num)
% if str_num==83
%     pause
% end
    if isempty(getstr)
        continue % skip empty string
    end
    if strcmpi(getstr(1),'#')
        continue % skip remark
    end
    
    if strcmpi(s1,getstr(1:length(s1)))
        index = index + 1;
        out_tmp(index,:) = str2num(getstr((length(s1)+1):end));
    end
end

% fclose(fin);
% sort
IV = sortrows(out_tmp(:,[end-1,end]),1);
IV(:,end) = -IV(:,end)/(width*1e-2*1e-3); % mA/cm^2
if size(IV,1)<2
    fprintf('\nFile %s is empty or there was problem with convergency in Silvaco.',infilename)
    return
end
while IV(1,end-1)==IV(2,end-1)
    IV = IV(2:end,:);
end
[CE,FF,J_SC,V_OC,J_maxpower,V_maxpower] = SCperformance0(IV(:,end-1),IV(:,end),isShow);
if isShow
    title(infilename)
end

out.IV = IV;
out.V_OC = V_OC;
out.J_SC = J_SC;
out.V_maxpower = V_maxpower;
out.J_maxpower = J_maxpower;
out.CE = CE;
out.FF = FF;

end
