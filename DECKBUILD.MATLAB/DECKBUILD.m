%{
   DECKBUILD.m -- This file is part of the treatment of
   Silvaco ATLAS files (Semiconductor Devices Simulator).

   Copyright (C) 2015-2015 Alexander V. Korovin
   <A.V.Korovin@rambler.ru>

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 3, or (at your option)
   any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <www.gnu.org/licenses/>.
%}

% ######################################################
% Created on 22.10.2015 by Alexander V. Korovin, Gif-sur-Yvette, France
% Last modifications : 
% ######################################################

% Load data from Silvaco deckbuild in-file and convert it for ATLAS
%% input variables
% infilename is the file name of Silvaco  deckbuild in-file
% outfilename is the output file name of Silvaco deckbuild in-file without
% SET, CONDITIONS and LOOPS commands
% externvars is the variables list that you need to change in your deckbuild in-file
% varcode is the code for displayed variable
% removestrings are strings that you need to remove
%% output variables
% out.stucture is the stucture name
% out.string 
% out.convertedFileName 
% out.inputFileName 
% out.variables 

function out = DECKBUILD(infilename,outfilename,externvars,removestrings)
removestrings = [];

if ~exist('infilename','var')
% only for test
%    infilename = 'd:\Korovin\work\CNRS\LGEP\IBC_HIT_full9.in';
    fprintf('\nInput filename is not specified! Use DECKBUILD(infilename,outfilename,externvars,isRun).')
    return
end
if ~exist('outfilename','var')
    outfilename = [infilename(1:end-3),'_(c).in'];
end
if isempty(outfilename)
    outfilename = [infilename(1:end-3),'_(c).in'];
end
% fin = fopen(infilename);
% if fin<=0
%     fprintf('\nFile %s is not found',infilename)
%     return
% end

% use this to change variable value in deckbuild structure file
%{
% only for test
externvars      = {'BSF_cSi_penetration_depth','0.1',0.100000000000000,'num'};
externvars{2,1} = 'cSiN_conc';
externvars{2,2} = '1e21';
externvars{2,3} = 1.00000000000000e+21;
externvars{2,4} ='num';
externvars{3,1} = 'isStructureShow';
externvars{3,2} = '1';
externvars{3,3} = 1;
externvars{3,4} ='num';
%}
if ~exist('externvars','var')
    variables = {};
    externvars = {};
else
    variables = {};
    for n_var=1:length(externvars)
        variables = addvar(externvars{n_var,1},externvars{n_var,2},variables,{});
    end
end

% use this to change variable value in deckbuild structure file
%{
% only for test
removestrings = {'/home/korovin_ale/silvaco_works/'};
%}
if ~exist('removestrings','var')
    removestrings = {};
end

%% main program
skipuntilELSE   = 0;

% load file to memory and split on strings
if ~exist(infilename,'file')
    fprintf('\nFile %s is not found',infilename)
    return
end
m = memmapfile(infilename);
instrings  = strsplit(char(m.Data.'),'\n','CollapseDelimiters',true).';
outstrings = {};
meshstrings = {};
isXMESHstart = false;
isXMESHfinish = false;
isYMESHstart = false;
isYMESHfinish = false;
xmin = 0;
ymin = 0;
for str_num=1:length(instrings)
    getstr = strtrim(instrings{str_num});
% str_num         = 0;
% while ~feof(fin)
%     getstr = strtrim(fgetl(fin));
%     str_num = str_num + 1;

%{
% only for test
fprintf('\nRead %i string',str_num)
if str_num==410
    pause
end
%}

    if isempty(getstr)
        continue % skip empty string
    end
    if strcmpi(getstr(1),'#')
        continue % skip remark
    end
    
    if ~isempty(regexpi(getstr,'if.end', 'once'))
        skipuntilELSE = skipuntilELSE(1:end-1);
        continue
    end
    % skipuntilELSE
    % 0 do all (not in IF CASE)
    % 1 (-1 out) skip until ELSE
    % 2 (-2 out) do until ELSE (skip after)
    % 3 skip all conditions of IF CASE
    switch skipuntilELSE(end)
        case 3
            % ignore all IF statements
            if ~isempty(regexpi(getstr,'if ', 'once')) && ~isempty(regexpi(getstr,'cond', 'once'))
                skipuntilELSE(end+1) = 3;
            end
            continue
        case 1
            % skip until ELSE
            if ~isempty(regexpi(getstr,'if ', 'once')) && ~isempty(regexpi(getstr,'cond', 'once'))
                skipuntilELSE(end+1) = 3;
                continue
            end
            if ~isempty(regexpi(getstr,'else', 'once'))
                skipuntilELSE(end) = -1;
            end
            continue
        case 2
            if ~isempty(regexpi(getstr,'else', 'once'))
                skipuntilELSE(end) = -2;
                continue
            end
        case -2
            if ~isempty(regexpi(getstr,'if ', 'once')) && ~isempty(regexpi(getstr,'cond', 'once'))
                skipuntilELSE(end+1) = 3;
            end
            continue
    end
    
    switch true
        case strcmpi(getstr(1:length('set')),'set')
            str_tmp = strtrim(getstr(length('set')+2:end));
            pos_tmp = regexp(str_tmp,'=', 'once');
            variables = addvar(strtrim(str_tmp(1:pos_tmp-1)),strtrim(str_tmp(pos_tmp+1:end)),variables,externvars);
        case strcmpi(getstr(1:length('if ')),'if ')
            bra = regexp(getstr,'(');
            cat = regexp(getstr,')');
            if analysecondition(multipledef(getstr(bra(1)+1:cat(end)-1),variables,1),variables)
                skipuntilELSE(end+1) = 2; % skip after ELSE
            else
                skipuntilELSE(end+1) = 1; % skip until ELSE
            end
        otherwise
            getstr = calcmathexpression(multipledef(getstr,variables,0));
            if ~isempty(removestrings)
                for n_remstr=1:length(removestrings)
                    getstr = strrep(getstr,removestrings{n_remstr},'');
                end
            end

            % form data for x.mesh
            if ~isXMESHfinish
                if ~isempty(regexpi(getstr,'x.mesh', 'once'));
                    isXMESHstart = true;
                    meshstrings{end+1} = getstr;
                    continue
                elseif isXMESHstart
                    isXMESHfinish = true;
                    % correct x.mesh data
                    xmax = findmax('x',variables);

                    xcoordinares = findmeshcoord(meshstrings);
                    xcoordinares = sortrows(xcoordinares(xcoordinares(:,1)>=xmin & xcoordinares(:,1)<=xmax,:),[1,2]);
                    [~,ind] = unique(xcoordinares(:,1));
                    xcoordinares = xcoordinares(ind,:);
                    for n=1:length(xcoordinares)
                        outstrings{end+1} = ['x.mesh location=',num2str(xcoordinares(n,1),'%20.15g'),' SPACING=',num2str(xcoordinares(n,2),'%20.15g')];
                    end
                    meshstrings = {};
                end
            end

            % form data for y.mesh
            if ~isYMESHfinish
                if ~isempty(regexpi(getstr,'y.mesh', 'once'));
                    isYMESHstart = true;
                    meshstrings{end+1} = getstr;
                    continue
                elseif isYMESHstart
                    isYMESHfinish = true;
                    % correct y.mesh data
                    ymax = findmax('y',variables);

                    ycoordinares = findmeshcoord(meshstrings);
                    ycoordinares = sortrows(ycoordinares(ycoordinares(:,1)>=ymin & ycoordinares(:,1)<=ymax,:),[1,2]);
                    [~,ind] = unique(ycoordinares(:,1));
                    ycoordinares = ycoordinares(ind,:);
                    for n=1:length(ycoordinares)
                        outstrings{end+1} = ['y.mesh location=',num2str(ycoordinares(n,1),'%20.15g'),' SPACING=',num2str(ycoordinares(n,2),'%20.15g')];
                    end
                    meshstrings = {};
                end
            end
            
            outstrings{end+1} = getstr;
%             fprintf(fout,'%s\n',getstr);
            continue
    end
end

% fclose(fin);
fout = fopen(outfilename,'w');
for n_str=1:length(outstrings)
    fprintf(fout,'%s\n',outstrings{n_str});
    outstrings{n_str} = [outstrings{n_str},'\n'];
end
fclose(fout);

fprintf('\nConversion is finished...')

out.stucture = outstrings;
out.string = getvar('OutputString',variables);
out.convertedFileName = outfilename;
out.inputFileName = infilename;
out.variables = variables;
end

%% Find max value for calculation region
function xmax = findmax(str,variables)
n = 1;
xmax = Inf;
while true
    tmp = getvar([str,int2str(n)],variables);
    n = n + 1;
    if isnan(tmp)
        break
    else
        xmax = tmp;
    end
end
end

function xcoordinares = findmeshcoord(meshstrings)
for n=1:length(meshstrings)
    tmp_str = strsplit(meshstrings{n},{' ','\t'},'CollapseDelimiters',true,'DelimiterType','Simple');
    for m=1:length(tmp_str)
        switch true
            case ~isempty(regexpi(tmp_str{m},'LOC', 'once'))
                tmp = strsplit(tmp_str{m},{'='},'CollapseDelimiters',true,'DelimiterType','Simple');
                xcoordinares(n,1) = str2num(tmp{2});
                continue
            case ~isempty(regexpi(tmp_str{m},'SPA', 'once'))
                tmp = strsplit(tmp_str{m},{'='},'CollapseDelimiters',true,'DelimiterType','Simple');
                xcoordinares(n,2) = str2num(tmp{2});
                continue
        end
    end    
end
end


