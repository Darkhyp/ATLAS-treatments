%{
   AtlasSimulations.m -- This file is part of the treatment of
   Silvaco ATLAS files (Semiconductor Devices Simulator) to analyse the Solar Cell performance.

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

% Run deckbuild in-file and calculate the solar cell performance
%% input variables
% InputData is the structure with parameters to run ATLAS simulations


function AtlasSimulations(InputData)

% addpath('d:\Korovin\work\CNRS\LGEP\DECKBUILD.MATLAB\')
% addpath('.\DECKBUILD.MATLAB\')
addpath('./DECKBUILD.MATLAB/')

if ~isfield(InputData,'UNIXserverName')
    InputData.UNIXserverName = 'silvacox2';
    fprintf('\n UNIX server ("InputData.UNIXserverName") is not defined. Use default (''silvacox2'')...')
end

if ~isfield(InputData,'ATLASparameter')
%     InputData.ATLASparameter = '-V 5.19.20.R -P 8';
%     InputData.ATLASparameter = '-V 5.20.2.R -P 8';
    InputData.ATLASparameter = '-V 5.21.10.C -P 8';
    fprintf('\n ATLAS parameters ("InputData.ATLASparameter") is not defined. Use default (%s)...',InputData.ATLASparameter)
end

InputData.OSstringEnd = '';
if ispc()
    InputData.UNIXserverName = 'dos';
    InputData.OSstring = ['C:\sedatools\exe\atlas.exe ',InputData.ATLASparameter,' '];
%     InputData.OSstring = 'C:\sedatools\exe\deckbuild.exe -run -n -noplot ';    
%     InputData.OSstring = 'C:\sedatools\etc\GuiAppStarter.exe -lib-dir-name deckbuild -exe-name deckbld -run ';
%     InputData.OSstring = 'C:\sedatools\lib\deckbuild\3.20.2.R\x86-nt\deckbuild.exe -run ';
%     InputData.OSstring = 'deckbuild -run -as -ascii ';
    if isfield(InputData,'isLog')
        InputData.OSstringEnd = [' -noplot -outfile log_tmp.out',InputData.OSstringEnd];
    else
        fprintf('\n Switch "InputData.isLog" (to create deckbuild log-file) is not defined. Use default (false)...')
    end
else
    switch InputData.UNIXserverName
        case 'silvacox'
%             InputData.OSstring = ['/usr/local/silvaco2014/lib/atlas/5.19.20.R/x86_64-linux/atlas.exe ',InputData.ATLASparameter,' '];
            InputData.OSstring = ['atlas ',InputData.ATLASparameter,' '];
        case 'silvacox2'
%             InputData.OSstring = ['silva2 "cd silvaco_works; atlas ',InputData.ATLASparameter,' '];
            InputData.OSstring = ['ssh korovin_ale@silvacox2 "cd silvaco_works; atlas ',InputData.ATLASparameter,' '];
            InputData.OSstringEnd = '"';
    end
end

nm = 1e-3; % microns

%%
if isfield(InputData,'isShow')
    isShowIV = InputData.isShow;
else
    isShowIV = false;
    fprintf('\n Switch "InputData.isShow" (Show IV curves) is not defined. Use default (false)...')
end
if isfield(InputData,'isReplace')
    isReplace = InputData.isReplace;
else
    isReplace = false;
    fprintf('\n Switch "InputData.isReplace" (Replace previously calculated data) is not defined. Use default (false)...')
end
if isfield(InputData,'isReplaceMat')
    isReplaceMat = InputData.isReplaceMat;
else
    isReplaceMat = false;
    fprintf('\n Switch "InputData.isReplaceMat" (Replace previously calculated mat-file) is not defined. Use default (false)...')
end
if isfield(InputData,'isSaveData')
    isSaveData = InputData.isSaveData;
else
    isSaveData = false;
    fprintf('\n Switch "InputData.isSaveData" (Save TonyPlot Data) is not defined. Use default (false)...')
end
if isfield(InputData,'isTony')
    isTony = InputData.isTony;
else
    isTony = false;
    fprintf('\n Switch "InputData.isTony" (extract TonyPlot Data) is not defined. Use default (false)...')
end
if isfield(InputData,'isReadOnly')
    isReadOnly = InputData.isReadOnly;
else
    isReadOnly = false;
    fprintf('\n Switch "InputData.isReadOnly" (no calculations, only read previously calculated data) is not defined. Use default (false)...')
end

if isfield(InputData,'externvars0')
    externvars0 = InputData.externvars0;
else
    externvars0 = [];
    fprintf('\n No external variables ("InputData.externvars0")...')
end
if isfield(InputData,'vars')
    vars = InputData.vars;
else
    fprintf('\n varying variables are not defined ("InputData.vars"). Exit...')
    return
end
if isfield(InputData,'outputDir')
    outputDir = InputData.outputDir;
else
    outputDir = '';
    fprintf('\n Output Directory is not defined ("InputData.outputDir")...')
end

if isfield(InputData,'infilename')
    infilename = InputData.infilename;
else
    fprintf('\n Input deckbuild file is not defined ("InputData.infilename")...')
    return
end
if isfield(InputData,'outfilename')
    outfilename = InputData.outfilename;
else
    fprintf('\n Temporary deckbuild file is not defined ("InputData.outfilename")...')
    fprintf('\n \tuse ["InputData.outfilename"+''tmp'']')
    strtmp = strsplit(InputData.infilename,'.');
    outfilename = [strtmp{1},'tmp','.',strtmp{2}];
end

if isfield(InputData,'DPM')
    fprintf('\nThe defect pool model is used in calculations.')
    isDPM = true;
    % default values for DPM
    if ~isfield(InputData.DPM,'regnum')
        InputData.DPM.regnum = 1;
    end
    if ~isfield(InputData.DPM,'Nmax_iter') % max number of iterations
        InputData.DPM.Nmax_iter = 12;
    end
    if ~isfield(InputData.DPM,'alpha') % convergence constant
        InputData.DPM.alpha = 0.26;
    end
    if ~isfield(InputData.DPM,'T')
        InputData.DPM.T = 300;
    end
    if ~isfield(InputData.DPM,'T0')
        InputData.DPM.T0 = 480;
    end
    if ~isfield(InputData.DPM,'Eg')
        InputData.DPM.Eg = 1.7;
    end
    if ~isfield(InputData.DPM,'sigma')
        InputData.DPM.sigma = 0.17;
    end
    if ~isfield(InputData.DPM,'U')
        InputData.DPM.U = 0.2; % the correlation energy, which is needed to place two electrons on the same defect
    end
    if ~isfield(InputData.DPM,'H')
        InputData.DPM.H = 5e21; % 1/cm3
    end
    if ~isfield(InputData.DPM,'NSiSi')
        InputData.DPM.NSiSi = 2e23; % 1/cm3
    end
    if ~isfield(InputData.DPM,'Ep')
        InputData.DPM.Ep = 1.27;
    end
    if ~isfield(InputData.DPM,'dx')
        InputData.DPM.dx = [];
    end
    if ~isfield(InputData.DPM,'dy')
        InputData.DPM.dy = [];
    end
else
    isDPM = false;
end

FinalDataFileName = InputData.FinalDataFileName;
if isfield(InputData,'originname')
    originname = InputData.originname;
else
    originname = FinalDataFileName;
    fprintf('\n "originname" is not defined. Use ther same as "FinalDataFileName".')
end


%% main procedure
out = {};

if ~exist(outputDir,'dir')
    mkdir(outputDir)
end

N0          = size(externvars0,1);
N1          = size(vars,1);
N_new       = [];
for i=N1:-1:1
    N(i) = length(vars{i,2});
    if N(i)>1
        N_new(end+1,1) = i;
        N_new(end,2) = N(i);
    end
end
if ~isempty(N_new)
    N_new           = sortrows(N_new,-2);
    N1_new          = size(N_new,1);
else
    N1_new          = 1;
end
n               = ones(size(N));
n_new           = ones(1,N1_new);
N_all           = prod(N);
V_OC(N_all)     = 0;
J_SC(N_all) 	= 0;
CE(N_all)       = 0;
FF(N_all)       = 0;
fnames{N_all}   = [];
X               = {};
X{N1_new}       = [];
ind             = 1;

nall            = 1;
if isReplaceMat || ~exist([FinalDataFileName,'.mat'],'file')
    while true
        Outfilename = '';
        externvars = externvars0;
        for i=1:N1
            externvars{N0+i,1} = vars{i,1};
            externvars{N0+i,2}   = num2str(vars{i,2}(n(i)),'%20.15g');
            if strcmpi(Outfilename,'')
                Outfilename = [vars{i,1},'_',num2str(vars{i,2}(n(i)),'%g')];
            else
                Outfilename = [Outfilename,', ',vars{i,1},'_',num2str(vars{i,2}(n(i)),'%g')];
            end
        end

        tic
        fprintf('\n\n%i of %i',ind,N_all)
        fprintf('\nRun deckbuild conversion...')
        outD = DECKBUILD(infilename,outfilename,externvars);
        toc

        filenameIV0 = ['IV_',outD.string,'.log'];
        filenameIV = [outputDir,'IV_',Outfilename,'.log'];
        filenameDark0 = [outD.string,'_dark.str'];
        filenameDark = [outputDir,Outfilename,'_dark.str'];
        filenameLight0V = [outD.string,'_light0V.str'];
        filenameLight0A = [outD.string,'_light0A.str'];
        outD.string = Outfilename;
        if length(filenameIV)>131
            filenameIV = filenameIV(1:131);
        end
        if isReplace || ~exist(filenameIV,'file')
            %% constract external file for photogeneration in the case of hybrid IBC (very robust)
            GR_FileName = getvar('GR_FileName',outD.variables);
            if ~isnan(GR_FileName)
                if ~exist(GR_FileName,'file')
                    Sweeptmp = getvar('SweepType',outD.variables);
                    switch Sweeptmp
                        case 3
                            fprintf('\n External photogeneration file preparation for BSF - het & emitter - homo')
                            isFull = strfind(getvar('GR_FileName',outD.variables),'half')==0;
                            if ~exist('G_het','var')
                                if isFull
                                    GR_FileName0 = strsplit(GR_FileName,{'G_AR_BSF_HET_EM_HOMO','ratio'});
                                else
                                    %half structure
                                    GR_FileName0 = strsplit(GR_FileName,{'Ghalf_AR_BSF_HET_EM_HOMO','ratio'});
                                end
                                GR_FileName0 = GR_FileName0{2};
                                % HET
                                G_het = load(['IBC_SC_air_AR',GR_FileName0,'.table']);
                                % HOMO
                                ind_tmp = regexpi(GR_FileName0,'_');
                                ind_Si = regexpi(GR_FileName0,'_Si');
                                ind_aSi = regexpi(GR_FileName0,'_aSi'); ind_aSi = ind_aSi(ind_aSi>ind_Si);
                                GR_FileName0 = GR_FileName0([1:(ind_aSi-1),ind_tmp(find(ind_aSi==ind_tmp)+1):end]);
                                G_homo = load(['IBC_SC_air_AR',GR_FileName0,'.table']);
                            end
                            % Hybrid
                            N = getvar('NX',outD.variables)+1;
                            ratio = getvar('BSF_area',outD.variables);
                            Ghybrid = zeros(max(size(G_het,1),size(G_homo,1)),N);
                            if isFull
                                % Full
                                N_r = round(ratio*N/2);
                                Ghybrid(1:size(G_het,1), [1:N_r,(N-N_r+1):N])   = G_het(:, ones(2*N_r,1));
                                Ghybrid(1:size(G_homo,1),(N_r+1):(N-N_r))       = G_homo(:,ones(N-2*N_r,1));
                            else
                                % half
                                N_r = round(ratio*N);
                                Ghybrid(1:size(G_het,1), 1:N_r)   = G_het(:, ones(N_r,1));
                                Ghybrid(1:size(G_homo,1),N_r+1:N) = G_homo(:,ones(N-N_r,1));
                            end
                            save(GR_FileName,'Ghybrid','-ascii','-double')
                        case 4
                            fprintf('\n External photogeneration file preparation for BSF - homo & emitter - het')
                            isFull = strfind(getvar('GR_FileName',outD.variables),'half')==0;
                            if ~exist('G_het','var')
                                if isFull
                                    GR_FileName0 = strsplit(GR_FileName,{'G_AR_BSF_HOMO_EM_HET','ratio'});
                                else
                                    %half structure
                                    GR_FileName0 = strsplit(GR_FileName,{'Ghalf_AR_BSF_HOMO_EM_HET','ratio'});
                                end
                                GR_FileName0 = GR_FileName0{2};
                                % HET
                                G_het = load(['IBC_SC_air_AR',GR_FileName0,'.table']);
                                % HOMO
                                ind_tmp = regexpi(GR_FileName0,'_');
                                ind_Si = regexpi(GR_FileName0,'_Si');
                                ind_aSi = regexpi(GR_FileName0,'_aSi'); ind_aSi = ind_aSi(ind_aSi>ind_Si);
                                GR_FileName0 = GR_FileName0([1:(ind_aSi-1),ind_tmp(find(ind_aSi==ind_tmp)+1):end]);
                                G_homo = load(['IBC_SC_air_AR',GR_FileName0,'.table']);
                            end
                            % Hybrid
                            N = getvar('NX',outD.variables)+1;
                            ratio = getvar('BSF_area',outD.variables);
                            Ghybrid = zeros(max(size(G_het,1),size(G_homo,1)),N);
                            if isFull
                                % Full
                                N_r = round(ratio*N/2);
                                Ghybrid(1:size(G_homo,1),[1:N_r,(N-N_r+1):N])   = G_homo(:,ones(2*N_r,1));
                                Ghybrid(1:size(G_het,1), (N_r+1):(N-N_r))       = G_het(:, ones(N-2*N_r,1));
                            else
                                % half
                                N_r = round(ratio*N);
                                Ghybrid(1:size(G_homo,1),1:N_r)   = G_homo(:,ones(N_r,1));
                                Ghybrid(1:size(G_het,1), N_r+1:N) = G_het(:, ones(N-N_r,1));
                            end
                            save(GR_FileName,'Ghybrid','-ascii','-double')
                    end
                    % copy .table file to another server
                    if strcmpi(InputData.UNIXserverName,'silvacox2')
                        system(['scp "',GR_FileName,'" korovin_ale@silvacox2:silvaco_works/']);
                    end
                    % figure, semilogy((0:(size(Ghybrid,1)-1))*0.5,Ghybrid,'.-'), hold on, semilogy((0:(size(G_homo,1)-1))*0.5,G_homo,'.-')
                    % figure, surf(log10(Ghybrid),'linestyle','none')
                end
            end
            
            fprintf('\nRun deckbuild...')
            if strcmpi(InputData.UNIXserverName,'silvacox2')
                system(['scp "',outfilename,'" korovin_ale@silvacox2:silvaco_works/'])
                InputData.OSstringEnd = ['; scp "',filenameIV0,'" korovin_ale@silvacox:silvaco_works/; scp "',filenameDark0,'" korovin_ale@silvacox:silvaco_works/; scp "',filenameLight0A,'" korovin_ale@silvacox:silvaco_works/"; scp "',filenameLight0V,'" korovin_ale@silvacox:silvaco_works/"'];
            end
%% the defect pool model correction
            if isDPM
                F0 = ATLASExtDefs2([],InputData.DPM.T0,InputData.DPM.T0,InputData.DPM.regnum,InputData,outfilename);
                if ~isempty(F0)
                    [F,outstrings] = ATLASExtDefs2(F0,InputData.DPM.T0,InputData.DPM.T,InputData.DPM.regnum,InputData,outfilename);

                    % save in-file
                    fout = fopen(outfilename,'w');
                    for n_str=1:length(outstrings)
                        fprintf(fout,'%s\n',outstrings{n_str});
                    end
                    fclose(fout);
                    
                    close(100)
                end
            end
%%            
            tic
            system([InputData.OSstring,outfilename,InputData.OSstringEnd]);
            toc

            if exist(filenameIV0,'file')
                movefile(filenameIV0,filenameIV)
            end
            if exist(filenameDark0,'file')
                movefile(filenameDark0,[outputDir,Outfilename,'_dark.str'])
            end
            if exist(filenameLight0V,'file')
                movefile(filenameLight0V,[outputDir,Outfilename,'_light0V.str'])
            end
            if exist(filenameLight0A,'file')
                movefile(filenameLight0A,[outputDir,Outfilename,'_light0A.str'])
            end
        end

        % SC performance extract
        out{end+1} = SCperformance(filenameIV,getvar('cellwidth_ill',outD.variables)*1e-6,isShowIV);
        out{end}.outstring = outD.string;
        fprintf('\n%s is finished.\n',outD.string)

        %% Show data (113 or 114 - Band bending)
        if isTony
            outT = TONYPLOTREAD(filenameDark,[]);
            if ~isempty(outT)
                out{end}.outT = outT;
                cond = find(outT.nodes_xyz(:,2)==0);
                y  = outT.nodes_xyz(cond,3);
                [~,indT] = unique(outT.calcvalues(:,1));
                unique_n = outT.calcvalues(indT,:);

                Ev = unique_n(cond,find(outT.valcodes==113)+1);
                Ec = unique_n(cond,find(outT.valcodes==114)+1);
                Ey = unique_n(cond,find(outT.valcodes==121)+1);
                Efull = unique_n(cond,find(outT.valcodes==103)+1);
                Charge = unique_n(cond,find(outT.valcodes==116)+1);

                if isSaveData
                    out_band_tmp = [y/nm,Ev,Ec,Ey,Efull,Charge];
                    if exist('out_band','var')
                        N_tmp = max(size(out_band_tmp,1),size(out_band,1));
                        if N_tmp>size(out_band,1), out_band(N_tmp,end) = 0; end
                        if N_tmp>size(out_band_tmp,1), out_band_tmp(N_tmp,end) = 0; end
                        out_band = [out_band,out_band_tmp];
                        out_band_names = [out_band_names,'\ty, nm\tE_v, eV (Ns=',num2str(getvar('Ndefsurf',outD.variables),'%g'),'cm^2)\tE_c, eV\tEx, V/cm\tEfull, V/cm\tCarge, C/cm'];
                    else
                        out_band = out_band_tmp;
                        out_band_names = ['y, nm\tE_v, eV (Ns=',num2str(getvar('Ndefsurf',outD.variables),'%g'),'cm^2)\tE_c, eV\tEx, V/cm\tEfull, V/cm\tCarge, C/cm'];
                    end
                end

                figure(1), 
                    subplot(121)
                    cond = y<(min(y)+max(y))/2;
                    semilogx((y(cond)+1e-5)/nm,[Ec(cond),Ev(cond)]), hold on
                    xlabel('from left side, nm')
                    h1 = gca;
                    subplot(122)
                    cond = y>(min(y)+max(y))/2;
                    semilogx((max(y)-y(cond)+1e-5)/nm,[Ec(cond),Ev(cond)]), hold on
                    xlabel('from right side, nm')
                    h2 = gca;
                    set(h2,'XDir','reverse')
                    set(h2,'YAxisLocation','right')
                    ylim(h1,[min([ylim(h1),ylim(h2)]),max([ylim(h1),ylim(h2)])])
                    ylim(h2,[min([ylim(h1),ylim(h2)]),max([ylim(h1),ylim(h2)])])
                    title([num2str(getvar('Ndefsurf',outD.variables),'%g'),' {cm^{-2}}'])
                figure(2), 
                    subplot(121)
                    cond = y<(min(y)+max(y))/2;
                    semilogx((y(cond)+1e-5)/nm,Efull(cond)), hold on
                    xlabel('from left side, nm')
                    ylabel('Electric field, V/cm')
                    h1 = gca;
                    subplot(122)
                    cond = y>(min(y)+max(y))/2;
                    semilogx((max(y)-y(cond)+1e-5)/nm,Efull(cond)), hold on
                    xlabel('from right side, nm')
                    h2 = gca;
                    set(h2,'XDir','reverse')
                    set(h2,'YAxisLocation','right')
                    ylim(h1,[min([ylim(h1),ylim(h2)]),max([ylim(h1),ylim(h2)])])
                    ylim(h2,[min([ylim(h1),ylim(h2)]),max([ylim(h1),ylim(h2)])])
                    title([num2str(getvar('Ndefsurf',outD.variables),'%g'),' {cm^{-2}}'])
                figure(3), 
                    subplot(121)
                    cond = y<(min(y)+max(y))/2;
                    semilogx((y(cond)+1e-5)/nm,Charge(cond)), hold on
                    xlabel('from left side, nm')
                    ylabel('Charge, C/cm')
                    h1 = gca;
                    subplot(122)
                    cond = y>(min(y)+max(y))/2;
                    semilogx((max(y)-y(cond)+1e-5)/nm,Charge(cond)), hold on
                    xlabel('from right side, nm')
                    h2 = gca;
                    set(h2,'XDir','reverse')
                    set(h2,'YAxisLocation','right')
                    ylim(h1,[min([ylim(h1),ylim(h2)]),max([ylim(h1),ylim(h2)])])
                    ylim(h2,[min([ylim(h1),ylim(h2)]),max([ylim(h1),ylim(h2)])])
                    title([num2str(getvar('Ndefsurf',outD.variables),'%g'),' {cm^{-2}}'])
            end
        end

        V_OC(ind) = out{ind}.V_OC;
        J_SC(ind) = out{ind}.J_SC;
        CE(ind) = out{ind}.CE;
        FF(ind) = out{ind}.FF;
        fnames{ind} = out{ind}.outstring;

        if N1_new==1 && N_new(1,2)==1
            X = vars{1,2};
            break;
        end
        for i=1:N1_new
            tnp = X{i};
            X{i} = [tnp,vars{N_new(i,1),2}(n_new(i))];
        end

        [n_new,sign] = set1(1,n_new,N_new(:,2));
        if sign, break; end
        n(N_new(:,1)) = n_new;
        ind = ind+1;
    end
    
    if isTony&&isSaveData
        fout = fopen([FinalDataFileName,' band.dat'],'w');
        fprintf(fout,[out_band_names,'\n']);
        fclose(fout);
        save([FinalDataFileName,' band.dat'],'out_band','-ascii','-double','-tabs','-append')
    end
    
    save([FinalDataFileName,'.mat'],'outD','out','vars','X','V_OC','J_SC','CE','FF','fnames','externvars0')
else
    load([FinalDataFileName,'.mat'])
    if size(externvars0,2)==1
        externvars0 = reshape(externvars0,2,[]).';
        save([FinalDataFileName,'.mat'])
    end
    if size(vars,2)==1
        vars = reshape(vars,2,[]).';
        save([FinalDataFileName,'.mat'])
    end
    %{
    while true
        for i=1:N1_new
            tnp = X{i};
            X{i} = [tnp,vars{N_new(i,1),2}(n_new(i))];
        end
        V_OC(ind) = out{ind}.V_OC;
        J_SC(ind) = out{ind}.J_SC;
        CE(ind) = out{ind}.CE;
        FF(ind) = out{ind}.FF;
        fnames{ind} = out{ind}.outstring;


        [n_new,sign] = set1(1,n_new,N_new(:,2));
        if sign, break; end
        n(N_new(:,1)) = n_new;
        ind = ind+1;
    end
  
    %}
end
% part 1 - save data for originLab (XYZ oputput format)
outXYZ = [];
str = '';
for n=length(X):-1:1
    outXYZ = [outXYZ,X{n}(:)];
    str = [str,vars{n,1},'	'];
end
outXYZ = [outXYZ,J_SC(:),V_OC(:)*1e3,CE(:),FF(:)];
str = [str,'J_SC	V_OC	CE	FF'];

ind = find(CE(:)==max(CE(:)));
fprintf('\nMaximum efficiency (#%i): ',ind)
fprintf('%g; ',outXYZ(ind,:))
fprintf('\t%s\t',originname)
fprintf('\n')

fout = fopen([originname,'XYZ.dat'],'w');
fprintf(fout,'%s\n',str);
fclose(fout);
save([originname,'XYZ.dat'],'outXYZ','-ascii','-double','-tabs','-append')

%% Show calculated data
MATLAB_version = version('-release');
MATLAB_version = num2str(MATLAB_version(1:end-1))+0.5*(MATLAB_version(end)=='b');

switch N1_new
    case 1
        if length(vars{1,2})==1
            return
        end
        name = vars{N_new(1),1};
        figure('Name',FinalDataFileName)
            subplot(221)
                plot(X{1},J_SC,'.-')
                xlabel(name,'Interpreter','none'), ylabel('J_{SC}, {mA/cm^2}')
            subplot(222)
                plot(X{1},V_OC,'.-')
                xlabel(name,'Interpreter','none'), ylabel('V_{OC}, mV')
            subplot(223)
                plot(X{1},CE,'.-')
                xlabel(name,'Interpreter','none'), ylabel('CE, %')
            subplot(224)
                plot(X{1},FF,'.-')
                xlabel(name,'Interpreter','none'), ylabel('FF, %')
savefig(gcf,[FinalDataFileName,'.fig'])

    case 2
        for i=N1_new:-1:1
            X{i} = reshape(X{i},N_new(:,2).');
        end
        V_OC = reshape(V_OC,N_new(:,2).');
        J_SC = reshape(J_SC,N_new(:,2).');
        CE = reshape(CE,N_new(:,2).');
        FF = reshape(FF,N_new(:,2).');
        fnames = reshape(fnames,N_new(:,2).');

        figure('Name',FinalDataFileName)
        Show1D(X{1}(:,1),vars{N_new(1),1},X{2}(1,:),vars{N_new(2),1},J_SC,  V_OC,  CE,  FF);

savefig(gcf,[FinalDataFileName,'.fig'])

        figure('Name',FinalDataFileName)
        Show1D(X{2}(1,:),vars{N_new(2),1},X{1}(:,1),vars{N_new(1),1},J_SC.',V_OC.',CE.',FF.');

        figure('Name',FinalDataFileName)
        Show2D(X,vars{N_new(2),1},vars{N_new(1),1},J_SC,V_OC,CE,FF);
    case 3
        for i=N1_new:-1:1
            X1D{i} = reshape(X{i},N_new(1,2),[]);
        end
        name1 = vars{N_new(1,1),1};
        name2 = vars{N_new(2,1),1};
        name3 = vars{N_new(3,1),1};
        x = X1D{1}(:,1);
        x2 = X1D{2}(1,:);
        x3 = X1D{3}(1,:);
        for n_t=numel(x2):-1:1
            dispname1{n_t} = [name2,'=',num2str(x2(n_t),'%g'),', ',name3,'=',num2str(x3(n_t),'%g')];
        end
        figure('Name',FinalDataFileName)
            subplot(221)
                if MATLAB_version>2014
                    for n_dm=1:length(dispname1)
                        semilogx(x,reshape(J_SC,N_new(1,2),[]),'.-','DisplayName',dispname1{n_dm})
                        hold on
                    end
                else
                    semilogx(x,reshape(J_SC,N_new(1,2),[]),'.-','DisplayName',dispname1)
                end
                xlabel(name1,'Interpreter','none'), ylabel('J_{SC}, {mA/cm^2}')
            subplot(222)
                if MATLAB_version>2014
                    for n_dm=1:length(dispname1)
                        semilogx(x,reshape(V_OC,N_new(1,2),[]),'.-','DisplayName',dispname1{n_dm})
                        hold on
                    end
                else
                    semilogx(x,reshape(V_OC,N_new(1,2),[]),'.-','DisplayName',dispname1)
                end
                xlabel(name1,'Interpreter','none'), ylabel('V_{OC}, mV')
            subplot(223)
                if MATLAB_version>2014
                    for n_dm=1:length(dispname1)
                        semilogx(x,reshape(CE,N_new(1,2),[]),'.-','DisplayName',dispname1{n_dm})
                        hold on
                    end
                else
                    semilogx(x,reshape(CE,N_new(1,2),[]),'.-','DisplayName',dispname1)
                end
                xlabel(name1,'Interpreter','none'), ylabel('CE, %')
            subplot(224)
                if MATLAB_version>2014
                    for n_dm=1:length(dispname1)
                        semilogx(x,reshape(FF,N_new(1,2),[]),'.-','DisplayName',dispname1{n_dm})
                        hold on
                    end
                else
                    semilogx(x,reshape(FF,N_new(1,2),[]),'.-','DisplayName',dispname1)
                end
                xlabel(name1,'Interpreter','none'), ylabel('FF, %')
                legend('show')
                
savefig(gcf,[FinalDataFileName,'.fig'])

        for i=N1_new:-1:1
            X{i} = reshape(X{i},N_new(:,2).');
        end
        V_OC = reshape(V_OC,N_new(:,2).');
        J_SC = reshape(J_SC,N_new(:,2).');
        CE = reshape(CE,N_new(:,2).');
        FF = reshape(FF,N_new(:,2).');
        fnames = reshape(fnames,N_new(:,2).');

        figure('Name',FinalDataFileName)
        for i=1:N_new(end)
            name = [vars{N_new(3),1},'=',num2str(X{3}(1,1,i),'%g')];
            subplot(221)
                if MATLAB_version>2014
                    for n_dm=1:length(dispname1)
                        scatter3(reshape(X{1}(:,:,i),[],1),reshape(X{2}(:,:,i),[],1),reshape(J_SC(:,:,i),[],1),'.','DisplayName',dispname1{n_dm})
                        hold on
                    end
                else
                    scatter3(reshape(X{1}(:,:,i),[],1),reshape(X{2}(:,:,i),[],1),reshape(J_SC(:,:,i),[],1),'.','DisplayName',dispname1)
                end
                hold on
            subplot(222)
                if MATLAB_version>2014
                    for n_dm=1:length(dispname1)
                        scatter3(reshape(X{1}(:,:,i),[],1),reshape(X{2}(:,:,i),[],1),reshape(V_OC(:,:,i),[],1),'.','DisplayName',dispname1{n_dm})
                        hold on
                    end
                else
                    scatter3(reshape(X{1}(:,:,i),[],1),reshape(X{2}(:,:,i),[],1),reshape(V_OC(:,:,i),[],1),'.','DisplayName',dispname1)
                end
                hold on
            subplot(223)
                if MATLAB_version>2014
                    for n_dm=1:length(dispname1)
                        scatter3(reshape(X{1}(:,:,i),[],1),reshape(X{2}(:,:,i),[],1),reshape(CE(:,:,i),[],1),'.','DisplayName',dispname1{n_dm})
                        hold on
                    end
                else
                    scatter3(reshape(X{1}(:,:,i),[],1),reshape(X{2}(:,:,i),[],1),reshape(CE(:,:,i),[],1),'.','DisplayName',dispname1)
                end
                hold on
            subplot(224)
                if MATLAB_version>2014
                    for n_dm=1:length(dispname1)
                        scatter3(reshape(X{1}(:,:,i),[],1),reshape(X{2}(:,:,i),[],1),reshape(FF(:,:,i),[],1),'.','DisplayName',dispname1{n_dm})
                        hold on
                    end
                else
                    scatter3(reshape(X{1}(:,:,i),[],1),reshape(X{2}(:,:,i),[],1),reshape(FF(:,:,i),[],1),'.','DisplayName',dispname1)
                end
                hold on
        end

        subplot(221)
            xlabel(vars{N_new(1),1},'Interpreter','none'), ylabel(vars{N_new(2),1},'Interpreter','none') 
            title('J_{SC}, {mA/cm^2}')
            hl = legend('show');
            set(hl,'Interpreter','none')
        subplot(222)
            xlabel(vars{N_new(1),1},'Interpreter','none'), ylabel(vars{N_new(2),1},'Interpreter','none') 
            title('V_{OC}, mV')
%             legend('show')
        subplot(223)
            xlabel(vars{N_new(1),1},'Interpreter','none'), ylabel(vars{N_new(2),1},'Interpreter','none') 
            title('CE, %')
%             legend('show')
        subplot(224)
            xlabel(vars{N_new(1),1},'Interpreter','none'), ylabel(vars{N_new(2),1},'Interpreter','none') 
            title('FF, %')
%             legend('show')
end

% part 2 - save data for originLab (1D format)
fn_tmp = [originname,'.dat'];
% Nx = N_new(1,2);
Nx = length(vars{1,2});
out_tmp = [reshape(X{1}(1:Nx),Nx,[]),reshape(J_SC,Nx,[]),reshape(V_OC,Nx,[]),reshape(CE,Nx,[]),reshape(FF,Nx,[])];
fnames0 = reshape(fnames,Nx,[]);
fout = fopen(fn_tmp,'w');
% fprintf(fout,'%s',vars{N_new(1),1});
fprintf(fout,'%s',vars{1,1});
if size(fnames0,2)==1
    fprintf(fout,'\t%s%g',fnames0{1});
else
    for n_val=1:size(fnames0,2)
        commas = regexp(fnames0{1,n_val},',');
        fprintf(fout,'\t%s%g',fnames0{1,n_val}(commas(1)+2:end));
    end
end
fprintf(fout,'\n');
fclose(fout);
save(fn_tmp,'out_tmp','-ascii','-double','-tabs','-append')
