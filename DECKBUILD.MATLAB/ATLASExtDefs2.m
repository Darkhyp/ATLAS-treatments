function [Fout,outstrings] = ATLASExtDefs2(F0,T0,T,aSiregnum,InputData,infilename)
isLinear = true;
isWorkfunc = false;
isWorkfunc = true;

Fout = []; 
% persistent defectXY def_ind    

q   = 1.602176565e-19;
k_B = 1.3806488e-23;

kT  = k_B/q*T;
kT0 = k_B/q*T0;

% Eg      = InputData.DPM.Eg;
% sigma   = InputData.DPM.sigma;
% U       = InputData.DPM.U;
% H       = InputData.DPM.H; % 1/cm3
% NSiSi   = InputData.DPM.NSiSi; % 1/cm3
% Ep      = InputData.DPM.Ep;

%%
ENDstr = {'solve init';
'output val.band charge';
'save outf=tmp.str'};
outfilename = 'TDPcalcs2.in';

% load file to memory and split on strings
if ~exist(infilename,'file')
    fprintf('\nFile %s is not found',infilename)
    return
end
m = memmapfile(infilename);
instrings  = strsplit(char(m.Data.'),'\n','CollapseDelimiters',true).';

%%
tmpstrings = {};
str_num=1;
while str_num<=length(instrings)
    getstr = instrings{str_num};
    str_num = str_num + 1;
    if isempty(getstr), continue, end

    tmpstrings{end+1} = getstr;
    while getstr(end)=='\' && str_num<length(instrings)
        getstr = instrings{str_num};
        str_num = str_num + 1;
        tmpstrings{end} = [tmpstrings{end}(1:end-1),' ',getstr];
    end    
end
instrings = tmpstrings;
checkstr = 'DEFECTS'; ind_DEFECTS = find(strncmpi(instrings,checkstr,length(checkstr)));
if isempty(ind_DEFECTS) % no defect statements
    outstrings = instrings;
    return
end
checkstr = 'solve'; ind_SOLVE = find(strncmpi(instrings,checkstr,length(checkstr)));
checkstr = 'mesh';  ind_MESH = find(strncmpi(instrings,checkstr,length(checkstr)));
checkstr = 'x.mesh';  ind_MESH2 = find(strncmpi(instrings,checkstr,length(checkstr)));
checkstr = 'y.mesh';  ind_MESH2 = [ind_MESH2,find(strncmpi(instrings,checkstr,length(checkstr)))];
checkstr = 'contact';  ind_CONTACT = find(strncmpi(instrings,checkstr,length(checkstr)));

% find workfunctions
for n_cont=length(ind_CONTACT):-1:1
    contact{n_cont} = strparser(instrings{ind_CONTACT(n_cont)});
end

% extract defect data
isREGnum = false;
for n_def=length(ind_DEFECTS):-1:1
    defectstr{n_def} = instrings{ind_DEFECTS(n_def)};
    defect{n_def} = strparser(instrings{ind_DEFECTS(n_def)});
    if defect{n_def}.number==aSiregnum
        isREGnum = true;
    end
end
if ~isREGnum % no aSi regions
    outstrings = instrings;
    return
end

% form deckbuild file (1)
outstrings1{1} = num2str(T,'MODELS TEMPERATURE=%g');
for str_num=ind_DEFECTS(end)+1:ind_SOLVE(1)-1
    outstrings1{end+1} = instrings{str_num};
end
for str_num=1:length(ENDstr)
    outstrings1{end+1} = ENDstr{str_num};
end

% prepare deckbuild run
if strcmpi(InputData.UNIXserverName,'silvacox2')
    InputData.OSstringEnd = ['; scp "',filenameIV0,'" korovin_ale@silvacox:silvaco_works/; scp "',filenameDark0,'" korovin_ale@silvacox:silvaco_works/; scp "',filenameLight,'" korovin_ale@silvacox:silvaco_works/"'];
end
    
n_iter = 2;
convergence_func = @(y0,y1) (1-InputData.DPM.alpha)*y0 + InputData.DPM.alpha*y1;            
while true
    % form deckbuild file (2)
    outstrings = {};
    for str_num=1:ind_DEFECTS(1)-1
        if (n_iter>=3)*0
            if str_num==ind_MESH
                outstrings{end+1} = 'MESH OUTFILE=tmp.str';
                continue
            end
            if sum(str_num==ind_MESH2)>0
                continue
            end
        end
        outstrings{end+1} = instrings{str_num};
    end
    
    % add defects
    for str_num=1:length(defectstr)
        outstrings{end+1} = defectstr{str_num};
    end
    % add final statements
    for str_num=1:length(outstrings1)
        outstrings{end+1} = outstrings1{str_num};
    end
    % save in-file
    fout = fopen(outfilename,'w');
    for n_str=1:length(outstrings)
        fprintf(fout,'%s\n',outstrings{n_str});
    end
    fclose(fout);

    if strcmpi(InputData.UNIXserverName,'silvacox2')
        system(['scp "',outfilename,'" korovin_ale@silvacox2:silvaco_works/'])
    end
    system([InputData.OSstring,outfilename,InputData.OSstringEnd]);

    % load Fermi level
    outT = TONYPLOTREAD('tmp.str',[]);
    if isempty(outT)
        break
    end

    if ~exist('x0','var')
        x = outT.nodes_xyz(:,2);
        y = outT.nodes_xyz(:,3);
        
        cond_aSi = outT.calcvalues(:,3)==aSiregnum;
        x0 = x( outT.calcvalues(cond_aSi,1)+1);
        y0 = y( outT.calcvalues(cond_aSi,1)+1);
        % boundary analyse
        for n_cont=1:length(contact)
            for n_cont2=1:length(outT.boundaries)
                if length(outT.boundaries{n_cont2}.w)<2, continue, end
                if strcmpi(outT.boundaries{n_cont2}.w{2},['"',contact{n_cont}.name,'"'])
                    [tmp1, tmp2] = ndgrid(outT.calcvalues(cond_aSi,1)+1,unique(outT.belements(outT.boundaries{n_cont2}.index,2:3)));
                    cond = sum(tmp1 == tmp2,2)==1;
                    contact{n_cont}.cond = cond;
                    % fig first closest
                    xx = x0(cond);
                    yy = y0(cond);
                end
            end
        end
    end
    % code = 111; % Electron QFL, eV
    % code = 113; % Ev, eV
    Fv = outT.calcvalues(cond_aSi,find(outT.valcodes==111)+1)-outT.calcvalues(cond_aSi,find(outT.valcodes==113)+1);

    % workfunction correction of Fv
    if isWorkfunc
        for n_cont=1:length(contact)
            if isfield(contact{n_cont},'workfun')
                Fv(contact{n_cont}.cond) = 3.87+1.7-contact{n_cont}.workfun;
            end
        end
    end

    Ev = outT.calcvalues(:,find(outT.valcodes==113)+1);
%     Ec = outT.calcvalues(:,find(outT.valcodes==114)+1);
%     Ey = outT.calcvalues(:,find(outT.valcodes==121)+1);
    Efull = outT.calcvalues(:,find(outT.valcodes==103)+1);
    Charge = outT.calcvalues(:,find(outT.valcodes==116)+1);

    
    % start
    Fout = [];
    xout = [];
    yout = [];
    for n_def=1:length(defect)
        if defect{n_def}.number==aSiregnum
            if n_iter==2
                defectXY{n_def}.cond = x0>=defect{n_def}.x_min & x0<=defect{n_def}.x_max & y0>=defect{n_def}.y_min & y0<=defect{n_def}.y_max;

                defectXY{n_def}.x  =  x0(defectXY{n_def}.cond);
                defectXY{n_def}.y  =  y0(defectXY{n_def}.cond);
                if n_def==1
                    defectXY{n_def}.N_xy = 1:length(defectXY{n_def}.x);
                else
                    defectXY{n_def}.N_xy = defectXY{n_def-1}.N_xy(end)+(1:length(defectXY{n_def}.x));
                end

                ngd_iter(defectXY{n_def}.N_xy,1) = defect{n_def}.ngd;
                egd_iter(defectXY{n_def}.N_xy,1) = defect{n_def}.egd;
                wgd_iter(defectXY{n_def}.N_xy,1) = defect{n_def}.wgd;
                nga_iter(defectXY{n_def}.N_xy,1) = defect{n_def}.nga;
                ega_iter(defectXY{n_def}.N_xy,1) = defect{n_def}.ega;
                wga_iter(defectXY{n_def}.N_xy,1) = defect{n_def}.wga;
            end
            defectXY{n_def}.Fv = Fv(defectXY{n_def}.cond);
Fout = [Fout;defectXY{n_def}.Fv];
xout = [xout;defectXY{n_def}.x];
yout = [yout;defectXY{n_def}.y];
        end
    end
Fv0(:,n_iter-1) = Fout;

    outtable = [];
    N_def = 0;
    for n_def=1:length(defect)
        N_def = N_def+1;
        if defect{n_def}.number==aSiregnum
%             fprintf('\n%i',n_def)
            % run the defect pool model
            if isempty(F0)
                [ngd,egd,wgd,nga,ega,wga] = func_defectPool([],defectXY{n_def}.Fv,kT0,kT,InputData.DPM.Eg,defect{n_def}.ntd,defect{n_def}.wtd,InputData.DPM.U,InputData.DPM.Ep,InputData.DPM.H,InputData.DPM.NSiSi,InputData.DPM.sigma);
            else
                [ngd,egd,wgd,nga,ega,wga] = func_defectPool(F0(defectXY{n_def}.N_xy),defectXY{n_def}.Fv,kT0,kT,InputData.DPM.Eg,defect{n_def}.ntd,defect{n_def}.wtd,InputData.DPM.U,InputData.DPM.Ep,InputData.DPM.H,InputData.DPM.NSiSi,InputData.DPM.sigma);
            end
            outtable = [outtable;[defectXY{n_def}.x,defectXY{n_def}.y,ngd(:),egd(:),wgd(:),nga(:),ega(:),wga(:),ones(size(ngd(:)))*[defect{n_def}.ntd,defect{n_def}.wtd,defect{n_def}.nta,defect{n_def}.wta]]];
            
            defectstr{N_def} = ['DEFECTS number=',num2str(defect{n_def}.number),...
                                ' X.MIN=',num2str(defect{n_def}.x_min), ...
                                ' X.MAX=',num2str(defect{n_def}.x_max), ...
                                ' Y.MIN=',num2str(defect{n_def}.y_min), ...
                                ' Y.MAX=',num2str(defect{n_def}.y_max), ...
                                ' F.DEFECTS=defects.c continuous numa=100 numd=100', ...
                                ' SIGTAE=',num2str(defect{n_def}.sigtae), ...
                                ' SIGTAH=',num2str(defect{n_def}.sigtah), ...
                                ' SIGTDE=',num2str(defect{n_def}.sigtde), ...
                                ' SIGTDH=',num2str(defect{n_def}.sigtdh), ...
                                ' SIGGAE=',num2str(defect{n_def}.siggae), ...
                                ' SIGGAH=',num2str(defect{n_def}.siggah), ...
                                ' SIGGDE=',num2str(defect{n_def}.siggde), ...
                                ' SIGGDH=',num2str(defect{n_def}.siggdh)];

        end
    end
    outtable(:,3) = convergence_func(ngd_iter(:,n_iter-1), outtable(:,3));
    outtable(:,4) = convergence_func(egd_iter(:,n_iter-1), outtable(:,4));
    outtable(:,5) = convergence_func(wgd_iter(:,n_iter-1), outtable(:,5));
    outtable(:,6) = convergence_func(nga_iter(:,n_iter-1), outtable(:,6));
    outtable(:,7) = convergence_func(ega_iter(:,n_iter-1), outtable(:,7));
    outtable(:,8) = convergence_func(wga_iter(:,n_iter-1), outtable(:,8));
    
    ngd_iter(:,n_iter) = outtable(:,3);
    egd_iter(:,n_iter) = outtable(:,4);
    wgd_iter(:,n_iter) = outtable(:,5);
    nga_iter(:,n_iter) = outtable(:,6);
    ega_iter(:,n_iter) = outtable(:,7);
    wga_iter(:,n_iter) = outtable(:,8);    
    
    test(n_iter-1) = norm(ngd_iter(:,n_iter-1)-ngd_iter(:,n_iter))+norm(nga_iter(:,n_iter-1)-nga_iter(:,n_iter));
    if test(end)<1e12 || n_iter>=InputData.DPM.Nmax_iter
        break
    end
    save('defs.dat','outtable', '-ascii','-double','-tabs');

    n_iter = n_iter+1;
    %%
    if InputData.isSaveData
        cond0 = x==0;
        out_band_tmp = [y(cond0),Ev(cond0),Efull(cond0),Charge(cond0)];
        if exist('out_band','var')
            N_tmp = max(size(out_band_tmp,1),size(out_band,1));
            if N_tmp>size(out_band,1), out_band(N_tmp,end) = 0; end
            if N_tmp>size(out_band_tmp,1), out_band_tmp(N_tmp,end) = 0; end
            out_band = [out_band,out_band_tmp];
            out_band_names = [out_band_names,'\ty, nm\tE_v, eV\tEful, V/cm\tCharge, C/cm'];
        else
            out_band = out_band_tmp;
            out_band_names = ['\tE_v, eV'];
        end
    end
    
nm = 1e-3;
    figure(100), 
        cond0 = x==0;
        subplot(321)
            cond1 = y<(min(y)+max(y))/2 & cond0;
            if isLinear
                plot((y(cond1))/nm,Ev(cond1))
                xlim([0,50])
            else
                semilogx((y(cond1)+1e-5)/nm,Ev(cond1))
            end
            hold on
            xlabel('from left side, nm')
            ylabel('E_v, eV')
            h1(1) = gca;
        subplot(323)
            if isLinear
                plot((y(cond1))/nm,Efull(cond1))
                xlim([0,50])
            else
                semilogx((y(cond1)+1e-5)/nm,Efull(cond1))
            end
            hold on
            xlabel('from left side, nm')
            ylabel('E_{full}, V/cm')
            h1(2) = gca;
        subplot(325)
            if isLinear
                plot((y(cond1))/nm,Charge(cond1))
                xlim([0,50])
            else
                semilogx((y(cond1)+1e-5)/nm,Charge(cond1))
            end
            hold on
            xlabel('from left side, nm')
            ylabel('Charge, C/cm')
            h1(3) = gca;

        subplot(322)
            cond1 = y>(min(y)+max(y))/2 & cond0;
            if isLinear
                plot((max(y)-y(cond1))/nm,Ev(cond1))
                xlim([0,50])
            else
                semilogx((max(y)-y(cond1)+1e-5)/nm,Ev(cond1))
            end
            hold on
            xlabel('from right side, nm')
            h2(1) = gca;
        subplot(324)
            cond1 = y>(min(y)+max(y))/2 & cond0;
            if isLinear
                plot((max(y)-y(cond1))/nm,Efull(cond1))
                xlim([0,50])
            else
                semilogx((max(y)-y(cond1)+1e-5)/nm,Efull(cond1))
            end
            hold on
            xlabel('from right side, nm')
            h2(2) = gca;
        subplot(326)
            cond1 = y>(min(y)+max(y))/2 & cond0;
            if isLinear
                plot((max(y)-y(cond1))/nm,Charge(cond1))
                xlim([0,50])
            else
                semilogx((max(y)-y(cond1)+1e-5),Charge(cond1))
            end
            hold on
            xlabel('from right side, nm')
            h2(3) = gca;

            set(h2,'XDir','reverse','YAxisLocation','right')
%             ylim(h1,[min([ylim(h1),ylim(h2)]),max([ylim(h1),ylim(h2)])])
%             ylim(h2,[min([ylim(h1),ylim(h2)]),max([ylim(h1),ylim(h2)])])
    
    
end

outstrings = {};
for str_num=1:ind_DEFECTS(1)-1
    outstrings{end+1} = instrings{str_num};
end
% add defects
for str_num=1:length(defectstr)
    outstrings{end+1} = defectstr{str_num};
end
% add old final statements
for str_num=ind_DEFECTS(end)+1:length(instrings)
    outstrings{end+1} = instrings{str_num};
end

FinalDataFileName='band';
    if InputData.isSaveData
        fout = fopen([FinalDataFileName,' band.dat'],'w');
        fprintf(fout,[out_band_names,'\n']);
        fclose(fout);
        save([FinalDataFileName,' band.dat'],'out_band','-ascii','-double','-tabs','-append')
    end

end


%% Additional functions

function obj = strparser(str)
    tmp = strsplit(str);
    for n=2:length(tmp)
        if ~isempty(regexpi(tmp{n},'=', 'once'))
            tmp2 = strsplit(tmp{n},'=','CollapseDelimiters',true);
            name = lower(regexprep(tmp2{1},'\.','_'));
            value = str2num(tmp2{2});
            if isempty(value)
                obj.(name) = tmp2{2};
            else
                obj.(name) = value;
            end
        else
            obj.(lower(regexprep(tmp{n},'\.','_'))) = 'true';
        end
    end
end
