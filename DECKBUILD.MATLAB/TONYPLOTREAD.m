%{
   TONYPLOTREAD.m -- This file is part of the treatment of
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

% Load data from Silvaco structure file for tonyplot
%% input variables
% infilename is the file name of Silvaco structure file
% regions is regions numers (all as default)
% varcode is the code for displayed variable
% code = 69; % Composition X
% code = -2; % abs Net doping, 1/cm3 ?
% code = 71; % donor conc, 1/cm3
% code = 72; % acceptor conc, 1/cm3
% code = 1035000; % donor Trap DOS #1, 1/cm3
% code = 1035001; % donor Trap DOS #2, 1/cm3
% code = 1036000; % donor ionized density #1, 1/cm3
% code = 1036001; % donor ionized density #2, 1/cm3
% code = 1038000; % acceptor Trap DOS #1, 1/cm3
% code = 1038001; % acceptor Trap DOS #2, 1/cm3
% code = 1039000; % acceptor ionized density #1, 1/cm3
% code = 1039001; % acceptor ionized density #2, 1/cm3
% code = 74; % ni, 1/cm3
% code = 75; % Nc, 1/cm3
% code = 76; % Nv, 1/cm3
% code = 430; % sqrt(n*p)(T), 1/cm3
% code = 431; % Nc(T), 1/cm3
% code = 432; % Nv(T), 1/cm3
% code = 115; % Net doping, 1/cm3
% code = 149; % Total doping (|Net doping|), 1/cm3
% code = 433; % Eg(T), eV
% code = 77; % Eg (band gap), eV
% code = 78; % electron affinity, eV
% code = 79; % optical intensity, W/cm2
% code = 208; % photogeneration rate, 1/cm3/s
% code = 209; % photon absorption rate, 1/cm3/s
% code = 118; % recombination rate, 1/cm3/s
% code = 616; % SRH recombination rate, 1/cm3/s
% code = 617; % Auger recombination rate, 1/cm3/s
% code = 100; % Potential, V
% code = 106; % electron conc, cm^(-3)
% code = 107; % hole conc, cm^(-3)
% code = 111; % Electron QFL, eV
% code = 112; % Hole QFL, eV
% code = 116; % Charge, C/cm2
% code = 113; % Ev, eV
% code = 435; % Ev(T), eV
% code = 114; % Ec, eV
% code = 434; % Ec(T), eV
% code = 120; % electric field Ex
% code = 121; % electric field Ey
% code = 103; % electric field E full
% code = 101; % e- mobility, cm2/V/s
% code = 143; % e- mobility X, cm2/V/s
% code = 144; % e- mobility Y, cm2/V/s
% code = 102; % h+ mobility, cm2/V/s
% code = 145; % h+ mobility X, cm2/V/s
% code = 146; % h+ mobility Y, cm2/V/s
% code = 210; % h+ current density, A/cm2
% code = 211; % e- current density, A/cm2
% code = 215; % Total current density, A/cm2
% code = 220; % Je- X, A/cm2
% code = 221; % Jh+ X, A/cm2
% code = 222; % Jtot X, A/cm2
% code = 224; % Je- Y, A/cm2
% code = 225; % Jh+ Y, A/cm2
% code = 226; % Jtot Y, A/cm2
% code = 524; % Cond current X, A/cm2
% code = 525; % Cond current Y, A/cm2
% code = 108; % Cond current density, A/cm2
% code = 103; % Eelectric field, V/cm
% code = 120; % E field X, V/cm
% code = 121; % E field Y, V/cm
%% output variables
% out.nodes_xyz     is the nodes coordinates [index, x, y, z]
% out.triangles	    is the elenets [index, region number, coordinate index1, coordinate index2, coordinate index3, coordinate index4,.??]
% out.boundaries	is the boundary defenitions
% out.belements	    is the boundary elements
% out.valcodes	    is the value codes
% out.calcvalues	is the calculated in Silvaco values with columns
% corresponding to out.valcodes

function out = TONYPLOTREAD(infilename,varcode,regions)


if ~exist('infilename','var')
    infilename = 'd:\Korovin\work\CNRS\LGEP\ExtGR0.158_AR_IBC_BSF1e+21_0.1homo2D_w30_c15_ratio0.05_emitter1e+21_0.1homo_w540_c30_half_Mesh1_100x300um_poly_nonperiodic_dark.str';
end
if ~exist('regions','var')
    regions = [];
end

% fin = fopen(infilename);
% if fin<=0
%     fprintf('\nFile %s is not found',infilename)
%     return
% end



% %{
tic
% load file to memory and split on strings
if ~exist(infilename,'file')
    fprintf('\nFile %s is not found',infilename)
    out = [];
    return
end
fprintf('\nLoad data from %s',infilename)
m = memmapfile(infilename);
s = char(m.Data.');
instrings = strsplit(s,'\n','CollapseDelimiters',true).';

% ind_calcvalues                = regexp(s,'\nn ');
% ind_calcvalues                = regexp(strings,'n ','once');
% ind_calcvalues                = isempty(ind_calcvalues{:});
% table(strings{ind_calcvalues})

N_nodes                     = length(regexp(s,'\nc '));
nodes_xyz(N_nodes,4)        = 0;
N_triangles                 = length(regexp(s,'\nt '));
triangles(N_triangles,8)    = 0;
N_belements                 = length(regexp(s,'\ne '));
belements(N_belements,4)    = 0;
N_calcvalues                = length(regexp(s,'\nn '));
calcvalues(N_calcvalues,1)  = 0; 
N_boundary                  = length(regexp(s,'\nr '));
boundaries{N_boundary,1}.r	= 0;
toc
%}

%{
nodes_xyz	= [];
triangles	= [];
belements	= [];
boundaries	= {};
calcvalues	= []; 
%}
% valcodes is column codes for calcvalues
% triangles(:,2) - region
% triangles(:,3:5) - all triangles
% triangles(:,6:8) - <0 (-1024;[87040,87040,87039])
% triangles(:,9) - ? (0,5)

unknown_j = [];
unknown_k = [];
valcodes = []; 
unknown_p = []; 
unknown_d = []; 
unknown_O = []; 

n_nodes         = 0;
n_triangles     = 0;
n_belements     = 0;
n_boundary      = 0;
n_calcvalues    = 0;


% str_num = 0;
% while ~feof(fin)
%     str_num = str_num + 1;
%     getstr = strtrim(fgetl(fin));
for str_num=1:length(instrings)
    getstr = strtrim(instrings{str_num});
%     getstr = strtrim(strings{str_num});
    %     fprintf('\nRead %i string',str_num)
    
    if isempty(getstr)
        continue % skip empty string
    end
    
    switch getstr(1)
        case 'c'
            n_nodes = n_nodes + 1;
            tmp = textscan(strtrim(getstr(2:end)),'%f','delimiter',{' ','\t'});
            nodes_xyz(n_nodes,:) = tmp{1}; 
        case 't'
            n_triangles = n_triangles + 1;
            tmp = textscan(strtrim(getstr(2:end)),'%f','delimiter',{' ','\t'});
            triangles(n_triangles,1:length(tmp{1})) = tmp{1}; 
        case 'e'
            n_belements = n_belements + 1;
            tmp = textscan(strtrim(getstr(2:end)),'%f','delimiter',{' ','\t'});
            belements(n_belements,:) = tmp{1}; 
        case 'j'
            tmp = textscan(strtrim(getstr(2:end)),'%f','delimiter',{' ','\t'});
            unknown_j(end+1,:) = tmp{1}; 
        case 'k'
            tmp = textscan(strtrim(getstr(2:end)),'%f','delimiter',{' ','\t'});
            unknown_k(end+1,:) = tmp{1}; 

        case 'r'
            n_boundary = n_boundary + 1;
            tmp = textscan(strtrim(getstr(2:end)),'%f','delimiter',{' ','\t'});
            boundaries{n_boundary}.r = tmp{1};
            boundaries{n_boundary}.index = [];
        case 'x'
            tmp = textscan(strtrim(getstr(2:end)),'%f','delimiter',{' ','\t'});
            boundaries{n_boundary}.x = tmp{1}; 
        case 'w'
            tmp = strsplit(strtrim(getstr(2:end)),{' ','\t'},'CollapseDelimiters',true);
            boundaries{n_boundary}.w = tmp; 
        case 'b'
            tmp = textscan(strtrim(getstr(2:end)),'%f','delimiter',{' ','\t'});
            boundaries{n_boundary}.index(end+1) = tmp{1}; 

        case 's'
            tmp = textscan(strtrim(getstr(2:end)),'%f','delimiter',{' ','\t'});
            valcodes(end+1,:) = tmp{1}; 
        case 'p'
            tmp = textscan(strtrim(getstr(2:end)),'%f','delimiter',{' ','\t'});
            unknown_p(end+1,:) = tmp{1}; 
        case 'd'
            tmp = textscan(strtrim(getstr(2:end)),'%f','delimiter',{' ','\t'});
            unknown_d(end+1,:) = tmp{1}; 
        case 'O'
            tmp = textscan(strtrim(getstr(2:end)),'%f','delimiter',{' ','\t'});
            unknown_O(end+1,:) = tmp{1}; 

        case 'n'
            n_calcvalues = n_calcvalues + 1;
            tmp = textscan(strtrim(getstr(2:end)),'%f','delimiter',{' ','\t'});
            calcvalues(n_calcvalues,1:length(tmp{1})) = tmp{1}; 
   end
    
end
toc

%% info
N_regions = max(triangles(:,2));
fprintf('\nNumber of nodes = %i',N_nodes);
fprintf('\n\ttriangles = %i',N_triangles);
fprintf('\n\tboundary elements = %i (%i sections: regions + contacts)',N_belements,N_boundary);
fprintf('\n\tcalculated nodes = %i',N_calcvalues);
fprintf('\n\tregions = %i.\n',N_regions);
out.nodes_xyz	= nodes_xyz;
out.triangles	= triangles;
out.boundaries	= boundaries;
out.belements	= belements;
out.valcodes	= valcodes;
out.calcvalues	= calcvalues;

if exist('varcode','var')
    if isempty(varcode)
       return
    end
    %% show data
    reg_colors = {[1,1,0.5];[1,0.5,1];[0.5,1,1];[0.5,0.5,1];[0.5,1,0.5]};
    % remove duplicated (why??)
    [~,ind] = unique(calcvalues(:,1));
    unique_n = calcvalues(ind,:);
    if varcode==113 || varcode==114
        %% band bending
        nc = find(valcodes==113)+1;
        nv = find(valcodes==114)+1;

        figure
        %{
        trisurf(triangles(:,3:5),nodes_xyz(:,2),nodes_xyz(:,3),calcvalues(ind,30),calcvalues(ind,30),'edgecolor','k','facecolor','interp');
        hold on
        trisurf(triangles(:,3:5),nodes_xyz(:,2),nodes_xyz(:,3),calcvalues(ind,31),calcvalues(ind,31),'edgecolor','k','facecolor','interp');
        %}
        for n_reg=1:N_regions
            if ~isempty(regions)
                if sum(regions==n_reg)==0
                    continue
                end
            end
        %     cond = unique_n(:,3)==n_reg;
            cond = triangles(:,2)==n_reg;
            ht1 = trisurf(triangles(cond,3:5),nodes_xyz(:,2),nodes_xyz(:,3),unique_n(:,nc),'FaceAlpha',0.7775,'FaceColor',reg_colors{n_reg},'edgecolor','k','displayname','Ec');
%             set(ht1,'edgecolor','none')
            hold on
            ht2 = trisurf(triangles(cond,3:5),nodes_xyz(:,2),nodes_xyz(:,3),unique_n(:,nv),'FaceAlpha',0.75,'FaceColor',reg_colors{n_reg},'edgecolor','k','displayname','Ev');
%             set(ht2,'edgecolor','none')
        end
        ht3 = trisurf(triangles(cond,3:5),nodes_xyz(:,2),nodes_xyz(:,3),0*unique_n(:,nv),'FaceAlpha',0.75,'FaceColor','w','edgecolor','k','displayname','Fermi level');
%         for n_reg=1:n_regions
%             cond = unique_n(:,3)==n_reg;
%             scatter3(nodes_xyz(cond,2),nodes_xyz(cond,3),unique_n(cond,30),'.','MarkerEdgeColor',reg_colors{n_reg});
%             hold on
%             scatter3(nodes_xyz(cond,2),nodes_xyz(cond,3),unique_n(cond,31),'.','MarkerEdgeColor',reg_colors{n_reg});
%         end
        zlabel('energy, eV')
    else
        n_val = find(valcodes==varcode)+1;
        figure
        ht = trisurf(triangles(:,3:5),nodes_xyz(:,2),nodes_xyz(:,3),unique_n(:,n_val),unique_n(:,n_val),'edgecolor','k','facecolor','interp');
        set(ht,'edgecolor','none')
    end
else
    return
end
Faces       = triangles(:,3:5);
Vertices	= nodes_xyz(:,2:3);

xlabel('x, {{\mu}m}')
ylabel('y, {{\mu}m}')
hax = gca;
xmin = min(nodes_xyz(:,2));
xmax = max(nodes_xyz(:,2));
ymin = min(nodes_xyz(:,3));
ymax = max(nodes_xyz(:,3));
% uicontrol('Style', 'popup','String', 'jet|hsv|hot|cool|gray','Position', [10 340 100 150],'Callback', @setmap);
% uicontrol('Style', 'pushbutton', 'String', 'Clear','Position', [10 10 50 20],'Callback', 'cla');
maxxy       = [xmax-xmin ymax-ymin];
pos         = [xmin+xmax ymin+ymax]/2;
dx          = maxxy;
hs_pos(1)	= uicontrol('Style','slider','Min',xmin,'Max',xmax,'Value',pos(1),'SliderStep',[1e-4,1e-2],'Position',[30,10,220,20],'Callback',{@shiftxy,hax,1});
hs_pos(2)	= uicontrol('Style','slider','Min',ymin,'Max',ymax,'Value',pos(2),'SliderStep',[1e-4,1e-2],'Position',[10,30,20,220],'Callback',{@shiftxy,hax,2});
    function shiftxy(source,event,ax,index) % #ok<INUSL>
        pos(index) = get(source,'Value');
        scaleTrisurf(ax);
    end

hs_zoom(1)	= uicontrol('Style','slider','Min',0,'Max',1e6,'Value',1e6,'SliderStep',[1e-4,1e-2],'Position',[430,470,220,20],'Callback',{@zoomxy,hax,1});
hs_zoom(2)	= uicontrol('Style','slider','Min',0,'Max',1e6,'Value',1e6,'SliderStep',[1e-4,1e-2],'Position',[650,250,20,220],'Callback',{@zoomxy,hax,2});
    function zoomxy(source,event,ax,index) % #ok<INUSL>
        dx(index) = maxxy(index)*get(source,'Value')/1e6/2;
        scaleTrisurf(ax);
    end

ht(1) = uicontrol('Style','text','Position',[260,10,50,20],'String','x pos');
ht(2) = uicontrol('Style','text','Position',[10,260,50,20],'String','y pos');
ht(3) = uicontrol('Style','text','Position',[360,470,60,20],'String','x zoom');

uicontrol('Style','checkbox','String', 'on/off','Position',[10,470,20,20],'Value',1,'Callback',@disableAll);
    function disableAll(source,~) % #ok<INUSL>
        if get(source,'Value')
            set(hs_pos,'Visible','on');
            set(hs_zoom,'Visible','on');
            set(ht,'Visible','on');
        else
            set(hs_pos,'Visible','off');
            set(hs_zoom,'Visible','off');
            set(ht,'Visible','off');
        end
    end


% ylim(gca,[99,100])
% set(hs_pos(1),'Position',[30,10,220,20])
% set(hs_pos(2),'Position',[10,30,20,220])
% 
% set(hsx_zoom,'Position',[430,470,220,20])
% set(hsy_zoom,'Position',[650,250,20,220])


return


figure
for n_reg=1:n_regions
    patch('vertices',nodes_xyz(:,2:3),'faces',triangles(triangles(:,2)==n_reg,3:5),'edgecol','k','facecol',reg_colors{n_reg})
    hold on
end
% or full
% patch('vertices',nodes_xyz(:,2:3),'faces',triangles(:,3:5),'edgecol','k','facecol',[1,1,.5])
hold on
%% regions boundary
N_BCs = length(boundaries); % region boundariey and contacts
for n_BCs=1:N_BCs
    BCind = belements(boundaries{n_BCs}.index,2:3);
    scatter3(nodes_xyz(BCind,2),nodes_xyz(BCind,3),nodes_xyz(BCind,4),'o','MarkerFaceColor',reg_colors{n_BCs})
end
xlabel('x, {{\mu}m}')
ylabel('y, {{\mu}m}')


n_val = 2; %(code = 29) ? const 143
n_val = 3; %(code = 513) region number




unknown_s(n_val-1)
figure
trisurf(triangles(:,3:5),nodes_xyz(:,2),nodes_xyz(:,3),calcvalues(ind,n_val),calcvalues(ind,n_val),'edgecolor','k','facecolor','interp');

set(gca,'zscale','log')




    function scaleTrisurf(ax)
        hc  = get(ax,'Children');

        xc = (Vertices(Faces(:,1),1)+Vertices(Faces(:,2),1)+Vertices(Faces(:,3),1))/3;
        yc = (Vertices(Faces(:,1),2)+Vertices(Faces(:,2),2)+Vertices(Faces(:,3),2))/3;
        cond = xc-pos(1)>=-dx(1) & xc-pos(1)<=dx(1) & yc-pos(2)>=-dx(2) & yc-pos(2)<=dx(2);
        for n_surf=1:length(hc)
        %     Faces = get(hc(n_surf),'Faces');
        %     Vertices = get(hc(n_surf),'Vertices');
        %     xc = (Vertices(Faces(:,1),index)+Vertices(Faces(:,2),index)+Vertices(Faces(:,3),index))/3;
        %     cond = xc-pos>=-val & xc-pos<=val;
            if sum(cond)~=0
                set(hc(n_surf),'Faces',Faces(cond,:));
            end
        end
    end


end

