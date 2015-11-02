% 21.10.15 were changed:
% cSiN_conc -> cSiN_concBSF
% cSiP_conc -> cSiP_concE

clc

%% switches
% InputData.isSeeResults = true;
% InputData.isLog = true;
% InputData.isReplace = true;
% InputData.isShowIV = true;
% InputData.isSaveData = true;
% InputData.isSeeResults = true;
% InputData.isTony = true;

%% structure definitions
% use this to change variable value in deckbuild structure file
%% (contacts width) simulations
%{
InputData.infilename = 'IBC_HIT_full9.in';
InputData.outfilename = 'IBC_HIT_full9_(b).in';
InputData.FinalDataFileName = 'd:\tempSC';
InputData.externvars0 = {...
    'cSiN_concBSF',                        '1e21';
    'cSiP_concE',                        '1e21';
    'np_x_gap_area',                    '0.05';
    'BSF_area',                         '0.05';
    'BSF_cSi_penetration_depth',        '0.1';
    'emitter_cSi_penetration_depth',    '0.1';
    'isStructureShow',                  '0';
    'DIR',                              '""';
    };

InputData.vars = {...
    'BSF_contact_area',                 [0.005,0.025,0.05];
    'emitter_contact_area',             [0.005,0.025,0.05,0.5];
    };
%}

%% (penetration depth) simulations
%{
InputData.infilename = 'IBC_HIT_full9.in';
InputData.outfilename = 'IBC_HIT_full9_(b).in';
InputData.FinalDataFileName = 'd:\tempSCpeneteraion';

InputData.externvars0 = {...
    'cSiN_concBSF',                        '1e21';
    'cSiP_concE',                        '1e21';
    'np_x_gap_area',                    '0.05';
    'BSF_contact_area',                 '0.005';
    'emitter_contact_area',             '0.005';
    'isStructureShow',                  '0';
    'DIR',                              '""';
    };

InputData.vars = {...
    'BSF_cSi_penetration_depth',        [0.05,0.1,0.2,0.5];
    'emitter_cSi_penetration_depth',    [0.05,0.1,0.2,0.5];
    'BSF_area',                         [0.01,0.05];
    };
%}

%{
InputData.infilename = 'IBC_HIT_full9.in';
InputData.outfilename = 'IBC_HIT_full9_(b).in';
InputData.FinalDataFileName = 'd:\tempSCpeneteraionnogap3';

% 'BSF_contact_area',                 '0.005';

InputData.externvars0 = {...
    'cSiN_concBSF',                        '1e21';
    'cSiP_concE',                        '1e21';
    'BSF_cSi_penetration_depth',        '0.05';
    'emitter_cSi_penetration_depth',    '0.05';
    'BSF_contact_area',                 '0.001';
    'emitter_contact_area',             '0.5';
    'np_x_gap_area',                    '0';
    'isStructureShow',                  '0';
    'DIR',                              '""';
    };

InputData.vars = {...
    'BSF_area',                         [0.001:0.001:0.006,0.008,0.01,0.015,0.02];
    };
%     'BSF_area',                         [0.005,0.0055,0.006,0.007,0.0075,0.0125,0.02:0.01:0.05];
%}

%{
InputData.infilename = 'IBC_HIT_full9.in';
InputData.outfilename = 'IBC_HIT_full9_(b).in';
InputData.FinalDataFileName = 'd:\tempSCpeneteraionnogap19';

InputData.externvars0 = {...
    'cSiN_concBSF',                        '1e19';
    'cSiP_concE',                        '1e19';
    'BSF_contact_area',                 '0.004';
    'BSF_area',                         '0.005';
    'emitter_contact_area',             '0.5';
    'np_x_gap_area',                    '0';
    'isStructureShow',                  '0';
    'DIR',                              '""';
    };

InputData.vars = {...
    'BSF_cSi_penetration_depth',        [0.01:0.02:0.06,0.1,0.2];
    'emitter_cSi_penetration_depth',    [0.01:0.02:0.07,0.1,0.15,0.2,0.3,0.5,1,2];
    };
%}

%{
InputData.infilename = 'IBC_HIT_full9.in';
InputData.outfilename = 'IBC_HIT_full9_(b).in';
InputData.FinalDataFileName = 'd:\tempSCpeneteraionnogapHET1';

InputData.externvars0 = {...
    'BSF_contact_area',                 '0.01';
    'BSF_area',                         '0.25';
    'emitter_contact_area',             '0.5';
    'np_x_gap_area',                    '0';
    'isStructureShow',                  '0';
    'DIR',                              '""';
    };

InputData.vars = {...
    'SweepType',                        2;
    };
%}

%{
InputData.infilename = 'IBC_HIT_full9.in';
InputData.outfilename = 'IBC_HIT_full9_(b).in';
%% heterojunctions
InputData.infilename = 'IBC_HIT_full10.in';
InputData.outfilename = 'IBC_HIT_full10_(b).in';

InputData.FinalDataFileName = 'd:\tempSCpeneteraionnogapHET';

InputData.externvars0 = {...
    'SweepType',                        '2';
    'cSi_thick',                        '100';
    'cSiN_concBSF_wafer',                  '1e15';
    'BSF_contact_area',                 '0.2';
    'emitter_contact_area',             '0.7';
    'np_x_gap_area',                    '0';
    'isStructureShow',                  '1';
    'DIR',                              '""';
    };

InputData.vars = {...
    'BSF_area',                         0.25;
    };
%}


%% new

%% heterojunctions + defects
%{
clear all
InputData.isReplaceMat = true;
InputData.isReplace = true;
    InputData.UNIXserverName = 'silvacox2';

InputData.infilename = 'IBC_HIT_full11.in';
InputData.outfilename = 'IBC_HIT_full11a.in';
% Ndef = 1e11;
% Ndef = 5e11;
Ndef = 1e12;
sigmadef = 1e-16;

InputData.outputDir = ['SC_hybrid_HETlayer',num2str(Ndef,'%g'),'/'];

InputData.FinalDataFileName = ['SC_hybrid_HET Ndeflayer',num2str(Ndef,'%g'),' s',num2str(sigmadef,'%g')];


InputData.externvars0 = {...
    'OutputString',                     'tmpHETdef';
    'cSi_thick',                        '100';
    'cSiN_concBSF_wafer',                  '1e15';
    'np_x_gap_area',                    '0';
    'DIR',                              '""';
    'DefectiveLayerThickness',          '0.001';
    'isInterfaceDefect',                '0';
    'Ndefsurf_BSF',                     num2str(Ndef);
    'SigmaDef_BSF',                     num2str(sigmadef);
    'Ndefsurf_emitter',                 num2str(Ndef);
    'SigmaDef_emitter',                 num2str(sigmadef);

    'aSiP_conc',                   '2.996e19';
    'NTAp',                        '1e22';
    'NTDp',                        '1e22';
    'WTAp',                        '0.03';
    'WTDp',                        '0.08';
    'NGAp',                        '1.76e19';
    'NGDp',                        '1.76e19';
    'EGAp',                        '1.7-1.3';
    'EGDp',                        '1.1';
    'WGAp',                        '0.24';
    'WGDp',                        '0.24';

    'aSiN_conc',                   '9.17e18';
    'NTAn',                        '1e22';
    'NTDn',                        '1e22';
    'WTAn',                        '0.03';
    'WTDn',                        '0.05';
    'NGAn',                        '1.76e19';
    'NGDn',                        '1.76e19';
    'EGAn',                        '1.7-0.7';
    'EGDn',                        '0.5';
    'WGAn',                        '0.24';
    'WGDn',                        '0.24';

    'i_aSi_thick',                      '0.00';

    'cSiP_concE',                        '1e21';
    'emitter_cSi_penetration_depth',    '0.01';
    'cSiN_concBSF',                        '1e21';
    'BSF_cSi_penetration_depth',	    '0.01';
    };
%     'isBSFdef',                         '1';
%     'isemitterdef',                     '1';
%     'SweepType',                        '2';

InputData.vars = {...
    'SweepType',                        2:4;
    'BSF_area',                         0.2:0.1:0.3;
    };
%     'SweepType',                        '3';
%     'isBSFdef',                         0:1;
%     'isemitterdef',                     0:1;

AtlasSimulations(InputData);
%}

%% heterojunctions (B,C,D types) + defects + intrinsic layer
% 21.10.15
%{
clear all
InputData.isReplaceMat = true;
InputData.isReplace = true;
    InputData.UNIXserverName = 'silvacox2';

InputData.infilename = 'IBC_HIT_full12.in';
InputData.outfilename = 'IBC_HIT_full12defs.in';
Ndef = 0;
% Ndef = 1e11;
% Ndef = 5e11;
sigmadef = 1e-16;

InputData.outputDir = ['SC_hybrid_HET_iaSi_Ndeflayer',num2str(Ndef,'%g'),'/'];

InputData.FinalDataFileName = ['SC_hybrid_HET_iaSi Ndeflayer',num2str(Ndef,'%g'),' s',num2str(sigmadef,'%g')];

InputData.externvars0 = {...
    'OutputString',                     'tmpHETdefint';
    'cSi_thick',                        '100';
    'cSiN_concBSF_wafer',                  '1e15';
    'np_x_gap_area',                    '0';
    'DIR',                              '""';
    'DefectiveLayerThickness',          '0.00';
    'isInterfaceDefect',                '0';
    'Ndefsurf_BSF',                     num2str(Ndef);
    'SigmaDef_BSF',                     num2str(sigmadef);
    'Ndefsurf_emitter',                 num2str(Ndef);
    'SigmaDef_emitter',                 num2str(sigmadef);
    'isBSFdef',                         '0';
    'isemitterdef',                     '0';
    'i_aSi_thick',                      '0.005';
    
    'aSiP_conc',                   '2.996e19';
    'NTAp',                        '1e22';
    'NTDp',                        '1e22';
    'WTAp',                        '0.03';
    'WTDp',                        '0.08';
    'NGAp',                        '1.76e19';
    'NGDp',                        '1.76e19';
    'EGAp',                        '1.7-1.3';
    'EGDp',                        '1.1';
    'WGAp',                        '0.24';
    'WGDp',                        '0.24';

    'aSiN_conc',                   '9.17e18';
    'NTAn',                        '1e22';
    'NTDn',                        '1e22';
    'WTAn',                        '0.03';
    'WTDn',                        '0.05';
    'NGAn',                        '1.76e19';
    'NGDn',                        '1.76e19';
    'EGAn',                        '1.7-0.7';
    'EGDn',                        '0.5';
    'WGAn',                        '0.24';
    'WGDn',                        '0.24';

    'aSiNi_concEmitter',            '1e15';
    'aSiPi_concEmitter',            '0';
    'NTAip',                        '1e22';
    'NTDip',                        '1e22';
    'WTAip',                        '0.03';
    'WTDip',                        '0.05';
    'NGAip',                        '1.054e17';
    'NGDip',                        '1.054e17';
    'EGAip',                        '1.7-1.027';
    'EGDip',                        '0.79';
    'WGAip',                        '0.525';
    'WGDip',                        '0.525';

    'aSiNi_concBSF',                '1e15';
    'aSiPi_concBSV',                '0';
    'NTAin',                        '1e22';
    'NTDin',                        '1e22';
    'WTAin',                        '0.03';
    'WTDin',                        '0.05';
    'NGAin',                        '1.054e17';
    'NGDin',                        '1.054e17';
    'EGAin',                        '1.7-1.027';
    'EGDin',                        '0.79';
    'WGAin',                        '0.525';
    'WGDin',                        '0.525';

    'cSiP_concE',                        '1e21';
    'emitter_cSi_penetration_depth',    '0.1';
    'cSiN_concBSF',                        '1e21';
    'BSF_cSi_penetration_depth',	    '0.1';
    };
%     'SweepType',                        '2';
%     'isBSFdef',                         '1';
%     'isemitterdef',                     '1';
% 

InputData.vars = {...
    'SweepType',                        2:4;
    'BSF_area',                         0.1:0.1:0.9;
    };
%     'isBSFdef',                         0:1;
%     'isemitterdef',                     0:1;
%     'SweepType',                        '3';

AtlasSimulations(InputData);
%}

%% heterojunctions + homojunction
% 22.10.15
%{
clear all
InputData.isReplaceMat = true;
InputData.isReplace = true;
    InputData.UNIXserverName = 'silvacox';
    InputData.ATLASparameter = ' -V 5.20.2.R';

InputData.infilename = 'IBC_HIT_full12.in';
InputData.outfilename = 'IBC_HIT_full12HHET.in';

Ndef = 0;
% Ndef = 1e11;
Ndef = 5e11;
sigmadef = 1e-16;

InputData.FinalDataFileName = ['SCdef',num2str(Ndef),'_s',num2str(sigmadef),'_HHET_IBC]'];

InputData.outputDir = ['SC_IBC/',InputData.FinalDataFileName,'/'];

InputData.externvars0 = {...
    'OutputString',                     'tmpHHET_cc';
    'cSi_thick',                        '100';
    'cSiN_concBSF_wafer',               '1e15';
    'np_x_gap_area',                    '0';
    'DIR',                              '""';
    'DefectiveLayerThickness',          '0.001';
    'isInterfaceDefect',                '0';
    'Ndefsurf_BSF',                     num2str(Ndef);
    'SigmaDef_BSF',                     num2str(sigmadef);
    'Ndefsurf_emitter',                 num2str(Ndef);
    'SigmaDef_emitter',                 num2str(sigmadef);
    'isBSFdef',                         '01';
    'isemitterdef',                     '01';
    'i_aSi_thick',                      '0.005*0';
    
    'aSiP_conc',                   '2.996e19';
    'NTAp',                        '1e22';
    'NTDp',                        '1e22';
    'WTAp',                        '0.03';
    'WTDp',                        '0.08';
    'NGAp',                        '1.76e19';
    'NGDp',                        '1.76e19';
    'EGAp',                        '1.7-1.3';
    'EGDp',                        '1.1';
    'WGAp',                        '0.24';
    'WGDp',                        '0.24';

    'aSiN_conc',                   '9.17e18';
    'NTAn',                        '1e22';
    'NTDn',                        '1e22';
    'WTAn',                        '0.03';
    'WTDn',                        '0.05';
    'NGAn',                        '1.76e19';
    'NGDn',                        '1.76e19';
    'EGAn',                        '1.7-0.7';
    'EGDn',                        '0.5';
    'WGAn',                        '0.24';
    'WGDn',                        '0.24';

    'aSiNi_concEmitter',            '1e15';
    'aSiPi_concEmitter',            '0';
    'NTAip',                        '1e22';
    'NTDip',                        '1e22';
    'WTAip',                        '0.03';
    'WTDip',                        '0.05';
    'NGAip',                        '1.054e17';
    'NGDip',                        '1.054e17';
    'EGAip',                        '1.7-1.027';
    'EGDip',                        '0.79';
    'WGAip',                        '0.525';
    'WGDip',                        '0.525';

    'aSiNi_concBSF',                '1e15';
    'aSiPi_concBSV',                '0';
    'NTAin',                        '1e22';
    'NTDin',                        '1e22';
    'WTAin',                        '0.03';
    'WTDin',                        '0.05';
    'NGAin',                        '1.054e17';
    'NGDin',                        '1.054e17';
    'EGAin',                        '1.7-1.027';
    'EGDin',                        '0.79';
    'WGAin',                        '0.525';
    'WGDin',                        '0.525';

    'emitter_cSi_penetration_depth',    '0.1';
    'BSF_cSi_penetration_depth',	    '0.1';
     
    'SweepType',                        '5';
    
    };
%     'cSiN_concBSF',                        '1e21';
%     'BSF_area',                         '0.1';

InputData.vars = {...
    'BSF_area',                         0.1:0.1:0.4;
    'cSiN_concBSF',                      logspace(16,19,7);
    'cSiP_concE',                        logspace(16,19,7);
    };
%     'SweepType',                        2:4;

%     'cellwidth',                        200:200:1000;
%     'cSiN_concBSF_wafer',                  [1e15,5e15,1e16];

AtlasSimulations(InputData);
%}

%% heterojunctions IBC-HET (type B)
%{
clear all
InputData.isReplaceMat = true;
InputData.isReplace = true;

InputData.infilename = 'IBC_HIT_full11.in';
InputData.outfilename = 'IBC_HIT_full11B.in';

InputData.outputDir = 'SC_IBC/HET_old SRV/';

InputData.FinalDataFileName = 'IBC-HET_old SRV';

% from Marie
InputData.externvars0 = {...
    'OutputString',                     'tmpHETB_old';
    'cSi_thick',                        '100';
    'np_x_gap_area',                    '0';
    'DIR',                              '""';

    'DefectiveLayerThickness',          '0.00';
    'isInterfaceDefect',                '0';
    'isBSFdef',                         '0';
    'isemitterdef',                     '0';
	'Ndefsurf_BSF',						'1e11';
	'SigmaDef_BSF'						'1e-16';
	'Ndefsurf_emitter'					'5e11';
	'SigmaDef_emitter'					'1e-16';

    'SweepType',                        '2';
    'cellwidth',                        '400';
    'cSiN_concBSF_wafer',                  '5e15';

    'i_aSi_thick',                      '0.00';
    'aSi_thick'                         '0.010';

    'BSF_area',                         '0.5';

	'aSiPi_concBSF'                     '0';
	'aSiNi_concEmitter'                 '0';
	
    };
%     'isBSFdef',                         '1';
%     'isemitterdef',                     '1';


%     'cSiN_concBSF_wafer',                  '1e15';
%     'Ndefsurf_BSF',                     num2str(Ndef);
%     'SigmaDef_BSF',                     num2str(sigmadef);
%     'Ndefsurf_emitter',                 num2str(Ndef);
%     'SigmaDef_emitter',                 num2str(sigmadef);

InputData.vars = {...
	'Sn_semsem'                         0:5:20;
	'Sp_semsem'                         0:5:20;
    };
%     'isBSFdef',                         0:1;
%     'isemitterdef',                     0:1;
%     'BSF_area',                         0.1:0.1:0.9;
%     'cellwidth',                        200:200:1000;
%     'cSiN_concBSF_wafer',                  [1e15,5e15,1e16];

%     'SweepType',                        2:4;
%     'isBSFdef',                         0:1;
%     'isemitterdef',                     0:1;
%     'SweepType',                        '3';


AtlasSimulations(InputData);
%}

%{
clear all
InputData.isReplaceMat = true;
InputData.isReplace = true;

InputData.infilename = 'IBC_HIT_full11.in';
InputData.outfilename = 'IBC_HIT_full11B.in';

InputData.outputDir = 'SC_IBC/HET_old check i/';

InputData.FinalDataFileName = 'IBC-HET_old check i';

% from Marie
InputData.externvars0 = {...
    'OutputString',                     'tmpHETB_old';
    'cSi_thick',                        '100';
    'np_x_gap_area',                    '0';
    'DIR',                              '""';

    'DefectiveLayerThickness',          '0.00';
    'isInterfaceDefect',                '0';
    'isBSFdef',                         '0';
    'isemitterdef',                     '0';
	'Ndefsurf_BSF',						'1e11';
	'SigmaDef_BSF'						'1e-16';
	'Ndefsurf_emitter'					'5e11';
	'SigmaDef_emitter'					'1e-16';

    'SweepType',                        '2';
    'cellwidth',                        '400';
    'cSiN_concBSF_wafer',                  '5e15';

	'aSi_thick'                         '0.010';
    'i_aSi_thick',                      '0.010';

    'BSF_area',                         '0.1';

	'aSiPi_concBSF'                     '0';
	'aSiNi_concEmitter'                 '0';
    };
%     'isBSFdef',                         '1';
%     'isemitterdef',                     '1';


%     'cSiN_concBSF_wafer',                  '1e15';
%     'Ndefsurf_BSF',                     num2str(Ndef);
%     'SigmaDef_BSF',                     num2str(sigmadef);
%     'Ndefsurf_emitter',                 num2str(Ndef);
%     'SigmaDef_emitter',                 num2str(sigmadef);

InputData.vars = {...
    'aSiNi_concBSF'                     [0 1e15];
	'aSiPi_concEmitter'                 [0 1e16];

    };
%     'isBSFdef',                         0:1;
%     'isemitterdef',                     0:1;
%     'BSF_area',                         0.1:0.1:0.9;
%     'cellwidth',                        200:200:1000;
%     'cSiN_concBSF_wafer',                  [1e15,5e15,1e16];

%     'SweepType',                        2:4;
%     'isBSFdef',                         0:1;
%     'isemitterdef',                     0:1;
%     'SweepType',                        '3';


AtlasSimulations(InputData);
%}

%{
clear all
InputData.isReplaceMat = true;
InputData.isReplace = true;

InputData.infilename = 'IBC_HIT_full11.in';
InputData.outfilename = 'IBC_HIT_full11B.in';

InputData.outputDir = 'SC_IBC/HET_old check 20/';

InputData.FinalDataFileName = 'IBC-HET_old check 20';

% from Marie
InputData.externvars0 = {...
    'OutputString',                     'tmpHETB_old';
    'cSi_thick',                        '100';
    'np_x_gap_area',                    '0';
    'DIR',                              '""';

    'DefectiveLayerThickness',          '0.00';
    'isInterfaceDefect',                '0';
    'isBSFdef',                         '0';
    'isemitterdef',                     '0';
	'Ndefsurf_BSF',						'1e11';
	'SigmaDef_BSF'						'1e-16';
	'Ndefsurf_emitter'					'5e11';
	'SigmaDef_emitter'					'1e-16';

    'SweepType',                        '2';
    'cellwidth',                        '400';
    'cSiN_concBSF_wafer',                  '5e15';

    'i_aSi_thick',                      '0.00';

    'BSF_area',                         '0.1';

	'aSiPi_concBSF'                     '0';
	'aSiNi_concEmitter'                 '0';
    };
%     'isBSFdef',                         '1';
%     'isemitterdef',                     '1';


%     'cSiN_concBSF_wafer',                  '1e15';
%     'Ndefsurf_BSF',                     num2str(Ndef);
%     'SigmaDef_BSF',                     num2str(sigmadef);
%     'Ndefsurf_emitter',                 num2str(Ndef);
%     'SigmaDef_emitter',                 num2str(sigmadef);

InputData.vars = {...
	'aSi_thick'                         [0.010,0.020];

    };
%     'isBSFdef',                         0:1;
%     'isemitterdef',                     0:1;
%     'BSF_area',                         0.1:0.1:0.9;
%     'cellwidth',                        200:200:1000;
%     'cSiN_concBSF_wafer',                  [1e15,5e15,1e16];

%     'SweepType',                        2:4;
%     'isBSFdef',                         0:1;
%     'isemitterdef',                     0:1;
%     'SweepType',                        '3';


AtlasSimulations(InputData);
%}

%{
clear all
InputData.isReplaceMat = true;
InputData.isReplace = true;

InputData.infilename = 'IBC_HIT_full11.in';
InputData.outfilename = 'IBC_HIT_full11B.in';

InputData.outputDir = 'SC_IBC/HET_old20nm/';

InputData.FinalDataFileName = 'IBC-HET_old20nm';

% from Marie
InputData.externvars0 = {...
    'OutputString',                     'tmpHETB_old';
    'cSi_thick',                        '100';
    'np_x_gap_area',                    '0';
    'DIR',                              '""';

    'DefectiveLayerThickness',          '0.00';
    'isInterfaceDefect',                '0';
    'isBSFdef',                         '0';
    'isemitterdef',                     '0';
	'Ndefsurf_BSF',						'1e11';
	'SigmaDef_BSF'						'1e-16';
	'Ndefsurf_emitter'					'5e11';
	'SigmaDef_emitter'					'1e-16';

    'SweepType',                        '2';
    'cellwidth',                        '400';

    'aSiP_conc',                   '2.996e19';
    'NTAp',                        '1e22';
    'NTDp',                        '1e22';
    'WTAp',                        '0.03';
    'WTDp',                        '0.08';
    'NGAp',                        '1.76e19';
    'NGDp',                        '1.76e19';
    'EGAp',                        '1.7-1.3';
    'EGDp',                        '1.1';
    'WGAp',                        '0.24';
    'WGDp',                        '0.24';

    'aSiN_conc',                   '9.17e18';
    'NTAn',                        '1e22';
    'NTDn',                        '1e22';
    'WTAn',                        '0.03';
    'WTDn',                        '0.05';
    'NGAn',                        '1.76e19';
    'NGDn',                        '1.76e19';
    'EGAn',                        '1.7-0.7';
    'EGDn',                        '0.5';
    'WGAn',                        '0.24';
    'WGDn',                        '0.24';

	'aSi_thick'                         '0.020';
    'i_aSi_thick',                      '0.00';
    };
%     'isBSFdef',                         '1';
%     'isemitterdef',                     '1';


%     'cSiN_concBSF_wafer',                  '1e15';
%     'Ndefsurf_BSF',                     num2str(Ndef);
%     'SigmaDef_BSF',                     num2str(sigmadef);
%     'Ndefsurf_emitter',                 num2str(Ndef);
%     'SigmaDef_emitter',                 num2str(sigmadef);

InputData.vars = {...
    'BSF_area',                         0.05:0.05:0.95;
    };
%     'isBSFdef',                         0:1;
%     'isemitterdef',                     0:1;
%     'BSF_area',                         0.1:0.1:0.9;
%     'cellwidth',                        200:200:1000;
%     'cSiN_concBSF_wafer',                  [1e15,5e15,1e16];

%     'SweepType',                        2:4;
%     'isBSFdef',                         0:1;
%     'isemitterdef',                     0:1;
%     'SweepType',                        '3';


AtlasSimulations(InputData);
%}

%{
clear all
InputData.isReplaceMat = true;
InputData.isReplace = true;

InputData.infilename = 'IBC_HIT_full11.in';
InputData.outfilename = 'IBC_HIT_full11B.in';

InputData.outputDir = 'SC_IBC/HETdef1e11_old/';

InputData.FinalDataFileName = 'IBC-HETdef1e11_old2';

% from Marie
InputData.externvars0 = {...
    'OutputString',                     'tmpHETBdef2_old';
    'cSi_thick',                        '100';
    'np_x_gap_area',                    '0';
    'DIR',                              '""';

    'DefectiveLayerThickness',          '0.001';
    'isInterfaceDefect',                '0';
	'Ndefsurf_BSF',						'1e11';
	'SigmaDef_BSF'						'1e-16';
	'Ndefsurf_emitter'					'5e11';
	'SigmaDef_emitter'					'1e-16';

    'SweepType',                        '2';

    'aSiP_conc',                   '2.996e19';
    'NTAp',                        '1e22';
    'NTDp',                        '1e22';
    'WTAp',                        '0.03';
    'WTDp',                        '0.08';
    'NGAp',                        '1.76e19';
    'NGDp',                        '1.76e19';
    'EGAp',                        '1.7-1.3';
    'EGDp',                        '1.1';
    'WGAp',                        '0.24';
    'WGDp',                        '0.24';

    'aSiN_conc',                   '9.17e18';
    'NTAn',                        '1e22';
    'NTDn',                        '1e22';
    'WTAn',                        '0.03';
    'WTDn',                        '0.05';
    'NGAn',                        '1.76e19';
    'NGDn',                        '1.76e19';
    'EGAn',                        '1.7-0.7';
    'EGDn',                        '0.5';
    'WGAn',                        '0.24';
    'WGDn',                        '0.24';

    'i_aSi_thick',                      '0.00';
    };
%     'isBSFdef',                         '1';
%     'isemitterdef',                     '1';


%     'cSiN_concBSF_wafer',                  '1e15';
%     'Ndefsurf_BSF',                     num2str(Ndef);
%     'SigmaDef_BSF',                     num2str(sigmadef);
%     'Ndefsurf_emitter',                 num2str(Ndef);
%     'SigmaDef_emitter',                 num2str(sigmadef);

InputData.vars = {...
    'BSF_area',                         0.05:0.05:0.95;
    'isBSFdef',                         0:1;
    'isemitterdef',                     0:1;
    };
%     'BSF_area',                         0.1:0.1:0.9;
%     'cellwidth',                        200:200:1000;
%     'cSiN_concBSF_wafer',                  [1e15,5e15,1e16];

%     'SweepType',                        2:4;
%     'isBSFdef',                         0:1;
%     'isemitterdef',                     0:1;
%     'SweepType',                        '3';


AtlasSimulations(InputData);
%}


%% homojunctions IBC-HOMO (type A1 - E changed)
%{
clear all
InputData.isReplaceMat = true;
InputData.isReplace = true;

InputData.infilename = 'IBC_HIT_full12.in';
InputData.outfilename = 'IBC_HIT_full12_tA1.in';

InputData.outputDir = 'SC_IBC/HOMO1_2/';

InputData.FinalDataFileName = 'IBC-HOMO1_2';

InputData.externvars0 = {...
    'OutputString',                     'tmpHOMOtA1';
    'cSi_thick',                        '100';
    'np_x_gap_area',                    '0';
    'DIR',                              '""';
    'DefectiveLayerThickness',          '0.00';
    'isInterfaceDefect',                '0';
    'isBSFdef',                         '0';
    'isemitterdef',                     '0';
    'SweepType',                        '1';
    'cSiN_concBSF',                        '1e21';
    };
%     'cSiP_concE',                        '1e21';

InputData.vars = {...
    'BSF_area',                         0.1:0.1:0.9;
    'cSiP_concE',                        logspace(18,21,10);
    };
%     'cSiN_concBSF',                        logspace(18,21,10);

%     'cellwidth',                        200:200:1000;
%     'cSiN_concBSF_wafer',                  [1e15,5e15,1e16];

AtlasSimulations(InputData);
%}
%% (penetration depth) simulations IBC-HOMO (type A1)
%{
clear all
InputData.isReplaceMat = true;
InputData.isReplace = true;

InputData.infilename = 'IBC_HIT_full12.in';
InputData.outfilename = 'IBC_HIT_full12_A1.in';

InputData.outputDir = 'SC_IBC/HOMO1_p/';

InputData.FinalDataFileName = 'IBC-HOMO1_p';

InputData.externvars0 = {...
    'OutputString',                     'tmpHOMOA1';
    'cSi_thick',                        '100';
    'np_x_gap_area',                    '0';
    'DIR',                              '""';
    'DefectiveLayerThickness',          '0.00';
    'isInterfaceDefect',                '0';
    'isBSFdef',                         '0';
    'isemitterdef',                     '0';
    'SweepType',                        '1';
    'cSiN_concBSF',                        '1e21';
    'cSiP_concE',                        '1e21';

    'BSF_area',                         '0.1';

    };
%     'cSiP_concE',                        '1e21';

InputData.vars = {...
    'BSF_cSi_penetration_depth',        [0.01,0.02,0.03,0.05,0.07,0.1:0.1:1];
    };
%     'cSiP_concE',                        logspace(18,21,10);
%     'BSF_area',                         0.1:0.1:0.9;
%     'cSiN_concBSF',                        logspace(18,21,10);

%     'cellwidth',                        200:200:1000;
%     'cSiN_concBSF_wafer',                  [1e15,5e15,1e16];

AtlasSimulations(InputData);
%}

%% homojunctions IBC-HOMO (type A2 - BSF changed)
%{
clear all
InputData.isReplaceMat = true;
InputData.isReplace = true;

InputData.infilename = 'IBC_HIT_full12.in';
InputData.outfilename = 'IBC_HIT_full12_tA2.in';

InputData.outputDir = 'SC_IBC/HOMO2_2/';

InputData.FinalDataFileName = 'IBC-HOMO2_2';

InputData.externvars0 = {...
    'OutputString',                     'tmpHOMOtA2';
    'cSi_thick',                        '100';
    'np_x_gap_area',                    '0';
    'DIR',                              '""';
    'DefectiveLayerThickness',          '0.00';
    'isInterfaceDefect',                '0';
    'isBSFdef',                         '0';
    'isemitterdef',                     '0';
    'SweepType',                        '1';
    'cSiP_concE',                        '1e21';
    };
%     'cSiN_concBSF',                        '1e21';

InputData.vars = {...
    'BSF_area',                         0.1:0.1:0.9;
    'cSiN_concBSF',                        logspace(18,21,10);
    };
%     'cSiP_concE',                        logspace(18,21,10);

%     'cellwidth',                        200:200:1000;
%     'cSiN_concBSF_wafer',                  [1e15,5e15,1e16];

AtlasSimulations(InputData);
%}
%% (penetration depth) simulations IBC-HOMO (type A2)
%{
clear all
InputData.isReplaceMat = true;
InputData.isReplace = true;

InputData.infilename = 'IBC_HIT_full12.in';
InputData.outfilename = 'IBC_HIT_full12_A2.in';

InputData.outputDir = 'SC_IBC/HOMO2_p/';

InputData.FinalDataFileName = 'IBC-HOMO2_p';

InputData.externvars0 = {...
    'OutputString',                     'tmpHOMOA2';
    'cSi_thick',                        '100';
    'np_x_gap_area',                    '0';
    'DIR',                              '""';
    'DefectiveLayerThickness',          '0.00';
    'isInterfaceDefect',                '0';
    'isBSFdef',                         '0';
    'isemitterdef',                     '0';
    'SweepType',                        '1';
    'cSiP_concE',                        '1e21';
    'cSiN_concBSF',                        '1e21';
    
    'BSF_area',                         '0.1';
    };
%     'cSiN_concBSF',                        '1e21';

InputData.vars = {...
    'emitter_cSi_penetration_depth',	[0.01,0.02,0.03,0.05,0.07,0.1:0.1:1];
    };
%     'BSF_area',                         0.1:0.1:0.9;
%     'cSiN_concBSF',                        logspace(18,21,10);
%     'cSiP_concE',                        logspace(18,21,10);

%     'cellwidth',                        200:200:1000;
%     'cSiN_concBSF_wafer',                  [1e15,5e15,1e16];

AtlasSimulations(InputData);
%}

%% (penetration depth = 0.1) simulations IBC-HOMO (type A as function of (conc_E, conc_BSF))
%{
%21.10.15
clear all
InputData.isReplaceMat = true;
InputData.isReplace = true;
    InputData.UNIXserverName = 'silvacox';

InputData.infilename = 'IBC_HIT_full12.in';
InputData.outfilename = 'IBC_HIT_full12_Acc.in';

depthE = [0.01,0.1];
depthBSF = [0.01,0.1];

for nE=1:length(depthE)
for nBSF=1:length(depthBSF)
penetration_depthE = depthE(nE);
penetration_depthBSF = depthBSF(nBSF);

InputData.FinalDataFileName = ['HOMO_cc[depthE',num2str(penetration_depthE),',depthBSF',num2str(penetration_depthBSF),']'];

InputData.outputDir = ['SC_IBC/',InputData.FinalDataFileName,'/'];

InputData.externvars0 = {...
    'OutputString',                     'tmpHOMO_cc';
    'cSi_thick',                        '100';
    'np_x_gap_area',                    '0';
    'DIR',                              '""';
    'DefectiveLayerThickness',          '0.00';
    'isInterfaceDefect',                '0';
    'isBSFdef',                         '0';
    'isemitterdef',                     '0';
    'SweepType',                        '1';
    'cSiP_concE',                       '1e21';
    'cSiN_concBSF',                     '1e21';
    
    'BSF_area',                         '0.1';
    'emitter_cSi_penetration_depth',	num2str(penetration_depthE);
    'BSF_cSi_penetration_depth',        num2str(penetration_depthBSF);
    };
%     'cSiN_concBSF',                        '1e21';

InputData.vars = {...
    'cSiN_concBSF',                        logspace(18,21,10);
    'cSiP_concE',                        logspace(18,21,10);
    };
%     'BSF_area',                         0.1:0.1:0.9;

%     'cellwidth',                        200:200:1000;
%     'cSiN_concBSF_wafer',                  [1e15,5e15,1e16];

AtlasSimulations(InputData);
end
end
%}


%% heterojunctions IBC HET-BSF/HOMO-E (type C)
% with defects (Ns = 1e11 1/cm3)
%{
clear all
InputData.isReplaceMat = true;
InputData.isReplace = true;


InputData.infilename = 'IBC_HIT_full11.in';
InputData.outfilename = 'IBC_HIT_full11C.in';

InputData.outputDir = 'SC_IBC/HYBRID_Cdef1e11_old2/';

InputData.FinalDataFileName = 'IBC-HYBRID_Cdef1e11_old2';

InputData.externvars0 = {...
    'OutputString',                     'tmpHYBRIDCdef2_old';
    'cSi_thick',                        '100';
    'np_x_gap_area',                    '0';
    'DIR',                              '""';

    'DefectiveLayerThickness',          '0.001';
    'isInterfaceDefect',                '0';
    'isBSFdef',                         '1';
    'isemitterdef',                     '1';
	'Ndefsurf_BSF',						'1e11';
	'SigmaDef_BSF'						'1e-16';
	'Ndefsurf_emitter'					'1e11';
	'SigmaDef_emitter'					'1e-16';

    'SweepType',                        '3';

    'aSiP_conc',                   '2.996e19';
    'NTAp',                        '1e22';
    'NTDp',                        '1e22';
    'WTAp',                        '0.03';
    'WTDp',                        '0.08';
    'NGAp',                        '1.76e19';
    'NGDp',                        '1.76e19';
    'EGAp',                        '1.7-1.3';
    'EGDp',                        '1.1';
    'WGAp',                        '0.24';
    'WGDp',                        '0.24';

    'aSiN_conc',                   '9.17e18';
    'NTAn',                        '1e22';
    'NTDn',                        '1e22';
    'WTAn',                        '0.03';
    'WTDn',                        '0.05';
    'NGAn',                        '1.76e19';
    'NGDn',                        '1.76e19';
    'EGAn',                        '1.7-0.7';
    'EGDn',                        '0.5';
    'WGAn',                        '0.24';
    'WGDn',                        '0.24';

    'i_aSi_thick',                      '0.00';

    'cSiP_concE',                        '1e21';
    'emitter_cSi_penetration_depth',    '0.01';
    'cSiN_concBSF',                        '1e21';
    'BSF_cSi_penetration_depth',	    '0.01';
    };

InputData.vars = {...
    'BSF_area',                         0.1:0.1:0.9;
    };
%     'cSiP_concE',                        logspace(18,21,10);
%     'cellwidth',                        200:200:1000;
%     'cSiN_concBSF_wafer',                  [1e15,5e15,1e16];

AtlasSimulations(InputData);
%}
% with defects (Ns = 5e11 1/cm3)
%{
clear all
InputData.isReplaceMat = true;
InputData.isReplace = true;


InputData.infilename = 'IBC_HIT_full11.in';
InputData.outfilename = 'IBC_HIT_full11C.in';

InputData.outputDir = 'SC_IBC/HYBRID_Cdef5e11_old2/';

InputData.FinalDataFileName = 'IBC-HYBRID_Cdef5e11_old2';

InputData.externvars0 = {...
    'OutputString',                     'tmpHYBRIDCdef2_old';
    'cSi_thick',                        '100';
    'np_x_gap_area',                    '0';
    'DIR',                              '""';

    'DefectiveLayerThickness',          '0.001';
    'isInterfaceDefect',                '0';
    'isBSFdef',                         '1';
    'isemitterdef',                     '1';
	'Ndefsurf_BSF',						'5e11';
	'SigmaDef_BSF'						'1e-16';
	'Ndefsurf_emitter'					'5e11';
	'SigmaDef_emitter'					'1e-16';

    'SweepType',                        '3';

    'aSiP_conc',                   '2.996e19';
    'NTAp',                        '1e22';
    'NTDp',                        '1e22';
    'WTAp',                        '0.03';
    'WTDp',                        '0.08';
    'NGAp',                        '1.76e19';
    'NGDp',                        '1.76e19';
    'EGAp',                        '1.7-1.3';
    'EGDp',                        '1.1';
    'WGAp',                        '0.24';
    'WGDp',                        '0.24';

    'aSiN_conc',                   '9.17e18';
    'NTAn',                        '1e22';
    'NTDn',                        '1e22';
    'WTAn',                        '0.03';
    'WTDn',                        '0.05';
    'NGAn',                        '1.76e19';
    'NGDn',                        '1.76e19';
    'EGAn',                        '1.7-0.7';
    'EGDn',                        '0.5';
    'WGAn',                        '0.24';
    'WGDn',                        '0.24';

    'i_aSi_thick',                      '0.00';

    'cSiP_concE',                        '1e21';
    'emitter_cSi_penetration_depth',    '0.01';
    'cSiN_concBSF',                        '1e21';
    'BSF_cSi_penetration_depth',	    '0.01';
    };

InputData.vars = {...
    'BSF_area',                         0.1:0.1:0.9;
    };
%     'cSiP_concE',                        logspace(18,21,10);
%     'cellwidth',                        200:200:1000;
%     'cSiN_concBSF_wafer',                  [1e15,5e15,1e16];

AtlasSimulations(InputData);
%}
% with defects (Ns = 1e12 1/cm3)
%{
clear all
InputData.isReplaceMat = true;
InputData.isReplace = true;


InputData.infilename = 'IBC_HIT_full11.in';
InputData.outfilename = 'IBC_HIT_full11C.in';

InputData.outputDir = 'SC_IBC/HYBRID_Cdef1e12_old2/';

InputData.FinalDataFileName = 'IBC-HYBRID_Cdef1e12_old2';

InputData.externvars0 = {...
    'OutputString',                     'tmpHYBRIDCdef2_old';
    'cSi_thick',                        '100';
    'np_x_gap_area',                    '0';
    'DIR',                              '""';

    'DefectiveLayerThickness',          '0.001';
    'isInterfaceDefect',                '0';
    'isBSFdef',                         '1';
    'isemitterdef',                     '1';
	'Ndefsurf_BSF',						'1e12';
	'SigmaDef_BSF'						'1e-16';
	'Ndefsurf_emitter'					'5e11';
	'SigmaDef_emitter'					'1e-16';

    'SweepType',                        '3';

    'aSiP_conc',                   '2.996e19';
    'NTAp',                        '1e22';
    'NTDp',                        '1e22';
    'WTAp',                        '0.03';
    'WTDp',                        '0.08';
    'NGAp',                        '1.76e19';
    'NGDp',                        '1.76e19';
    'EGAp',                        '1.7-1.3';
    'EGDp',                        '1.1';
    'WGAp',                        '0.24';
    'WGDp',                        '0.24';

    'aSiN_conc',                   '9.17e18';
    'NTAn',                        '1e22';
    'NTDn',                        '1e22';
    'WTAn',                        '0.03';
    'WTDn',                        '0.05';
    'NGAn',                        '1.76e19';
    'NGDn',                        '1.76e19';
    'EGAn',                        '1.7-0.7';
    'EGDn',                        '0.5';
    'WGAn',                        '0.24';
    'WGDn',                        '0.24';

    'i_aSi_thick',                      '0.00';

    'cSiP_concE',                        '1e21';
    'emitter_cSi_penetration_depth',    '0.01';
    'cSiN_concBSF',                        '1e21';
    'BSF_cSi_penetration_depth',	    '0.01';
    };

InputData.vars = {...
    'BSF_area',                         0.1:0.1:0.9;
    };
%     'cSiP_concE',                        logspace(18,21,10);
%     'cellwidth',                        200:200:1000;
%     'cSiN_concBSF_wafer',                  [1e15,5e15,1e16];

AtlasSimulations(InputData);
%}

% without defects
%{
clear all
InputData.isReplaceMat = true;
InputData.isReplace = true;


InputData.infilename = 'IBC_HIT_full12.in';
InputData.outfilename = 'IBC_HIT_full12_tC.in';

InputData.outputDir = 'SC_IBC/HYBRID_C_old_2/';

InputData.FinalDataFileName = 'IBC-HYBRID_C_old_2';

InputData.externvars0 = {...
    'OutputString',                     'tmpHYBRID_C_old_2';
    'cSi_thick',                        '100';
    'np_x_gap_area',                    '0';
    'DIR',                              '""';

    'DefectiveLayerThickness',          '0.00';
    'isInterfaceDefect',                '0';
    'isBSFdef',                         '0';
    'isemitterdef',                     '0';

    'SweepType',                        '3';

    'aSiP_conc',                   '2.996e19';
    'NTAp',                        '1e22';
    'NTDp',                        '1e22';
    'WTAp',                        '0.03';
    'WTDp',                        '0.08';
    'NGAp',                        '1.76e19';
    'NGDp',                        '1.76e19';
    'EGAp',                        '1.7-1.3';
    'EGDp',                        '1.1';
    'WGAp',                        '0.24';
    'WGDp',                        '0.24';

    'aSiN_conc',                   '9.17e18';
    'NTAn',                        '1e22';
    'NTDn',                        '1e22';
    'WTAn',                        '0.03';
    'WTDn',                        '0.05';
    'NGAn',                        '1.76e19';
    'NGDn',                        '1.76e19';
    'EGAn',                        '1.7-0.7';
    'EGDn',                        '0.5';
    'WGAn',                        '0.24';
    'WGDn',                        '0.24';

    'i_aSi_thick',                      '0.00';

    'cSiP_concE',                        '1e21';
    'cSiN_concBSF',                        '1e21';
     };

InputData.vars = {...
    'BSF_area',                         0.1:0.1:0.9;
    'cSiP_concE',                        logspace(18,21,10);
    };
%     'cellwidth',                        200:200:1000;
%     'cSiN_concBSF_wafer',                  [1e15,5e15,1e16];

AtlasSimulations(InputData);
%}
% various penetration depth of dopant in homojunction (HET with defects Ns = 5e11 1/cm3)
%{
clear all
InputData.isReplaceMat = true;
InputData.isReplace = true;


InputData.infilename = 'IBC_HIT_full12.in';
InputData.outfilename = 'IBC_HIT_full12C.in';

InputData.outputDir = 'SC_IBC/HYBRID_C_old_p/';

InputData.FinalDataFileName = 'IBC-HYBRID_C_old_p';

InputData.externvars0 = {...
    'OutputString',                     'tmpHYBRIDC_old';
    'cSi_thick',                        '100';
    'np_x_gap_area',                    '0';
    'DIR',                              '""';

    'DefectiveLayerThickness',          '0.001';
    'isInterfaceDefect',                '0';
    'isBSFdef',                         '1';
    'isemitterdef',                     '1';
	'Ndefsurf_BSF',						'5e11';
	'SigmaDef_BSF'						'1e-16';
	'Ndefsurf_emitter'					'5e11';
	'SigmaDef_emitter'					'1e-16';

    'SweepType',                        '3';

    'aSiP_conc',                   '2.996e19';
    'NTAp',                        '1e22';
    'NTDp',                        '1e22';
    'WTAp',                        '0.03';
    'WTDp',                        '0.08';
    'NGAp',                        '1.76e19';
    'NGDp',                        '1.76e19';
    'EGAp',                        '1.7-1.3';
    'EGDp',                        '1.1';
    'WGAp',                        '0.24';
    'WGDp',                        '0.24';

    'aSiN_conc',                   '9.17e18';
    'NTAn',                        '1e22';
    'NTDn',                        '1e22';
    'WTAn',                        '0.03';
    'WTDn',                        '0.05';
    'NGAn',                        '1.76e19';
    'NGDn',                        '1.76e19';
    'EGAn',                        '1.7-0.7';
    'EGDn',                        '0.5';
    'WGAn',                        '0.24';
    'WGDn',                        '0.24';

    'i_aSi_thick',                      '0.00';

    'cSiP_concE',                        '1e21';
    'cSiN_concBSF',                        '1e21';
    'BSF_area',                         '0.9';
     };

InputData.vars = {...
    'emitter_cSi_penetration_depth',	[0.01,0.02,0.03,0.05,0.07,0.1:0.1:1];
    };
%     'cSiP_concE',                        logspace(18,21,10);
%     'cellwidth',                        200:200:1000;
%     'cSiN_concBSF_wafer',                  [1e15,5e15,1e16];

AtlasSimulations(InputData);
%}

%% heterojunctions IBC HOMO-BSF/HET-E (type D)
% with defects (Ns = 1e11 1/cm3)
%{
clear all
InputData.isReplaceMat = true;
InputData.isReplace = true;


InputData.infilename = 'IBC_HIT_full11.in';
InputData.outfilename = 'IBC_HIT_full11D.in';

InputData.outputDir = 'SC_IBC/HYBRID_Ddef1e11_old2/';

InputData.FinalDataFileName = 'IBC-HYBRID_Ddef1e11_old2';

InputData.externvars0 = {...
    'OutputString',                     'tmpHYBRIDDdef2_old';
    'cSi_thick',                        '100';
    'np_x_gap_area',                    '0';
    'DIR',                              '""';

    'DefectiveLayerThickness',          '0.001';
    'isInterfaceDefect',                '0';
    'isBSFdef',                         '1';
    'isemitterdef',                     '1';
	'Ndefsurf_BSF',						'1e11';
	'SigmaDef_BSF'						'1e-16';
	'Ndefsurf_emitter'					'1e11';
	'SigmaDef_emitter'					'1e-16';

    'SweepType',                        '4';

    'aSiP_conc',                   '2.996e19';
    'NTAp',                        '1e22';
    'NTDp',                        '1e22';
    'WTAp',                        '0.03';
    'WTDp',                        '0.08';
    'NGAp',                        '1.76e19';
    'NGDp',                        '1.76e19';
    'EGAp',                        '1.7-1.3';
    'EGDp',                        '1.1';
    'WGAp',                        '0.24';
    'WGDp',                        '0.24';

    'aSiN_conc',                   '9.17e18';
    'NTAn',                        '1e22';
    'NTDn',                        '1e22';
    'WTAn',                        '0.03';
    'WTDn',                        '0.05';
    'NGAn',                        '1.76e19';
    'NGDn',                        '1.76e19';
    'EGAn',                        '1.7-0.7';
    'EGDn',                        '0.5';
    'WGAn',                        '0.24';
    'WGDn',                        '0.24';

    'i_aSi_thick',                      '0.00';

    'cSiP_concE',                        '1e21';
    'emitter_cSi_penetration_depth',    '0.01';
    'cSiN_concBSF',                        '1e21';
    'BSF_cSi_penetration_depth',	    '0.01';
     };

InputData.vars = {...
    'BSF_area',                         0.1:0.1:0.9;
    };
%     'cSiN_concBSF',                        logspace(18,21,10);
%     'cellwidth',                        200:200:1000;
%     'cSiN_concBSF_wafer',                  [1e15,5e15,1e16];

AtlasSimulations(InputData);
%}
% with defects (Ns = 5e11 1/cm3)
%{
clear all
InputData.isReplaceMat = true;
InputData.isReplace = true;


InputData.infilename = 'IBC_HIT_full11.in';
InputData.outfilename = 'IBC_HIT_full11D.in';

InputData.outputDir = 'SC_IBC/HYBRID_Ddef5e11_old2/';

InputData.FinalDataFileName = 'IBC-HYBRID_Ddef5e11_old2';

InputData.externvars0 = {...
    'OutputString',                     'tmpHYBRIDDdef2_old';
    'cSi_thick',                        '100';
    'np_x_gap_area',                    '0';
    'DIR',                              '""';

    'DefectiveLayerThickness',          '0.001';
    'isInterfaceDefect',                '0';
    'isBSFdef',                         '1';
    'isemitterdef',                     '1';
	'Ndefsurf_BSF',						'5e11';
	'SigmaDef_BSF'						'1e-16';
	'Ndefsurf_emitter'					'5e11';
	'SigmaDef_emitter'					'1e-16';

    'SweepType',                        '4';

    'aSiP_conc',                   '2.996e19';
    'NTAp',                        '1e22';
    'NTDp',                        '1e22';
    'WTAp',                        '0.03';
    'WTDp',                        '0.08';
    'NGAp',                        '1.76e19';
    'NGDp',                        '1.76e19';
    'EGAp',                        '1.7-1.3';
    'EGDp',                        '1.1';
    'WGAp',                        '0.24';
    'WGDp',                        '0.24';

    'aSiN_conc',                   '9.17e18';
    'NTAn',                        '1e22';
    'NTDn',                        '1e22';
    'WTAn',                        '0.03';
    'WTDn',                        '0.05';
    'NGAn',                        '1.76e19';
    'NGDn',                        '1.76e19';
    'EGAn',                        '1.7-0.7';
    'EGDn',                        '0.5';
    'WGAn',                        '0.24';
    'WGDn',                        '0.24';

    'i_aSi_thick',                      '0.00';

    'cSiP_concE',                        '1e21';
    'emitter_cSi_penetration_depth',    '0.01';
    'cSiN_concBSF',                        '1e21';
    'BSF_cSi_penetration_depth',	    '0.01';
     };

InputData.vars = {...
    'BSF_area',                         0.1:0.1:0.9;
    };
%     'cSiN_concBSF',                        logspace(18,21,10);
%     'cellwidth',                        200:200:1000;
%     'cSiN_concBSF_wafer',                  [1e15,5e15,1e16];

AtlasSimulations(InputData);
%}
% with defects (Ns = 1e12 1/cm3)
%{
clear all
InputData.isReplaceMat = true;
InputData.isReplace = true;


InputData.infilename = 'IBC_HIT_full11.in';
InputData.outfilename = 'IBC_HIT_full11D.in';

InputData.outputDir = 'SC_IBC/HYBRID_Ddef1e12_old2/';

InputData.FinalDataFileName = 'IBC-HYBRID_Ddef1e12_old2';

InputData.externvars0 = {...
    'OutputString',                     'tmpHYBRIDDdef2_old';
    'cSi_thick',                        '100';
    'np_x_gap_area',                    '0';
    'DIR',                              '""';

    'DefectiveLayerThickness',          '0.001';
    'isInterfaceDefect',                '0';
    'isBSFdef',                         '1';
    'isemitterdef',                     '1';
	'Ndefsurf_BSF',						'1e12';
	'SigmaDef_BSF'						'1e-16';
	'Ndefsurf_emitter'					'1e11';
	'SigmaDef_emitter'					'1e-16';

    'SweepType',                        '4';

    'aSiP_conc',                   '2.996e19';
    'NTAp',                        '1e22';
    'NTDp',                        '1e22';
    'WTAp',                        '0.03';
    'WTDp',                        '0.08';
    'NGAp',                        '1.76e19';
    'NGDp',                        '1.76e19';
    'EGAp',                        '1.7-1.3';
    'EGDp',                        '1.1';
    'WGAp',                        '0.24';
    'WGDp',                        '0.24';

    'aSiN_conc',                   '9.17e18';
    'NTAn',                        '1e22';
    'NTDn',                        '1e22';
    'WTAn',                        '0.03';
    'WTDn',                        '0.05';
    'NGAn',                        '1.76e19';
    'NGDn',                        '1.76e19';
    'EGAn',                        '1.7-0.7';
    'EGDn',                        '0.5';
    'WGAn',                        '0.24';
    'WGDn',                        '0.24';

    'i_aSi_thick',                      '0.00';

    'cSiP_concE',                        '1e21';
    'emitter_cSi_penetration_depth',    '0.01';
    'cSiN_concBSF',                        '1e21';
    'BSF_cSi_penetration_depth',	    '0.01';
     };

InputData.vars = {...
    'BSF_area',                         0.1:0.1:0.9;
    };
%     'cSiN_concBSF',                        logspace(18,21,10);
%     'cellwidth',                        200:200:1000;
%     'cSiN_concBSF_wafer',                  [1e15,5e15,1e16];

AtlasSimulations(InputData);
%}

% without defects
%{
clear all
InputData.isReplaceMat = true;
InputData.isReplace = true;


InputData.infilename = 'IBC_HIT_full12.in';
InputData.outfilename = 'IBC_HIT_full12_tD.in';

InputData.outputDir = 'SC_IBC/HYBRID_D_old_2/';

InputData.FinalDataFileName = 'IBC-HYBRID_D_old_2';

InputData.externvars0 = {...
    'OutputString',                     'tmpHYBRID_D_old_2';
    'cSi_thick',                        '100';
    'np_x_gap_area',                    '0';
    'DIR',                              '""';

    'DefectiveLayerThickness',          '0.00';
    'isInterfaceDefect',                '0';
    'isBSFdef',                         '0';
    'isemitterdef',                     '0';

    'SweepType',                        '4';

    'aSiP_conc',                   '2.996e19';
    'NTAp',                        '1e22';
    'NTDp',                        '1e22';
    'WTAp',                        '0.03';
    'WTDp',                        '0.08';
    'NGAp',                        '1.76e19';
    'NGDp',                        '1.76e19';
    'EGAp',                        '1.7-1.3';
    'EGDp',                        '1.1';
    'WGAp',                        '0.24';
    'WGDp',                        '0.24';

    'aSiN_conc',                   '9.17e18';
    'NTAn',                        '1e22';
    'NTDn',                        '1e22';
    'WTAn',                        '0.03';
    'WTDn',                        '0.05';
    'NGAn',                        '1.76e19';
    'NGDn',                        '1.76e19';
    'EGAn',                        '1.7-0.7';
    'EGDn',                        '0.5';
    'WGAn',                        '0.24';
    'WGDn',                        '0.24';

    'i_aSi_thick',                      '0.00';

    'cSiP_concE',                        '1e21';
    'cSiN_concBSF',                        '1e21';
     };

InputData.vars = {...
    'BSF_area',                         0.1:0.1:0.9;
    'cSiN_concBSF',                        logspace(18,21,10);
    };
%     'cellwidth',                        200:200:1000;
%     'cSiN_concBSF_wafer',                  [1e15,5e15,1e16];

AtlasSimulations(InputData);
%}
% various penetration depth of dopant in homojunction (HET with defects Ns = 5e11 1/cm3)
%{
clear all
InputData.isReplaceMat = true;
InputData.isReplace = true;


InputData.infilename = 'IBC_HIT_full12.in';
InputData.outfilename = 'IBC_HIT_full12D.in';

InputData.outputDir = 'SC_IBC/HYBRID_D_old_p/';

InputData.FinalDataFileName = 'IBC-HYBRID_D_old_p';

InputData.externvars0 = {...
    'OutputString',                     'tmpHYBRIDD_old';
    'cSi_thick',                        '100';
    'np_x_gap_area',                    '0';
    'DIR',                              '""';

    'DefectiveLayerThickness',          '0.001';
    'isInterfaceDefect',                '0';
    'isBSFdef',                         '1';
    'isemitterdef',                     '1';
	'Ndefsurf_BSF',						'5e11';
	'SigmaDef_BSF'						'1e-16';
	'Ndefsurf_emitter'					'5e11';
	'SigmaDef_emitter'					'1e-16';

    'SweepType',                        '4';

    'aSiP_conc',                   '2.996e19';
    'NTAp',                        '1e22';
    'NTDp',                        '1e22';
    'WTAp',                        '0.03';
    'WTDp',                        '0.08';
    'NGAp',                        '1.76e19';
    'NGDp',                        '1.76e19';
    'EGAp',                        '1.7-1.3';
    'EGDp',                        '1.1';
    'WGAp',                        '0.24';
    'WGDp',                        '0.24';

    'aSiN_conc',                   '9.17e18';
    'NTAn',                        '1e22';
    'NTDn',                        '1e22';
    'WTAn',                        '0.03';
    'WTDn',                        '0.05';
    'NGAn',                        '1.76e19';
    'NGDn',                        '1.76e19';
    'EGAn',                        '1.7-0.7';
    'EGDn',                        '0.5';
    'WGAn',                        '0.24';
    'WGDn',                        '0.24';

    'i_aSi_thick',                      '0.00';

    'cSiP_concE',                        '1e21';
    'cSiN_concBSF',                        '1e21';
    'BSF_area',                         '0.1';
     };

InputData.vars = {...
    'BSF_cSi_penetration_depth',        [0.01,0.02,0.03,0.05,0.07,0.1:0.1:1];
    };
%     'cSiN_concBSF',                        logspace(18,21,10);
%     'cellwidth',                        200:200:1000;
%     'cSiN_concBSF_wafer',                  [1e15,5e15,1e16];

AtlasSimulations(InputData);
%}



%% 1D heterojunctions+homojunctions
%24.10.15
% %{
% +intrinsic aSi:H layer
clear all
InputData.isReplaceMat = true;
InputData.isReplace = true;
    InputData.UNIXserverName = 'silvacox2';
    InputData.UNIXserverName = 'silvacox';
    InputData.ATLASparameter = ' -V 5.19.20.R';
    InputData.ATLASparameter = ' -V 5.20.2.R';

% defect pool model params
InputData.DPM.regnum = 1;
InputData.DPM.T = 300; % operating (room) temperature
InputData.DPM.T0 = 480; % production temperature
InputData.DPM.Eg = 1.7;
InputData.DPM.sigma = 0.19;
InputData.DPM.U = 0.2; % the correlation energy, which is needed to place two electrons on the same defect
InputData.DPM.H = 5e21; % 1/cm3
InputData.DPM.NSiSi = 2e23; % 1/cm3
InputData.DPM.Ep = 1.27;
InputData.DPM.dy = 1e-3; 

InputData.isTony = true;
InputData.isSaveData = true;
InputData.isReadOnly = true;

isPNN = 0;
% isPNN = 1;
% SigmaDef = '1e-17';
SigmaDef = '1e-16';

InputData.infilename  = '1D_HET_defectiveL.in';
InputData.outfilename = '1D_HET_defectiveL_a.in';

% name = ['1D_HET_defectiveL_iaSi_s',SigmaDef,'_'];
name = ['1D[25.10.15]_HETSiO2_defectiveL_iaSi_s',SigmaDef,'_'];
if isPNN
    InputData.outputDir = [name,'pnn/'];
    InputData.FinalDataFileName = [name,'pnn'];
    InputData.originname = [name,'pnn'];
else
    InputData.outputDir = [name,'nnp/'];
    InputData.FinalDataFileName = [name,'nnp'];
    InputData.originname = [name,'nnp'];
end

InputData.externvars0 = {...
    'OutputString',                     'tmpa';
    'isPNN',                            num2str(isPNN);
    'cSi_thick',                        '100';
    'IsStructureShow',                  '1';
	'isInterfaceDefect',                '0';
    'DIR',                              '""';
    'SigmaDef',                         SigmaDef;
    
    % fitted with DPM
    % F-Ev = 0.299984eV
    'aSiP_conc',                   '2.708e19';
    'NTAp',                        '2e21';
    'NTDp',                        '2e21';
    'WTAp',                        '0.03';
    'WTDp',                        '0.085';
    'NGAp',                        '4.514e+19';
    'NGDp',                        '4.514e+19';
    'EGAp',                        '1.7-1.4879';
    'EGDp',                        '1.2521';
    'WGAp',                        '0.2687';
    'WGDp',                        '0.2687';

    % F-Ev = 1.50005eV
    'aSiN_conc',                   '6.81e17';
    'NTAn',                        '2e21';
    'NTDn',                        '2e21';
    'WTAn',                        '0.03';
    'WTDn',                        '0.065';
    'NGAn',                        '6.5898e+17';
    'NGDn',                        '6.5898e+17';
    'EGAn',                        '1.7-0.9325';
    'EGDn',                        '0.6967';
    'WGAn',                        '0.2687';
    'WGDn',                        '0.2687';

    % F-Ev = 1.09516eV
    'aSiNi_concEmitter',            '1e15';
    'aSiPi_concEmitter',            '0';
    'NTAip',                        '2e21';
    'NTDip',                        '2e21';
    'WTAip',                        '0.03';
    'WTDip',                        '0.045';
    'NGAip',                        '5.0943e+16';
    'NGDip',                        '5.0943e+16';
    'EGAip',                        '1.7-1.2116';
    'EGDip',                        '0.9758';
    'WGAip',                        '0.5530';
    'WGDip',                        '0.5530';

    'aSiNi_concBSF',                '1e15';
    'aSiPi_concBSF',                '0';
    'NTAin',                        '2e21';
    'NTDin',                        '2e21';
    'WTAin',                        '0.03';
    'WTDin',                        '0.045';
    'NGAin',                        '5.0943e+16';
    'NGDin',                        '5.0943e+16';
    'EGAin',                        '1.7-1.2116';
    'EGDin',                        '0.9758';
    'WGAin',                        '0.5530';
    'WGDin',                        '0.5530';

    'emitter_cSi_penetration_depth',    '0.1';
    'BSF_cSi_penetration_depth',	    '0.1';
     
    'SweepType',                        '2';
    'DefectiveLayerThicknessFront',     '0';
    'DefectiveLayerThickness',          '0';
    'i_aSi_thick_front',                '0.00';
    'i_aSi_thick',                      '0.010';
    };

InputData.vars = {...
    'Ndefsurf',                         5e11;
    };
%     'Ndefsurf',                         0*logspace(10,12,3+2*4);
%     'DefectiveLayerThicknessFront',     0*[0,1e-3];
%     'DefectiveLayerThickness',          0*[0,1e-3];
%     'i_aSi_thick_front',                0*[0,5e-3];
%     'i_aSi_thick',                      0*[0,5e-3];
%     'SweepType',                        [1:7];
%     'emitterSiO2_thick',                [0,10e-3];
%     'BSFSiO2_thick',                    [0,10e-3];

AtlasSimulations(InputData);
%}

%{
% +intrinsic aSi:H layer
clear all
InputData.isReplaceMat = true;
InputData.isReplace = true;

InputData.isTony = true;
InputData.isSaveData = true;
InputData.isReadOnly = true;

% isPNN = 0;
isPNN = 1;
% SigmaDef = '1e-17';
SigmaDef = '1e-16';


InputData.infilename = '1D_HET_defectiveL.in';
InputData.outfilename = '1D_HET_defectiveL_b.in';

% name = ['1D_HET_defectiveL_iaSi_s',SigmaDef,'_'];
name = ['1D_HETSiO2_defectiveL_iaSi_s',SigmaDef,'_'];
if isPNN
    InputData.outputDir = [name,'pnn/'];
    InputData.FinalDataFileName = [name,'pnn'];
    InputData.originname = [name,'pnn'];
else
    InputData.outputDir = [name,'nnp/'];
    InputData.FinalDataFileName = [name,'nnp'];
    InputData.originname = [name,'nnp'];
end

InputData.externvars0 = {...
    'OutputString',                     'tmpb';
    'isPNN',                            num2str(isPNN);
    'cSi_thick',                        '100';
    'IsStructureShow',                  '1';
	'isInterfaceDefect',                '0';
    'DIR',                              '""';
    'SigmaDef',                         SigmaDef;
    };
%     'DefectiveLayerThicknessFront',     '0';
%     'DefectiveLayerThickness',          '0';
%     'i_aSi_thick_front',                '0.00';
%     'i_aSi_thick',                      '0.010';

InputData.vars = {...
    'Ndefsurf',                         logspace(10,12,31);
    'DefectiveLayerThicknessFront',     [0,1e-3];
    'DefectiveLayerThickness',          [0,1e-3];
	'i_aSi_thick_front',                [0,1e-2];
	'i_aSi_thick',                      [0,1e-2];
    };
%     'emitterSiO2_thick',                [0,10e-3];
%     'BSFSiO2_thick',                    [0,10e-3];

AtlasSimulations(InputData);
%}

%{
clear all
InputData.isReplaceMat = true;
InputData.isReplace = true;

isPNN = 0;
isPNN = 1;
% SigmaDef = '1e-17';
SigmaDef = '1e-16';

InputData.isTony = true;
InputData.isSaveData = true;
InputData.isReadOnly = true;

InputData.infilename = '1D_HET_defectiveL.in';
InputData.outfilename = '1D_HET_defectiveLa.in';


% name = ['1D_HET_defectiveL_s',SigmaDef,'_'];
name = ['1D_HETSiO2_defectiveL_s',SigmaDef,'_'];
if isPNN
    InputData.outputDir = [name,'pnn/'];
    InputData.FinalDataFileName = [name,'pnn'];
    InputData.originname = [name,'pnn'];
else
    InputData.outputDir = [name,'nnp/'];
    InputData.FinalDataFileName = [name,'nnp'];
    InputData.originname = [name,'nnp'];
end

InputData.externvars0 = {...
    'OutputString',                     'tmp';
    'isPNN',                            num2str(isPNN);
    'cSi_thick',                        '100';
    'IsStructureShow',                  '1';
	'isInterfaceDefect',                '0';
    'DIR',                              '""';
    'SigmaDef',                         SigmaDef;
	'i_aSi_thick_front',                '0.00';
	'i_aSi_thick',                      '0.00';
    'DefectiveLayerThicknessFront',     '0';
    'DefectiveLayerThickness',          '0';
    };

InputData.vars = {...
    'emitterSiO2_thick',                [0,1e-3];
    'BSFSiO2_thick',                    [0,1e-3];
    };
%     'Ndefsurf',                         logspace(10,12,31);
%     'DefectiveLayerThicknessFront',     [0,1e-3];
%     'DefectiveLayerThickness',          [0,1e-3];

AtlasSimulations(InputData);
%}


