function sessList = SessionList()
% Description 
%Enter the list of session descriptions

sessList = {};
dataDir = 'D:\Data\4jason\jason\cells_reclassified_jason\';


%------------------------------------------------------------
s.animal = '6318';
s.name = '06';
s.group = 'CPP';
s.split = 1;
s.tfiledir = [dataDir 'CPP\6318_06\tfiles'];
s.tfileglob = '';
tfileprefix = '6318_9_6_tt';
s.sameTTfmt = [tfileprefix '%%'];
s.epochs.sleep1 = [tfileprefix '*_tsleep1.*'];
s.epochs.maze1 = [tfileprefix '*_tmaze1.*'];
s.epochs.sleep2 = [tfileprefix '*_tsleep2.*'];
s.other = '';
sessList{end+1} = session(s);

%------------------------------------------------------------
s.animal = '6318';
s.name = '16';
s.group = 'CPP';
s.split = 1;
s.tfiledir = [dataDir 'CPP\6318_16\tfiles'];
s.tfileglob = '';
tfileprefix = '6318_9_16_tt';
s.sameTTfmt = [tfileprefix '%%'];
s.epochs.sleep1 = [tfileprefix '*_tsleep1.*'];
s.epochs.maze1 = [tfileprefix '*_tmaze1.*'];
s.epochs.sleep2 = [tfileprefix '*_tsleep2.*'];
s.other = '';
sessList{end+1} = session(s);

%------------------------------------------------------------
s.animal = '6401';
s.name = '05';
s.group = 'CPP';
s.split = 1;
s.tfiledir = [dataDir 'CPP\6401_05\tfiles'];
s.tfileglob = '';
tfileprefix = '6401_12_5_tt';
s.sameTTfmt = [tfileprefix '%%'];
s.epochs.sleep1 = [tfileprefix '*_tsleep1.*'];
s.epochs.maze1 = [tfileprefix '*_tmaze1.*'];
s.epochs.sleep2 = [tfileprefix '*_tsleep2.*'];
s.other = '';
sessList{end+1} = session(s);

%------------------------------------------------------------
s.animal = '6401';
s.name = '08';
s.group = 'CPP';
s.split = 1;
s.tfiledir = [dataDir 'CPP\6401_08\tfiles'];
s.tfileglob = '';
tfileprefix = '6401_Dec_8_tt';
s.sameTTfmt = [tfileprefix '%%'];
s.epochs.sleep1 = [tfileprefix '*_tsleep1.*'];
s.epochs.maze1 = [tfileprefix '*_tmaze1.*'];
s.epochs.sleep2 = [tfileprefix '*_tsleep2.*'];
s.other = '';
sessList{end+1} = session(s);

%------------------------------------------------------------
s.animal = '6401';
s.name = '10';
s.group = 'CPP';
s.split = 1;
s.tfiledir = [dataDir 'CPP\6401_10\tfiles'];
s.tfileglob = '';
tfileprefix = '6401_12_10_tt';
s.sameTTfmt = [tfileprefix '%%'];
s.epochs.sleep1 = [tfileprefix '*_tsleep1.*'];
s.epochs.maze1 = [tfileprefix '*_tmaze1.*'];
s.epochs.sleep2 = [tfileprefix '*_tsleep2.*'];
s.other = '';
sessList{end+1} = session(s);

%------------------------------------------------------------
s.animal = '6481';
s.name = '02';
s.group = 'CPP';
s.split = 1;
s.tfiledir = [dataDir 'CPP\6481_02\tfiles'];
s.tfileglob = '';
tfileprefix = '6481_02_tt';
s.sameTTfmt = [tfileprefix '%%'];
s.epochs.sleep1 = [tfileprefix '*_tsleep1.*'];
s.epochs.maze1 = [tfileprefix '*_tmaze1.*'];
s.epochs.sleep2 = [tfileprefix '*_tsleep2.*'];
s.other = '';
sessList{end+1} = session(s);

%------------------------------------------------------------
s.animal = '6495';
s.name = '16';
s.group = 'CPP';
s.split = 1;
s.tfiledir = [dataDir 'CPP\6495_16\tfiles'];
s.tfileglob = '';
tfileprefix = '6495_12_16_tt';
s.sameTTfmt = [tfileprefix '%%'];
s.epochs.sleep1 = [tfileprefix '*_tsleep1.*'];
s.epochs.maze1 = [tfileprefix '*_tmaze1.*'];
s.epochs.sleep2 = [tfileprefix '*_tsleep2.*'];
s.other = '';
sessList{end+1} = session(s);

%------------------------------------------------------------
s.animal = '6495';
s.name = '19';
s.group = 'CPP';
s.split = 1;
s.tfiledir = [dataDir 'CPP\6495_19\tfiles'];
s.tfileglob = '';
tfileprefix = '6495_12_19_tt';
s.sameTTfmt = [tfileprefix '%%'];
s.epochs.sleep1 = [tfileprefix '*_tsleep1.*'];
s.epochs.maze1 = [tfileprefix '*_tmaze1.*'];
s.epochs.sleep2 = [tfileprefix '*_tsleep2.*'];
s.other = '';
sessList{end+1} = session(s);

%------------------------------------------------------------
s.animal = '6495';
s.name = '22';
s.group = 'CPP';
s.split = 1;
s.tfiledir = [dataDir 'CPP\6495_22\tfiles'];
s.tfileglob = '';
tfileprefix = '6495_12_22_tt';
s.sameTTfmt = [tfileprefix '%%'];
s.epochs.sleep1 = [tfileprefix '*_tsleep1.*'];
s.epochs.maze1 = [tfileprefix '*_tmaze1.*'];
s.epochs.sleep2 = [tfileprefix '*_tsleep2.*'];
s.other = '';
sessList{end+1} = session(s);

%------------------------------------------------------------
s.animal = '6496';
s.name = '02';
s.group = 'CPP';
s.split = 1;
s.tfiledir = [dataDir 'CPP\6496_02\tfiles'];
s.tfileglob = '';
tfileprefix = '6496_02_tt';
s.sameTTfmt = [tfileprefix '%%'];
s.epochs.sleep1 = [tfileprefix '*_tsleep1.*'];
s.epochs.maze1 = [tfileprefix '*_tmaze1.*'];
s.epochs.sleep2 = [tfileprefix '*_tsleep2.*'];
s.other = '';
sessList{end+1} = session(s);

%------------------------------------------------------------
s.animal = '6496';
s.name = '05';
s.group = 'CPP';
s.split = 1;
s.tfiledir = [dataDir 'CPP\6496_05\tfiles'];
s.tfileglob = '';
tfileprefix = '6496_05_TET_';
s.sameTTfmt = [tfileprefix '%%'];
s.epochs.sleep1 = [tfileprefix '*_tsleep1.*'];
s.epochs.maze1 = [tfileprefix '*_tmaze1.*'];
s.epochs.sleep2 = [tfileprefix '*_tsleep2.*'];
s.other = '';
sessList{end+1} = session(s);

%------------------------------------------------------------
s.animal = '6496';
s.name = '08';
s.group = 'CPP';
s.split = 1;
s.tfiledir = [dataDir 'CPP\6496_08\tfiles'];
s.tfileglob = '';
tfileprefix = '6496_08_TET_';
s.sameTTfmt = [tfileprefix '%%'];
s.epochs.sleep1 = [tfileprefix '*_tsleep1.*'];
s.epochs.maze1 = [tfileprefix '*_tmaze1.*'];
s.epochs.sleep2 = [tfileprefix '*_tsleep2.*'];
s.other = '';
sessList{end+1} = session(s);

%------------------------------------------------------------
s.animal = '6496';
s.name = '13';
s.group = 'CPP';
s.split = 1;
s.tfiledir = [dataDir 'CPP\6496_13\tfiles'];
s.tfileglob = '';
tfileprefix = '6496_13_TET_';
s.sameTTfmt = [tfileprefix '%%'];
s.epochs.sleep1 = [tfileprefix '*_tsleep1.*'];
s.epochs.maze1 = [tfileprefix '*_tmaze1.*'];
s.epochs.sleep2 = [tfileprefix '*_tsleep2.*'];
s.other = '';
sessList{end+1} = session(s);

%------------------------------------------------------------
s.animal = '6591';
s.name = '02';
s.group = 'CPP';
s.split = 1;
s.tfiledir = [dataDir 'CPP\6591_02\tfiles'];
s.tfileglob = '';
tfileprefix = '6591_02_TET_';
s.sameTTfmt = [tfileprefix '%%'];
s.epochs.sleep1 = [tfileprefix '*_tsleep1.*'];
s.epochs.maze1 = [tfileprefix '*_tmaze1.*'];
s.epochs.sleep2 = [tfileprefix '*_tsleep2.*'];
s.other = '';
sessList{end+1} = session(s);

%------------------------------------------------------------
s.animal = '6591';
s.name = '04';
s.group = 'CPP';
s.split = 1;
s.tfiledir = [dataDir 'CPP\6591_04\tfiles'];
s.tfileglob = '';
tfileprefix = '6591_04_TET_';
s.sameTTfmt = [tfileprefix '%%'];
s.epochs.sleep1 = [tfileprefix '*_tsleep1.*'];
s.epochs.maze1 = [tfileprefix '*_tmaze1.*'];
s.epochs.sleep2 = [tfileprefix '*_tsleep2.*'];
s.other = '';
sessList{end+1} = session(s);

%------------------------------------------------------------
s.animal = '6568';
s.name = '02';
s.group = 'CPP';
s.split = 1;
s.tfiledir = [dataDir 'CPP\6568_02\tfiles'];
s.tfileglob = '';
tfileprefix = '6568_02_tt';
s.sameTTfmt = [tfileprefix '%%'];
s.epochs.sleep1 = [tfileprefix '*_tsleep1.*'];
s.epochs.maze1 = [tfileprefix '*_tmaze1.*'];
s.epochs.sleep2 = [tfileprefix '*_tsleep2.*'];
s.other = '';
sessList{end+1} = session(s);

%------------------------------------------------------------
s.animal = '6568';
s.name = '04';
s.group = 'CPP';
s.split = 1;
s.tfiledir = [dataDir 'CPP\6568_04\tfiles'];
s.tfileglob = '';
tfileprefix = '6568_04_tt';
s.sameTTfmt = [tfileprefix '%%'];
s.epochs.sleep1 = [tfileprefix '*_tsleep1.*'];
s.epochs.maze1 = [tfileprefix '*_tmaze1.*'];
s.epochs.sleep2 = [tfileprefix '*_tsleep2.*'];
s.other = '';
sessList{end+1} = session(s);

%------------------------------------------------------------
s.animal = '6568';
s.name = '07';
s.group = 'CPP';
s.split = 1;
s.tfiledir = [dataDir 'CPP\6568_07\tfiles'];
s.tfileglob = '';
tfileprefix = '6568_07_tt';
s.sameTTfmt = [tfileprefix '%%'];
s.epochs.sleep1 = [tfileprefix '*_tsleep1.*'];
s.epochs.maze1 = [tfileprefix '*_tmaze1.*'];
s.epochs.sleep2 = [tfileprefix '*_tsleep2.*'];
s.other = '';
sessList{end+1} = session(s);
%------------------------------------------------------------
s.animal = '6318';
s.name = '23';
s.group = 'Saline';
s.split = 1;
s.tfiledir = [dataDir 'Saline\6318_23\tfiles'];
s.tfileglob = '';
tfileprefix = '6318_9_23_tt';
s.sameTTfmt = [tfileprefix '%%'];
s.epochs.sleep1 = [tfileprefix '*_tsleep1.*'];
s.epochs.maze1 = [tfileprefix '*_tmaze1.*'];
s.epochs.sleep2 = [tfileprefix '*_tsleep2.*'];
s.other = '';
sessList{end+1} = session(s);

%------------------------------------------------------------
s.animal = '6318';
s.name = '28';
s.group = 'Saline';
s.split = 1;
s.tfiledir = [dataDir 'Saline\6318_28\tfiles'];
s.tfileglob = '';
tfileprefix = '6318_9_28_tt';
s.sameTTfmt = [tfileprefix '%%'];
s.epochs.sleep1 = [tfileprefix '*_tsleep1.*'];
s.epochs.maze1 = [tfileprefix '*_tmaze1.*'];
s.epochs.sleep2 = [tfileprefix '*_tsleep2.*'];
s.other = '';
sessList{end+1} = session(s);

%------------------------------------------------------------
s.animal = '6401';
s.name = '04';
s.group = 'Saline';
s.split = 1;
s.tfiledir = [dataDir 'Saline\6401_04\tfiles'];
s.tfileglob = '';
tfileprefix = '6401_12_4_tt';
s.sameTTfmt = [tfileprefix '%%'];
s.epochs.sleep1 = [tfileprefix '*_tsleep1.*'];
s.epochs.maze1 = [tfileprefix '*_tmaze1.*'];
s.epochs.sleep2 = [tfileprefix '*_tsleep2.*'];
s.other = '';
sessList{end+1} = session(s);

%------------------------------------------------------------
s.animal = '6401';
s.name = '07';
s.group = 'Saline';
s.split = 1;
s.tfiledir = [dataDir 'Saline\6401_07\tfiles'];
s.tfileglob = '';
tfileprefix = '6401_12_7_TET_';
s.sameTTfmt = [tfileprefix '%%'];
s.epochs.sleep1 = [tfileprefix '*_tsleep1.*'];
s.epochs.maze1 = [tfileprefix '*_tmaze1.*'];
s.epochs.sleep2 = [tfileprefix '*_tsleep2.*'];
s.other = '';
sessList{end+1} = session(s);

%------------------------------------------------------------
s.animal = '6401';
s.name = '09';
s.group = 'Saline';
s.split = 1;
s.tfiledir = [dataDir 'Saline\6401_09\tfiles'];
s.tfileglob = '';
tfileprefix = '6401_12_9_TET_';
s.sameTTfmt = [tfileprefix '%%'];
s.epochs.sleep1 = [tfileprefix '*_tsleep1.*'];
s.epochs.maze1 = [tfileprefix '*_tmaze1.*'];
s.epochs.sleep2 = [tfileprefix '*_tsleep2.*'];
s.other = '';
sessList{end+1} = session(s);

%------------------------------------------------------------
s.animal = '6481';
s.name = '03';
s.group = 'Saline';
s.split = 1;
s.tfiledir = [dataDir 'Saline\6481_03\tfiles'];
s.tfileglob = '';
tfileprefix = '6481_03_tt';
s.sameTTfmt = [tfileprefix '%%'];
s.epochs.sleep1 = [tfileprefix '*_tsleep1.*'];
s.epochs.maze1 = [tfileprefix '*_tmaze1.*'];
s.epochs.sleep2 = [tfileprefix '*_tsleep2.*'];
s.other = '';
sessList{end+1} = session(s);

%------------------------------------------------------------
s.animal = '6495';
s.name = '15';
s.group = 'Saline';
s.split = 1;
s.tfiledir = [dataDir 'Saline\6495_15\tfiles'];
s.tfileglob = '';
tfileprefix = '6495_12_15_TET_';
s.sameTTfmt = [tfileprefix '%%'];
s.epochs.sleep1 = [tfileprefix '*_tsleep1.*'];
s.epochs.maze1 = [tfileprefix '*_tmaze1.*'];
s.epochs.sleep2 = [tfileprefix '*_tsleep2.*'];
s.other = '';
sessList{end+1} = session(s);

%------------------------------------------------------------
s.animal = '6495';
s.name = '17';
s.group = 'Saline';
s.split = 1;
s.tfiledir = [dataDir 'Saline\6495_17\tfiles'];
s.tfileglob = '';
tfileprefix = '6495_12_17_TET_';
s.sameTTfmt = [tfileprefix '%%'];
s.epochs.sleep1 = [tfileprefix '*_tsleep1.*'];
s.epochs.maze1 = [tfileprefix '*_tmaze1.*'];
s.epochs.sleep2 = [tfileprefix '*_tsleep2.*'];
s.other = '';
sessList{end+1} = session(s);

%------------------------------------------------------------
s.animal = '6495';
s.name = '18';
s.group = 'Saline';
s.split = 1;
s.tfiledir = [dataDir 'Saline\6495_18\tfiles'];
s.tfileglob = '';
tfileprefix = '6495_12_18_TET_';
s.sameTTfmt = [tfileprefix '%%'];
s.epochs.sleep1 = [tfileprefix '*_tsleep1.*'];
s.epochs.maze1 = [tfileprefix '*_tmaze1.*'];
s.epochs.sleep2 = [tfileprefix '*_tsleep2.*'];
s.other = '';
sessList{end+1} = session(s);

%------------------------------------------------------------
s.animal = '6495';
s.name = '21';
s.group = 'Saline';
s.split = 1;
s.tfiledir = [dataDir 'Saline\6495_21\tfiles'];
s.tfileglob = '';
tfileprefix = '6495_12_21_TET_';
s.sameTTfmt = [tfileprefix '%%'];
s.epochs.sleep1 = [tfileprefix '*_tsleep1.*'];
s.epochs.maze1 = [tfileprefix '*_tmaze1.*'];
s.epochs.sleep2 = [tfileprefix '*_tsleep2.*'];
s.other = '';
sessList{end+1} = session(s);

%------------------------------------------------------------
s.animal = '6495';
s.name = '23';
s.group = 'Saline';
s.split = 1;
s.tfiledir = [dataDir 'Saline\6495_23\tfiles'];
s.tfileglob = '';
tfileprefix = '6495_12_23_tt';
s.sameTTfmt = [tfileprefix '%%'];
s.epochs.sleep1 = [tfileprefix '*_tsleep1.*'];
s.epochs.maze1 = [tfileprefix '*_tmaze1.*'];
s.epochs.sleep2 = [tfileprefix '*_tsleep2.*'];
s.other = '';
sessList{end+1} = session(s);

%------------------------------------------------------------
s.animal = '6496';
s.name = '01';
s.group = 'Saline';
s.split = 1;
s.tfiledir = [dataDir 'Saline\6496_01\tfiles'];
s.tfileglob = '';
tfileprefix = '6496_01_tt';
s.sameTTfmt = [tfileprefix '%%'];
s.epochs.sleep1 = [tfileprefix '*_tsleep1.*'];
s.epochs.maze1 = [tfileprefix '*_tmaze1.*'];
s.epochs.sleep2 = [tfileprefix '*_tsleep2.*'];
s.other = '';
sessList{end+1} = session(s);

%------------------------------------------------------------
s.animal = '6496';
s.name = '04';
s.group = 'Saline';
s.split = 1;
s.tfiledir = [dataDir 'Saline\6496_04\tfiles'];
s.tfileglob = '';
tfileprefix = '6496_04_TET_';
s.sameTTfmt = [tfileprefix '%%'];
s.epochs.sleep1 = [tfileprefix '*_tsleep1.*'];
s.epochs.maze1 = [tfileprefix '*_tmaze1.*'];
s.epochs.sleep2 = [tfileprefix '*_tsleep2.*'];
s.other = '';
sessList{end+1} = session(s);

%------------------------------------------------------------
s.animal = '6496';
s.name = '06';
s.group = 'Saline';
s.split = 1;
s.tfiledir = [dataDir 'Saline\6496_06\tfiles'];
s.tfileglob = '';
tfileprefix = '6496_06_TET_';
s.sameTTfmt = [tfileprefix '%%'];
s.epochs.sleep1 = [tfileprefix '*_tsleep1.*'];
s.epochs.maze1 = [tfileprefix '*_tmaze1.*'];
s.epochs.sleep2 = [tfileprefix '*_tsleep2.*'];
s.other = '';
sessList{end+1} = session(s);

%------------------------------------------------------------
s.animal = '6496';
s.name = '07';
s.group = 'Saline';
s.split = 1;
s.tfiledir = [dataDir 'Saline\6496_07\tfiles'];
s.tfileglob = '';
tfileprefix = '6496_07_TET_';
s.sameTTfmt = [tfileprefix '%%'];
s.epochs.sleep1 = [tfileprefix '*_tsleep1.*'];
s.epochs.maze1 = [tfileprefix '*_tmaze1.*'];
s.epochs.sleep2 = [tfileprefix '*_tsleep2.*'];
s.other = '';
sessList{end+1} = session(s);

%------------------------------------------------------------
s.animal = '6496';
s.name = '09';
s.group = 'Saline';
s.split = 1;
s.tfiledir = [dataDir 'Saline\6496_09\tfiles'];
s.tfileglob = '';
tfileprefix = '6496_09_TET_';
s.sameTTfmt = [tfileprefix '%%'];
s.epochs.sleep1 = [tfileprefix '*_tsleep1.*'];
s.epochs.maze1 = [tfileprefix '*_tmaze1.*'];
s.epochs.sleep2 = [tfileprefix '*_tsleep2.*'];
s.other = '';
sessList{end+1} = session(s);

%------------------------------------------------------------
s.animal = '6496';
s.name = '12';
s.group = 'Saline';
s.split = 1;
s.tfiledir = [dataDir 'Saline\6496_12\tfiles'];
s.tfileglob = '';
tfileprefix = '6496_12_TET_';
s.sameTTfmt = [tfileprefix '%%'];
s.epochs.sleep1 = [tfileprefix '*_tsleep1.*'];
s.epochs.maze1 = [tfileprefix '*_tmaze1.*'];
s.epochs.sleep2 = [tfileprefix '*_tsleep2.*'];
s.other = '';
sessList{end+1} = session(s);

%------------------------------------------------------------
s.animal = '6591';
s.name = '01';
s.group = 'Saline';
s.split = 1;
s.tfiledir = [dataDir 'Saline\6591_01\tfiles'];
s.tfileglob = '';
tfileprefix = '6591_01_TET_';
s.sameTTfmt = [tfileprefix '%%'];
s.epochs.sleep1 = [tfileprefix '*_tsleep1.*'];
s.epochs.maze1 = [tfileprefix '*_tmaze1.*'];
s.epochs.sleep2 = [tfileprefix '*_tsleep2.*'];
s.other = '';
sessList{end+1} = session(s);

%------------------------------------------------------------
s.animal = '6591';
s.name = '03';
s.group = 'Saline';
s.split = 1;
s.tfiledir = [dataDir 'Saline\6591_03\tfiles'];
s.tfileglob = '';
tfileprefix = '6591_03_TET_';
s.sameTTfmt = [tfileprefix '%%'];
s.epochs.sleep1 = [tfileprefix '*_tsleep1.*'];
s.epochs.maze1 = [tfileprefix '*_tmaze1.*'];
s.epochs.sleep2 = [tfileprefix '*_tsleep2.*'];
s.other = '';
sessList{end+1} = session(s);

%------------------------------------------------------------
s.animal = '6568';
s.name = '01';
s.group = 'Saline';
s.split = 1;
s.tfiledir = [dataDir 'Saline\6568_01\tfiles'];
s.tfileglob = '';
tfileprefix = '6568_01_tt';
s.sameTTfmt = [tfileprefix '%%'];
s.epochs.sleep1 = [tfileprefix '*_tsleep1.*'];
s.epochs.maze1 = [tfileprefix '*_tmaze1.*'];
s.epochs.sleep2 = [tfileprefix '*_tsleep2.*'];
s.other = '';
sessList{end+1} = session(s);

%------------------------------------------------------------
s.animal = '6568';
s.name = '05';
s.group = 'Saline';
s.split = 1;
s.tfiledir = [dataDir 'Saline\6568_05\tfiles'];
s.tfileglob = '';
tfileprefix = '6568_05_tt';
s.sameTTfmt = [tfileprefix '%%'];
s.epochs.sleep1 = [tfileprefix '*_tsleep1.*'];
s.epochs.maze1 = [tfileprefix '*_tmaze1.*'];
s.epochs.sleep2 = [tfileprefix '*_tsleep2.*'];
s.other = '';
sessList{end+1} = session(s);

%------------------------------------------------------------
s.animal = '6568';
s.name = '08';
s.group = 'Saline';
s.split = 1;
s.tfiledir = [dataDir 'Saline\6568_08\tfiles'];
s.tfileglob = '';
tfileprefix = '6568_08_tt';
s.sameTTfmt = [tfileprefix '%%'];
s.epochs.sleep1 = [tfileprefix '*_tsleep1.*'];
s.epochs.maze1 = [tfileprefix '*_tmaze1.*'];
s.epochs.sleep2 = [tfileprefix '*_tsleep2.*'];
s.other = '';
sessList{end+1} = session(s);

%------------------------------------------------------------