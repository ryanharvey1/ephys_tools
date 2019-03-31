function sessList = SessionList()
% Description 
%Enter the list of session descriptions

sessList = {};
dataDir = 'T:\CPP1\JessicasData\';                  % don't forget the '\' at the end

%------------------------------------------------------------
s.animal = '7508';
s.name = '01';
s.group = 'novelCPP';
s.split = 0;
s.tfiledir = [dataDir 'CPP_Data\7508\7508_01\tfiles'];
s.tfileglob = 'TT*.t';
tfileprefix = 'TT_';
s.sameTTfmt = [tfileprefix '##'];
s.epochs.sleep1 = [14969140.37,20969140.37];
s.epochs.maze1  = [21569140.37,39569140.37];
s.epochs.sleep2 = [41952304.4,47952304.4];
s.other = '';
sessList{end+1} = session(s);

%------------------------------------------------------------
s.animal = '7508';
s.name = '04';
s.group = 'CPP';
s.split = 0;
s.tfiledir = [dataDir 'CPP_Data\7508\7508_04\tfiles'];
s.tfileglob = 'TT*.t';
tfileprefix = 'TT_';
s.sameTTfmt = [tfileprefix '##'];
s.epochs.sleep1 = [29691048.5,35691048.5];
s.epochs.maze1  = [36291048.5,54291048.5];
s.epochs.sleep2 = [58412324.1,64412324.1];
s.other = '';
sessList{end+1} = session(s);

%------------------------------------------------------------
s.animal = '7508';
s.name = '07';
s.group = 'CPP';
s.split = 0;
s.tfiledir = [dataDir 'CPP_Data\7508\7508_07\tfiles'];
s.tfileglob = 'TT*.t';
tfileprefix = 'TT_';
s.sameTTfmt = [tfileprefix '##'];
s.epochs.sleep1 = [21929091.99,27929091.99];
s.epochs.maze1  = [28529091.99,46529091.99];
s.epochs.sleep2 = [48858541.37,54858541.37];
s.other = '';
sessList{end+1} = session(s);

%------------------------------------------------------------
s.animal = '7508';
s.name = '10';
s.group = 'CPP';
s.split = 0;
s.tfiledir = [dataDir 'CPP_Data\7508\7508_10\tfiles'];
s.tfileglob = 'TT*.t';
tfileprefix = 'TT_';
s.sameTTfmt = [tfileprefix '##'];
s.epochs.sleep1 = [20229866.3,26229866.3];
s.epochs.maze1  = [26829866.3,44829866.3];
s.epochs.sleep2 = [46411849.17,52411849.17];
s.other = '';
sessList{end+1} = session(s);

%------------------------------------------------------------
s.animal = '7508';
s.name = '13';
s.group = 'CPP';
s.split = 0;
s.tfiledir = [dataDir 'CPP_Data\7508\7508_13\tfiles'];
s.tfileglob = 'TT*.t';
tfileprefix = 'TT_';
s.sameTTfmt = [tfileprefix '##'];
s.epochs.sleep1 = [31380189.21,37380189.21];
s.epochs.maze1  = [37980189.21,55980189.21];
s.epochs.sleep2 = [59559886.22,65559886.22];
s.other = '';
sessList{end+1} = session(s);

%------------------------------------------------------------
s.animal = '7508';
s.name = '17';
s.group = 'emptyCPP';
s.split = 0;
s.tfiledir = [dataDir 'CPP_Data\7508\7508_17\tfiles'];
s.tfileglob = 'TT*.t';
tfileprefix = 'TT_';
s.sameTTfmt = [tfileprefix '##'];
s.epochs.sleep1 = [28197463.66,34197463.66];
s.epochs.maze1  = [34797463.66,52797463.66];
s.epochs.sleep2 = [51121043.48,57121043.48];
s.other = '';
sessList{end+1} = session(s);

%------------------------------------------------------------
s.animal = '7651';
s.name = '01';
s.group = 'novelCPP';
s.split = 0;
s.tfiledir = [dataDir 'CPP_Data\7651\7651_01\tfiles'];
s.tfileglob = 'TT*.t';
tfileprefix = 'TT_';
s.sameTTfmt = [tfileprefix '##'];
s.epochs.sleep1 = [17985833.02,23985833.02];
s.epochs.maze1  = [24585833.02,42585833.02];
s.epochs.sleep2 = [44558711.44,50558711.44];
s.other = '';
sessList{end+1} = session(s);

%------------------------------------------------------------
s.animal = '7651';
s.name = '04';
s.group = 'CPP';
s.split = 0;
s.tfiledir = [dataDir 'CPP_Data\7651\7651_04\tfiles'];
s.tfileglob = 'TT*.t';
tfileprefix = 'TT_';
s.sameTTfmt = [tfileprefix '##'];
s.epochs.sleep1 = [31622590.58,37622590.58];
s.epochs.maze1  = [38222590.58,56222590.58];
s.epochs.sleep2 = [58951877.91,64951877.91];
s.other = '';
sessList{end+1} = session(s);

%------------------------------------------------------------

s.animal = '7651';
s.name = '07';
s.group = 'CPP';
s.split = 0;
s.tfiledir = [dataDir 'CPP_Data\7651\7651_07\tfiles'];
s.tfileglob = 'TT*.t';
tfileprefix = 'TT_';
s.sameTTfmt = [tfileprefix '##'];
s.epochs.sleep1 = [24926973.36,30926973.36];
s.epochs.maze1  = [31526973.36,49526973.36];
s.epochs.sleep2 = [51723638.37,57723638.37];
s.other = '';
sessList{end+1} = session(s);

%------------------------------------------------------------
s.animal = '7651';
s.name = '10';
s.group = 'CPP';
s.split = 0;
s.tfiledir = [dataDir 'CPP_Data\7651\7651_10\tfiles'];
s.tfileglob = 'TT*.t';
tfileprefix = 'TT_';
s.sameTTfmt = [tfileprefix '##'];
s.epochs.sleep1 = [29071468.66,5071468.66];
s.epochs.maze1  = [35671468.66,53671468.66];
s.epochs.sleep2 = [54686859.33,60686859.33];
s.other = '';
sessList{end+1} = session(s);

%------------------------------------------------------------
s.animal = '7651';
s.name = '13';
s.group = 'CPP';
s.split = 0;
s.tfiledir = [dataDir 'CPP_Data\7651\7651_13\tfiles'];
s.tfileglob = 'TT*.t';
tfileprefix = 'TT_';
s.sameTTfmt = [tfileprefix '##'];
s.epochs.sleep1 = [29969081.89,35969081.89];
s.epochs.maze1  = [36569081.89,54569081.89];
s.epochs.sleep2 = [56343527.14,62343527.14];
s.other = '';
sessList{end+1} = session(s);

%------------------------------------------------------------
s.animal = '7651';
s.name = '17';
s.group = 'Unknown_Bug_CPP';
s.split = 0;
s.tfiledir = [dataDir 'CPP_Data\7651\7651_17\tfiles'];
s.tfileglob = 'TT*.t';
tfileprefix = 'TT_';
s.sameTTfmt = [tfileprefix '##'];
s.epochs.sleep1 = [13671987.07,19671987.07];
s.epochs.maze1  = [20271987.07,38271987.07];
s.epochs.sleep2 = [35918229.01,41918229.01];
s.other = '';
sessList{end+1} = session(s);

%------------------------------------------------------------

s.animal = '7653';
s.name = '01';
s.group = 'novelCPP';
s.split = 0;
s.tfiledir = [dataDir 'CPP_Data\7653\7653_01\tfiles'];
s.tfileglob = 'TT*.t';
tfileprefix = 'TT_';
s.sameTTfmt = [tfileprefix '##'];
s.epochs.sleep1 = [26833857.61,32833857.61];
s.epochs.maze1  = [33433857.61,51433857.61];
s.epochs.sleep2 = [55708800.33,61708800.33];
s.other = '';
sessList{end+1} = session(s);

%------------------------------------------------------------

s.animal = '7653';
s.name = '04';
s.group = 'CPP';
s.split = 0;
s.tfiledir = [dataDir 'CPP_Data\7653\7653_04\tfiles'];
s.tfileglob = 'TT*.t';
tfileprefix = 'TT_';
s.sameTTfmt = [tfileprefix '##'];
s.epochs.sleep1 = [35997314.75,41997314.75];
s.epochs.maze1  = [42597314.75,60597314.75];
s.epochs.sleep2 = [63196799.39,69196799.39];
s.other = '';
sessList{end+1} = session(s);

%------------------------------------------------------------

s.animal = '7653';
s.name = '07';
s.group = 'CPP';
s.split = 0;
s.tfiledir = [dataDir 'CPP_Data\7653\7653_07\tfiles'];
s.tfileglob = 'TT*.t';
tfileprefix = 'TT_';
s.sameTTfmt = [tfileprefix '##'];
s.epochs.sleep1 = [25158403.33,31158403.33];
s.epochs.maze1  = [31758403.33,49758403.33];
s.epochs.sleep2 = [48748883.9,54748883.9];
s.other = '';
sessList{end+1} = session(s);

%------------------------------------------------------------

s.animal = '7653';
s.name = '10';
s.group = 'CPP';
s.split = 0;
s.tfiledir = [dataDir 'CPP_Data\7653\7653_10\tfiles'];
s.tfileglob = 'TT*.t';
tfileprefix = 'TT_';
s.sameTTfmt = [tfileprefix '##'];
s.epochs.sleep1 = [24341158.1,30341158.1];
s.epochs.maze1  = [30941158.1,48941158.1];
s.epochs.sleep2 = [48447852.43,54447852.43];
s.other = '';
sessList{end+1} = session(s);

%------------------------------------------------------------

s.animal = '7653';
s.name = '13';
s.group = 'CPP';
s.split = 0;
s.tfiledir = [dataDir 'CPP_Data\7653\7653_13\tfiles'];
s.tfileglob = 'TT*.t';
tfileprefix = 'TT_';
s.sameTTfmt = [tfileprefix '##'];
s.epochs.sleep1 = [20610630.36,26610630.36];
s.epochs.maze1  = [27210630.36,45210630.36];
s.epochs.sleep2 = [43968368.15,49968368.15];
s.other = '';
sessList{end+1} = session(s);

%------------------------------------------------------------

s.animal = '7653';
s.name = '16';
s.group = 'CPP';
s.split = 0;
s.tfiledir = [dataDir 'CPP_Data\7653\7653_16\tfiles'];
s.tfileglob = 'TT*.t';
tfileprefix = 'TT_';
s.sameTTfmt = [tfileprefix '##'];
s.epochs.sleep1 = [21419228.07,27419228.07];
s.epochs.maze1  = [28019228.07,46019228.07];
s.epochs.sleep2 = [44766526.84,50766526.84];
s.other = '';
sessList{end+1} = session(s);

%------------------------------------------------------------

s.animal = '7658';
s.name = '01';
s.group = 'novelCPP';
s.split = 0;
s.tfiledir = [dataDir 'CPP_Data\7658\7658_01\tfiles'];
s.tfileglob = 'TT*.t';
tfileprefix = 'TT_';
s.sameTTfmt = [tfileprefix '##'];
s.epochs.sleep1 = [43977800.12,49977800.12];
s.epochs.maze1  = [50577800.12,68577800.12];
s.epochs.sleep2 = [71142769.57,77142769.57];
s.other = '';
sessList{end+1} = session(s);

%------------------------------------------------------------

s.animal = '7658';
s.name = '06';
s.group = 'CPP';
s.split = 0;
s.tfiledir = [dataDir 'CPP_Data\7658\7658_06\tfiles'];
s.tfileglob = 'TT*.t';
tfileprefix = 'TT_';
s.sameTTfmt = [tfileprefix '##'];
s.epochs.sleep1 = [25750579.55,31750579.55];
s.epochs.maze1  = [32350579.55,50350579.55];
s.epochs.sleep2 = [53127399.62,59127399.62];
s.other = '';
sessList{end+1} = session(s);

%------------------------------------------------------------

s.animal = '7658';
s.name = '09';
s.group = 'CPP';
s.split = 0;
s.tfiledir = [dataDir 'CPP_Data\7658\7658_09\tfiles'];
s.tfileglob = 'TT*.t';
tfileprefix = 'TT_';
s.sameTTfmt = [tfileprefix '##'];
s.epochs.sleep1 = [23963209.79,29963209.79];
s.epochs.maze1  = [30563209.79,48563209.79];
s.epochs.sleep2 = [50974489.32,56974489.32];
s.other = '';
sessList{end+1} = session(s);

%------------------------------------------------------------

s.animal = '7658';
s.name = '12';
s.group = 'CPP';
s.split = 0;
s.tfiledir = [dataDir 'CPP_Data\7658\7658_12\tfiles'];
s.tfileglob = 'TT*.t';
tfileprefix = 'TT_';
s.sameTTfmt = [tfileprefix '##'];
s.epochs.sleep1 = [19779682.65,25779682.65];
s.epochs.maze1  = [26379682.65,44379682.65];
s.epochs.sleep2 = [47371122.73,53371122.73];
s.other = '';
sessList{end+1} = session(s);

%------------------------------------------------------------

s.animal = '7658';
s.name = '15';
s.group = 'CPP';
s.split = 0;
s.tfiledir = [dataDir 'CPP_Data\7658\7658_15\tfiles'];
s.tfileglob = 'TT*.t';
tfileprefix = 'TT_';
s.sameTTfmt = [tfileprefix '##'];
s.epochs.sleep1 = [26105427.28,32105427.28];
s.epochs.maze1  = [32705427.28,50705427.28];
s.epochs.sleep2 = [52731777.78,58731777.78];
s.other = '';
sessList{end+1} = session(s);

%------------------------------------------------------------

s.animal = '7658';
s.name = '18';
s.group = 'CPP';
s.split = 0;
s.tfiledir = [dataDir 'CPP_Data\7658\7658_18\tfiles'];
s.tfileglob = 'TT*.t';
tfileprefix = 'TT_';
s.sameTTfmt = [tfileprefix '##'];
s.epochs.sleep1 = [19494439.42,25494439.42];
s.epochs.maze1  = [26094439.42,44094439.42];
s.epochs.sleep2 = [46154712.01,52154712.01];
s.other = '';
sessList{end+1} = session(s);

%------------------------------------------------------------

%------------------------------------------------------------
s.animal = '7508';
s.name = '03';
s.group = 'novelSaline';
s.split = 0;
s.tfiledir = [dataDir 'Saline_Data\7508\7508_03\tfiles'];
s.tfileglob = 'TT*.t';
tfileprefix = 'TT_';
s.sameTTfmt = [tfileprefix '##'];
s.epochs.sleep1 = [30745628.15,36745628.15];
s.epochs.maze1  = [37345628.15,55345628.15];
s.epochs.sleep2 = [57470684.78,63470684.78];
s.other = '';
sessList{end+1} = session(s);

%------------------------------------------------------------
s.animal = '7508';
s.name = '05';
s.group = 'Saline';
s.split = 0;
s.tfiledir = [dataDir 'Saline_Data\7508\7508_05\tfiles'];
s.tfileglob = 'TT*.t';
tfileprefix = 'TT_';
s.sameTTfmt = [tfileprefix '##'];
s.epochs.sleep1 = [22245312.33,28245312.33];
s.epochs.maze1  = [28845312.33,46845312.33];
s.epochs.sleep2 = [48997223.81,54997223.81];
s.other = '';
sessList{end+1} = session(s);

%------------------------------------------------------------
s.animal = '7508';
s.name = '08';
s.group = 'Saline';
s.split = 0;
s.tfiledir = [dataDir 'Saline_Data\7508\7508_08\tfiles'];
s.tfileglob = 'TT*.t';
tfileprefix = 'TT_';
s.sameTTfmt = [tfileprefix '##'];
s.epochs.sleep1 = [24964549.57,30964549.57];
s.epochs.maze1  = [31564549.57,49564549.57];
s.epochs.sleep2 = [49596404.54,55596404.54];
s.other = '';
sessList{end+1} = session(s);

%------------------------------------------------------------
s.animal = '7508';
s.name = '12';
s.group = 'Saline';
s.split = 0;
s.tfiledir = [dataDir 'Saline_Data\7508\7508_12\tfiles'];
s.tfileglob = 'TT*.t';
tfileprefix = 'TT_';
s.sameTTfmt = [tfileprefix '##'];
s.epochs.sleep1 = [23904201.63,29904201.63];
s.epochs.maze1  = [30504201.63,48504201.63];
s.epochs.sleep2 = [48473396.87,54473396.87];
s.other = '';
sessList{end+1} = session(s);

%------------------------------------------------------------
s.animal = '7508';
s.name = '15';
s.group = 'Saline';
s.split = 0;
s.tfiledir = [dataDir 'Saline_Data\7508\7508_15\tfiles'];
s.tfileglob = 'TT*.t';
tfileprefix = 'TT_';
s.sameTTfmt = [tfileprefix '##'];
s.epochs.sleep1 = [31747391.09,37747391.09];
s.epochs.maze1  = [38347391.09,56347391.09];
s.epochs.sleep2 = [54571013.14,60571013.14];
s.other = '';
sessList{end+1} = session(s);

%------------------------------------------------------------
s.animal = '7508';
s.name = '16';
s.group = 'Saline';
s.split = 0;
s.tfiledir = [dataDir 'Saline_Data\7508\7508_16\tfiles'];
s.tfileglob = 'TT*.t';
tfileprefix = 'TT_';
s.sameTTfmt = [tfileprefix '##'];
s.epochs.sleep1 = [24252414.25,30252414.25];
s.epochs.maze1  = [30852414.25,48852414.25];
s.epochs.sleep2 = [47683259.56,53683259.56];
s.other = '';
sessList{end+1} = session(s);

%------------------------------------------------------------
%------------------------------------------------------------

s.animal = '7651';
s.name = '02';
s.group = 'novelSaline';
s.split = 0;
s.tfiledir = [dataDir 'Saline_Data\7651\7651_02\tfiles'];
s.tfileglob = 'TT*.t';
tfileprefix = 'TT_';
s.sameTTfmt = [tfileprefix '##'];
s.epochs.sleep1 = [-3332596.95,2667403.05];
s.epochs.maze1  = [3267403.05,21267403.05];
s.epochs.sleep2 = [19934131.78,25934131.78];
s.other = '';
sessList{end+1} = session(s);

%------------------------------------------------------------
s.animal = '7651';
s.name = '05';
s.group = 'Saline';
s.split = 0;
s.tfiledir = [dataDir 'Saline_Data\7651\7651_05\tfiles'];
s.tfileglob = 'TT*.t';
tfileprefix = 'TT_';
s.sameTTfmt = [tfileprefix '##'];
s.epochs.sleep1 = [29206299.09,35206299.09];
s.epochs.maze1  = [35806299.09,53806299.09];
s.epochs.sleep2 = [53357750.28,59357750.28];
s.other = '';
sessList{end+1} = session(s);

%------------------------------------------------------------

s.animal = '7651';
s.name = '08';
s.group = 'Saline';
s.split = 0;
s.tfiledir = [dataDir 'Saline_Data\7651\7651_08\tfiles'];
s.tfileglob = 'TT*.t';
tfileprefix = 'TT_';
s.sameTTfmt = [tfileprefix '##'];
s.epochs.sleep1 = [28851275.74,34851275.74];
s.epochs.maze1  = [35451275.74,53451275.74];
s.epochs.sleep2 = [52143367.6,58143367.6];
s.other = '';
sessList{end+1} = session(s);

%------------------------------------------------------------
s.animal = '7651';
s.name = '12';
s.group = 'Saline';
s.split = 0;
s.tfiledir = [dataDir 'Saline_Data\7651\7651_12\tfiles'];
s.tfileglob = 'TT*.t';
tfileprefix = 'TT_';
s.sameTTfmt = [tfileprefix '##'];
s.epochs.sleep1 = [19686454.73,25686454.73];
s.epochs.maze1  = [26286454.73,44286454.73];
s.epochs.sleep2 = [40631722.82,46631722.82];
s.other = '';
sessList{end+1} = session(s);

%------------------------------------------------------------
s.animal = '7651';
s.name = '14';
s.group = 'Saline';
s.split = 0;
s.tfiledir = [dataDir 'Saline_Data\7651\7651_14\tfiles'];
s.tfileglob = 'TT*.t';
tfileprefix = 'TT_';
s.sameTTfmt = [tfileprefix '##'];
s.epochs.sleep1 = [23924426.1,29924426.1];
s.epochs.maze1  = [30524426.1,48524426.1];
s.epochs.sleep2 = [45174049.2,51174049.2];
s.other = '';
sessList{end+1} = session(s);

%------------------------------------------------------------
s.animal = '7651';
s.name = '16';
s.group = 'Saline';
s.split = 0;
s.tfiledir = [dataDir 'Saline_Data\7651\7651_16\tfiles'];
s.tfileglob = 'TT*.t';
tfileprefix = 'TT_';
s.sameTTfmt = [tfileprefix '##'];
s.epochs.sleep1 = [22993079.1,28993079.1];
s.epochs.maze1  = [29593079.1,47593079.1];
s.epochs.sleep2 = [44517218.79,50517218.79];
s.other = '';
sessList{end+1} = session(s);

%------------------------------------------------------------
%------------------------------------------------------------

s.animal = '7653';
s.name = '02';
s.group = 'novelSaline';
s.split = 0;
s.tfiledir = [dataDir 'Saline_Data\7653\7653_02\tfiles'];
s.tfileglob = 'TT*.t';
tfileprefix = 'TT_';
s.sameTTfmt = [tfileprefix '##'];
s.epochs.sleep1 = [31932065.8,37932065.8];
s.epochs.maze1  = [38532065.8,56532065.8];
s.epochs.sleep2 = [59052487.03,65052487.03];
s.other = '';
sessList{end+1} = session(s);

%------------------------------------------------------------

s.animal = '7653';
s.name = '05';
s.group = 'Saline';
s.split = 0;
s.tfiledir = [dataDir 'Saline_Data\7653\7653_05\tfiles'];
s.tfileglob = 'TT*.t';
tfileprefix = 'TT_';
s.sameTTfmt = [tfileprefix '##'];
s.epochs.sleep1 = [20019459.25,26019459.25];
s.epochs.maze1  = [26619459.25,44619459.25];
s.epochs.sleep2 = [44612529.38,50612529.38];
s.other = '';
sessList{end+1} = session(s);

%------------------------------------------------------------

s.animal = '7653';
s.name = '08';
s.group = 'Saline';
s.split = 0;
s.tfiledir = [dataDir 'Saline_Data\7653\7653_08\tfiles'];
s.tfileglob = 'TT*.t';
tfileprefix = 'TT_';
s.sameTTfmt = [tfileprefix '##'];
s.epochs.sleep1 = [18140975.96,24140975.96];
s.epochs.maze1  = [24740975.96,42740975.96];
s.epochs.sleep2 = [40757219.23,46757219.23];
s.other = '';
sessList{end+1} = session(s);

%------------------------------------------------------------

s.animal = '7653';
s.name = '12';
s.group = 'Saline';
s.split = 0;
s.tfiledir = [dataDir 'Saline_Data\7653\7653_12\tfiles'];
s.tfileglob = 'TT*.t';
tfileprefix = 'TT_';
s.sameTTfmt = [tfileprefix '##'];
s.epochs.sleep1 = [21111527.83,27111527.83];
s.epochs.maze1  = [27711527.83,45711527.83];
s.epochs.sleep2 = [44272328.49,50272328.49];
s.other = '';
sessList{end+1} = session(s);

%------------------------------------------------------------

s.animal = '7653';
s.name = '15';
s.group = 'Saline';
s.split = 0;
s.tfiledir = [dataDir 'Saline_Data\7653\7653_15\tfiles'];
s.tfileglob = 'TT*.t';
tfileprefix = 'TT_';
s.sameTTfmt = [tfileprefix '##'];
s.epochs.sleep1 = [15774117.25,21774117.25];
s.epochs.maze1  = [22374117.25,40374117.25];
s.epochs.sleep2 = [38590279.27,44590279.27];
s.other = '';
sessList{end+1} = session(s);

%------------------------------------------------------------

s.animal = '7653';
s.name = '17';
s.group = 'Saline';
s.split = 0;
s.tfiledir = [dataDir 'Saline_Data\7653\7653_17\tfiles'];
s.tfileglob = 'TT*.t';
tfileprefix = 'TT_';
s.sameTTfmt = [tfileprefix '##'];
s.epochs.sleep1 = [17895522.81,23895522.81];
s.epochs.maze1  = [24495522.81,42495522.81];
s.epochs.sleep2 = [40668423.54,46668423.54];
s.other = '';
sessList{end+1} = session(s);

%------------------------------------------------------------

s.animal = '7658';
s.name = '02';
s.group = 'novelSaline';
s.split = 0;
s.tfiledir = [dataDir 'Saline_Data\7658\7658_02\tfiles'];
s.tfileglob = 'TT*.t';
tfileprefix = 'TT_';
s.sameTTfmt = [tfileprefix '##'];
s.epochs.sleep1 = [20295488.93,26295488.93];
s.epochs.maze1  = [26895488.93,44895488.93];
s.epochs.sleep2 = [47757457.06,53757457.06];
s.other = '';
sessList{end+1} = session(s);

%------------------------------------------------------------

s.animal = '7658';
s.name = '07';
s.group = 'Saline';
s.split = 0;
s.tfiledir = [dataDir 'Saline_Data\7658\7658_07\tfiles'];
s.tfileglob = 'TT*.t';
tfileprefix = 'TT_';
s.sameTTfmt = [tfileprefix '##'];
s.epochs.sleep1 = [21276877.77,27276877.77];
s.epochs.maze1  = [27876877.77,45876877.77];
s.epochs.sleep2 = [48852443.78,54852443.78];
s.other = '';
sessList{end+1} = session(s);

%------------------------------------------------------------

s.animal = '7658';
s.name = '10';
s.group = 'Saline';
s.split = 0;
s.tfiledir = [dataDir 'Saline_Data\7658\7658_10\tfiles'];
s.tfileglob = 'TT*.t';
tfileprefix = 'TT_';
s.sameTTfmt = [tfileprefix '##'];
s.epochs.sleep1 = [23055786.11,29055786.11];
s.epochs.maze1  = [29655786.11,47655786.11];
s.epochs.sleep2 = [48464359.22,54464359.22];
s.other = '';
sessList{end+1} = session(s);

%------------------------------------------------------------

s.animal = '7658';
s.name = '13';
s.group = 'Saline';
s.split = 0;
s.tfiledir = [dataDir 'Saline_Data\7658\7658_13\tfiles'];
s.tfileglob = 'TT*.t';
tfileprefix = 'TT_';
s.sameTTfmt = [tfileprefix '##'];
s.epochs.sleep1 = [20471785.62,26471785.62];
s.epochs.maze1  = [27071785.62,45071785.62];
s.epochs.sleep2 = [45385354.56,51385354.56];
s.other = '';
sessList{end+1} = session(s);

%------------------------------------------------------------

s.animal = '7658';
s.name = '16';
s.group = 'Saline';
s.split = 0;
s.tfiledir = [dataDir 'Saline_Data\7658\7658_16\tfiles'];
s.tfileglob = 'TT*.t';
tfileprefix = 'TT_';
s.sameTTfmt = [tfileprefix '##'];
s.epochs.sleep1 = [21849269.59,27849269.59];
s.epochs.maze1  = [28449269.59,46449269.59];
s.epochs.sleep2 = [45983979.47,51983979.47];
s.other = '';
sessList{end+1} = session(s);

%------------------------------------------------------------

s.animal = '7658';
s.name = '20';
s.group = 'Saline';
s.split = 0;
s.tfiledir = [dataDir 'Saline_Data\7658\7658_20\tfiles'];
s.tfileglob = 'TT*.t';
tfileprefix = 'TT_';
s.sameTTfmt = [tfileprefix '##'];
s.epochs.sleep1 = [24787801.64,30787801.64];
s.epochs.maze1  = [31387801.64,49387801.64];
s.epochs.sleep2 = [52226928.54,58226928.54];
s.other = '';
sessList{end+1} = session(s);

%------------------------------------------------------------

%------------------------------------------------------------


