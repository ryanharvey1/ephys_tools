% HDLethCrcl
% Generate HD barplots or polar plots (as specified) binned by the specified amount.
% Both raster types as saved to .jpg files in current directory.
% Written by AW for Lethbridge Neuralynx Configs and any .pvd file with
% direction data in the d column
% Modified to Run on data collected on original Hendrik Never Quite Worked
% Program or Data Collected on Valery Version

%Test Data Set
% sess_info.rat = 'TestData';
% datapath = '\\huxley.resrch.uleth.ca\EP\ppc-seq-lrng\Head Direction Test 2011-04-22_11-38-27\';
% tfilespath = [];
% sess_info.sess_num = 1000000;
% spike_data = [];
% EventFileDur.North = (939806060-930396483);
% EventFileDur.East = (947613803-940306660);
% EventFileDur.South = (957723418-948113511);
% EventFileDur.West = (962687168-957723419);

% P0017
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Work Computer I think:
% sess_info.rat = 'P0017';
% datapath = 'C:\Data\PPCseqLrng\2011-03-28_07-18-02PPC_p\';
% tfilespath = 'C:\Data\PPCseqLrng\2011-03-28_07-18-02PPC_p\TT\TT16\tfiles';
% sess_info.sess_num = 1; %for no barriers alteration

%Work Computer:
% 4\2\2011
% sess_info.rat = 'P0017';
% datapath = '\\huxley.resrch.uleth.ca\EP\ppc-seq-lrng\P0017\RandomLights\2011-04-02_17-14-54PPCp\';
% tfilespath = '\\huxley.resrch.uleth.ca\EP\ppc-seq-lrng\P0017\RandomLights\2011-04-02_17-14-54PPCp\tfiles';
% sess_info.sess_num = 3.012;

% 4/3/2011
% sess_info.rat = 'P0017';
% datapath = '\\huxley.resrch.uleth.ca\EP\ppc-seq-lrng\P0017\RandomLights\2011-04-03_14-43-40_PPCp\';
% tfilespath = '\\huxley.resrch.uleth.ca\EP\ppc-seq-lrng\P0017\RandomLights\2011-04-03_14-43-40_PPCp\tfiles';
% sess_info.sess_num = 3.021; %for sequence1

% 4\4\2011 PPC
% sess_info.rat = 'P0017';
% datapath = '\\huxley.resrch.uleth.ca\EP\ppc-seq-lrng\P0017\RandomLights\2011-04-04_07-13-56pPPC\';
% datapath = '\\huxley.resrch.uleth.ca\EP\ppc-seq-lrng\P0017\RandomLights\2011-04-04_07-13-56pPPC\';
%  % A-C ranked principle cells:
% % % tfilespath = '\\huxley.resrch.uleth.ca\EP\ppc-seq-lrng\P0017\RandomLights\2011-04-04_07-13-56pPPC\tfiles';
% % % interneuron
% tfilespath = '\\huxley.resrch.uleth.ca\EP\ppc-seq-lrng\P0017\RandomLights\2011-04-04_07-13-56pPPC\TT\interneuron';
% sess_info.sess_num = 3.031;

% 4\4\2011 HPC
% sess_info.rat = 'P0017';
% 
% sess_info.sess_num = 3.031;

% 4\6\2011
% sess_info.rat = 'P0017';
% datapath = '\\huxley.resrch.uleth.ca\EP\ppc-seq-lrng\P0017\RandomLights\2011-04-06_15-29-41PPCp\';
% % tfilespath = '\\huxley.resrch.uleth.ca\EP\ppc-seq-lrng\P0017\RandomLights\2011-04-06_15-29-41PPCp\tfiles';
% sess_info.sess_num = 3.052;
% % % all
% % % A-C ranked principle cells:
% % % tfilespath = '\\huxley.resrch.uleth.ca\EP\ppc-seq-lrng\P0017\RandomLights\2011-04-06_15-29-41PPCp\TT\tfiles';
% % % interneuron
% tfilespath = '\\huxley.resrch.uleth.ca\EP\ppc-seq-lrng\P0017\RandomLights\2011-04-06_15-29-41PPCp\TT\interneuron_tfiles';

% Sequence Task (P0017)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4/15/2011
% sess_info.rat = 'P0017';
% datapath = '\\huxley.resrch.uleth.ca\EP\ppc-seq-lrng\P0017\Sequence1\2011-04-15_10-06-50_PPCp\';
% tfilespath = '\\huxley.resrch.uleth.ca\EP\ppc-seq-lrng\P0017\Sequence1\2011-04-15_10-06-50_PPCp\tfiles';
% sess_info.sess_num = 4.031; %for sequence1

% 4/17/2011
% sess_info.rat = 'P0017';
% datapath = '\\huxley.resrch.uleth.ca\EP\ppc-seq-lrng\P0017\Sequence1\2011-04-17_15-44-53_PPCp\';
% tfilespath = '\\huxley.resrch.uleth.ca\EP\ppc-seq-lrng\P0017\Sequence1\2011-04-17_15-44-53_PPCp\tfiles';
% sess_info.sess_num = 4.052; %for sequence1

% 4/18/2011
% sess_info.rat = 'P0017';
% datapath = '\\huxley.resrch.uleth.ca\EP\ppc-seq-lrng\P0017\Sequence1\2011-04-18_06-53-55_PPCp\';
% tfilespath = '\\huxley.resrch.uleth.ca\EP\ppc-seq-lrng\P0017\Sequence1\2011-04-18_06-53-55_PPCp\tfiles';
% sess_info.sess_num = 4.061; %for sequence1

% 4/18/2011
% sess_info.rat = 'P0017';
% datapath = '\\huxley.resrch.uleth.ca\EP\ppc-seq-lrng\P0017\Sequence1\2011-04-18_17-07-00_PPCp\';
% tfilespath = '\\huxley.resrch.uleth.ca\EP\ppc-seq-lrng\P0017\Sequence1\2011-04-18_17-07-00_PPCp\tfiles';
% sess_info.sess_num = 4.062; %for sequence1

% 4/19/2011
% sess_info.rat = 'P0017';
% datapath = '\\huxley.resrch.uleth.ca\EP\ppc-seq-lrng\P0017\Sequence1\2011-04-19_07-10-19_PPCp\';
% tfilespath = '\\huxley.resrch.uleth.ca\EP\ppc-seq-lrng\P0017\Sequence1\2011-04-19_07-10-19_PPCp\tfiles';
% sess_info.sess_num = 4.071; %for sequence1

% 4/20/2011
% sess_info.rat = 'P0017';
% datapath = '\\huxley.resrch.uleth.ca\EP\ppc-seq-lrng\P0017\Sequence1\2011-04-20_07-15-47_PPCp\';
% tfilespath = '\\huxley.resrch.uleth.ca\EP\ppc-seq-lrng\P0017\Sequence1\2011-04-20_07-15-47_PPCp\tfiles';
% sess_info.sess_num = 4.081; %for sequence1

% 4/20/2011 m2
% sess_info.rat = 'P0017';
% datapath = '\\huxley.resrch.uleth.ca\EP\ppc-seq-lrng\P0017\Sequence1\2011-04-20_15-21-37_PPCp\';
% tfilespath = '\\huxley.resrch.uleth.ca\EP\ppc-seq-lrng\P0017\Sequence1\2011-04-20_15-21-37_PPCp\tfiles';
% sess_info.sess_num = 4.082; %for sequence1

% % 4/21/2011
% sess_info.rat = 'P0017';
% datapath = '\\huxley.resrch.uleth.ca\EP\ppc-seq-lrng\P0017\Sequence1\2011-04-21_07-03-28_PPCp\';
% tfilespath = '\\huxley.resrch.uleth.ca\EP\ppc-seq-lrng\P0017\Sequence1\2011-04-21_07-03-28_PPCp\tfiles';
% sess_info.sess_num = 4.091; %for sequence1

% 4/22/2011
% sess_info.rat = 'P0017';
% datapath = '\\huxley.resrch.uleth.ca\EP\ppc-seq-lrng\P0017\Sequence1\2011-04-22_12-10-37_PPCp\';
% tfilespath = '\\huxley.resrch.uleth.ca\EP\ppc-seq-lrng\P0017\Sequence1\2011-04-22_12-10-37_PPCp\tfiles';
% sess_info.sess_num = 4.102; %for sequence1

% 4/24/2011
sess_info.rat = 'P0017';
datapath = '\\huxley.resrch.uleth.ca\EP\ppc-seq-lrng\P0017\Sequence1\2011-04-24_15-50-58_PPCp\';
tfilespath = '\\huxley.resrch.uleth.ca\EP\ppc-seq-lrng\P0017\Sequence1\2011-04-24_15-50-58_PPCp\tfiles';
sess_info.sess_num = 4.122; %for sequence1

% 4/26/2011
% sess_info.rat = 'P0017';
% datapath = '\\huxley.resrch.uleth.ca\EP\ppc-seq-lrng\P0017\Sequence1\2011-04-26_07-03-28_PPCp\';
% tfilespath = '\\huxley.resrch.uleth.ca\EP\ppc-seq-lrng\P0017\Sequence1\2011-04-26_07-03-28_PPCp\tfiles';
% sess_info.sess_num = 4.14; %for sequence1

% %5/1/2011
% sess_info.rat = 'P0017';
% datapath = '\\huxley.resrch.uleth.ca\EP\ppc-seq-lrng\P0017\Sequence1\2011-05-01_11-49-45_PPC_p\';
% tfilespath = '\\huxley.resrch.uleth.ca\EP\ppc-seq-lrng\P0017\Sequence1\2011-05-01_11-49-45_PPC_p\tfiles';
% sess_info.sess_num = 4.19; %for sequence1

% 5/2/2011
% sess_info.rat = 'P0017';
% datapath = '\\huxley.resrch.uleth.ca\EP\ppc-seq-lrng\P0017\Sequence1\2011-05-02_09-43-38_PPC_p\';
% tfilespath = '\\huxley.resrch.uleth.ca\EP\ppc-seq-lrng\P0017\Sequence1\2011-05-02_09-43-38_PPC_p\tfiles';
% sess_info.sess_num = 4.20; %for sequence1

%5/3/2011
% sess_info.rat = 'P0017';
% datapath = '\\huxley.resrch.uleth.ca\EP\ppc-seq-lrng\P0017\Sequence1\2011-05-03_08-35-19_PPC_p\';
% tfilespath = '\\huxley.resrch.uleth.ca\EP\ppc-seq-lrng\P0017\Sequence1\2011-05-03_08-35-19_PPC_p\TT\tfiles';
% sess_info.sess_num = 4.21;

% P0100
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 9/2/2011
% sess_info.rat = 'P0100';
% datapath = '\\huxley.resrch.uleth.ca\EP\ppc-seq-lrng\P0100\RandomLights\2011-09-02_10-13-59_PPCp\';
% % % A-C ranked principle cells:
% % % tfilespath = '\\huxley.resrch.uleth.ca\EP\ppc-seq-lrng\P0100\RandomLights\2011-09-02_10-13-59_PPCp\tfiles';
% % % interneuron
% tfilespath = '\\huxley.resrch.uleth.ca\EP\ppc-seq-lrng\P0100\RandomLights\2011-09-02_10-13-59_PPCp\TT\interneuron_tfiles';
% sess_info.sess_num = 3.02; %for sequence1

% 9/4/2011 PPC
% sess_info.rat = 'P0100';
% datapath = '\\huxley.resrch.uleth.ca\EP\ppc-seq-lrng\P0100\RandomLights\2011-09-04_12-04-21_PPCp\';
% tfilespath = '\\huxley.resrch.uleth.ca\EP\ppc-seq-lrng\P0100\RandomLights\2011-09-04_12-04-21_PPCp\tfiles';
% sess_info.sess_num = 3.03; %for sequence1

% 9/4/2011 HPC
% sess_info.rat = 'P0100';
% datapath = '\\huxley.resrch.uleth.ca\EP\ppc-seq-lrng\P0100\RandomLights\2011-09-04_12-04-21_PPCp\';
% tfilespath = '\\huxley.resrch.uleth.ca\EP\ppc-seq-lrng\P0100\RandomLights\2011-04-04_07-13-56_HPCp\tfiles';
% sess_info.sess_num = 3.03; %for sequence1

% 9/5/2011 PPC
% sess_info.rat = 'P0100';
% datapath = '\\huxley.resrch.uleth.ca\EP\ppc-seq-lrng\P0100\RandomLights\2011-09-04_12-04-21_PPCp\';
% sess_info.sess_num = 3.04; %for sequence1

% 9/5/2011 HPC
% sess_info.rat = 'P0100';
% datapath = '\\huxley.resrch.uleth.ca\EP\ppc-seq-lrng\P0100\RandomLights\2011-09-04_12-04-21_PPCp\';
% tfilespath = '\\huxley.resrch.uleth.ca\EP\ppc-seq-lrng\P0100\RandomLights\2011-09-05_12-25-51_HPCp\tfiles';
% sess_info.sess_num = 3.04; %for sequence1

% 9/21/2011
% sess_info.rat = 'P0100';
% datapath = '\\huxley.resrch.uleth.ca\EP\ppc-seq-lrng\P0100\Sequence1_Disambiguation\2011-09-21_09-57-40_PPCp\';
% tfilespath = '\\huxley.resrch.uleth.ca\EP\ppc-seq-lrng\P0100\Sequence1_Disambiguation\2011-09-21_09-57-40_PPCp\TT\tfiles';
% sess_info.sess_num = 4.14; %for sequence1 

% 9/28/2011
% sess_info.rat = 'P0100';
% datapath = '\\huxley.resrch.uleth.ca\EP\ppc-seq-lrng\P0100\Sequence1_Disambiguation\2011-09-28_09-48-57_PPC_p\';
% tfilespath = '\\huxley.resrch.uleth.ca\EP\ppc-seq-lrng\P0100\Sequence1_Disambiguation\2011-09-28_09-48-57_PPC_p\TT\tfiles';
% sess_info.sess_num = 4.21; %for sequence1 

% P0662
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 8/12/2012
% sess_info.rat = 'P0662';
% datapath = '\\huxley.resrch.uleth.ca\EP\ppc-seq-lrng\P0662\RandomLights_Sess3_xx\Data\2012-08-12_12-48-04_PPCp\';
% % A-C ranked principle cells:
% % tfilespath = '\\huxley.resrch.uleth.ca\EP\ppc-seq-lrng\P0662\RandomLights_Sess3_xx\Data\2012-08-12_12-48-04_PPCp\tfiles';
% % interneuron
% tfilespath = '\\huxley.resrch.uleth.ca\EP\ppc-seq-lrng\P0662\RandomLights_Sess3_xx\Data\2012-08-12_12-48-04_PPCp\TT\interneuron';
% sess_info.sess_num = 3.01; %for sequence1

% 8/14/2012
% sess_info.rat = 'P0662';
% datapath = '\\huxley.resrch.uleth.ca\EP\ppc-seq-lrng\P0662\RandomLights_Sess3_xx\Data\2012-08-14_07-59-59_PPCp\';
% tfilespath = '\\huxley.resrch.uleth.ca\EP\ppc-seq-lrng\P0662\RandomLights_Sess3_xx\Data\2012-08-14_07-59-59_PPCp\tfiles';
% sess_info.sess_num = 3.03; %for sequence1

% 8/16/2012
% sess_info.rat = 'P0662';
% datapath = '\\huxley.resrch.uleth.ca\EP\ppc-seq-lrng\P0662\RandomLights_Sess3_xx\Data\2012-08-16_08-19-07_PPCp\';
% tfilespath = '\\huxley.resrch.uleth.ca\EP\ppc-seq-lrng\P0662\RandomLights_Sess3_xx\Data\2012-08-16_08-19-07_PPCp\tfiles';
% sess_info.sess_num = 3.05; %for sequence1

% 8/17/2012
% sess_info.rat = 'P0662';
% datapath = '\\huxley.resrch.uleth.ca\EP\ppc-seq-lrng\P0662\RandomLights_Sess3_xx\Data\2012-08-17_08-32-42_PPCp\';
% % A-C ranked principle cells:
% % tfilespath = '\\huxley.resrch.uleth.ca\EP\ppc-seq-lrng\P0662\RandomLights_Sess3_xx\Data\2012-08-17_08-32-42_PPCp\tfiles';
% % A-C ranked interneurons cells:
% tfilespath = '\\huxley.resrch.uleth.ca\EP\ppc-seq-lrng\P0662\RandomLights_Sess3_xx\Data\2012-08-18_09-23-41_M1_Onward_PPCp\TT\interneuron_tfiles';
% sess_info.sess_num = 3.08; %for sequence1

% 8/18/2012
% sess_info.rat = 'P0662';
% datapath = '\\huxley.resrch.uleth.ca\EP\ppc-seq-lrng\P0662\RandomLights_Sess3_xx\Data\2012-08-18_09-23-41_M1_Onward_PPCp\';
% % A-C ranked principle cells:
% % tfilespath = '\\huxley.resrch.uleth.ca\EP\ppc-seq-lrng\P0662\RandomLights_Sess3_xx\Data\2012-08-17_08-32-42_PPCp\tfiles';
% % A-C ranked interneurons cells:
% tfilespath = '\\huxley.resrch.uleth.ca\EP\ppc-seq-lrng\P0662\RandomLights_Sess3_xx\Data\2012-08-18_09-23-41_M1_Onward_PPCp\TT\interneuron_tfiles';
% sess_info.sess_num = 3.07; %for sequence1

% 8/19/2012
% sess_info.rat = 'P0662';
% datapath = '\\huxley.resrch.uleth.ca\EP\ppc-seq-lrng\P0662\RandomLights_Sess3_xx\Data\2012-08-19_14-11-14_PPCp\';
% tfilespath = '\\huxley.resrch.uleth.ca\EP\ppc-seq-lrng\P0662\RandomLights_Sess3_xx\Data\2012-08-19_14-11-14_PPCp\tfiles';
% sess_info.sess_num = 3.08; %for sequence1

% 8/20/2012
% sess_info.rat = 'P0662';
% datapath = '\\huxley.resrch.uleth.ca\EP\ppc-seq-lrng\P0662\RandomLights_Sess3_xx\Data\2012-08-20_07-35-37_PPCp\';
% tfilespath = '\\huxley.resrch.uleth.ca\EP\ppc-seq-lrng\P0662\RandomLights_Sess3_xx\Data\2012-08-20_07-35-37_PPCp\tfiles\';
% sess_info.sess_num = 3.09; %for sequence1

% 8/21/2012
% sess_info.rat = 'P0662';
% datapath = '\\huxley.resrch.uleth.ca\EP\ppc-seq-lrng\P0662\RandomLights_Sess3_xx\Data\2012-08-21_08-49-01_PPCp\';
% % A-C ranked principle cells:
% % tfilespath = '\\huxley.resrch.uleth.ca\EP\ppc-seq-lrng\P0662\RandomLights_Sess3_xx\Data\2012-08-21_08-49-01_PPCp\tfiles\';
% % A-C ranked interneurons cells:
% tfilespath = '\\huxley.resrch.uleth.ca\EP\ppc-seq-lrng\P0662\RandomLights_Sess3_xx\Data\2012-08-21_08-49-01_PPCp\TT\interneuron_tfiles';
% sess_info.sess_num = 3.10; %for sequence1

% 8/23/2012
% sess_info.rat = 'P0662';
% datapath = '\\huxley.resrch.uleth.ca\EP\ppc-seq-lrng\P0662\RandomLights_Sess3_xx\Data\2012-08-23_07-51-24_PPCp\';
% tfilespath = '\\huxley.resrch.uleth.ca\EP\ppc-seq-lrng\P0662\RandomLights_Sess3_xx\Data\2012-08-23_07-51-24_PPCp\tfiles\';
% sess_info.sess_num = 3.12; %for sequence1

% 8/24/2012
% sess_info.rat = 'P0662';
% datapath = '\\huxley.resrch.uleth.ca\EP\ppc-seq-lrng\P0662\RandomLights_Sess3_xx\Data\2012-08-24_08-04-48_PPCp\';
% tfilespath = '\\huxley.resrch.uleth.ca\EP\ppc-seq-lrng\P0662\RandomLights_Sess3_xx\Data\2012-08-24_08-04-48_PPCp\tfiles\';
% sess_info.sess_num = 3.13; %for sequence1

% 8/25/2012
% sess_info.rat = 'P0662';
% datapath = '\\huxley.resrch.uleth.ca\EP\ppc-seq-lrng\P0662\RandomLights_Sess3_xx\Data\2012-08-25_07-36-40_PPCp\';
% % A-C ranked principle cells:
% % tfilespath = '\\huxley.resrch.uleth.ca\EP\ppc-seq-lrng\P0662\RandomLights_Sess3_xx\Data\2012-08-25_07-36-40_PPCp\tfiles\';
% % A-C ranked interneurons cells:
% tfilespath = '\\huxley.resrch.uleth.ca\EP\ppc-seq-lrng\P0662\RandomLights_Sess3_xx\Data\2012-08-25_07-36-40_PPCp\TT\interneuron tfiles';
% sess_info.sess_num = 3.14; %for sequence1

% 8/26/2012
% sess_info.rat = 'P0662';
% datapath = '\\huxley.resrch.uleth.ca\EP\ppc-seq-lrng\P0662\RandomLights_Sess3_xx\Data\2012-08-26_09-09-41_PPCp\';
% tfilespath = '\\huxley.resrch.uleth.ca\EP\ppc-seq-lrng\P0662\RandomLights_Sess3_xx\Data\2012-08-26_09-09-41_PPCp\tfiles\';
% sess_info.sess_num = 3.15; %for sequence1

% 8/27/2012
% sess_info.rat = 'P0662';
% datapath = '\\huxley.resrch.uleth.ca\EP\ppc-seq-lrng\P0662\RandomLights_Sess3_xx\Data\2012-08-27_07-20-54_PPCp\';
% tfilespath = '\\huxley.resrch.uleth.ca\EP\ppc-seq-lrng\P0662\RandomLights_Sess3_xx\Data\2012-08-27_07-20-54_PPCp\tfiles\';
% sess_info.sess_num = 3.16; %for sequence1

% 8/29/2012
% sess_info.rat = 'P0662';
% datapath = '\\huxley.resrch.uleth.ca\EP\ppc-seq-lrng\P0662\RandomLights_Sess3_xx\Data\2012-08-29_07-24-22_PPCp\';
% % % A-C ranked principle cells:
% % % tfilespath = '\\huxley.resrch.uleth.ca\EP\ppc-seq-lrng\P0662\RandomLights_Sess3_xx\Data\2012-08-29_07-24-22_PPCp\tfiles\';
% % % A-C ranked interneurons cells:
% tfilespath = '\\huxley.resrch.uleth.ca\EP\ppc-seq-lrng\P0662\RandomLights_Sess3_xx\Data\2012-08-29_07-24-22_PPCp\TT\interneuron\';
% sess_info.sess_num = 3.18; %for sequence1

% 8/30/2012
% sess_info.rat = 'P0662';
% datapath = '\\huxley.resrch.uleth.ca\EP\ppc-seq-lrng\P0662\RandomLights_Sess3_xx\Data\2012-08-30_07-44-27_PPCp\';
% tfilespath = '\\huxley.resrch.uleth.ca\EP\ppc-seq-lrng\P0662\RandomLights_Sess3_xx\Data\2012-08-30_07-44-27_PPCp\tfiles\';
% sess_info.sess_num = 3.19; %for sequence1

% 9/1/2012
% sess_info.rat = 'P0662';
% datapath = '\\huxley.resrch.uleth.ca\EP\ppc-seq-lrng\P0662\RandomLights_Sess3_xx\Data\2012-09-01_07-23-45_PPCp\';
% %  A-C ranked principle cells:
% % % tfilespath = '\\huxley.resrch.uleth.ca\EP\ppc-seq-lrng\P0662\RandomLights_Sess3_xx\Data\2012-09-01_07-23-45_PPCp\tfiles\';
% % % A-C ranked interneurons cells:
% tfilespath = '\\huxley.resrch.uleth.ca\EP\ppc-seq-lrng\P0662\RandomLights_Sess3_xx\Data\2012-09-01_07-23-45_PPCp\TT\interneuron';
% sess_info.sess_num = 3.21; %for sequence1

% 9/3/2012
% sess_info.rat = 'P0662';
% datapath = '\\huxley.resrch.uleth.ca\EP\ppc-seq-lrng\P0662\RandomLights_Sess3_xx\Data\2012-09-03_08-16-51_PPCp\';
%  % % A-C ranked principle cells:
% % % % tfilespath = '\\huxley.resrch.uleth.ca\EP\ppc-seq-lrng\P0662\RandomLights_Sess3_xx\Data\2012-09-03_08-16-51_PPCp\tfiles\';
% % % % A-C ranked interneurons cells:
% tfilespath = '\\huxley.resrch.uleth.ca\EP\ppc-seq-lrng\P0662\RandomLights_Sess3_xx\Data\2012-09-03_08-16-51_PPCp\TT\interneuron';
% sess_info.sess_num = 3.23; %for sequence1

% 9/4/2012
% sess_info.rat = 'P0662';
% datapath = '\\huxley.resrch.uleth.ca\EP\ppc-seq-lrng\P0662\RandomLights_Sess3_xx\Data\2012-09-04_07-27-36_PPCp\';
% %  A-C ranked principle cells:
% % tfilespath = '\\huxley.resrch.uleth.ca\EP\ppc-seq-lrng\P0662\RandomLights_Sess3_xx\Data\2012-09-04_07-27-36_PPCp\tfiles\';
% % A-C ranked interneurons cells:
% tfilespath = '\\huxley.resrch.uleth.ca\EP\ppc-seq-lrng\P0662\RandomLights_Sess3_xx\Data\2012-09-04_07-27-36_PPCp\TT\interneuron\';
% sess_info.sess_num = 3.24; %for sequence1

% 9/5/2012
% sess_info.rat = 'P0662';
% datapath = '\\huxley.resrch.uleth.ca\EP\ppc-seq-lrng\P0662\RandomLights_Sess3_xx\Data\2012-09-05_07-26-10_PPCp\';
% %  A-C ranked principle cells:
% % % tfilespath = '\\huxley.resrch.uleth.ca\EP\ppc-seq-lrng\P0662\RandomLights_Sess3_xx\Data\2012-09-05_07-26-10_PPCp\tfiles\';
% % % A-C ranked interneurons cells:
% tfilespath = '\\huxley.resrch.uleth.ca\EP\ppc-seq-lrng\P0662\RandomLights_Sess3_xx\Data\2012-09-05_07-26-10_PPCp\TT\interneuron\';
% sess_info.sess_num = 3.25; %for sequence1

% P1268
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 5/27/2013
% sess_info.rat = 'P1268';
% datapath = '\\huxley.resrch.uleth.ca\EP\ppc-seq-lrng\P1268\3_RandomLights\2013-05-27_07-51-46_PPCp\';
% tfilespath = '\\huxley.resrch.uleth.ca\EP\ppc-seq-lrng\P1268\3_RandomLights\2013-05-27_07-51-46_PPCp\tfiles';
% sess_info.sess_num = 3.10; %for sequence1

% 5/29/2013
% sess_info.rat = 'P1268';
% datapath = '\\huxley.resrch.uleth.ca\EP\ppc-seq-lrng\P1268\3_RandomLights\2013-05-29_07-34-17_PPCp\';
% tfilespath = '\\huxley.resrch.uleth.ca\EP\ppc-seq-lrng\P1268\3_RandomLights\2013-05-29_07-34-17_PPCp\tfiles';
% sess_info.sess_num = 3.12; %for sequence1

% 6/5/2013
% sess_info.rat = 'P1268';
% datapath = '\\huxley.resrch.uleth.ca\EP\ppc-seq-lrng\P1268\3_RandomLights\2013-06-05_08-37-34_PPCp\';
% tfilespath = '\\huxley.resrch.uleth.ca\EP\ppc-seq-lrng\P1268\3_RandomLights\2013-06-05_08-37-34_PPCp\tfiles';
% sess_info.sess_num = 3.18; %for sequence1

% 6/8/2013
% sess_info.rat = 'P1268';
% datapath = '\\huxley.resrch.uleth.ca\EP\ppc-seq-lrng\P1268\3_RandomLights\2013-06-08_08-25-00_PPCp\';
% tfilespath = '\\huxley.resrch.uleth.ca\EP\ppc-seq-lrng\P1268\3_RandomLights\2013-06-08_08-25-00_PPCp\tfiles';
% sess_info.sess_num = 3.21; %for sequence1

% 6/10/2013
% sess_info.rat = 'P1268';
% datapath = '\\huxley.resrch.uleth.ca\EP\ppc-seq-lrng\P1268\3_RandomLights\2013-06-10_07-52-19_PPCp\';
% %  A-C ranked principle cells:
% % % tfilespath = '\\huxley.resrch.uleth.ca\EP\ppc-seq-lrng\P1268\3_RandomLights\2013-06-10_07-52-19_PPCp\tfiles';
% % % A-C ranked interneurons cells:
% tfilespath = '\\huxley.resrch.uleth.ca\EP\ppc-seq-lrng\P1268\3_RandomLights\2013-06-10_07-52-19_PPCp\TT\interneuron';
% sess_info.sess_num = 3.23; %for sequence1

%Work Computer:
addpath 'C:\Users\aaron.wilber\Documents\MATLAB\Circle Maze Analyses'
addpath 'C:\Users\aaron.wilber\Documents\MATLAB\CircStat2011f'

%Laptop:
% addpath 'C:\Users\Aaron\Documents\MATLAB\circ_stat'
% addpath 'C:\Users\Aaron\Documents\MATLAB\Neural Analysis Scripts\Circle Maze Analyses'

diary off; % in case it is on from a bombed previous run
diary([datapath sess_info.rat 'sess' num2str(sess_info.sess_num), 'run_', strrep(datestr(now,0),':','-'),'.txt']);

%HD OPTIONS:                
HDBinSize = 6;
disp(['Using HD Bin Size of ' int2str(HDBinSize) ' degrees.']);

%Video Tracker .pvd file OPTIONS:
vt_flag = 0;
% 0 = use standard vtm_.pvd files
% 1 = use new vtPDm_.pvd files

%PLOT OPTIONS:
plot_flag = 0;
% 0 = polar plot
% 1 = bar plot

%MAZE LIMITING OPTIONS: 
maze_flag = 1; 
disp(['maze flag = ' int2str(maze_flag)])     
disp('key:')
disp('maze_flag = 0 - display all trials')
disp('maze_flag = 1 - display only maze 1 trials')
disp('maze_flag = 2 - display only maze 2 trials')

split_maze_flag = 2; 
disp(['split maze flag = ' int2str(split_maze_flag)])     
disp('key:')
disp('maze_flag = 0 - display all trials')
disp('maze_flag = 1 - display only first 1/2 of maze session')
disp('maze_flag = 2 - display only second 1/2 of maze session')

min_spikes_cutoff = 500; %500; 
if maze_flag > 0
    min_spikes_cutoff = min_spikes_cutoff/2;
end
disp(['Minimum Spikes Cutoff = ' int2str(min_spikes_cutoff)])     

vel_thresh_flag = 1;
min_still_dur = 60;
disp(['still_maze flag = ' int2str(vel_thresh_flag)])
disp(['minimum still duration (s) = ' int2str(min_still_dur)])
disp('key:')
disp('maze_flag = 0 - display all trials')
disp('maze_flag = 1 - apply velocity threshold')

% Automatically Assign Reader Type Flag
TF0 = strcmp(sess_info.rat, 'P0100');
TF = strcmp(sess_info.rat, 'P0017');
if TF0 == 0 && TF == 0
    ValReaderFlag = 1;
else
    ValReaderFlag = 0;
end
clear TF TF0
disp(['ValReaderFlag = ' int2str(ValReaderFlag)])     
disp('key:')
disp('ValReaderFlag = 0 - use Hendrik software version events decoder')
disp('ValReaderFlag = 1 - use Valery software version events decoder')

if maze_flag == 0
    %because criterion is applied to each maze session separately and goal
    %is at least X, e.g., 500, spikes total:
    min_spikes_cutoff = min_spikes_cutoff/2;
end   

stim_flag = 1;
if stim_flag == 1
    poststim_blackout = 20;  % time in ms after stim ends to continue blocking stim artifact
    disp(['Recoding pv hds for stim duration & using post-stim blackout of ' int2str(poststim_blackout) ' ms.']);
end

%Still working on this:
% if exist([datapath 'TT\Analyses\HD' '*.mat'], 'file')
%     Load([datapath 'TT\Analyses\HD' '*.mat'])
%     disp('Data Cache File Found and Loaded')
%     existingCacheFlag = 1;
% else
%     existingCacheFlag = 0;
% end

if vt_flag == 0
    PVDFile1 = [datapath 'vtm1.pvd'];
    PVDFile2 = [datapath 'vtm2.pvd'];
else
    PVDFile1 = [datapath 'vtPDm1.pvd'];
    PVDFile2 = [datapath 'vtPDm2.pvd'];
end

EventFile = [datapath 'events.txt'];
%Load Event File Data, i.e., maze/sleep epochs
if ~exist('event_data', 'var') || isempty(event_data)
      disp(['loading Events file: ' EventFile]);  
      if sess_info.sess_num == 1000000 
            epoch_times = [813944833 841867903; 930396483 962687168; 962687170 962687570];   
      else
            if exist([datapath 'Events.txt'], 'file')
                if ValReaderFlag == 0;
                    [~, epoch_times] = ReadEventsLethCrcl(datapath); 
                    [seqevents, ~] = read_events4ul([datapath '\Events.txt']);
                else
                    [~, epoch_times] = ReadEventsLethCrclVal(datapath); 
                    [seqevents, ~] = read_events4ulVal([datapath '\Events.txt']);
                end
                   % VT froze so no pre sleep for this session, need to adjust pvd
                   % naming and epoch times list accordingly:
                   TF = strcmp(sess_info.rat, 'P0662');
                   if TF == 1 && sess_info.sess_num == 3.07
                   epoch_times(2:5,1:2) = epoch_times(1:4,1:2);
                   epoch_times(1,1:2) = 0;
                   end
                   clear TF
            end
      end % if sess_info.sess_num == 1000000 
      clear p EventFile
end;

if split_maze_flag > 0
    if maze_flag == 1
    %     half midpoint
        m1_midPoint = epoch_times(2,1)+((epoch_times(2,2) - epoch_times(2,1))/2);
        % first quarter midpoint
    %     m1_midPoint = epoch_times(2,1)+((epoch_times(2,2) - epoch_times(2,1))/4);
%         m1_midPoint = 7291516300;
        if split_maze_flag == 1
            epoch_times(2,2) = m1_midPoint;
        else
            epoch_times(2,1) = m1_midPoint;
        end
    else %assuming it is 2 for now can't think if why to do this for 1 maze session
            %     half midpoint
        m1_midPoint = epoch_times(4,1)+((epoch_times(4,2) - epoch_times(4,1))/2);
        % first quarter midpoint
    %     m1_midPoint = epoch_times(2,1)+((epoch_times(2,2) - epoch_times(2,1))/4);
        if split_maze_flag == 1
            epoch_times(4,2) = m1_midPoint;
        else
            epoch_times(4,1) = m1_midPoint;
        end
    end
end

% assign m1/m2 times variables
m1_times = epoch_times(2,:); 
%only do if maze 2 session exists
if length(epoch_times) > 3;
   m2_times = epoch_times(4,:);
else
   m2_times = [];
end;

if vel_thresh_flag == 1
    disp('Calculating still intervals');
    [prop_sleep, max_int, med_int, sleep_status, sleep_intervals, raw_vel] = no_mov_stats2(PVDFile1, m1_times, min_still_dur);
    if length(epoch_times) > 3 && ~(maze_flag == 1)
       [prop_sleep_m2, max_int_m2, med_int_m2, sleep_status_m2, sleep_intervals_m2, raw_vel_m2] = no_mov_stats2(PVDFile2, m2_times, min_still_dur);
    end
end

clear epoch_times
disp('done');

%Load PVD File Data
if exist(PVDFile1, 'file')
	if ~exist('pv_data', 'var') || isempty(pv_data)
      disp(['loading PVD files: ' PVDFile1 ' & ' PVDFile2]);
      pv_data1 = load(PVDFile1);
      %To fix a bug in Neuralynx HD code related to offsets see emails with
      %Tate
      angles = pv_data1(:,5);
      pv_data1(:,5) = mod(angles,360);
      clear angles
      
      if vel_thresh_flag == 1
          disp('Recoding still intervals as HD450');
          [rows ~] = size(sleep_intervals);
          for i = 1:rows
              exclude = find(pv_data1(:,1) > sleep_intervals(i,1) & pv_data1(:,1) < sleep_intervals(i,2));
              % eventually need to recode this to something else since HD
              % 450 used to be for something else
              pv_data1(exclude,5) = 450;
              clear exclude
          end %for i = 1:rows
          clear rows i
      end % if vel_thresh_flat ==1
      
       if exist(PVDFile2, 'file') && ~(maze_flag == 1)
            pv_data2 = load(PVDFile2);
            %To fix a bug in Neuralynx HD code related to offsets see emails with
            %Tate
            angles = pv_data2(:,5);
            pv_data2(:,5) = mod(angles,360);
            clear angles
            
            if vel_thresh_flag == 1
                  disp('Recoding still intervals as HD450');
                  [rows ~] = size(sleep_intervals_m2);
                  for i = 1:rows
                      exclude = find(pv_data2(:,1) > sleep_intervals_m2(i,1) & pv_data2(:,1) < sleep_intervals_m2(i,2));
                      % eventually need to recode this to something else since HD
                      % 450 used to be for something else
                      pv_data2(exclude,5) = 450;
                      clear exclude
                  end %for i = 1:rows
                  clear rows i
            end % if vel_thresh_flat ==1
      
      else
            disp('No maze two PVD file Found');
      end; %if exist nested
      
      if maze_flag == 1 || maze_flag == 0
           pv_data1 = pv_data1(pv_data1(:,1)>m1_times(1) & pv_data1(:,1) < m1_times(2),:);
           DataChecks.PVsessiontimemin = ((pv_data1(end,1) - pv_data1(1,1))/1e6)/60;
           if pv_data1(1,1) > m1_times(1)
                m1_times(1) = pv_data1(1,1);
                disp('adjusted m1 time 1 to match pv1 time 1');
           end
           if maze_flag == 1
                pv_data = pv_data1;
                disp('limiting to maze 1 pv_data.'); 
           end
      else %must be 2
          pv_data2 = pv_data2(pv_data2(:,1)>m2_times(1) & pv_data2(:,1) < m2_times(2),:);  
          pv_data = pv_data2;
          DataChecks.PVsessiontimemin = ((pv_data(end,1) - pv_data(1,1))/1e6)/60;
          disp('limiting to maze 2 pv_data.');
          if pv_data2(1,1) > m2_times(1)
                m2_times(1) = pv_data2(1,1);
                disp('adjusted m2 time 1 to match pv2 time 1');
          end
      end %if maze flag =1   
      
      if  maze_flag == 0 
          pv_data2 = pv_data2(pv_data2(:,1)>m2_times(1) & pv_data2(:,1) < m2_times(2),:);
          DataChecks.PVsessiontimemin2 = ((pv_data2(end,1) - pv_data2(1,1))/1e6)/60;
          DataChecks.PVsessiontimemin1 = DataChecks.PVsessiontimemin;
          DataChecks.PVsessiontimemin = DataChecks.PVsessiontimemin2 + DataChecks.PVsessiontimemin1;
          pv_data = [pv_data1; pv_data2];
          disp('Using all pv_data.');
          if pv_data2(1,1) > m2_times(1)
                m2_times(1) = pv_data2(1,1);
                disp('adjusted m2 time 1 to match pv2 time 1');
          end
      end % if maze_flag ==0 
      clear PVDFile1 PVDFile2;  
      
      %Remove time in pot and transition from pot to table for this session: 
      % Need to add a rat limiter before turning back on, I think it's
      % P0017
%       if sess_info.sess_num == 2.52
%            min_pv_time = 3970738234;
%            pv_data = pv_data(pv_data(:,1)>min_pv_time, :); 
%            disp(['limiting to after placed on maze for session ' str2num(sess_info.sess_num)])
%       end %if sess_info.sess_num 
      
      disp('done');	
    end;% if~exist
 else
	error(['Unable to find PVD file: ' PVDFile]);
 end; %PVDFile1

% STIM Blockout Time Recoded as HD 1050 in PVD file
if stim_flag == 1
    disp('Recoding HD as 1050 for stim blockout times...')
    % determine stim duration
    d = dir([datapath '\*.INI']);
    if ~isempty(d)
       stim_dur = str2num(read_ini([datapath '\' d(1).name], 'Stim_Param', 'Stimdur'));
    else
       disp('INI file not found, querying user for stim duration');
       stim_dur_cell = inputdlg('Enter Stim Dur (ms):', 'Stim Duration Query', 1, {'350'});
       if isempty(stim_dur_cell)
          disp('No stim duration specified.  Aborting...');
          return;
       end;
       stim_dur = str2num(stim_dur_cell{1});
    end;

    disp(['Stim Duration is: ' num2str(stim_dur) ' ms']);

    % Deterimine the time of the stims
    stim_times= seqevents(seqevents(:,4)==1, 1);
    stim_times(:,2) = stim_times + (stim_dur*1e3) + (poststim_blackout*1e3);

    %set up loop to repeat if both maze sessions should be combined:
    if maze_flag == 0
        sessi_end = 2;
    else
        sessi_end = 1;
    end %if maze flag == 0 

    for sessi = 1:sessi_end
        %load appropriate pv_data
        if sessi == 1 && maze_flag < 2
            pv_data = pv_data1;
        else
            pv_data = pv_data2;
        end %if sessi = 1    
        if maze_flag == 2 
            pv_data = pv_data2;
        end  
        for i = 1:length(stim_times)
            %find pvdata with this matching head direction 
            %Save time by skipping to first appropriate value
            startj = find(pv_data(:,1) >= stim_times(i,1),1);
            if isempty(startj)
                startj = length(pv_data);
            end
            for j = startj:length(pv_data)
                if pv_data(j,1) >= stim_times(i,1) && pv_data(j,1) < stim_times(i,2)
                    pv_data(j,5) = 1050;
                end %if pv_data
            end %for j
        end % for i
        if sessi == 1
            pv_data1 = pv_data;
        else %must be 2
            pv_data2 = pv_data;
        end %if sessi
    end % for sessi
end % if stim_flag ==1
      
%RandNumTest 
% rand_list = randi(360,(length(pv_data1)),1);
% if sess_info.sess_num == 1000000
%     pv_data1(1:end,5) = rand_list;
% end

%Single Pos(s) Test
% if sess_info.sess_num == 1000000
%     pv_data1(1:100,5) = 72;
%     pv_data1(101:200,5) = 210;
%     pv_data1(201:400,5) = 45;
%     pv_data1(401:900,5) = 170;
% end

%LOAD SPIKES
if ~(sess_info.sess_num == 1000000) %b/c test data set doesn't have spikes   
    if ~exist('spike_data', 'var') || isempty(spike_data)
        spike_data = {};
        disp('loadings spikes');
        tfile_list = FindFiles('*.t', 'StartingDirectory', tfilespath);
        ts_data = LoadSpikes(tfile_list);
        keep_celli = true(1,length(ts_data));
        spike_totals.spikes = zeros(1,length(ts_data));
        for i = 1:length(ts_data)
            times = Data(ts_data{i}); % extract data from the silly ts object
            %tsdata times in microseconds presumably because uses 32 bit ts
            %encoding so would run out of space is stored in microseconds
            times = times*100;  % convert to microseconds
            
              if maze_flag == 2 || maze_flag == 0
                times2 = times(times>pv_data2(1,1) & times<pv_data2(end,1));  
                n_m2 = sum(times2>m2_times(1) & times2<m2_times(2));
                time.m2timemin = ((m2_times(2) - m2_times(1))/1e6)/60;
                DataChecks.sessiontimemin = time.m2timemin;
                
                if maze_flag == 0
                    times1 = times(times>pv_data1(1,1) & times<pv_data1(end,1));
                    n_m1 = sum(times1>m1_times(1) & times1<m1_times(2));
                    time.m1timemin = ((m1_times(2) - m1_times(1))/1e6)/60;
                    spike_totals.spikes(i) = n_m1 + n_m2;
                    DataChecks.sessiontimemin = time.m1timemin+time.m2timemin;
                    times = [times1; times2];
                else%maze flag must be 2 so:
                    spike_totals.spikes(i) = n_m2;
                    times = times2;
                    n_m1 = n_m2;
                end %if maze_flag == 0
              else % maze flag must be 1, eventually do something better her:
                times1 = times(times>pv_data1(1,1) & times<pv_data1(end,1));
                n_m1 = sum(times1>m1_times(1) & times1<m1_times(2));
                spike_totals.spikes(i) = n_m1; %need to collapse this with n_m1 at some point
                time.m1timemin = ((m1_times(2) - m1_times(1))/1e6)/60;
                DataChecks.sessiontimemin = time.m1timemin;
                spike_totals.spikes(i) = n_m1;   
                n_m2 = n_m1;
                times = times1;
              end %if maze_flag == 2 || 0

              %Excludes spikes that don't fire enough on the maze
              if n_m1<min_spikes_cutoff && n_m2<min_spikes_cutoff % 
                 keep_celli(i) = 0;  
                 disp(['Exclude: ' tfile_list{i} ...
                            '  ' int2str(spike_totals.spikes(i)) '-' ' spikes (m1 && || m2)']);%eventually make print: int2str(n_m2) when m2 exists, commented out for now
              else
                 spike_data{i} = times;
                 disp(['Load: ' tfile_list{i} ...
                    '  ' int2str(spike_totals.spikes(i)) '-' ' spikes (m1 && || m2)']); %eventually make print: int2str(n_m2) when m2 exists, commented out for now
              end %if n_m1<min_spikes_cutoff && ...
       end; %for i
       spike_data = spike_data(keep_celli);
       tfile_list = tfile_list(keep_celli);
       spike_totals.spikes = spike_totals.spikes(keep_celli);
       clear ts_data n_m1 n_m2 i tfilespath times times1 times2 keep_celli
    end; %if ~exist
end %if ~(sess_info.sess_num == 1000000)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        HEAD DIRECTION ANALYSES                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

disp('Computing head direction occupancy...')
%parse head direction duration
% hdlist - [degrees, number spikes, duration of time head was pointed in this direction]
n_cells = length(spike_data);

%TEST A FAKE CELL THAT "FIRES" ONLY WHEN THE RAT IS FACING ONE DIRECTION
%Designed for data set P0017 4/6:
% disp('Loading fake cell 27...')
% spike_data{27}(:,1) = pv_data(15838:37394,1);
% n_cells = n_cells+1;
% would be nice to calc spike count for fake cell and compare to HD counts

hdlist = zeros(363,(2+n_cells));
hdlist(:,1) = [0:359 450 527 1050];

hdONoffList{363} = zeros;
for i = 1:363
hdONoffList{i}(1,1:2) = zeros;
end %for i
clear i

%set up loop to repeat if both maze sessions should be combined:
if maze_flag == 0
    sessi_end = 2;
else% must be 1 or 2 (i.e., limit 2 1 session)
    sessi_end = 1;
end %if maze flag == 0  
            
%calculate duration spent at each degree of head direction
for i = 1:length(hdlist)   
        % do this for each degree for each session
        for sessi = 1:sessi_end         
            %make sure last data point from first pv_data doesn't get
            %grabbed as an hdonset for the second pv data (thus counting
            %all spike data between maze sessions:
            hdonset = -1;
            hdofftime = -1;
            
            %load appropriate pv_data
            if sessi == 1 && maze_flag < 2
                pv_data = pv_data1;
            else
                pv_data = pv_data2;
            end %if sessi = 1    
            if maze_flag == 2 %if maze 2 data only being run (i.e., sess i 1, maze flag 2):
                pv_data = pv_data2;
            end  
             
            %find pvdata with this matching head direction 
            %Save time by skipping to first appropriate value
            startj = find(pv_data(:,5) == hdlist(i,1),1);
            if isempty(startj)
                startj = length(pv_data);
            end
            
            for j = startj:length(pv_data)  
                if j > startj %don't test for a change until second pv_data hd value
                    if pv_data(j,5) == hdlist(i,1) && ~(pv_data(j-1,5) == hdlist(i,1)) && hdofftime == -1
                        hdonset = pv_data(j,1);
                    end %if
                else
                    % the first data point is unique so grab it because it should 
                    % match the current hd angle, but we'll test to be safe!:
                    if pv_data(j,5) == hdlist(i,1)
                        hdonset = pv_data(j,1);
                    else
                        %There is only an error if we didn't intentionally skip to
                        %the end of the data set because it didn't contain any of
                        %the current heading values
                        if startj < length(pv_data)
                        disp('Error skipped to end of data set but shouldn not have done so!!!')
                        end
                    end %if pv_data...
                end %if j> startj

                % grab time this head direction ended and calculate duration at
                % this head direction
                if j+1 <= length(pv_data)
                    if ~(pv_data(j+1,5) == hdlist(i,1)) && hdonset > -1 ...
                            || j+1 == length(pv_data) && hdonset > -1; 
                        hdofftime = pv_data(j+1,1);
                        hdlist(i,2) = hdlist(i,2) + (hdofftime - hdonset);
                        if hdONoffList{i}(1,1) > 0
                            hdONoffList{i}(length(hdONoffList{i}(:,1))+1,1) = hdonset;
                            hdONoffList{i}((length(hdONoffList{i}(:,2))),2) = hdofftime;
                        else %first row is zeros to start and these need replaced:
                            hdONoffList{i}(1,1) = hdonset;
                            hdONoffList{i}(1,2) = hdofftime;                                
                        end
                        hdonset = -1;
                        hdofftime = -1;
                    end%if~...
                end % j+1 < length(pv_data)...
            end %for j
        end %for sessi = 1:2        
end % for i
clear i j startj hdonset hdofftime sessi sessi_end

disp('Counting spikes for each head direction...')
%Create Histogram of FR by HD
for cell_num = 1:n_cells
%     if cell_num < 27
        [p, name] = fileparts(tfile_list{cell_num});
        newname = name;
        newname(newname=='_') = '-';
%     else
%         newname = 'fakeCell27';
%     end
    disp(['CELL NUMBER ' newname]);
    cur_spike_data = spike_data{cell_num}(:,1);
    
    % Something would like this would be much faster, but I need to think about
    % what the third term should be:
    % cur_spike_dist = interp1(pv_data(:,1), cur_spike_data, missing third term);
    % cur_spike_dist = spike_dist{cell_num};

      %This loop should calculate the number of spikes per degree of head
      %direction occupied then divide by the amount of time the rat occupied
      %that head direction
      for j = 1:length(cur_spike_data)
          for k = 1:length(hdONoffList)
              for l = 1:length(hdONoffList{k}(:,1))
                    if cur_spike_data(j) >= hdONoffList{k}(l,1) && hdONoffList{k}(l,1) > 0 ...
                            && cur_spike_data(j) < hdONoffList{k}(l,2) && hdONoffList{k}(l,2) > 0
                        % if timestamps in the range then look up the hd that
                        % corresponds to that timestamp (i.e., k), and increase 
                        % the number of spikes for that head direction by 1,
                        % one column for each cell, first cell in column 3
                        hdlist(k,(2 + cell_num)) = hdlist(k,(2 + cell_num))+1;
                    end %if  
              end %l
          end % for k
      end % for j
    spike_totals.HD(cell_num) = sum(hdlist(:,cell_num+2));
    disp(['Make sure total spikes for HD counts: ' int2str(spike_totals.HD(cell_num)) ...
' matches the total number spikes for that cell: ' int2str(spike_totals.spikes(cell_num))]);
end %for n_cells
clear cell_num j k l num_spikes

disp('Binning Data...')
%Bin HD Data:
% -1 then plus 1 for length because last row is for HD 450 Neuralynx error
% and this value must be removed for the # rows calculation then put back in
hdBinned = zeros(((length(hdlist)-3)/HDBinSize)+3,(n_cells+2));
% length-1 because last row needs to be treated differently b/c HD450 coding, see below:
hdBinned(1:end-3,1) = HDBinSize:HDBinSize:360;
hdBinned(end-2,1) = 450;
hdBinned(end-1,1) = 527;
hdBinned(end,1) = 1050;
[m ~] = size(hdBinned);
for i = 1:(m-3)
    if i > 1
        % sum across binwidth and convert to microseconds to seconds
        hdBinned(i,2) = sum(hdlist((hdBinned(i-1,1)+1):hdBinned(i,1),2))/1e6;
        for j = 1:n_cells
            %for each cell count spikes across bin and divide by seconds to
            %convert to firing rate in spikes/s
             hdBinned(i,2+j) = sum(hdlist((hdBinned(i-1,1)+1):hdBinned(i,1),2+j))/hdBinned(i,2);
        end % for j  
        clear j
    else % for first bin the formula is different, need to skip zero bin which 
        % is contaminated by lost HD, note first bin is a Binsize-1!
        % should be ok b/c FR but make sure I keep this in mind when using
        % this data for additional analyses:
        hdBinned(1,2) = ((sum(hdlist(2:HDBinSize,2)))/1e6);
        for j = 1:n_cells
            hdBinned(1,2+j) = sum(hdlist(2:HDBinSize,2+j))/hdBinned(1,2);
        end % for j  
        clear j
    end %if i    
end    

% Because last bin is coding HD450, i.e., a single "degree" and this needs to be treated separately:  
hdBinned(i+1:i+3,1) = hdlist(end-2:end,1); % copy label
hdBinned(i+1,2) = hdlist(end-2,2)/1e6; %convert duration to sec
hdBinned(i+2,2) = hdlist(end-1,2)/1e6; %convert duration to sec
hdBinned(i+3,2) = hdlist(end,2)/1e6; %convert duration to sec
HDBinnedFR = zeros(n_cells,1);
HDListFR = zeros(n_cells,1);
for j = 1:n_cells
    %divide # spikes in this category by seconds to convert to firing rate in spikes/s
     hdBinned(i+1,2+j) = hdlist(end-2,2+j)/hdBinned(i+1,2);
     hdBinned(i+2,2+j) = hdlist(end-1,2+j)/hdBinned(i+2,2);
     hdBinned(i+3,2+j) = hdlist(end,2+j)/hdBinned(i+3,2);
     
     %computing average firing rate for cells for HDlist and HDBinned and make sure they agree:
     %don't include 450:
     HDBinnedSpikes_temp = hdBinned(1:end-3,2+j).*hdBinned(1:end-3,2);
     HDBinnedFR(j) = sum(HDBinnedSpikes_temp(1:end));
     clear HDBinnedSpikes_temp
     %don't include 0 or 450 AND skip 1st bin (binned) :
     HDListFR(j) = sum(hdlist(2:end-3,2+j));
     disp(['Make sure #spikes matches for HDlist (by degree): ' int2str( HDListFR(j)) ...
     ' vs HDbinned data: ' int2str(HDBinnedFR(j))]);
end % for j  
clear i j

%Display Some Data Checks
if sess_info.sess_num < 1000000
    DataChecks.hdONoffList.CumTimeMicroSec = 0;
    for i=1:length(hdONoffList)
        [row ~] = size(hdONoffList{1,i});
        for j = 1:row
            DataChecks.hdONoffList.CumTimeMicroSec = DataChecks.hdONoffList.CumTimeMicroSec ...
                + (hdONoffList{1,i}(j,2) - hdONoffList{1,i}(j,1));
        end
        clear row
    end
    DataChecks.hdONoffList.CumTimeMin = (DataChecks.hdONoffList.CumTimeMicroSec/1e6)/60;
    DataChecks.hdlist.CumOccupancySec = (sum(hdlist(2:end,2)))/1e6;
    DataChecks.hdlist.CumOccupancySecInclZero = sum(hdlist(1:end,2))/1e6;
    DataChecks.hdlist.CumOccupancyMin = DataChecks.hdlist.CumOccupancySec/60;
    %no zero to skip for hdbinned
    DataChecks.hdbinned.CumOccupancySec = sum(hdBinned(1:end,2));
    DataChecks.hdbinned.CumOccupancyMin = DataChecks.hdbinned.CumOccupancySec/60;
    DataChecks.binVsList.CumOccupDifSec = DataChecks.hdlist.CumOccupancySec - DataChecks.hdbinned.CumOccupancySec;
    %NEED TO ADD 1050 HERE
    DataChecks.UnaccountedSec.FiveTwentySeven = (hdlist(end-1,2)/1e6);
    DataChecks.UnaccountedSec.fourFifty = (hdlist(end-2,2)/1e6);
    DataChecks.UnaccountedSec.Zero = hdlist(1,2)/1e6;
    DataChecks.UnaccountedSec.Total = DataChecks.UnaccountedSec.Zero + DataChecks.UnaccountedSec.fourFifty + DataChecks.UnaccountedSec.FiveTwentySeven;
    DataChecks.PVsessiontimeminMinZero = DataChecks.PVsessiontimemin - (DataChecks.UnaccountedSec.Zero/60);
    DataChecks.UnaccountedMin.PropHDlost = (DataChecks.UnaccountedSec.Zero/DataChecks.hdlist.CumOccupancySecInclZero)*100;
    disp(['Make sure total session times including possible HD 450s & 527s (but not HD 0s) match: ' int2str(DataChecks.hdlist.CumOccupancyMin) ...
    ' for HD by degree list vs ' int2str(DataChecks.hdbinned.CumOccupancyMin) ...
    ' for HDbinned vs ' int2str(DataChecks.PVsessiontimeminMinZero) ' in pv_data.']);
    disp(['Head Direction Lost: ' int2str(DataChecks.UnaccountedMin.PropHDlost) ' Percent of time.']);
    disp(['Head Direction coded as 450 for ' int2str(DataChecks.UnaccountedSec.fourFifty) ' seconds.']);    
end %sess_info.sess_num < 1000000

if sess_info.sess_num == 1000000
    Extracted.North = sum(hdlist(225:314,2));
    Extracted.East = sum(hdlist(1:44,2))+ sum(hdlist(63:72,2));
    Extracted.South = sum(hdlist(45:135,2));
    Extracted.West = sum(hdlist(136:224,2));
    
    DifferenceSec.North = (EventFileDur.North - Extracted.North)/1e6;
    DifferenceSec.East = (EventFileDur.East- Extracted.East)/1e6;
    DifferenceSec.South = (EventFileDur.South - Extracted.South)/1e6;
    DifferenceSec.West = (EventFileDur.West - Extracted.West)/1e6;

    time450sec = hdlist(361,2)/1e6;
    EventFileDur.Total = (962687169-930396482);
    EventFileDur.Unaccounted = EventFileDur.Total - (EventFileDur.East+EventFileDur.North+EventFileDur.South+EventFileDur.West);
    UnaccountedSec.EventFileDur = (EventFileDur.Total - (EventFileDur.East+EventFileDur.North+EventFileDur.South+EventFileDur.West))/1e6;
    UnaccountedSec.Extracted = (DifferenceSec.North + DifferenceSec.South + DifferenceSec.East + DifferenceSec.West);

    %check binned too, I didn't before and there was a mistake there!!! 
    %NOTE 0 doesn't get plotted on graph so this bin looks empty on the
    %plot!  Looks like hd was lost/0 during much of West Time.
    ExtractedBin.North = sum(hdBinned(45:63,2));
    ExtractedBin.East = sum(hdBinned(1:9,2))+ sum(hdBinned(63:72,2));
    ExtractedBin.South = sum(hdBinned(10:27,2));
    ExtractedBin.West = sum(hdBinned(28:44,2));
    ExtractedBin.Zero = hdlist(1,2)/1e6;
end %sess_info.sess_num == 1000000

% close any figures open from last run        
close all

%Print Figures to .jpeg file:
fh1 = gcf;
bar(hdBinned(1:end-3,1),hdBinned(1:end-3,2))
box off
% Create xlabel
xlabel('Heading Degrees','FontSize',30);
% Create ylabel
ylabel('Seconds','FontSize',30);
title('Heading Occupancy','FontSize',30);
xlim([0 365]);
print(fh1, '-djpeg', 'hdOccupancy');

if plot_flag == 1 % Create bar HD figure
    fh2 = figure;
    for i = 1:n_cells 
        %get neuron name
    %     if cell_num < 27
            [~, name] = fileparts(tfile_list{n_cells});
            newname = [name 'bar' HDBinSize];
            newname(newname=='_') = '-';
    %     else
    %         newname = 'fakeCell27';
    %     end
        bar(hdBinned(1:end-3,1),hdBinned(1:end-3,i+2))
        % Create xlabel
        xlabel('Head Direction (Degrees)','FontSize',30);
        % Create ylabel
        ylabel('Firing Rate (Spikes/s)','FontSize',30);
        title(newname,'FontSize',30);
        xlim([0 365]);
        print(fh2, '-djpeg', newname);

    end 
    clear i newname fh2 name p
end %if plot flag == 1

if plot_flag == 0% Create polar HD figure
    fh2 = figure;
    for i = 1:n_cells 
    %     get neuron name
        [~, name] = fileparts(tfile_list{i});
        newname = [name 'polar' num2str(HDBinSize)];
        newname(newname=='_') = '-';  
        theta = circ_ang2rad(hdBinned(1:end-3,1));
        theta(end+1) = hdBinned(1,1);
        rho = hdBinned(1:end-3,i+2);
        rho(end+1) = hdBinned(1,i+2);
        polar(theta,rho)
        clear theta rho
        print(fh2, '-djpeg', newname);    
    end    
    clear i newname fh2 name
end %if plot_flat == 0

if sess_info.sess_num < 1000000
    filename = [datapath 'TT\Analyses\HD\M2ndHalf_CacheIntrnrnPCsBin6' num2str(sess_info.rat) 'sess' num2str(sess_info.sess_num, '%03d') '.mat'];
    disp(['saving data to file: ' filename]);
    eval(['save ' filename ' tfile_list sess_info DataChecks HDBinSize hdlist hdBinned spike_data time m1_times m2_times pv_data']);
end %sess_info.sess_num < 1000000

HDVec

diary off; 