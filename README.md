# ephys_tools

## Contributors

<!-- ALL-CONTRIBUTORS-LIST:START - Do not remove or modify this section -->
<!-- prettier-ignore-start -->
<!-- markdownlint-disable -->
<table>
<tr>     
     <td align="center"><a href="https://github.com/ryanharvey1"><img src="https://avatars3.githubusercontent.com/u/16714674?s=400&u=08b04d641912e453b02db42f08b6e5e2bb353166&v=4" width="100px;" alt=""/><br /><sub><b>Ryan Harvey</b></sub></a><br /><a href="https://github.com/ryanharvey1"</td>
    <td align="center"><a href="https://github.com/lolaBerkowitz"><img src="https://avatars2.githubusercontent.com/u/31447162?s=400&u=cf4002e1ac6d7642fde99e19b4abbbd56463d1c6&v=4" width="100px;" alt=""/><br /><sub><b>Laura Berkowitz</b></sub></a><br /><a href="https://github.com/lolaBerkowitz"</td>
<tr> 
</table>    
       
University of New Mexico, Psychology Department

This repo contains code used to process and analyze our in-vivo electrophysiology tetrode data as well as code to analyze data from a few behavioral tasks. Everything is under development and many functions are not general or polished, so if you find something that needs changing, please open an issue.  

***

#### Getting started with new project

1. Open a terminal and clone or download the repository  

   ```
   git clone https://github.com/ryanharvey1/ephys_tools.git
   ```

2. Open matlab and navigate to the repository 

3. Run the start up script to add necessary scripts and functions to your path

   ```matlab
   startup 
   ```

4. In matlab, cd to where you want your project [raw and processed data] to live. 

5. run `create_project_folder.m`, *located within "preprocessing"*, to create all the required directories. replace 'your project name' with what you want your project to be called.  

   ```matlab
   create_project_folder('your project name')
   ```

6. Populate the newly created folders with any data you have already collected. Notes on how to do this are in create_project_folder.m's documentation at the top of the function. 

Currently, ephys_tools supports Neuralynx data and Open Ephys formats (more formats can be integrated in the future). However, it does support several spike sorters ([Spike Sort 3D](https://neuralynx.com/software/spikesort-3d), [MClust](https://github.com/adredish/MClust-Spike-Sorting-Toolbox), & [Kilosort2](https://github.com/MouseLand/Kilosort2) / [Phy](https://github.com/cortex-lab/phy))

***

#### Spike sorting & processing your data

ephys_tools has functions to convert your spike sorted data into a single .mat format that lives within a "Sorted" subfolder in your raw data folder. You can spike sort your Neuralynx data in Spike Sort 3D, MClust, or Kilosort2 / Phy. 

##### Kilosort2 / Phy

Dependencies:

* [Kilosort2](https://github.com/MouseLand/Kilosort2) 
* [sortingQuality](https://github.com/cortex-lab/sortingQuality)
* [spikes](https://github.com/cortex-lab/spikes)
* [npy-matlab](https://github.com/kwikteam/npy-matlab)

1. To convert your Neuralynx .ncs channels into a single .dat file to run through Kilosort, cd to your current session folder (it will look something like 2019-08-15_19-25-40 and should be located within subfolders: project name > data > animal ID) and run `ncs2dat`

2. At the same time, it's a good idea to run `filter_raw_dat.m` which high pass filters your channels for later average waveform extraction. 

3. Run `get_lfp.m` to make .xml and .lfp files. Open up your .lfp or .dat file in neuroscope to mark bad channels. 

4. Run Kilosort and manually adjust your units with Phy. (If you created the .xml file from step 3, you can run `KS2Wrapper.m` to run kilosort2)

5. Now that you have saved and are happy with your Phy results, cd to your current session folder and run `after_spikesort_cleanup.main`. This will create a 'Sorted' folder that contains your your spike time stamps and quality metrics. 

6. Now is a good time to set up animal metadata. run `handleAnimalMetaData` and answer a couple questions in the command window. Animal metadata stores info like strain, birth date, surgery info, manipulation info, recording sites, electrode turn record, recording session info, etc.

7. cd back to your current session folder and run `postprocess`. This will compile spikes, lfp, and behavior, calculate several features, and save them to a processed data .m file in the 'ProcessedData' folder. These processed data files can then be used for much or all future analyses. Below shows the contents of a processed data file. 

   ```matlab
   data = load('F:\Projects\HPCatn\ProcessedData\HPCatn05_S20181229143332.mat')
   
   data = 
     struct with fields:
                BasicLoco: [1×1 struct]
                   Spikes: {53×1 cell}
                ThPrecess: {53×3 cell}
                  avgwave: {53×1 cell}
               binned_vel: {[3758×1 double]  [6319×1 double]}
           date_processed: '13-Feb-2020'
                   events: [2×2 double]
                   frames: [70891×5 double]
                 hdTuning: {53×3 cell}
       hdTuning_corrected: {53×3 cell}
                      ifr: {53×2 cell}
                      lfp: [1×1 struct]
             linear_track: {[1×1 struct]}
             maze_size_cm: [360 100]
                mazetypes: {'circ track'  'box'}
                 measures: [53×66×3 double]
                openfield: {[]  [1×1 struct]}
                      rat: 'HPCatn05'
                  ratemap: {53×3 cell}
               samplerate: 30
                sessionID: 'S20181229143332'
         session_duration: [12.5156 21.0444]
             session_path: 'F:\Projects\HPCatn\data\HPCatn05\2018-12-29_14-33-32'
                 spikesID: [1×1 struct]
            thetaautocorr: {53×3 cell}
             ts_timescale: 'sec'
                 varnames: {1×66 cell}
   ```

    

##### Spike Sort 3D & MClust

The Spike sort 3d & MClust work flows are very much the same as the steps above starting at step 4. For specific differences in processing spike sort 3d, see documentation in `after_spikesort_cleanup.m` 

ephys_tools contains a forked version of MClust-4.4 with dark mode and a few other added features.



### Acknowledgements

ephys_tools structure and workflow was largely inspired by [CMBHOME]( https://github.com/hasselmonians/CMBHOME) & [buzcode](https://github.com/buzsakilab/buzcode) and relies on many other packages located within the 'external_packages' folder. 

