function data=targetDensityPlot(varargin)
% Scrapes target search results from Allen Brain Atlas website, creates
% data structure and plots projection density. 
%
%
% By L.Berkowitz & R.Harvey 2018

%Inputs: 
%       injection_structure: Choose abbreviation as provided by the Allen Brain Atlas naming
%       structure. default('grey') will search all structures.
%       transgenic_lines: default('0') is WT
%       target: structure of interest
%       num_inputs: maximum number of top projection regions by density
%Output: 
%       data: data structure from allen brain webstite :) 
%
%
p = inputParser;
p.addParameter('injection_structure','grey');
p.addParameter('transgenic_lines', '0');
p.addParameter('target', 'POST');
p.addParameter('num_inputs',20);

p.parse(varargin{:});

injection_structure = p.Results.injection_structure;
transgenic_lines = p.Results.transgenic_lines;
target = p.Results.target;
num_inputs=p.Results.num_inputs;


%% TARGET SEARCH ALL GREY MATTER INPUTS TO POST

options = weboptions('Timeout',60);

query=['http://api.brain-map.org/api/v2/data/query.json?criteria=service::mouse_connectivity_injection_structure[injection_structures$eq',...
    injection_structure,'][transgenic_lines$eq',transgenic_lines,'][primary_structure_only$eqtrue][target_domain$eq',target,']'];
data=webread(query,options);

if data.success==0
    error('No viable experiments found. Check abbreviation against abbreviation provided by Allen Brain Atlas')
end

abv=unique({data.msg.structure_abbrev},'stable'); %find all unique inputs and keeps order (i.e. highest to lowest density)
abv=abv(1:num_inputs); %keep just the top number of inputs

%Create figure
fig=figure;fig.Color=[1 1 1];
subplot(2,1,1)
for i=1:length(abv)
   loc=strfind({data.msg.structure_abbrev},abv{i});
   loc=find(~cellfun(@isempty,loc));
   
   plot(zeros(length(loc),1)+i,[data.msg(loc).sum],'.k','MarkerSize',20);hold on
   
   mean_tgt_vol(i)=mean([data.msg(loc).sum]);
end
set(gca,'XTick',1:size(abv,2), 'XTickLabel',abv)
ylabel('Avg Tgt Vol mm^3')
xlabel('Injection Sites')
grid on
title(['Projections to ',target])
xlim([.5 num_inputs+.5])

subplot(2,1,2)
imagesc(mean_tgt_vol)
colormap hot
axis equal tight
set(gca,'YTick',1, 'YTickLabel',[],'XTick',1:size(abv,2), 'XTickLabel',abv)
title('Avg Tgt Vol mm^3')
xlabel('Injection Sites')
