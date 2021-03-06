%The purpose of this script is to pull connectivity data from the Allen
%Brain Institute website. 

%Dependencies: 
%   jsonlab � 1.5 toolbox
%   getAllenStructureList.m
%   getProjectionDataFromExperiment.m


%Pull List of Structure Names
S=getAllenStructureList; %loads table of structure name data

%Find Abbreviations of structures of interest
NameIdx=contains(S{:,4},'POST'); %String is structure of interest
abv=S{NameIdx,3:4};

%Display areas identified
if ~isempty(abv)
    disp('Structures identified :')
    for i=1:size(abv,1)
%         str = sprintf(" %s, ",abv{:,1})
        disp([num2str(i),' ',abv{i,:}])
    end
    x=input('Omit regions. Enter the number corresponding to the region(s) to omit in a column vector (i.e. [3;4;6]) ');
    abv(x,:)=[];
else
    error('Error. Name not identified on structure list. Please peruse S and try another keyword')
end

clear NameIdx 

%Loop to pull data for each structure of interest. 

for i=1:size(abv,1)
    
%Pull IDs of target Experiments
ID = findAllenExperiments('injection',abv{i,2}); %output is row vector of Experiment IDs

if isempty(ID)
    continue %go to next iteration if there are no experiments for keyword
end
%Pull Data from ID list
    for ii=1:size(ID',1)
        %Pull Data and store in structure
        Data.(abv{i,2}).(['exp_',num2str(ID(1,ii))]) = getProjectionDataFromExperiment(ID(1,ii));   
    end

end


