
%%Variables needing pre-set values
fr=30; %frame rate 
binsize=3; %to use for HB detection


load('params.rawcoords.mat') %%loads table containing the raw coords and xMax/Min and yMax/Min

params.transcoords{i}=[];
%%Converting the Raw Coordinates to import into the table
for i = 1:size(params.PathName,1)-1 %%size gives 2 values, # of rows and columns, by subtracting 1 we keep it on the interval 1:296
 
        dia=params.dia(i); %%establishes the variables for the function "transformed coordinates" while looping through all 296 files
        xMax= params.Xmax(i);
        xMin= params.Xmin(i);
        yMax= params.Ymax(i);
        yMin= params.Ymin(i);
        trackData= params.rawcoords{i};
        trackData_x= trackData(:,1);
        trackData_y= trackData(:,2);
        
   
 [params.transcoords{i}(:,1),params.transcoords{i}(:,2)] = transformCoordinates(dia,xMax,xMin,yMax,yMin,trackData_x,trackData_y);
         
end


%%ComputeOF smooths out raw coordinates
param_idx=params.PathName; %serves as index for getParam
params.smoothcoords{j}=[];
for j=1:size(param_idx,1)
    close all
  %% Maze Parameters;
     %Velocity Filter
        sqrXDiffVel=(diff(params.transcoords{j,1}(:,1))).^2; %replace dsDAta(:,2) with ray X
        sqrYDiffVel=(diff(params.transcoords{j,1}(:,2))).^2; % """" raw Y
        pathDist=sqrt((sqrXDiffVel + sqrYDiffVel));
        pathVel=(pathDist)*fr;
        dsData=params.transcoords{j}(pathVel<100,:); %remove data points greater than 100 cm/s outliers due to tracking error
        
        %Smooth Data
        %pad ends with last data point
        keepLength=length(dsData(:,1));
        padsize=15;
        endPadding=repmat([dsData(end,1) dsData(end,2)],padsize,1);
        frontPadding=repmat([dsData(1,1) dsData(1,2)],padsize,1);
        padDsData=[frontPadding; dsData(:,1) dsData(:,2); endPadding];
        
        %smooth Data
        addpath(genpath('F:\Users\BClarkLab\GoogleDrive\MatlabDir\chronux_2_11'))
        x_smooth=runline(padDsData(1:end,1),5,1); %Smooth data
        y_smooth=runline(padDsData(1:end,2),5,1); %Smooth data
        
        %remove padding
        params.smoothcoords{j}(:,1)=x_smooth(padsize+1:keepLength+(padsize-1),1); %size of padding
        params.smoothcoords{j}(:,2)=y_smooth(padsize+1:keepLength+(padsize-1),1);
        
        
        % compute total distance and instantaneous velocity
        [ params.pathL{j},params.pathIV{j}] = compute_pathCalc(x_smooth,y_smooth,fr);
  
        
        %calculate verticies for quadrants
        
        quadrants=createZones([0,0],params.dia{j},'numQuad',16);
        
        %Calculate verticies for wall area (%80)
        outsideWall=createZones([0,0],params.dia{j},'type','annulus');

    
    
        %Calculate dwell time per quadrant
        for i=1:size(quadrants,2)-1
            tempXin=[0 quadrants(i,1) quadrants(i+1,1)];
            tempYin=[0 quadrants(i,2) quadrants(i+1,2)];
            [in,~]=inpolygon(params.smoothcoords{j,1}(:,1),params.smoothcoords{j,1}(:,2),tempXin,tempYin);
            params.dwellQuad{j}=sum(in)/fr;
        end
        
        %Bin coordinates for determining home base & making heat maps
        xMax=params.Xmax(j); xMin=params.Xmin(j); yMax=params.Ymax(j); yMin=params.Ymin(j);
        
        xedge=linspace(xMin,xMax,round(params.dia(j)/binsize));
        yedge=linspace(yMin,yMax,round(params.dia(j)/binsize));
        
        [map] = histcounts2(params.smoothcoords{j,1}(:,1),params.smoothcoords{j,1}(:,2),xedge,yedge);
        map=flipud(imrotate(map,90));
        map=map/fr;
        
        addpath(genpath('F:\Users\BClarkLab\GoogleDrive\MatlabDir\BClarkToolbox\Analysis\Utils'))
        [BW,maskedImage,x,y,fieldarea,X] = segmentImage('maps',map);
        
        %% save fig
        saveas(gcf,param_idx{j},'jpeg'); 
end


% %Index out of Params to plot transformed coords and analyze
%%Loops through each time point for graphing and analysis
% for i = 1:size(file,1)-1
%     if regexp(params.PathName(i),'4[m]+o');
% %         Identifiers(i,1)= params.PathName(i);
% figure; plot(cos(0:2*pi/1000:2*pi)*39,sin(0:2*pi/1000:2*pi)*39,'-k');
%         hold on; plot(params.transcoords{i}(:,1),params.transcoords{i}(:,2),'-r');
%         axis image
%         pause(0.5)
%         close
%     elseif  regexp(params.PathName(i),'7[m]+o');
% %         Identifiers(i,2)= params.PathName(i);
%     elseif  regexp(params.PathName(i),'2[m]+o');
% %         Identifiers(i,3)= params.PathName(i);
%     elseif regexp(params.PathName(i),'\lgOF*');;
% %         Identifiers(i,4)= params.PathName(i);
%         figure; plot(cos(0:2*pi/1000:2*pi)*101,sin(0:2*pi/1000:2*pi)*101,'-k');
%         hold on; plot(params.transcoords{i}(:,1),params.transcoords{i}(:,2),'-r');
%         axis image
%         pause(0.5)
%         close
%     end
% end



%% Creates the zones and annuli for the dwell time
% for i = 1:size(file,1)-1
%   
%     third_radius = (params.dia(i)/2)/3; %%calculates a third of the radius to use to plot the snaller circle zones in line__
%     x0=0;%%origin
%     y0=0;%%origin
%     r=params.dia(i)/2;%%radius
%     teta=-pi:0.01:pi;
%     x=r*cos(teta)+x0;
%     y=r*sin(teta)+y0;
%     figure
%     plot(x,y)
%     hold on
%     scatter(x0,y0,'or')
%     axis square 
% %     ----------------------------------------
% %     divide your circle to n sectors
%     n=16;
%     tet=linspace(-pi,pi,n+1); %%allows for even space between each line??
%     xi=r*cos(tet)+x0;
%     yi=r*sin(tet)+y0;
%     
%         for k=1:numel(xi)
%         plot([x0 xi(k)],[y0 yi(k)])
%         hold on
%         end
%     hold on 
%     plot((sin(0:2*pi/1000:2*pi)*((params.dia(i)/2)-third_radius)),cos(0:2*pi/1000:2*pi)*((params.dia(i)/2)-third_radius))
%     plot((sin(0:2*pi/1000:2*pi)*((params.dia(i)/2)-(third_radius*2))),cos(0:2*pi/1000:2*pi)*((params.dia(i)/2)-(third_radius*2)))
%     hold on
%     plot((params.transcoords{i,1}(:,1)), params.transcoords{i,1}(:,2))
%     axis image
%     
% end 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%CODE GRAVEYARD%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%         
%     elseif isequal('7mo',(params.PathName))
%         dia=params(:,3);
%         xMax=params(:,4);
%         xMin=params(:,5);
%         yMax=params(:,6);
%         yMin=params(:,7);
%         trackData_x= params.rawcoords{:,1};
%         trackData_y= params.rawcoords{:,2};
%         
%     elseif isequal('10mo',(params.PathName))
%         dia=params(:,3);
%         xMax=params(:,4);
%         xMin=params(:,5);
%         yMax=params(:,6);
%         yMin=params(:,7);
%         trackData_x= params.rawcoords{:,1};
%         trackData_y= params.rawcoords{:,2};
%         
%     elseif isequal('lgOF',(params.PathName))
%         dia=params(:,3);
%         xMax=params(:,4);
%         xMin=params(:,5);
%         yMax=params(:,6);
%         yMin=params(:,7);
%         trackData_x= params.rawcoords{:,1};
%         trackData_y= params.rawcoords{:,2};
%         
% 
% strfind(params.PathName(i),'4mo','ForceCellOutput',true);

% 
%     if regexpi(params.PathName(i),'\4*');
%         Identifiers(i,1)= params.PathName(i);
%     elseif regexp(params.PathName(i),'\7*');;
%         Identifiers(i,2)= params.PathName(i);
%     elseif regexp(params.PathName(i),'\10*');
%         Identifiers(i,3)= params.PathName(i);
%     elseif regexp(params.PathName(i),'\lgOF*');;
%         Identifiers(i,4)= params.PathName(i);
%         
%     end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%KEEP%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % %%Checks to see if the script is working, puts each timepoint into
%               a different column 
% % % % if regexp(params.PathName(i),'4[m]+o');
% % % %         Identifiers(i,1)= params.PathName(i);
% % % %     elseif  regexp(params.PathName(i),'7[m]+o');
% % % %         Identifiers(i,2)= params.PathName(i);
% % % %     elseif  regexp(params.PathName(i),'2[m]+o');
% % % %         Identifiers(i,3)= params.PathName(i);
% % % %     elseif regexp(params.PathName(i),'\lgOF*');;
% % % %         Identifiers(i,4)= params.PathName(i);
% % % %         
% % % %     end





%%%%%%%%%%%%Renames OF_Excel paths to match lab computer 1
% % % % % %          temp=strcat(dataDir,names{i});  
% % % % % %     
% % % % % %          params.PathName{i}=temp;
