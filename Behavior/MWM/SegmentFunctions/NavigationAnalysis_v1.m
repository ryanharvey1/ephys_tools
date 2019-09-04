%Navigational Strategy Analysis 

%Within Trajectory iteration
%Analyze Looping
%Analyze Whole Trajectory 
    %Search Error
    %Heading Error
    %Time in goal corridor
        %Determine if Direct Trajectory or Directed Path 
        %if Search Error = <30m*s average search error & <20deg absolute
        %heading error 
            %trajectory is direct path; 
            %run segmentation measures
            %continue
        %elseif time in goal corridor is >80% 
            %trajectory is directed path; 
            %run segmentation measures
            %continue
        %elseif segment path through velocity thresholding
            %for segment=1:size(segData.sbj.trial.segments,2)
                %time in wall zone
                %time in wider wall zone
                %mean distance to plat
                %path centroid to plat
                %eccentricity
                %time in plat radius (6x platform radius diameter)
                %Detect loop 
                %Determine Strategy 
                    %Thigmotaxis >35%time in closer wallzone and/or
                    %>65%time in wider wall zone
                    %Chaining >80% time in annulus zone
                    %Focal Search <.35 mean dist. to swim path centroid &
                    %or <.3 mean distand to present goal 
                    %Self orienting: Segment contains a loop 
                    %Perseverance <.35 mean dist. to swim path centroid &
                    %or <.3 mean distand to prior goal (matching to place
                    %only) 