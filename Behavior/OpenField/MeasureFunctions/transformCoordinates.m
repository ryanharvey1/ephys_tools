function [ trackData_x,trackData_y] = transformCoordinates(dia,x_max,x_min,y_max,y_min,trackData_x,trackData_y) 
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here

    %Find Origin
%     x_Origin= x_max-((x_max-x_min)/2);
%     y_Origin= y_max-((y_max-y_min)/2);

    x_Origin=median(x_min:x_max); 
    y_Origin=median(y_min:y_max);
      
    %% Convert tracking data to cm
    trackData_x=((trackData_x-(x_Origin))*(dia/(x_max-x_min)));
    trackData_y=((trackData_y-(y_Origin))*(dia/(y_max-y_min)));
end

