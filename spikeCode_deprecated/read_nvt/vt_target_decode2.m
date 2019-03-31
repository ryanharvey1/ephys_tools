%Decodes a Neuralynx VT record target value
%Input:
%	target_value:  the 32 bit integer representation of the encoded target.
%Output:
%	red: 1 if target is a red color target
%	green: 1 if target is a green color target
%	x: x coordinate of the center of mass of the target
%	y: y coordinate of the center of mass of the target
function [red, green, x, y] = vt_target_decode2(target_value)  

binary_target_string = dec2bin(target_value, 32);
red = bin2dec(binary_target_string(2));
green = bin2dec(binary_target_string(3));
 
x = bin2dec(binary_target_string(21:32));
y = bin2dec(binary_target_string(5:16));
end
