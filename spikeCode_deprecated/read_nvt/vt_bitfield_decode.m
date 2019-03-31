%VT_BITFIELD_DECODE 		Decodes a Neuralynx VT bitfield
%	[red, green, blue, intensity, x, y] = vt_bitfield_decode(bitfield_value)
%	The values used for both target and point data in video tracker records are actually bitfields. 
%	Since Matlab does not easily handle bitfileds, this function will decode Neuralynx VT bitfields
%	for use in Matlab. 
%
%	Input:
%		bitfield_value:  the 32 bit bitfiled value.
%	Output:
%		red: 1 if bitfield has a red color component, 0 if not.
%		green:  1 if bitfield has a green color component, 0 if not.
%		blue:  1 if bitfield has a blue color component, 0 if not.
%		intensity:  1 if bitfield has a intensity color component, 0 if not.
%		x: extracted X coordinate for the bitfield
%		y: extracted X coordinate for the bitfield
%
% v1.0.0
function [red, green, blue, intensity, x, y] = vt_bitfield_decode(bitfield_value)  

	binary_bitfield_string = dec2bin(bitfield_value, 32);
	red = bin2dec(binary_bitfield_string(2));
	green = bin2dec(binary_bitfield_string(3));
	blue = bin2dec(binary_bitfield_string(4)); 
	intensity = bin2dec(binary_bitfield_string(17));
	 
	x = bin2dec(binary_bitfield_string(21:32));
	y = bin2dec(binary_bitfield_string(5:16));
end