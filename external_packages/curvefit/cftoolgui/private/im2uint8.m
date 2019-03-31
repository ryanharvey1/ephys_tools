function u = im2uint8(img, typestr)
%IM2UINT8 Convert image to eight-bit unsigned integers.
%   IM2UINT8 takes an image as input, and returns an image of
%   class uint8.  If the input image is of class uint8, the
%   output image is identical to it.  If the input image is of
%   class double or uint16 , im2uint8 returns the equivalent
%   image of class uint8, rescaling or offsetting the data as
%   necessary.
%
%   I2 = IM2UINT8(I1) converts the intensity image I1 to uint8,
%   rescaling the data if necessary.
%
%   RGB2 = IM2UINT8(RGB1) converts the truecolor image RGB1 to
%   uint8, rescaling the data if necessary.
%
%   BW2 = IM2UINT8(BW1) converts the binary image BW1 to uint8.
%
%   X2 = IM2UINT8(X1,'indexed') converts the indexed image X1 to
%   uint8, offsetting the data if necessary. Note that it is
%   not always possible to convert an indexed image to
%   uint8. If X1 is double, then the maximum value of X1 must
%   be 256 or less.  If X1 is uint16, the maximum value of X1
%   must be 255 or less.
%
%   See also DOUBLE, IM2DOUBLE, IM2UINT16, UINT8, UINT16, IMAPPROX.

%   Copyright 2000-2004 The MathWorks, Inc.
%   $Revision: 1.4.2.4 $  $Date: 2010/10/08 16:36:45 $

if isa(img, 'uint8')
    u = img; 
   
elseif isa(img, 'double') || isa(img, 'uint16')
    if nargin==1
        if islogical(img)      
            % binary image
            u = uint8(img~=0);
        else
            % intensity image; call MEX-file
            u = grayto8(img);
        end
        
    elseif nargin==2
    
        if ~ischar(typestr) || (typestr(1) ~= 'i')
            error(message('curvefit:im2uint8:invalidInArgs'));
        else 
            if (isa(img, 'uint16'))
                if (max(img(:)) > 255)
                    error(message('curvefit:im2uint8:tooManyColors'));
                else
                    u = uint8(img);
                end
                
            else
                % img is double
                if max(img(:))>=257 
                    error(message('curvefit:im2uint8:tooManyColors'));
                elseif min(img(:))<1
                    error(message('curvefit:im2uint8:invalidIndexedImage'));
                else
                    u = uint8(img-1);
                end
            end
        end
    else
        error(message('curvefit:im2uint8:invalidInArgs'));
    end
else
    error(message('curvefit:im2uint8:unsupportedClass'));
end
