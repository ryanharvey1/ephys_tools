% parseSchmitt          detect segments of buffered threshold crossings
%
% y = parseSchmitt( x, hi, lo )
%
% to quality as a segment by standard threshold crossing, a segment has to be
% continuously above TH
% 
% in this routine, a segment has to 
% -be continuously about lo, 
% -begin with a neg->pos hi crossing following a neg->pos lo crossing
% -end with a pos->neg hi crossing followed by a pos->neg lo crossing
% 
% For instance, for lo/hi=5/10, the sequence 
% 
%       [ 5 7 12 13 7 15 3 ] 
% 
% yields the segment [ 3 6 ],
% whereas a standard TH-crossing yields [ 3 4 ] and [ 6 ]
%
% calls         SchmittTrigger

% 26-jan-13 ES

% revisions
% 20-aug-13 bug fix

% note that an alternative definition, to be continuously above hi
% (with all the rest the same) is equivalent to thresholding by hi alone;
% and the alterbative definition of being continusouly above hi with the
% beginning/end marked by the flanking lo crossing is equivalent to
% hi-thresholding and edge expansion
 
function y = parseSchmitt( x, hi, lo )

y = [];
x = x( : );
n = length( x );
if length( hi ) == 2
    lo = hi( 2 );
    hi = hi( 1 );
end
if ~n || isempty( hi ) || isempty( lo )
    return
end

% detect buffered hi-crossings
cr1 = SchmittTrigger( x, hi, lo );
cr2 = SchmittTrigger( flipud( x ), hi, lo );
if isempty( cr1 ) || isempty( cr2 ) 
    return
end
% remove single-sample segments
cr2 = flipud( n - cr2 + 1 );
[ ix i1 i2 ] = intersect( cr1, cr2 ); 
cr1( i1 ) = []; 
cr2( i2 ) = [];
if isempty( cr1 ) || isempty( cr2 ) 
    return
end
% combine the rest
mat = sortrows( [ [ cr1 ones( length( cr1 ), 1 ) ]; [ cr2 -ones( length( cr2 ), 1 ) ] ], 1 );
ridx = find( abs( diff( mat( :, 2 ) ) ) ~= 2 ); % 1: lo->hi; -1: hi->lo
mat( ridx, : ) = [];
if mat( 1, 2 ) == -1
    mat( 1, : ) = [];
end
if mat( end, 2 ) == 1
    mat( end, : ) = [];
end
if ~iseven( size( mat, 1 ) )
    error( 'bug!' )
end
m = size( mat, 1 );
y = reshape( mat( :, 1 )', [ 2 m / 2 ] )';

return

% EOF
