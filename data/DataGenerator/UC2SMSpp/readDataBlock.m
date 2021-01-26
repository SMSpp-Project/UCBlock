function res = readDataBlock(fid, endFIELD, varargin)
               %
               % read data from a block in the opened file identified with the FID
% endFIELD is a string character indicating the end of the field to read
%

if ( nargin >= 3 ),
nbRows = varargin{1};
else
nbRows = inf;
end

        res = [];
tline = fgetl(fid);
iRow = 0;
while ( ischar(tline) && isempty(strfind(tline,endFIELD)) ),
res = [res str2num(tline)];
iRow = iRow + 1;
if ( iRow > nbRows );
break;
end
        tline = fgetl(fid);
end