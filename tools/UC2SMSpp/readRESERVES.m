function rr = readRESERVES( ficName, rr )

fid = fopen(ficName,'r');

tline = fgetl(fid);
%if ( isempty(strfind(tline, 'RESERVES' )) ),
if ( isempty(strfind(tline, 'RESERVOIRS' )) ),
error('inappropriate file');
end
        quoteMark='''';

tline = fgetl(fid);
iR = 1;
while (ischar(tline)),
%Uncover the reservoir name
rr(iR).name = strrep(strtrim(tline),quoteMark,'');

%read field on type
tline = fgetl(fid);
%
ss = strsplit(tline);
rr(iR).v0 = str2num(ss{2});
rr(iR).vfinmin = str2num(ss{3});
rr(iR).vfinmax = str2num(ss{4});

if ( ~isempty(strfind(tline,'-1000.0')) ),
% new VU format
        tline = fgetl(fid);
if ( isempty(strfind(tline,'VU_EAU')) ),
error(strcat('Format error while reading water values for ', rr(iR).name));
end
%
tline = fgetl(fid);
nbP = str2num(tline);
% read partition
tline = fgetl(fid);
rr(iR).wpart = str2num(tline);
% read vu
tline = fgetl(fid);
rr(iR).wvalues = str2num(tline);
if ( ( length(rr(iR).wpart) ~= nbP ) || ( length(rr(iR).wvalues) ~= nbP ) )
error(strcat('Water value information incorrectly specified: dimension error ', rr(iR).name));
end
else
% old VU format
rr(iR).wvalues = str2num(ss{5});
end

%read apports
tline = fgetl(fid);
if ( isempty(strfind(tline,'APPORTS')) ),
error(strcat('Format error while reading water values for ', rr(iR).name));
end
rr(iR).inflows = readDataBlock(fid, 'MAX OPTIMISATION');

%read block max opti
readDataBlock(fid, 'MIN OPTIMISATION');
%read block max opti
readDataBlock(fid, 'MAX LISSAGE');

rr(iR).vmax = readDataBlock(fid, 'MIN LISSAGE');
% number of rows to read
        nbRows = length(rr(iR).vmax) / 6;
rr(iR).vmin = readDataBlock(fid, 'UNDEFINED', nbRows );

% go to the next reservoir
        tline = fgetl(fid);

iR = iR + 1;
if ( ~isempty(strfind(tline,'FIN')) )
break;
end
        end

fclose(fid);