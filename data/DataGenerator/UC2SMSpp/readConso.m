function [load, prim, tele] = readConso( ficName ),

        fid = fopen(ficName,'r');

tline = fgetl(fid);
if ( isempty(strfind(tline, 'CONSO' )) ),
error('inappropriate file');
end
        quoteMark='''';
% -- second line contains junk
tline = fgetl(fid);

% -- third line is junk too
        tline = fgetl(fid);

% -- likewise fourth
tline = fgetl(fid);

load = [];
for k=1:16,
tline = fgetl(fid);
load=[load str2num(tline)];
end

% -- read 2 junk lines
tline = fgetl(fid);
tline = fgetl(fid);

prim = [];
for k=1:16,
tline = fgetl(fid);
prim=[prim str2num(tline)];
end

% -- read 2 junk lines
tline = fgetl(fid);
tline = fgetl(fid);

tele = [];
for k=1:16,
tline = fgetl(fid);
tele=[tele str2num(tline)];
end

fclose(fid);

