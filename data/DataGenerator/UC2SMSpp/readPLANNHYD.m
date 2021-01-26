function [totGen, totPri, totTel, uu] = readPLANNHYD( ficName, uu ),

fid = fopen(ficName,'r');

tline = fgetl(fid);
if ( isempty(strfind(tline, 'PLANN_HYD' )) ),
error('inappropriate file');
end
        quoteMark='''';
% -- second line contains junk
tline = fgetl(fid);

totGen = [];
totPri = [];
totTel = [];

% -- now for each valley we need to do things
% 3rd line contains usefull data
tline = fgetl(fid);
while (ischar(tline)),
ss = strsplit(strtrim(tline));
nbRes = str2num(ss{3});
nbUsi = str2num(ss{2});

% Read Reservoirs
for i=1:nbRes,
% Nom du Reservoir
        tline = fgetl(fid);
%
% Nom destocke ou volume
tline = fgetl(fid);
if ( ~isempty(strfind(tline, 'DESTOCKE' )) ),
tline = fgetl(fid);
end
if ( ~isempty(strfind(tline, 'VOLUME' )) ),
for k=1:17,
tline = fgetl(fid);
end
        end
end

% Read Usines
for i=1:nbUsi,
% Nom du Usine
        tline = fgetl(fid);

ss = strsplit(tline, quoteMark);
usiName = strrep(ss{2},quoteMark,'');

iU  = find(strcmp({uu.name}, usiName)==1);
%read 18 lines for debit
        tline = fgetl(fid);
if ( isempty(strfind(tline, 'DEBIT' )) ),
error('inappropriate file');
end
        tline = fgetl(fid); %initial flow rate
if ( ~isempty(iU) ),
uu(iU).d0 = str2num(tline);
end
        flowR=[];
for k=1:16,
tline = fgetl(fid);
flowR=[flowR str2num(tline)];
end
if ( ~isempty(iU) ),
uu(iU).initflow = flowR;
end

%read line puissance
        tline = fgetl(fid);
if ( isempty(strfind(tline, 'PUISSANCE' )) ),
error('inappropriate file');
end
%read line for initial power output
        tline = fgetl(fid);
%
pwr=[];
for k=1:16,
tline = fgetl(fid);
pwr=[pwr str2num(tline)];
end
%

%read line primaire
        tline = fgetl(fid);
if ( isempty(strfind(tline, 'PRIMAIRE' )) ),
error('inappropriate file');
end
%read line for initial power output
        tline = fgetl(fid);
%
prim=[];
for k=1:16,
tline = fgetl(fid);
prim=[prim str2num(tline)];
end

%read line primaire
        tline = fgetl(fid);
if ( isempty(strfind(tline, 'TELEREGLAGE' )) ),
error('inappropriate file');
end
%read line for initial power output
        tline = fgetl(fid);
%
tele=[];
for k=1:16,
tline = fgetl(fid);
tele=[tele str2num(tline)];
end

%
if ( ~isempty(iU) ),
if (isempty(totGen) ),
totGen=pwr;

totPri = prim;
totTel = tele;
else
if (length(totGen) ~= length(pwr)),
error('aie');
end
        totGen = totGen + pwr;

if (length(totPri) ~= length(prim)),
error('aie');
end
        totPri = totPri + prim;

if (length(totTel) ~= length(tele)),
error('aie');
end
        totTel = totTel + tele;
end
        end
end

% read new line for new valley
tline = fgetl(fid);

if ( ~isempty(strfind(tline,'FIN')) )
break;
end
        end





fclose(fid);
