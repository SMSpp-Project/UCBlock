function grp = readGroupes( ficName )

grp = struct('name',{}, 'Paux', {}, 'cfix', {}, 'cprop', {}, 'cquad', {}, 'stcost', {}, 'tup', {}, 'tdn', {}, 'minup', {}, 'mindn', {}, 'MSG', {}, 'FL', {}, 'P0', {}, 'sDur', {}, 'Ppri', {}, 'Ptel', {});

fid = fopen(ficName,'r');

tline = fgetl(fid);
if ( isempty(strfind(tline, 'GROUPES' )) ),
error('inappropriate file');
end
        quoteMark='''';

tline = fgetl(fid);
iG = 1;
while (ischar(tline)),
%Uncover the Groupe name
grp(iG).name = strrep(strtrim(tline),quoteMark,'');

%Read puissance à l'arret : Pau
tline = fgetl(fid);
%
ss = strsplit(strtrim(tline));
grp(iG).Paux = str2num(ss{2});

% les dernier trois champs sont la marge primaire
grp(iG).Ppri = (str2num(ss{10}) + str2num(ss{11}) + str2num(ss{12}))/3;
grp(iG).Ptel = (str2num(ss{7}) + str2num(ss{8}) + str2num(ss{9}))/3;

%Read coût : cout fixe, cout prop, coeff quadratique
tline = fgetl(fid);
%
ss =  strsplit(strtrim(tline));
grp(iG).cfix = str2num(ss{1});
grp(iG).cprop = str2num(ss{2});
grp(iG).cquad = str2num(ss{3});

% 4th line : Start up Costs
        tline = fgetl(fid);
ss = strsplit(strtrim(tline));

type=strrep(strtrim(ss{1}),quoteMark,'');;
C1  =str2num(ss{2});
C2 = str2num(ss{3});
C3 = str2num(ss{4});
mxC= str2num(ss{5});

stcost = 0.0;
% The following is based on a standard off duration of 4 hours
if ( ~isempty(strfind(type,'L')) )
stcost = min(C1 + C2*log(4*60.0 / C3), mxC);
end
if ( ~isempty(strfind(type,'E')) )
stcost = max(min(C1 + C2*(1.0 - exp(4*60.0 / C3)), mxC),0.0);
end
if ( ( isempty(strfind(type,'L')) ) && ( isempty(strfind(type,'E')) ) )
error('Unknown type');
end
grp(iG).stcost = stcost;

% 5th line : Ramps and min up / min down
tline = fgetl(fid);
ss = strsplit(strtrim(tline));

grp(iG).tup = str2num(ss{1})*60.0; % en MW/h
grp(iG).tdn = str2num(ss{2})*60.0; % en MW/h
grp(iG).mindn=str2num(ss{3})/60.0; % en heures
grp(iG).minup=str2num(ss{4})/60.0; % en heures

% 6th line : Nombre de Paliers
%
tline = fgetl(fid);

% 7th line : Les paliers
tline = fgetl(fid);
ss = strsplit(strtrim(tline));

grp(iG).MSG = str2num(ss{1});
grp(iG).FL  = str2num(ss{end});

% go to the next Groupe
        tline = fgetl(fid);

iG = iG + 1;
if ( ~isempty(strfind(tline,'FIN')) )
break;
end
        end

fclose(fid);
