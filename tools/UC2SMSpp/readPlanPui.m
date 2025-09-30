function grp = readPlanPui( ficName, grp ),

        fid = fopen(ficName,'r');

tline = fgetl(fid);
if ( isempty(strfind(tline, 'PLANPUI' )) ),
error('inappropriate file');
end
        quoteMark='''';
tline = fgetl(fid);

% 3rd line contains usefull data
tline = fgetl(fid);
while (ischar(tline)),
% Read name of this Groupe
%Uncover the Groupe name
thfName = strrep(strtrim(tline),quoteMark,'');

%find corresponding groupe in grp
        I = find( strcmp( thfName, {grp.name} ) );
if ( length(I) ~= 1 ),
error('oops');
end
        iG0 = I(1);

% this contains stateDuration
        tline = fgetl(fid);
ss = strsplit(strtrim(tline));

grp(iG0).sDur = str2num( ss{3} )/60.0 ; % en h

% read line with 'P'
tline = fgetl(fid);
if ( isempty(strfind(tline, 'P' )) ),
error('something went wrong');
end
% line with initial power
tline = fgetl(fid);
ss = strsplit(strtrim(tline));
grp(iG0).P0 = str2num( ss{1} );

% read 70 blanc lines
for i=1:70,
tline = fgetl(fid);
end

% read new line for new unit
tline = fgetl(fid);

if ( ~isempty(strfind(tline,'FIN')) )
break;
end
        end

fclose(fid);

