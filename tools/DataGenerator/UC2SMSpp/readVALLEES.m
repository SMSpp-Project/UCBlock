function vv = readVALLEES( ficName, vv, aa, rr, uu )
              %
              % Read hydro Valleys, but take as input reservoirs information and usine
% information
%

fid = fopen(ficName,'r');

tline = fgetl(fid);
if ( isempty(strfind(tline, 'VALLEES' )) ),
error('inappropriate file');
end
        quoteMark='''';

tline = fgetl(fid);
iV = 1;
while (ischar(tline)),
%Uncover the valley name
ss = strsplit(tline);
vv(iV).name = strrep(ss{1},quoteMark,'');

%read some size data
tline = fgetl(fid);
ss = strsplit(tline);
nbNode = str2num(ss{1});
nbRes = str2num(ss{2});

vv(iV).reservoirs = [];
vv(iV).flow = zeros(nbRes,nbRes+1); %the nbRes + 1st column is the final puit node
%vv(iV).updelay = zeros(nbRes,nbRes+1);
%vv(iV).dndelay = zeros(nbRes,nbRes+1);

vv(iV).arc = aa;
for i=1:nbNode,
tline = fgetl(fid);
ss = strsplit(tline,quoteMark);
ss2 = strsplit(ss{7});

usiName = ss{2}; %strrep(ss{1},quoteMark,'');
rAmont  = ss{4}; %strrep(ss{2},quoteMark,'');
rAval   = strtrim(ss{6}); %strtrim(strrep(ss{3},quoteMark,''));

iAm = find(strcmp({rr.name}, rAmont)==1);
iAv = find(strcmp({rr.name}, rAval)==1);

iU  = find(strcmp({uu.name}, usiName)==1);

if ( ~ismember(iAm, vv(iV).reservoirs) )
vv(iV).reservoirs = [vv(iV).reservoirs iAm];
end
if ( ~ismember(iAv, vv(iV).reservoirs) )
vv(iV).reservoirs = [vv(iV).reservoirs iAv];
end

        ilAm = find(vv(iV).reservoirs == iAm );
if ( isempty(rAval) )
ilAv = nbRes + 1;
else
ilAv = find(vv(iV).reservoirs == iAv );
end

vv(iV).flow( ilAm, ilAv ) = iU;

%vv(iV).arc(i) = aa; %struct('amont',{},'aval',{},'updelay',{},'dndelay',{},'usiIdx',{});

vv(iV).arc(i).amont = ilAm;
vv(iV).arc(i).aval  = ilAv;
vv(iV).arc(i).updelay = str2num(ss2{2});
vv(iV).arc(i).dndelay = str2num(ss2{3});
vv(iV).arc(i).usiIdx  = iU;

%vv(iV).updelay( ilAm, ilAv ) = str2num(ss2{2});
%vv(iV).dndelay( ilAm, ilAv ) = str2num(ss2{3});
end

% go to the next usine
        tline = fgetl(fid);

iV = iV + 1;
if ( ~isempty(strfind(tline,'FIN')) )
break;
end
        end

fclose(fid);

