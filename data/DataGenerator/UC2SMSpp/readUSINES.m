function uu = readUSINES( ficName, uu )

pdfs = struct('start_time',{},'flowc',{},'puisc',{},'primc',{},'telc',{},'flowd',{},'puisd',{},'primd',{},'teld',{});

fid = fopen(ficName,'r');

tline = fgetl(fid);
if ( isempty(strfind(tline, 'USINES' )) ),
error('inappropriate file');
end
        quoteMark='''';

tline = fgetl(fid);
iU = 1;
while (ischar(tline)),
%Uncover the usine name
uu(iU).name = strrep(strtrim(tline),quoteMark,'');

% default type is 0
uu(iU).type = 0;

%Read pimary pourcentage
        tline = fgetl(fid);
%
ss = strsplit(tline);
uu(iU).prim_pour = str2num(ss{5});
%
% -- pick up the type if it is a pump
if ( strcmp(strrep(ss{2},quoteMark,''),'P') ),
%disp(strcat('The usine:', uu(iU).name, ' is a pump'));
uu(iU).type = 1;
end

%
%Read gradients
tline = fgetl(fid);
ss = strsplit(strtrim(tline));
uu(iU).gradup = str2num(ss{1});
uu(iU).graddn = str2num(ss{2});

nbPdf = str2num(ss{5});
for i=1:nbPdf
uu(iU).pdfset(i) = readPDF(fid,pdfs, uu(iU).prim_pour);
end

% go to the next usine
        tline = fgetl(fid);

iU = iU + 1;
if ( ~isempty(strfind(tline,'FIN')) )
break;
end
        end

fclose(fid);