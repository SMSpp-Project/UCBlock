function pp = readPDF(fid, pp, ind)
              % read a single set of pdfs from fid
% the indicator states if primary is in pourcentage or not

tline = fgetl(fid);
ss = strsplit(tline);

pp(1).start_time = str2num(ss{1});

nbC = str2num(ss{2});
nbD = str2num(ss{3});

for i=1:nbC
        tline = fgetl(fid);
xD = str2num(tline);
pp(1).flowc(i) = xD(1);
pp(1).puisc(i) = xD(2);
if ( ind < 0 ),
pp(1).primc(i) = xD(3);
pp(1).telc(i) = xD(4);
else
pp(1).telc(i) = xD(3);
end
        end

for i=1:nbD,
tline = fgetl(fid);
xD = str2num(tline);

pp(1).flowd(i) = xD(1);
pp(1).puisd(i) = xD(2);
if ( ind < 0 ),
pp(1).primd(i) = xD(3);
pp(1).teld(i) = xD(4);
else
pp(1).teld(i) = xD(3);
end
        end
