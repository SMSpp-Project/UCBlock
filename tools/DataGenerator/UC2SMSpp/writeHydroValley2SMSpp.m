function res = writeHydroValley2SMSpp( fic, T, dT, idxBlock, idxV, vallees, usines, reservoirs ),
%
% fic is the file handler : it is not closed
%
% idxV is the valley index
%

totPieces=0;
iU=[vallees(idxV).arc(:).usiIdx]'; %setdiff(unique(vallees(idxV).flow),0);

%for each usine, pick the last pdfset with starttime 30
%
iPdf = ones(length(iU),1);
for i=1:length(iU),
        IdxPdf = find([usines(iU(i)).pdfset(:).start_time]==30);
if ( isempty(IdxPdf) ),
iPdf(i) = 1;
else
iPdf(i) = IdxPdf(end);
end
        end

% the number of pieces is directly given by the set of points
% unless the unit is a pump, in which case we have to put 1
for i=1:length(iU),
if (usines(iU(i)).type == 0),
totPieces = totPieces + length(usines(iU(i)).pdfset(iPdf(i)).flowc);
else
totPieces = totPieces + 1;
end
        end

nbR=length(vallees(idxV).reservoirs);
nbA=length(vallees(idxV).arc);

fprintf(fic,'\t group: UnitBlock_%d {\n',idxBlock);
fprintf(fic,'\t\t // The %s valley\n', vallees(idxV).name);
fprintf(fic,'\t\t dimensions:\n');
fprintf(fic,'\t\t\t NumberReservoirs=%d;\n',length(vallees(idxV).reservoirs));
%fprintf(fic,'\t\t\t NumberArcs=%d;\n',length(find(unique(vallees(idxV).flow)>0)));
fprintf(fic,'\t\t\t NumberArcs=%d;\n',length(vallees(idxV).arc));
fprintf(fic,'\t\t\t NumberIntervals = %d;\n', T );
fprintf(fic,'\t\t\t TotalNumberPieces=%d;\n',totPieces);
fprintf(fic,'\n');
fprintf(fic,'\t\t variables:\n');
fprintf(fic,'\t\t\t int StartArc(NumberArcs);\n');
fprintf(fic,'\t\t\t int EndArc(NumberArcs);\n');
fprintf(fic,'\n');
fprintf(fic,'\t\t\t double Inflows(NumberReservoirs,NumberIntervals);\n');
fprintf(fic,'\t\t\t double InitialVolumetric(NumberReservoirs);\n');
fprintf(fic,'\t\t\t double MinVolumetric(NumberReservoirs,NumberIntervals);\n');
fprintf(fic,'\t\t\t double MaxVolumetric(NumberReservoirs,NumberIntervals);\n');
fprintf(fic,'\n');
fprintf(fic,'\t\t\t int UphillFlow(NumberArcs);\n');
fprintf(fic,'\t\t\t int DownhillFlow(NumberArcs);\n');
fprintf(fic,'\t\t\t int InitialFlowRate(NumberArcs);\n');
%
fprintf(fic,'\t\t\t double DeltaRampUp(NumberIntervals,NumberArcs);\n');
fprintf(fic,'\t\t\t double DeltaRampDown(NumberIntervals,NumberArcs);\n');
%
fprintf(fic,'\t\t\t double MinFlow(NumberIntervals,NumberArcs);\n');
fprintf(fic,'\t\t\t double MaxFlow(NumberIntervals,NumberArcs);\n');
fprintf(fic,'\t\t\t double MinPower(NumberIntervals,NumberArcs);\n');
fprintf(fic,'\t\t\t double MaxPower(NumberIntervals,NumberArcs);\n');
fprintf(fic,'\t\t\t int NumberPieces(NumberArcs);\n');
%
fprintf(fic,'\n');
fprintf(fic,'\t\t\t double LinearTerm(TotalNumberPieces);\n');
fprintf(fic,'\t\t\t double ConstantTerm(TotalNumberPieces);\n');

% --other stuff

fprintf(fic, '\t\t // group attributes:\n');
fprintf(fic, '\t\t\t\t :type = "HydroUnitBlock" ;\n');

fprintf(fic,'\t\t data:\n');
fprintf(fic,'\t\t\t StartArc = %s ;\n',regexprep(int2str([vallees(idxV).arc(:).amont]-1),'\s+',',') );
fprintf(fic,'\t\t\t EndArc = %s ;\n',regexprep(int2str([vallees(idxV).arc(:).aval]-1),'\s+',',') );
fprintf(fic,'\n');
% Reservoir Stuff

fprintf(fic,'\t\t\t Inflows=\n');
for i=1:nbR-1,
iR = vallees(idxV).reservoirs(i);
% inflows are in m3 / s -- we need to convert in m3 / time step
iF = regexprep( num2str(reservoirs(iR).inflows*dT*3600 ),'\s+',',');
fprintf(fic,'\t\t\t\t%s,\n',iF);
end
        iR = vallees(idxV).reservoirs(nbR);
iF = regexprep(num2str(reservoirs(iR).inflows*dT*3600 ),'\s+',',');
fprintf(fic,'\t\t\t\t%s;\n',iF);

fprintf(fic,'\t\t\t InitialVolumetric=');
for i=1:nbR,
if ( i==nbR )
fprintf(fic,'%f;\n',reservoirs(vallees(idxV).reservoirs(i)).v0);
else
fprintf(fic,'%f,',reservoirs(vallees(idxV).reservoirs(i)).v0);
end
        end

fprintf(fic,'\t\t\t MinVolumetric=\n');
for i=1:nbR-1,
iR = vallees(idxV).reservoirs(i);

iF = regexprep(num2str(reservoirs(iR).vmin),'\s+',',');
fprintf(fic,'\t\t\t\t%s,\n',iF);
end
        iR = vallees(idxV).reservoirs(nbR);
iF = regexprep(num2str(reservoirs(iR).vmin),'\s+',',');
fprintf(fic,'\t\t\t\t%s;\n',iF);

fprintf(fic,'\t\t\t MaxVolumetric=\n');
for i=1:nbR-1,
iR = vallees(idxV).reservoirs(i);

iF = regexprep(num2str(reservoirs(iR).vmax),'\s+',',');
fprintf(fic,'\t\t\t\t%s,\n',iF);
end
        iR = vallees(idxV).reservoirs(nbR);
iF = regexprep(num2str(reservoirs(iR).vmax),'\s+',',');
fprintf(fic,'\t\t\t\t%s;\n',iF);
%
fprintf(fic,'\n');

% Arc Stuff
fprintf(fic,'\t\t\t UphillFlow = %s ;\n',regexprep(int2str(ceil([vallees(idxV).arc(:).updelay]/(dT*60.0))),'\s+',',') );
fprintf(fic,'\t\t\t DownhillFlow = %s ;\n',regexprep(int2str(ceil([vallees(idxV).arc(:).dndelay]/(dT*60.0))),'\s+',',') );
fprintf(fic,'\t\t\t InitialFlowRate = %s ;\n',regexprep(num2str( [usines([vallees(idxV).arc(:).usiIdx]).d0]*3600*dT ),'\s+',',') );

%Gradients
        iStr = regexprep(num2str( ([usines([vallees(idxV).arc(:).usiIdx]).gradup]*(60*dT))*3600*dT),'\s+',',');
fprintf(fic,'\t\t\t DeltaRampUp = \n' );
for it=1:T-1,
fprintf(fic,'\t\t\t\t%s,\n',iStr);
end
fprintf(fic,'\t\t\t\t%s ;\n',iStr);

iStr = regexprep(num2str( ([usines([vallees(idxV).arc(:).usiIdx]).graddn]*(60*dT))*3600*dT),'\s+',',');
fprintf(fic,'\t\t\t DeltaRampDown = \n'  );
for it=1:T-1,
fprintf(fic,'\t\t\t\t%s,\n',iStr);
end
fprintf(fic,'\t\t\t\t%s ;\n',iStr);

%Minflow, maxflow stuff
fprintf(fic,'\t\t\t MinFlow=');
for it=1:T
for iA=1:nbA,
iUsi=vallees(idxV).arc(iA).usiIdx;

ii=find(iU==iUsi);

% the real minimum flow rate is given in the discrete set of points
% ...
if ( usines(iUsi).type == 0 ),
minF = min(usines(iUsi).pdfset(iPdf(ii)).flowd)*dT*3600;

% -- check if fine with the real initial flow rate
        ominF = minF;
minF = min(min(round(1e6*usines(iUsi).initflow)/1e6)*dT*3600,minF);
if ( (minF < ominF) && (it==1) ),
warning(strcat('Changes made to min flow rate for usine:', usines(iUsi).name));
usines(iUsi).minfchanged=1;
end
else
minF = -1.0*max(usines(iUsi).pdfset(iPdf(ii)).flowd)*dT*3600;
end
if ( (iA==nbA) && (it==T)),
fprintf(fic,'%f',minF);
else
fprintf(fic,'%f,',minF);
end
        end
end
fprintf(fic,';\n');
%
fprintf(fic,'\t\t\t MaxFlow=');
for it=1:T
for iA=1:nbA,
iUsi=vallees(idxV).arc(iA).usiIdx;

ii=find(iU==iUsi);

if ( usines(iUsi).type == 0 ),
maxF = sum(usines(iUsi).pdfset(iPdf(ii)).flowc)*dT*3600;

% -- check if fine with the real initial flow rate
        omaxF = maxF;
maxF = max(max(usines(iUsi).initflow)*dT*3600,maxF);
if ( (maxF > omaxF) && (it==1) ),
warning(strcat('Changes made to max flow rate for usine:', usines(iUsi).name));
end
else
maxF = -1.0*min(usines(iUsi).pdfset(iPdf(ii)).flowd)*dT*3600;
end
if ( (iA==nbA) && (it==T)),
fprintf(fic,'%f',maxF);
else
fprintf(fic,'%f,',maxF);
end

        end
end
fprintf(fic,';\n');
%
fprintf(fic,'\t\t\t MinPower=');
for it=1:T
for iA=1:nbA,
iUsi=vallees(idxV).arc(iA).usiIdx;

ii=find(iU==iUsi);

%min power is also related to discrete case
minP = min(abs(usines(iUsi).pdfset(iPdf(ii)).puisd));
%if ( ~isempty(usines(iUsi).minfchanged) ),
%   minP=0.0;
%end
if ( (iA==nbA) && (it==T)),
fprintf(fic,'%f',minP);
else
fprintf(fic,'%f,',minP);
end
        end
end
fprintf(fic,';\n');
%
fprintf(fic,'\t\t\t MaxPower=');
for it=1:T
for iA=1:nbA,
iUsi=vallees(idxV).arc(iA).usiIdx;

ii=find(iU==iUsi);

maxP = abs(sum(usines(iUsi).pdfset(iPdf(ii)).puisc));
if ( (iA==nbA) && (it==T)),
fprintf(fic,'%f',maxP);
else
fprintf(fic,'%f,',maxP);
end
        end
end
fprintf(fic,';\n');
%
fprintf(fic,'\t\t\t NumberPieces=');
for iA=1:nbA,
iUsi=vallees(idxV).arc(iA).usiIdx;

ii=find(iU==iUsi);

if ( usines(iUsi).type == 0 ),
nbP = length(usines(iUsi).pdfset(iPdf(ii)).puisc);
else
nbP = 1; % a pump only has one piece
end
if ( iA==nbA ),
fprintf(fic,'%d',nbP);
else
fprintf(fic,'%d,',nbP);
end
        end
fprintf(fic,';\n');
%
fprintf(fic,'\n');
% Generation stuff
fprintf(fic,'\t\t\t LinearTerm=');
for iA=1:nbA,
iUsi=vallees(idxV).arc(iA).usiIdx;

ii=find(iU==iUsi);

if ( usines(iUsi).type == 0),
slope=usines(iUsi).pdfset(iPdf(ii)).puisc ./ max(usines(iUsi).pdfset(iPdf(ii)).flowc*dT*3600,1);
else
slope=-sum(usines(iUsi).pdfset(iPdf(ii)).puisc)./ max(sum(usines(iUsi).pdfset(iPdf(ii)).flowc*dT*3600),1);
end

        linT = regexprep(num2str(slope),'\s+',',');
if ( iA==nbA ),
fprintf(fic,'%s',linT);
else
fprintf(fic,'%s,',linT);
end
        end
fprintf(fic,';\n');
%
fprintf(fic,'\t\t\t ConstantTerm=');
for iA=1:nbA,
iUsi=vallees(idxV).arc(iA).usiIdx;

ii=find(iU==iUsi);

slope=usines(iUsi).pdfset(iPdf(ii)).puisc ./ max(usines(iUsi).pdfset(iPdf(ii)).flowc,1);
intercept=cumsum(usines(iUsi).pdfset(iPdf(ii)).puisc);
fls=cumsum(usines(iUsi).pdfset(iPdf(ii)).flowc);

cT = intercept - slope.*fls;
if ( usines(iUsi).type == 1 ),
cT = 0.0;
end
        conT = regexprep(num2str(cT),'\s+',',');

if ( iA==nbA ),
fprintf(fic,'%s',conT);
else
fprintf(fic,'%s,',conT);
end
        end
fprintf(fic,';\n');

%

fprintf(fic,'\n');
fprintf(fic,'\t\t } // group UnitBlock_%d\n',idxBlock);



