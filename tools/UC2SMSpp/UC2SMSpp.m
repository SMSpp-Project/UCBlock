%%
%
% -- UC to SMS++
%

% Variants to look into
%
% Spinning reserves:
%   None
%   Prim
%   PplusS
        SpinRes = 'none';
%SpinRes = 'prim';
%SpinRes = 'PplusS';

% Thermal:
%   noHydro
%   wHydro
%   pHydro -- some subset of cascading systems -- to use with Caution
%
%ThermOpt = 'noHydro';
%ThermOpt = 'wHydro';
ThermOpt = 'pHydro';

VnamesToExclude = struct('name',{});
VnamesToExclude(end+1).name = 'DORDOGNE';
VnamesToExclude(end+1).name = 'TRUY-SUP';
VnamesToExclude(end+1).name = 'TRUY-INF';
VnamesToExclude(end+1).name = 'AIN';
VnamesToExclude(end+1).name = 'ARC';
VnamesToExclude(end+1).name = 'ROMANCHE';
VnamesToExclude(end+1).name = 'DRAC1';
VnamesToExclude(end+1).name = 'DRAC2';
VnamesToExclude(end+1).name = 'ROSELEND';
VnamesToExclude(end+1).name = 'H-ISERE';
VnamesToExclude(end+1).name = 'PRAGNERE';
VnamesToExclude(end+1).name = 'ARIEGE';
VnamesToExclude(end+1).name = 'ORLU';
VnamesToExclude(end+1).name = 'MONTAHUT';
VnamesToExclude(end+1).name = 'DURANCE';
VnamesToExclude(end+1).name = 'CHASSEZA';
VnamesToExclude(end+1).name = 'LARDIT';
VnamesToExclude(end+1).name = 'BEYSSAC';
%VnamesToExclude(end+1).name = 'ARRENS';
VnamesToExclude(end+1).name = 'VICDESSO';
VnamesToExclude(end+1).name = 'ARAING';
VnamesToExclude(end+1).name = 'CERE';
VnamesToExclude(end+1).name = 'MAULDE';
VnamesToExclude(end+1).name = 'PIQUES';
VnamesToExclude(end+1).name = 'AUDE';
VnamesToExclude(end+1).name = 'MARCIL';
%VnamesToExclude(2).name = 'AIN';

dataSetsUC = struct('basePth', {}, 'name',{});

dataSetsUC(1).basePth = '/Users/ali/Desktop/EDF/20100623/20090907';
dataSetsUC(1).name = '20090907';
dataSetsUC(2).basePth = '/Users/ali/Desktop/EDF/20100623/20090908';
dataSetsUC(2).name = '20090908';
dataSetsUC(3).basePth = '/Users/ali/Desktop/EDF/20100623/20091005';
dataSetsUC(3).name = '20091005';
dataSetsUC(4).basePth = '/Users/ali/Desktop/EDF/20100623/20100311';
dataSetsUC(4).name = '20100311';
dataSetsUC(5).basePth = '/Users/ali/Desktop/EDF/20100623/20100323';
dataSetsUC(5).name = '20100323';
dataSetsUC(6).basePth = '/Users/ali/Desktop/EDF/20100623/20100623';
dataSetsUC(6).name = '20100623';
dataSetsUC(7).basePth = '/Users/ali/Desktop/EDF/20100623/20101231';
dataSetsUC(7).name = '20101231';


for iSets=1:7 %length(dataSetsUC),
        basePath = dataSetsUC(iSets).basePth; %'/Users/ali/Desktop/EDF/20100623/20101231';

GrpFile = '/GROUPES';
PlanFile= '/PLANPUI';

HydFile = '/PLANN_HYD';

ConsoFile='/CONSO';

T = 96;
dT=0.5;

% -- Recover Thermal Data
        ficGROUPES = [basePath GrpFile];
ficPLANPUI = [basePath PlanFile];

grp = readGroupes( ficGROUPES );
grp = readPlanPui( ficPLANPUI, grp );

% -- Recover Hydro Data
[vallees, usines, reservoirs] = readHydroValleys( basePath );

% -- set the name
if ( strcmp( ThermOpt,'pHydro') ),
UCName = strcat(dataSetsUC(iSets).name,'_',ThermOpt,'_', int2str( length(vallees) - length(VnamesToExclude) ), 'A_' ,SpinRes);
else
UCName = strcat(dataSetsUC(iSets).name,'_',ThermOpt,'_',SpinRes);
end

% -- Figure out if some indexes need to be excluded:
VexcIdx = [];
subUsiIdx=[];
if ( strcmp( ThermOpt,'pHydro') ),
for i=1:length(VnamesToExclude),
        i0 = find(strcmp({vallees(:).name},VnamesToExclude(i).name));
VexcIdx = [VexcIdx; i0];
end

for i=1:length(VexcIdx),
        subUsiIdx = [subUsiIdx; [vallees(VexcIdx(i)).arc(:).usiIdx]'];
end
        end

% -- Get Hydro Base Program
ficPLANHYD = [basePath HydFile];
[totHydroGen, totHydroGenPr, totHydroGenTel, usines]= readPLANNHYD( ficPLANHYD, usines );

% total load to recover from the to be removed usines
uusub = usines( subUsiIdx );
[totHydroGen, totHydroGenPr, totHydroGenTel ]= readPLANNHYD( ficPLANHYD, uusub );

% -- Get Load
ficConso = [basePath ConsoFile];
[Load, PriLoad, SecLoad] = readConso( ficConso );

% -- Figure out if Hydro data is consistent in relation to the initial
% schedule
for idxV=1:length(vallees),
        reservoirs = determineHydroConsistency( idxV, dT, T, vallees, usines, reservoirs );
end

%% -- Convert to SMSpp structure
%
%
%
%
if ( strcmp( ThermOpt,'noHydro') ),
disp('No Hydro option: recomputing load');
lload   = Load - totHydroGen; %net load
priload = max(PriLoad - totHydroGenPr,0);
telload = max(SecLoad - totHydroGenTel,0);
elseif ( strcmp( ThermOpt,'pHydro') ),
disp('Partial Hydro option: recomputing load');
lload   = Load - totHydroGen; %net load
priload = max(PriLoad - totHydroGenPr,0);
telload = max(SecLoad - totHydroGenTel,0);
else
lload   = Load;
priload = PriLoad;
telload = SecLoad;
end

%fic = fopen(strcat('C:\LocalDriveD\LocLinux\SMS++\UC\', UCName, '.txt'),'w');
%C:\LocalDriveD\Tools\plan4res-new\p4r-env\scripts\add-ons\install\sms++\examples
%fic = fopen(strcat('C:\LocalDriveD\Tools\plan4res-new\p4r-env\scripts\add-ons\install\sms++\examples\', UCName, '.txt'),'w');
fic = fopen(strcat('/Users/ali/Desktop/EDF/', UCName, '.txt'),'w');
fprintf(fic, 'netcdf \\%s {\n',UCName);
fprintf(fic, '\n');
fprintf(fic, '// group attributes:\n');
fprintf(fic, '\t\t :SMS++_file_type = 1 ;\n');
fprintf(fic, '\n');
fprintf(fic, 'group: Block_0 {\n');
fprintf(fic, '\t dimensions:\n');
fprintf(fic, '\t\t TimeHorizon = %d ;\n', T );
if ( strcmp( ThermOpt,'noHydro') ),
fprintf(fic, '\t\t NumberUnits = %d ;\n', length(grp));
else
fprintf(fic, '\t\t NumberUnits = %d ;\n', length(grp)+length(vallees)-length(VexcIdx));
end
%fprintf(fic, '\t\t NumberIntervals = 1 ;\n');
%
if ( strcmp(SpinRes,'prim') || strcmp(SpinRes,'PplusS') )
fprintf(fic, '\t\t NumberPrimaryZones = 1 ;\n');
end
if ( strcmp(SpinRes,'PplusS') )
fprintf(fic, '\t\t NumberSecondaryZones = 1 ;\n');
end
fprintf(fic, '\t variables:\n');
fprintf(fic, '\t\t double ActivePowerDemand(TimeHorizon) ;\n');
if ( strcmp(SpinRes,'prim') || strcmp(SpinRes,'PplusS') )
fprintf(fic, '\t\t double PrimaryDemand(TimeHorizon) ;\n');
end
if ( strcmp(SpinRes,'PplusS') )
fprintf(fic, '\t\t double SecondaryDemand(TimeHorizon) ;\n');
end
fprintf(fic, '\n');
fprintf(fic, '\t // group attributes:\n');
fprintf(fic, '\t\t :type = "UCBlock" ;\n');
fprintf(fic, '\t data:\n');
fprintf(fic, '\n');
fprintf(fic,'\t\t ActivePowerDemand=');
iF = regexprep(num2str(lload),'\s+',',');
fprintf(fic,'\t\t\t%s;\n',iF);
if ( strcmp(SpinRes,'prim') || strcmp(SpinRes,'PplusS') )
fprintf(fic,'\t\t PrimaryDemand=');
iF = regexprep(num2str(priload),'\s+',',');
fprintf(fic,'\t\t\t%s;\n',iF);
end
if ( strcmp(SpinRes,'PplusS') )
fprintf(fic,'\t\t SecondaryDemand=');
iF = regexprep(num2str(telload),'\s+',',');
fprintf(fic,'\t\t\t%s;\n',iF);
end
fprintf(fic, '\n');
fprintf(fic,'\t group: NetworkData {\n');
fprintf(fic,'\t\t dimensions:\n');
fprintf(fic,'\t\t\t	NumberNodes = 1 ;\n');
fprintf(fic,'\t } // group NetworkData\n');
fprintf(fic, '\n');

% for each Thermal block add a group
for iG=1:length(grp),
        fprintf(fic,'\t group: UnitBlock_%d {\n',iG-1);
fprintf(fic,'\t\t variables:\n');
fprintf(fic,'\t\t\t	double MinPower ;\n');
fprintf(fic,'\t\t\t	double MaxPower ;\n');
fprintf(fic,'\t\t\t	double FixedConsumption;\n');
fprintf(fic,'\t\t\t	double DeltaRampUp ;\n');
fprintf(fic,'\t\t\t	double DeltaRampDown ;\n');
if ( strcmp(SpinRes,'prim') || strcmp(SpinRes,'PplusS') )
fprintf(fic,'\t\t\t	double PrimaryRho ;\n');
end
if ( strcmp(SpinRes,'PplusS') )
fprintf(fic,'\t\t\t	double SecondaryRho ;\n');
end
fprintf(fic,'\t\t\t	double QuadTerm ;\n');
fprintf(fic,'\t\t\t	double StartUpCost ;\n');
fprintf(fic,'\t\t\t	double LinearTerm ;\n');
fprintf(fic,'\t\t\t	double ConstTerm ;\n');
fprintf(fic,'\t\t\t	double InitialPower ;\n');
fprintf(fic,'\t\t\t	int64 InitUpDownTime ;\n');
fprintf(fic,'\t\t\t	uint64 MinUpTime ;\n');
fprintf(fic,'\t\t\t	uint64 MinDownTime ;\n');
fprintf(fic, '\n');
fprintf(fic, '\t\t // group attributes:\n');
fprintf(fic,'\t\t\t	:type = "ThermalUnitBlock" ;\n');
fprintf(fic,'\t\t data:\n');
fprintf(fic,'\t\t\t MinPower = %10.6f ;\n', grp(iG).MSG );
fprintf(fic,'\t\t\t MaxPower = %10.6f ;\n', grp(iG).FL );
fprintf(fic,'\t\t\t FixedConsumption = %10.6f ;\n', grp(iG).Paux );

fprintf(fic,'\t\t\t DeltaRampUp = %10.6f ;\n', grp(iG).tup*dT );
fprintf(fic,'\t\t\t DeltaRampDown = %10.6f ;\n', grp(iG).tdn*dT );

if ( strcmp(SpinRes,'prim') || strcmp(SpinRes,'PplusS') )
fprintf(fic,'\t\t\t PrimaryRho = %10.6f ;\n', (1.0*grp(iG).Ppri) / (max(grp(iG).FL,1))  );
end
if ( strcmp(SpinRes,'PplusS') )
fprintf(fic,'\t\t\t SecondaryRho = %10.6f ;\n', (1.0*grp(iG).Ptel) / (max(grp(iG).FL,1)) );
end

fprintf(fic,'\t\t\t ConstTerm = %10.6f ;\n', grp(iG).cfix*dT );
fprintf(fic,'\t\t\t LinearTerm = %10.6f ;\n', grp(iG).cprop*dT );
fprintf(fic,'\t\t\t QuadTerm = %10.6f ;\n', grp(iG).cquad*dT );

fprintf(fic,'\t\t\t StartUpCost = %10.6f ;\n', grp(iG).stcost );
if ( (grp(iG).P0 < grp(iG).MSG) || (grp(iG).P0 > grp(iG).FL) )
warning(strcat('Initial power out of bounds :', grp(iG).name));
end
fprintf(fic,'\t\t\t InitialPower = %10.6f ;\n', max(min(grp(iG).P0,grp(iG).FL),grp(iG).MSG) );
sign = 1 - 2*(grp(iG).P0 <= 0);
fprintf(fic,'\t\t\t InitUpDownTime = %d ;\n', sign*floor(grp(iG).sDur/dT) );

fprintf(fic,'\t\t\t MinUpTime = %10.6f ;\n', floor(grp(iG).minup/dT) );
fprintf(fic,'\t\t\t MinDownTime = %10.6f ;\n', floor(grp(iG).mindn/dT) );

fprintf(fic,'\t\t } // group UnitBlock_%d\n',iG-1);
fprintf(fic, '\n');
end
if ( strcmp( ThermOpt,'pHydro') ),
idxV0 = length(grp);
% -- for each valley add the relevant information
for idxV=1:length(vallees),
if (~ismember( idxV, VexcIdx ) ),

writeHydroValley2SMSpp(fic, T, dT, idxV0, idxV, vallees, usines, reservoirs);

idxV0 = idxV0 + 1;
fprintf(fic, '\n');
end
        end
end
if ( strcmp( ThermOpt,'wHydro') ),
% -- for each valley add the relevant information
for idxV=1:length(vallees),
        writeHydroValley2SMSpp(fic, T, dT, length(grp)+idxV-1, idxV, vallees, usines, reservoirs);
fprintf(fic, '\n');
end
        end
% -- terminate UC block

fprintf(fic, '} // group Block_0\n');
fprintf(fic, '}\n');
fprintf(fic, '\n');
fclose(fic);

end


