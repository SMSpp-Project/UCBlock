function reservoirs = determineHydroConsistency( idxV, dT, T, vallees, usines, reservoirs )
                      %
                      %
                      %
                      deltat = dT*3600;

% for each reservoir compute the reservoir level
nbR = length(vallees(idxV).reservoirs);
nbA = length(vallees(idxV).arc);

volumes = zeros(nbR, T);
for iR = 1: nbR,
idxR = vallees(idxV).reservoirs(iR);

volumes(iR,:) = reservoirs(idxR).v0 + cumsum(reservoirs(idxR).inflows*deltat);
end

%now go through all arcs
for iA=1:nbA,
iUsi=vallees(idxV).arc(iA).usiIdx;

% -- round initflow data
        usiFlow = round( usines(iUsi).initflow*1e6 ) / 1e6;

iAm =vallees(idxV).arc(iA).amont;
iAv =vallees(idxV).arc(iA).aval;

tauup = ceil(vallees(idxV).arc(iA).updelay / (dT*60.0) );
taudn = ceil(vallees(idxV).arc(iA).dndelay / (dT*60.0) );

%
% -- this usine iUsi flows (unless it is a pump) -- from iR to I(j)
%
if (usines(iUsi).type > 0),
% then this is a pump
volumes(iAm,:) = volumes(iAm,:) + cumsum(usiFlow*deltat);
if ( iAv <= nbR ),
volumes(iAv,:) = volumes(iAv,:) - cumsum(usiFlow*deltat);
end
else
% this is a turbine
flux = cumsum([zeros(1,tauup) usiFlow*deltat]);
volumes(iAm,:) = volumes(iAm,:) - flux(1:T);
if ( iAv <= nbR ),
flux = cumsum([zeros(1,taudn) usiFlow*deltat]);
volumes(iAv,:) = volumes(iAv,:) + flux(1:T);
end
        end
end

% figure out if the volumes work out with the bounds
for iR = 1: nbR,
idxR = vallees(idxV).reservoirs(iR);

Iup = find( reservoirs(idxR).vmax - volumes(iR,:) < -1e-8);
Idn = find( volumes(iR,:) - reservoirs(idxR).vmin < -1e-8);

if ( ~isempty(Iup) || ~isempty(Idn) ),
%plot([1:T],reservoirs(idxR).vmin, 'k', [1:T],volumes(iR,:), 'b', [1:T],reservoirs(idxR).vmax, 'k');
% -- update minimal reservoir level first:
reservoirs(idxR).vmin = max( min(volumes(iR,:),reservoirs(idxR).vmin), zeros(1,T) );

Idn = find( volumes(iR,:) - reservoirs(idxR).vmin < -1e-8);
% -- if still not empty
if ( ~isempty( Idn ) ),
% -- move initial level up
vlack = min(volumes(iR,:) - reservoirs(idxR).vmin);

reservoirs(idxR).v0 = reservoirs(idxR).v0 - vlack;

volumes(iR,:) = volumes(iR,:) + abs(vlack);
end

% -- check upper volumes
        Iup = find( reservoirs(idxR).vmax - volumes(iR,:) < -1e-8);
if ( ~isempty(Iup) ),
reservoirs(idxR).vmax = max( reservoirs(idxR).vmax, volumes(iR,:) ); %max(1.1*max( reservoirs(idxR).vmax, volumes(iR,:) ), 0.1*ones(1,T));
end

% -- check final volumes
if ( reservoirs(idxR).vfinmin > volumes(iR,end) )
reservoirs(idxR).vfinmin = min( volumes(iR,end), reservoirs(idxR).vmin(end) );
end
if ( reservoirs(idxR).vfinmax < volumes(iR,end) )
reservoirs(idxR).vfinmax = max( volumes(iR,end), reservoirs(idxR).vmax(end) );
end

%figure;
%plot([1:T],reservoirs(idxR).vmin, 'k', [1:T],volumes(iR,:), 'b', [1:T],reservoirs(idxR).vmax, 'k');
disp(strcat('Probleme sur Reservoir', reservoirs(idxR).name));
%close all;
end
        end

%
