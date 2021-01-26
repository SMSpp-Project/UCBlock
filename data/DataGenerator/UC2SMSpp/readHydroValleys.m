%%%%%%
%
% script to read the data from APOGEE related to Hydro generation
%
function [vallees, usines, reservoirs] = readHydroValleys( pathName ),

        ficRESERVES = [pathName '/RESERVES'];
ficUSINE    = [pathName '/USINES'];
ficVALLEES  = [pathName '/VALLEES'];


usines = struct('name',{},'type',{},'prim_pour',{}, 'gradup', {}, 'graddn',{}, 'pdfset',{}, 'initflow',{}, 'd0',{} );
reservoirs = struct('name',{},'inflows',{},'vmax',{},'vmin',{},'v0',{},'vfinmin',{},'vfinmax',{},'wpart',{},'wvalues',{});
arcs = struct('amont',{},'aval',{},'updelay',{},'dndelay',{},'usiIdx',{});

vallees = struct('name',{},'reservoirs',{},'flow',{},'arc',{});

% read information from the reservoirs into the data structures
        reservoirs = readRESERVES( ficRESERVES, reservoirs );

%
usines = readUSINES( ficUSINE, usines );

%
vallees = readVALLEES( ficVALLEES, vallees, arcs, reservoirs, usines );


