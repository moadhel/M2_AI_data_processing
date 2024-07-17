function display_position(file, llh, color, sat_llh)

% Create a kml trajectory file and open it in GoogleEarth
%
% simple_kml(file, t, llh);
% OR
% simple_kml(file, t, llh, color);
% simple_kml(file, t, llh, color, ephemeris);
% simple_kml(file, t, llh, color, ephemeris, mes_gps);
%
%
% Inputs :
%	file : Name of the kml file
%
%   t :    [1,N] vector or scalar ("t" and "llh" sizes must match)
%          Time of the positions
%          (use GPS TOW if you want to display the satellites positions)
%
%   llh :  [3,N] matrix or [3,1] vector ("t" and "llh" sizes must match)
%          LLH coordinates of the trajectory
%          (Lat, Lon in degrees; Height in m)
%
%   color : color of the trajectory (char string)
%
%   ephemeris : [1,M] struct containing the ephemeris or the almanac
%               of the satellites (must be in the ISAE standard format)
%
%   mes_gps : struct containing the raw GPS measurements
%             (must be in the ISAE standard format)
%
%
% Droits réservés ISAE.
% Contact : Benoit Priot
%             benoit.priot@isae.fr
%             (+33)5 61 33 83 66
%
% version 2013_11
%



% creation de l'objet avec le nom du fichier kml
traj = kml_creator(file);






%--------------------------------------------------------------------------
% options du kml
if (nargin < 3)
    traj.color = 'white';    % 'white' par défaut
else
    traj.color = color;    % 'white' par défaut
end

% référence de l'altitude (sea, ground, 2D)
traj.reference = '2D';     % 'sea' par défaut

% correction des coordonnées pour l'affichage
traj.lat_corr = 0;       % (en °)   0 par défaut
traj.lon_corr = 0;       % (en °)   0 par défaut
traj.h_corr = -40;       % (en m) -40 par défaut
%--------------------------------------------------------------------------


% coordonnées WGS84 LLH (satellite fixe)
if ((nargin >= 4) && ~isempty(sat_llh))
for m = 1:4
    sat.SVID = m;             % id du satellite
    sat.llh = sat_llh(:,m);    % lat, lon en deg, h en m
    traj = add_satellite(traj, sat);
end
end

%traj = add_satellite(traj, ephemeris);


% ajout des positions de la trajectoire
if (size(llh, 1) ~= 3)
    llh = llh';
end
if (size(llh, 1) ~= 3)
    error('Invalid input ''llh'' size.');
end


for k = 1:size(llh, 2)
	traj = add_position(traj, k, llh(:,k), [1:4]);
end


% création du fichier
try
create_and_view(traj);
end

end


