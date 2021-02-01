clear all, clc

% My Code
% ____________________________________________________________________

fid = fopen('TLE.txt');
line1=fgetl(fid);
line2=fgetl(fid);
fclose(fid);

SatNum=line1(3:7);
SatClass=line1(8);
LaunchYear=line1(10:11);
LaunchNum=line1(12:14);
LaunchPiece=line1(15:17);
EpochYear=line1(19:20);
EpochDay=line1(21:32);
MeanMotionFirstDeriv=line1(34:43);
MeanMotionSecDeriv=line1(45:52);
BSTARdrag=line1(54:61);
EphType=line1(63);
ElementNum=line1(65:68);
CheckSum=line1(69);
Inc=line2(9:16);
RAAN=line2(18:25);
Ecc=line2(27:33);
AP=line2(35:42);
MA=line2(44:51);
MM=line2(53:63);
RevNum=line2(64:68);
% ____________________________________________________________________








% CITATION: Code taken from https://www.mathworks.com/matlabcentral/fileexchange/72098-keplerian-orbit
% ____________________________________________________________________
% Input Orbit parameters

i = str2double(Inc);          % Inclination [°]
O = str2double(RAAN);          % Right Ascension of Ascending Node [°]
e = str2double(Ecc)*10^(-7);  % Eccentricity [-]
o = str2double(AP);          % Argument of perigee [°] 
M0 = str2double(MA);         % Mean anomaly [°]
n = str2double(MM);          % [revolutions/day]

% Input Constants

mu =398600.4418;    % Geocentric gravitational constant [km^3*s^(-2)] 
Rz = 6378;          % Equatorial radius of Earth [km]

% Convert degrees to radians

O = O*pi/180;      %[rad]
i = i*pi/180;      %[rad]
M0 = M0*pi/180;    %[rad]
o = o*pi/180;      %[rad]

% Calculations
n2=n*2*pi/86164;                   % n to [rad/s]
a=mu^(1/3)/n2^(2/3);               % Major semi-axis [km]
p  = a*(1-e^2);                    % Semi-latus rectus [km]
rp = a*(1-e);                      % Radius of perigee [km]
ra = a*(1+e);                      % Radius of apogee [km]
vp = sqrt(mu*(2/rp-1/a));          % Velocity at perigee km/s]
va = sqrt(mu*(2/ra-1/a));          % Velocity at apogee [km/s]
T  = 2*pi/n2;                      % period [s]
hodiny   = floor(T/3600);                   % hours        
minuty = floor((T-hodiny*3600)/60);         % minutes        
sekundy = floor(T-hodiny*3600-minuty*60);   % seconds  
% ____________________________________________________________________













% My code:
% ____________________________________________________________________

[ea, ta]=kepler1(M0,e);

[r,v] = orb2rv(p,e,i,O,o,ta);

fileID = fopen('KERV.txt','w');
fprintf(fileID,'Position in ECI (km):\n');
fprintf(fileID,'[%i\t%i\t%i]\n\n',r(1),r(2),r(3));
fprintf(fileID,'Velocity in ECI (km/s):\n');
fprintf(fileID,'[%i\t%i\t%i]\n',v(1),v(2),v(3));
fclose(fileID);


% ____________________________________________________________________










% CREDIT: CODE TAKEN FROM https://www.mathworks.com/matlabcentral/fileexchange/35455-convert-keplerian-orbital-elements-to-a-state-vector
% ____________________________________________________________________



%----------------------- Begin Code Sequence -----------------------------%
% Purpose:                                                                %
% Convert a given set of Keplerian Orbital Elements to state vectors in   %
% in the ECI frame of reference                                           %
%                                                                         %
%-------------------------------------------------------------------------%
%                                                                         %
% Inputs:                                                                 %
%--------                                                                  
%p                      [1 x N]                         Semilatus Rectum
%                                                       (km)
%
%e                      [1 x N]                         Eccentricity Magnitude
%                                                       (unitless)
%
%i                      [1 x N]                         inclination
%                                                       (radians)
%
%O                      [1 x N]                         Right Ascention of
%                                                       the ascending node
%                                                       (radians)
%
%o                      [1 x N]                         Argument of perigee
%                                                       (radians)
%
%nu                     [1 x N]                         True Anomaly
%                                                       (radians)
%
%truLon                 [1 x N]                         True Longitude
%                                                       (radians)
%
%argLat                 [1 x N]                         Argument of Latitude
%                                                       (radians)
%
%lonPer                 [1 x N]                         Longitude of Periapse
%                                                       (radians)
%
%mu                     double                          Gravitational Constant
%                                                       Defaults to Earth if
%                                                       not specified
%
% Outputs:
%---------                                                                %
%r_ECI                  [3 x N]                         Position Vector in
%                                                       ECI coordinate
%                                                       frame of reference
%
%v_ECI                  [3 x N]                         Velocity vector in
%                                                       ECI coordinate
%                                                       frame of reference
%
% References:
%-------------
%Vallado,D. Fundamentals of Astrodynamics and Applications. 2007.
%
% Function Dependencies:
%------------------
%None
%------------------------------------------------------------------       %
% Programed by Darin Koblick  03-04-2012                                  %
%------------------------------------------------------------------       %
function [r,v] = orb2rv(p,e,i,O,o,nu,truLon,argLat,lonPer,mu)
if ~exist('mu','var');  mu = 398600.4418; end
%Make all inputs consistent w/ dimensions
p = p(:); e = e(:); i = i(:); O = O(:); o = o(:); nu = nu(:);
if exist('truLon','var')
truLon = truLon(:); argLat = argLat(:); lonPer = lonPer(:);
end
ietol = 1e-8;
idx = e < ietol & mod(i,pi) < ietol;
if any(idx); o(idx) = 0; O(idx) = 0; nu(idx) = truLon(idx); end
idx = e < ietol & mod(i,pi) > ietol;
if any(idx); o(idx) = 0; nu(idx) = argLat(idx); end
idx = e > ietol & mod(i,pi) < ietol;
if any(idx); O(idx) = 0; o(idx) = lonPer(idx); end
%Find rPQW and vPQW
rPQW = cat(2,p.*cos(nu)./(1 +e.*cos(nu)),p.*sin(nu)./(1+e.*cos(nu)),zeros(size(nu)));
vPQW = cat(2,-sqrt(mu./p).*sin(nu),sqrt(mu./p).*(e+cos(nu)),zeros(size(nu)));
%Create Transformation Matrix
PQW2IJK = NaN(3,3,size(p,1));
cO = cos(O); sO = sin(O); co = cos(o); so = sin(o); ci = cos(i); si = sin(i);
PQW2IJK(1,1,:) = cO.*co-sO.*so.*ci; PQW2IJK(1,2,:) = -cO.*so-sO.*co.*ci; PQW2IJK(1,3,:) = sO.*si;
PQW2IJK(2,1,:) = sO.*co+cO.*so.*ci; PQW2IJK(2,2,:) = -sO.*so+cO.*co.*ci; PQW2IJK(2,3,:) = -cO.*si;
PQW2IJK(3,1,:) = so.*si;            PQW2IJK(3,2,:) = co.*si;             PQW2IJK(3,3,:) = ci;
%Transform rPQW and vPQW to rECI and vECI
r = multiDimMatrixMultiply(PQW2IJK,rPQW)';  
v = multiDimMatrixMultiply(PQW2IJK,vPQW)';
end
function c = multiDimMatrixMultiply(a,b)
c = NaN(size(b));
c(:,1) = sum(bsxfun(@times,a(1,:,:),permute(b,[3 2 1])),2);
c(:,2) = sum(bsxfun(@times,a(2,:,:),permute(b,[3 2 1])),2);
c(:,3) = sum(bsxfun(@times,a(3,:,:),permute(b,[3 2 1])),2);
end
% ____________________________________________________________________










% CREDIT: CODE TAKEN FROM https://www.mathworks.com/matlabcentral/fileexchange/39031-kepler-s-equation
% ____________________________________________________________________
function [eanom, tanom] = kepler1 (manom, ecc)
% solve Kepler's equation for circular,
% elliptic and hyperbolic orbits
% Danby's method
% input
%  manom = mean anomaly (radians)
%  ecc   = orbital eccentricity (non-dimensional)
% output
%  eanom = eccentric anomaly (radians)
%  tanom = true anomaly (radians)
% Orbital Mechanics with Matlab
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global niter
% define convergence criterion
ktol = 1.0e-10;
pi2 = 2 * pi;
xma = manom - pi2 * fix(manom/pi2);
% initial guess
if (ecc == 0)
    
    % circular orbit
    tanom = xma;
    
    eanom = xma;
    
    return;
    
elseif (ecc < 1)
    
    % elliptic orbit
    eanom = xma + 0.85 * sign(sin(xma)) * ecc;
else
    
    % hyperbolic orbit
    eanom = log(2 * xma / ecc + 1.8);
    
end
% perform iterations
niter = 0;
while(1)
    
    if (ecc < 1)
        % elliptic orbit
        s = ecc * sin(eanom);
        c = ecc * cos(eanom);
        f = eanom - s - xma;
        fp = 1 - c;
        fpp = s;
        fppp = c;
        
    else
        % hyperbolic orbit
        s = ecc * sinh(eanom);
        c = ecc * cosh(eanom);
        f = s - eanom - xma;
        fp = c - 1;
        fpp = s;
        fppp = c;
        
    end
    niter = niter + 1;
    % check for convergence
    if (abs(f) <= ktol || niter > 20)
        
        break;
        
    end
    % update eccentric anomaly
    delta = -f / fp;
    deltastar = -f / (fp + 0.5 * delta * fpp);
    deltak = -f / (fp + 0.5 * deltastar * fpp ...
        + deltastar * deltastar * fppp / 6);
    eanom = eanom + deltak;
    
end
if (niter > 20)
    
    clc; home;
    
    fprintf('\n\n   more than 20 iterations in kepler1 \n\n');
    
    keycheck;
    
end
% compute true anomaly
if (ecc < 1)
    
    % elliptic orbit
    sta = sqrt(1 - ecc * ecc) * sin(eanom);
    
    cta = cos(eanom) - ecc;
    
else
    
    % hyperbolic orbit
    sta = sqrt(ecc * ecc - 1) * sinh(eanom);
    
    cta = ecc - cosh(eanom);
    
end
tanom = atan3(sta, cta);
end

function y = atan3 (a, b)
% four quadrant inverse tangent
% input
%  a = sine of angle
%  b = cosine of angle
% output
%  y = angle (radians; 0 =< c <= 2 * pi)
% Orbital Mechanics with Matlab
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
epsilon = 0.0000000001;
pidiv2 = 0.5 * pi;
if (abs(a) < epsilon)
    y = (1 - sign(b)) * pidiv2;
    return;
else
    c = (2 - sign(a)) * pidiv2;
end
if (abs(b) < epsilon)
    y = c;
    return;
else
    y = c + sign(a) * sign(b) * (abs(atan(a / b)) - pidiv2);
end
end
% ____________________________________________________________________