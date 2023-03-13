function [Pq] = powerCalc()
%This function calculates the power generated over the mission lifetime by 
%utilizing the powerVsAngle and the faceOrientation functions
clear all
close all
clc

%Load all of the mission trajectory data
S = load('BallisticLunarTraj.mat');
%Define the time array
tarray = S.tarray;
%tarray = S.tarray2;
%Define the instantaneous satelite position vectors array
spacecraftpos = S.spacecraftpos;
%spacecraftpos = S.spacecraftpos2;
%Define the instantaneous satelite velocity vectors array
spacecraftvel = S.spacecraftvel;
%spacecraftvel = S.spacecraftvel2;
%Define the instantaneous sun position vectors array
sunpos_ECI = S.sunpos_ECI;
%sunpos_ECI = S.sunpos2;

%Redefine the satellite position, satellite velocity, and sun position
%arrays and transpose them
r = spacecraftpos';
v = spacecraftvel';
r_sun = sunpos_ECI';

G = 6.674*10^-11;  %Earth Gravitational constant [N(m/kg)^2]
M = 5.972*10^24;   %Earth Mass [kg]

mu = G*M;

%Spin rate of the satellite
w_s = 2*pi;

%Define number of time steps in mission life
n = length(tarray);

%Convert position and velocity arrays into vectors
rs = zeros(3*n,1);
vs = zeros(3*n,1);
for j = 1:n
    rs(3*j-2:3*j) = r(j,:);
    vs(3*j-2:3*j) = v(j,:);
end

%Call the function that takes velocity and position vectors and returns the
%orbital elements
[a,e,E,i,w,Om,P,tau,A,B] = vec2orbElem(rs,vs,mu);
nu = 2.*atan2(sqrt(1+e).*cos(E/2),sqrt(1-e).*cos(E/2));

%Defined the vector from the sun to the satellite
r_sunsat = r_sun - r;

%Given the orbital elements, find the orientation of the satellite in ECI
%reference frame
[N,R] = body2ECI(Om,nu,w,i,w_s,tarray,n);

%Call the powerVsAngle function to get power as a function of the angle of
%incidence
[th,P] = powerVsAngle();

%Find the actual angle of incidence at every time step by finding the
%cosine and sine of the angle and then using the atan2 function
cosq = dot(r_sunsat,N.N2,2)./(magn(N.N2).*magn(r_sunsat));
sinq = magn(cross(r_sunsat,N.N2,2))./(magn(N.N2).*magn(r_sunsat));
thq = atan2(sinq,cosq);

%Check angle step to make sure it is no more than 7 degrees...if it is you
%need to get better data that has a smaller change in angle between the
%time steps
diffthq = [];
for j = 1:n-1
    diffthq = [ diffthq 180/pi*(abs(thq(j) - thq(j+1)))];
    if diffthq(j) > 7
        diffthq(j)
    end
end

%Use MATLABs interpolation function to interpolate a power value that
%corresponds to the angle at each time step using the relation derived in
%powerVsAngle function
Pq = interp1(th,P,thq);

%Plot the power produced over the mission life
figure
plot(tarray,Pq)
title('Power Production over Mission Life')
xlabel('Time (s)')
ylabel('Power (W)')
axis([0 tarray(end) 0 8])
whitebg('k')

end
%%
function [N,R] = body2ECI(Om,nu,w,i,w_s,tarray,n)
%This function takes the orbital elements of the trajectory and the spin of
%the satellite and determines the orientation of the faces in the ECI
%coordinate system

for j = 1:n
    %Rotation matrix from the orbital elements
    B_Q_N(:,:,j) = [cos(Om(j)-pi).*cos(w(j)+nu(j))-cos(i(j)).*sin(Om(j)-pi).*sin(w(j)+nu(j)), ...
            -cos(Om(j)-pi).*sin(w(j)+nu(j))-cos(i(j)).*cos(w(j)+nu(j)).*sin(Om(j)-pi), ...
            sin(Om(j)-pi).*sin(i(j)); ...
            cos(w(j)+nu(j)).*sin(Om(j)-pi)+cos(Om(j)-pi).*cos(i(j)).*sin(w(j)+nu(j)),...
            cos(Om(j)-pi).*cos(i(j)).*cos(w(j)+nu(j))-sin(Om(j)-pi).*sin(w(j)+nu(j)),...
            -cos(Om(j)-pi).*sin(i(j));...
            sin(i(j)).*sin(w(j)+nu(j)), cos(w(j)+nu(j)).*sin(i(j)), cos(i(j))];
    %Spin of the satellite as a function of time
    spin = mod(w_s.*tarray(j),2*pi);
    %Rotation matrix from the satellite spin
    orien_B(:,:,j) = [cos(spin) 0 sin(spin);0 1 0; -sin(spin) 0 cos(spin)];
end

%The normal vectors to each face 
N1 = zeros(n,3);
N2 = zeros(n,3);
N3 = zeros(n,3);
for j = 1:n
    %Total rotation matrix at each time step
    R(:,:,j) = orien_B(:,:,j)*B_Q_N(:,:,j);
    %Rotate the face normal vectors
    N1(j,:) = R(1,:,j);
    N2(j,:) = R(2,:,j);
    N3(j,:) = R(3,:,j);
end

%Define the face vectors that are opposite
N4 = -N1;
N5 = -N2;
N6 = -N3;

%Output normal vectors in a structure array
N = struct('N1',N1,'N2',N2,'N3',N3,'N4',N4,'N5',N5,'N6',N6);

end
%%
function [mag] = magn(mat)
%Find the magnitude of each of n vectors from a matrix of nx3 vectors
mag = sum(mat.^2,2).^0.5;

end