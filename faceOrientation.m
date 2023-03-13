function [P_avg] = faceOrientation(r_sunsat)
clc

%The rate at which the satelite is spinning about it's spin axis
w_s = 2*pi;

%Number of steps
n = 1413;

%Step size
dt = 0.01;

%Set orbital elements and define time vector
Om = zeros(n,1);
w = zeros(n,1);
nu = zeros(n,1);
i = zeros(n,1);
tarray = linspace(0,dt*n,n)';

%% CHANGE THIS SECTION %%
%Uncomment lines 25 and 26 to run code separately.  Uncomment line 29 to
%call from powerVsAngle code

%Code to run separately
% r_sunear = 1; %[m]
% r_sunsat = r_sunear*[ones(n,1) ones(n,1) zeros(n,1)];

%Code to run from powerVsAngle code
r_sunsat = repmat(r_sunsat,n,1);

%% Create rotation matrix from orbital elements and spin of satellite

%Initialize rotation matrix at each step for the orbital elements
B_Q_N = zeros(3,3,n);
%Initialize rotation matrix at each step for the satelite spin
orien_B = zeros(3,3,n);

%Create rotation matrices for the orbital rotatiin and the satelite spin
for j = 1:n
    B_Q_N(:,:,j) = [cos(Om(j)-pi).*cos(w(j)+nu(j))-cos(i(j)).*sin(Om(j)-pi).*sin(w(j)+nu(j)), ...
            -cos(Om(j)-pi).*sin(w(j)+nu(j))-cos(i(j)).*cos(w(j)+nu(j)).*sin(Om(j)-pi), ...
            sin(Om(j)-pi).*sin(i(j)); ...
            cos(w(j)+nu(j)).*sin(Om(j)-pi)+cos(Om(j)-pi).*cos(i(j)).*sin(w(j)+nu(j)),...
            cos(Om(j)-pi).*cos(i(j)).*cos(w(j)+nu(j))-sin(Om(j)-pi).*sin(w(j)+nu(j)),...
            -cos(Om(j)-pi).*sin(i(j));...
            sin(i(j)).*sin(w(j)+nu(j)), cos(w(j)+nu(j)).*sin(i(j)), cos(i(j))];
    spin = mod(w_s.*tarray(j),2*pi);
    orien_B(:,:,j) = [cos(spin) 0 sin(spin);0 1 0; -sin(spin) 0 cos(spin)];
end

%Initialize the normal vectors to each face
N1 = zeros(n,3);
N2 = zeros(n,3);
N3 = zeros(n,3);

%Initialize the entire rotation matrix
R = zeros(3,3,n);

%Rotate the normal vectors by the rotation matrix at each step
for j = 1:n
    R(:,:,j) = orien_B(:,:,j)*B_Q_N(:,:,j);
    N1(j,:) = R(1,:,j);
    N2(j,:) = R(2,:,j);
    N3(j,:) = R(3,:,j);
end

%Define the opposite faces' normal vectors
N4 = -N1;
N5 = -N2;
N6 = -N3;

%The body coordinate system in ECI coordinates...X and Z are unused but
%defined for clarity
% X = N1./repmat(magn(N1),1,3);
Y = N2./repmat(magn(N2),1,3);
% Z = N3./repmat(magn(N3),1,3);

%Dot product of vector from satelite to sun and normal vector of each
%surface of the satelite
dot1 = dot(r_sunsat,N1,2);
dot2 = dot(r_sunsat,N2,2);
dot3 = dot(r_sunsat,N3,2);
dot4 = dot(r_sunsat,N4,2);
dot5 = dot(r_sunsat,N5,2);
dot6 = dot(r_sunsat,N6,2);

%Cross product of vector from satelite to sun and normal vector of each
%surface of the satelite
cross1 = cross(r_sunsat,N1,2);
cross2 = cross(r_sunsat,N2,2);
cross3 = cross(r_sunsat,N3,2);
cross4 = cross(r_sunsat,N4,2);
cross5 = cross(r_sunsat,N5,2);
cross6 = cross(r_sunsat,N6,2);

%Magitude of the sun-satelite vector times the normal vector magnitude
mag1 = magn(r_sunsat).*magn(N1);
mag2 = magn(r_sunsat).*magn(N2);
mag3 = magn(r_sunsat).*magn(N3);
mag4 = magn(r_sunsat).*magn(N4);
mag5 = magn(r_sunsat).*magn(N5);
mag6 = magn(r_sunsat).*magn(N6);

%Sin(theta) where theta is the angle between the sun-satelite vector and
%the normal face vector
sinth1 = magn(cross1)./mag1;
sinth2 = magn(cross2)./mag2;
sinth3 = magn(cross3)./mag3;
sinth4 = magn(cross4)./mag4;
sinth5 = magn(cross5)./mag5;
sinth6 = magn(cross6)./mag6;

%Cos(theta) where theta is the angle between the sun-satelite vector and
%the normal face vector
costh1 = dot1./mag1;
costh2 = dot2./mag2;
costh3 = dot3./mag3;
costh4 = dot4./mag4;
costh5 = dot5./mag5;
costh6 = dot6./mag6;

%Solve for the angle between the sun-sat vector and normal face vector
th1 = atan2(sinth1,costh1);
th2 = atan2(sinth2,costh2);
th3 = atan2(sinth3,costh3);
th4 = atan2(sinth4,costh4);
th5 = atan2(sinth5,costh5);
th6 = atan2(sinth6,costh6);

%Component of the sun-sat vector along the normal to the face direction
L1 = repmat(dot1./magn(N1),1,3).*N1;
L3 = repmat(dot3./magn(N3),1,3).*N3;

crossphi1 = cross((r_sunsat-L1),Y,2);
crossphi3 = cross((r_sunsat-L3),-Y,2);

dotphi1 = dot((r_sunsat-L1),Y,2);
dotphi3 = dot((r_sunsat-L3),-Y,2);

magphi1 = magn(r_sunsat-L1).*magn(Y);
magphi3 = magn(r_sunsat-L3).*magn(-Y);

%The sine function of the two face angles
sinphi1 = magn(crossphi1)./magphi1;
sinphi3 = magn(crossphi3)./magphi3;

%The cosine function of the two face angles
cosphi1 = dotphi1./magphi1;
cosphi3 = dotphi3./magphi3;

%Angle on the plane of the face from the faux x coordinate
phi1 = atan2(sinphi1,cosphi1);
phi3 = atan2(sinphi3,cosphi3);

phi1(isnan(phi1)) = 0;
phi3(isnan(phi3)) = 0;

%% Now given the angles of incidence and the clockwise angles, determine
%% when faces are shadowed or exposed to the sun

%Initialize percent area on each face exposed to the sun
F1a = zeros(n,1);
F1b = zeros(n,1);
F2 = zeros(n,1);
F3a = zeros(n,1);
F3b = zeros(n,1);
F4 = zeros(n,1);
F5 = zeros(n,1);
F6 = zeros(n,1);
for j = 1:n
    %If angle of incidence is less than 90
    if th1(j) < pi/2
        %Unshadowed side is shown
        F1b(j) = 1;
        %Find height of trapezoidal shadow
        h = 1.5*tan(th1(j));
        %Adjust height if it goes off the edge
        if h > 1
            h = 1;
        end
        %If clockwise angle is casting a shadow
        if pi < phi1(j) < 2*pi
            %Modify angle to be less than pi
            phi1(j) = mod(phi1(j),pi);
            %If angle is close to 90 so that it will cause the tan function
            %to give unreasonable results
            if abs(phi1(j)-pi/2) < 0.1;
                %If so, approximate area as a rectangle
                A_s = h;
            elseif tan(phi1(j))<h
                A_s = 0.5*tan(phi1(j));
            else
                %Area of trapezoidal shadow
                a = 1 - h/tan(phi1(j));
                A_s = 0.5*(a+1)*h;
            end
            %Percent Area that is exposed
            F1a(j) = 1 - A_s;
        else
            %if the clockwise angle does not throw a shadow, face is
            %exposed
            F1a(j) = 1;
        end
    end
    %If angle of incidence is less than 90
    if th3(j) < pi/2
        %Unshadowed side is shown
        F3b(j) = 1;
        %Find height of trapezoidal shadow
        h = tan(th3(j));
        %Adjust height if it goes off the edge
        if h > 1.5
            h = 1.5;
        end
        %If clockwise angle is casting a shadow
        if pi < phi3(j) < 2*pi
            %Modify angle to be less than pi
            phi3(j) = mod(phi3(j),pi);
            %If angle is close to 90 so that it will cause the tan function
            %to give unreasonable results
            if abs(phi3(j)-pi/2) < 0.1;
                %If so, approximate area as a rectangle
                A_s = h;
            elseif tan(phi3(j))<h
                A_s = 0.5*tan(phi3(j));
            else
                %Area of trapezoidal shadow
                A_s = h/2*(2-h/tan(phi3(j)));
            end
            %Percent Area that is exposed
            F3a(j) = 1 - A_s/1.5;        
        else
            %if the clockwise angle does not throw a shadow, face is
            %exposed
            F3a(j) = 1;
        end
    end
    %If angle of incidence is less than 90
    if th2(j) < pi/2
        %Face is exposed
        F2(j) = 1;
    end
    %If angle of incidence is less than 90
    if th4(j) < pi/2
        %Face is exposed
        F4(j) = 1;
    end
    %If angle of incidence is less than 90
    if th5(j) < pi/2
        %Face is exposed
        F5(j) = 1;
    end
    %If angle of incidence is less than 90
    if th6(j) < pi/2
        %Face is exposed
        F6(j) = 1;
    end  
end

%% THIS IS THE SECTION TO BE CHANGED %%

%Once the solar panels have been designed, edit the numCells variables
%below to correspond to the correct number of TASC cells that are able to
%be fit on each face

%Area of one TASC cell
A_cell = 2.277*10^(-4);  %m^2

%Number of cells per face (CHANGE THESE)
numCells1a = 15;
numCells1b = 15;
numCells2  = 15;
numCells3a = 15;
numCells3b = 15;
numCells4  = 15;
numCells5  = 15;
numCells6  = 15;

%Area of each face that is covered in solar panels
A1a = A_cell*numCells1a;
A1b = A_cell*numCells1b;
A2  = A_cell*numCells2;
A3a = A_cell*numCells3a;
A3b = A_cell*numCells3b;
A4  = A_cell*numCells4;
A5  = A_cell*numCells5;
A6  = A_cell*numCells6;

%Unit length of the cubesat (DELETE THIS ONCE SOLAR PANELS HAVE BEEN
%DESIGNED)
U = 0.1;  %m

%Surface area of each face of the cubesat (DELETE THIS ONCE SOLAR PANELS
%HAVE BEEN DESIGNED)
A1a = 1*U^2;
A1b = 0.5*U^2;
A2  = 3*U^2;
A3a = 1.5*U^2;
A3b = 1.5*U^2;
A4  = 1.5*U^2;
A5  = 3*U^2;
A6  = 3*U^2;

%Percentage of the faces that are covered by solar panels (DELETE THIS ONCE
%SOLAR PANELS HAVE BEEN DESIGNED)
PerFace = 0.6;

%% Now consider how much power is generated for the given angles

%Maximum power produced by one of the TASC cells
P_pan = 0.027*100^2; %W/cm^2

%Power generated per face considering the angle of inceidence and the
%clockwise angle
P1a = P_pan.*F1a.*A1a.*costh1*PerFace;
P1b = P_pan.*F1b.*A1b.*costh1*PerFace;
P2  = P_pan.*F2.*A2.*costh2*PerFace;
P3a = P_pan.*F3a.*A3a.*costh3*PerFace;
P3b = P_pan.*F3b.*A3b.*costh3*PerFace;
P4  = P_pan.*F4.*A4.*costh4*PerFace;
P5  = P_pan.*F5.*A5.*costh5*PerFace;
P6  = P_pan.*F6.*A6.*costh6*PerFace;

%Initialize all variables for the average power per face
P1a_avg = 0;
P1b_avg = 0;
P2_avg  = 0;
P3a_avg = 0;
P3b_avg = 0;
P4_avg  = 0;
P5_avg  = 0;
P6_avg  = 0;

%Find the average power generated by averaging the power generated between
%two steps
for j = 2:n
    P1a_avg = P1a_avg + (P1a(j) + P1a(j-1))/2*dt;
    P1b_avg = P1b_avg + (P1b(j) + P1b(j-1))/2*dt;
    P2_avg  = P2_avg + (P2(j) + P2(j-1))/2*dt;
    P3a_avg = P3a_avg + (P3a(j) + P3a(j-1))/2*dt;
    P3b_avg = P3b_avg + (P3b(j) + P3b(j-1))/2*dt;
    P4_avg  = P4_avg + (P4(j) + P4(j-1))/2*dt;
    P5_avg  = P5_avg + (P5(j) + P5(j-1))/2*dt;
    P6_avg  = P6_avg + (P6(j) + P6(j-1))/2*dt;
end

%Sum the power over all faces
P_avg = (P1a_avg + P1b_avg + P2_avg + P3a_avg + P3b_avg + P4_avg + P5_avg + P6_avg)/(dt*n);

%Change show plot variable if you wish to see the power generated on each
%face plotted individually
showplot = 0;
if showplot == 1
    figure(1)
    subplot(3,3,1),plot(tarray,P1a)
    title('Face 1a')
    subplot(3,3,4),plot(tarray,P1b)
    title('Face 1b')
    subplot(3,3,7),plot(tarray,P2)
    title('Face 2')
    subplot(3,3,2),plot(tarray,P3a)
    title('Face 3a')
    subplot(3,3,5),plot(tarray,P3b)
    title('Face 3b')
    subplot(3,3,8),plot(tarray,P4)
    title('Face 4')
    subplot(3,3,3),plot(tarray,P5)
    title('Face 5')
    subplot(3,3,6),plot(tarray,P6)
    title('Face 6')
    subplot(3,3,9),plot(tarray,repmat(P_avg,n,1))
    title('Average Power')
end

end

function [mag] = magn(mat)
%Find the magnitude of each of n vectors from a matrix of nx3 vectors
mag = sum(mat.^2,2).^0.5;

end