function [th,P] = powerVsAngle()
%This function creates a relation between the power produced from the solar
%cells and the angle of incidence that the sun has on the spinning face.
clc

%Create a vector for the distance from the sun to the earth in ECI
%coordinates
r_sunear = 149597870700; %[m]
r = r_sunear*[0 1 0];

%Run the initial power calculation for this vector
[P] = faceOrientation(r);

%Define the end point of the function (only need to flip it around 180
%degrees since it is symmetric) and the step size for the change in angle.
thend = pi;
dth = 0.005;

%Calculate average power at each delta angle step
for th = 0:dth:thend
    %Rotate the vector to the sun by the new angle
    r_sunsat1 = r*[cos(th) sin(th) 0; -sin(th) cos(th) 0; 0 0 1];
    %Call function that calculates the average power at that angle
    [P_avg] = faceOrientation(r_sunsat1);
    %Add calculated power to the array
    P = [P P_avg];
end

%Define corresponding angle vector
th = linspace(0,thend,thend/dth + 2);

%Plot the function of average power versus angle of incidence on the
%spining face that has just been calculated
figure
plot(th,P)
title('Power Average over Satellite Spin')
xlabel('Angle of Incidence from the spinning Axis (rad)')
ylabel('Average Power (W)')
axis([0 pi 0 7])
whitebg('k')

end

