function axes = satelliteAxes(z_input)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%normalize z
z = z_input / norm(z_input);

% Choose a random vector v that is not parallel to z
if z(1) ~= 0 || z(2) ~= 0
    v = [z(2), -z(1), 0];
else
    v = [0, z(3), -z(2)];
end

% Make v orthogonal to z
v = v - dot(v, z) * z / dot(z, z);
v = v / norm(v);

%find x : take cross product of z with a random (non-parallel) vector
x = cross(z, v);
%find y : cross z with x (guarantees right-handed coordinate system)
y = cross(z, x);

%final coordinate axes: each axis is a column of the resulting matrix
axes = [x' y' z'];

end
