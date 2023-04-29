function r_sunsat = satelliteCoords(z_input, r_sunearth)
 %satelliteCoords - Outputs a vector from the sun to the satellite based on
 %the satellite's coordinate system
 %      satelliteCoords(z_input, r_sunearth) uses the 3 x 1 (?*) normal 
 %      vector to the satellite's top face (z_input) and the 3 x 1 (?) 
 %      vector from the sun to the earth (r_sunearth) to calculate the 
 %      vector from the sun to the satellite in satellite coordinates. The
 %      satellite coordinate reference frame is determined by the normal
 %      vector to the satellite's top face (z_input), where z_input
 %      corresponds to the z-axis, and the x- and y- axes are chosen
 %      arbitrarily such that they form a right-handed coordinate system
 %      with z-input. (Note: they are chosen arbitrarily because the satellite 
 %      will be spinning, and thus the initial position of the x- and y- 
 %      will have minimal impact on the total power calculation, although 
 %      this function does not yet account for spin). 
 %     
 %      *see note at bottom
 %
 %      INPUTS:
 %      z_input    3 x 1 normal vector to satellite's top face (top face = 
 %                 face opposite to propulsion firing). ECI coordinates.
 %      r_sunearth Direction vector from the earth to the sun in ECI
 %                 coordinates.
 %                  
 %      
 %      OUTPUTS:
 %      r_sunsat   Direction vector from the sun to the satellite in
 %                 satellite coordinates. (note: satellite to sun? doesn't
 %                 matter much anyways)

 %      Author: Rahma Abdullah (ra567)
 %      Date: 4/29/23

%normalize z
z = z_input / norm(z_input);

% Choose a random vector v that is not parallel to z
if z(1) ~= 0 || z(2) ~= 0
    v = [z(2), -z(1), 0];
else
    v = [0, z(3), -z(2)];
end

% Make v orthogonal to z and normalize v
v = v - dot(v, z) * z / dot(z, z);
v = v / norm(v);

%find x : take cross product of z with a random (non-parallel) vector
x = cross(v,z);
%find y : cross z with x (guarantees right-handed coordinate system)
y = cross(z,x);

%final coordinate axes: each axis is a column of the resulting matrix
axes = [x' y' z'];
disp(axes);

%make r_sunearth a column vector if necessary
r_sunearth = reshape(r_sunearth, [], 1);

%solve for b-coordinates of r_sunearth in terms of axes
r_sunsat = axes \ r_sunearth;
r_sunsat = r_sunsat / norm(r_sunsat);

%note: idk how row and column vectors/matrices work in matlab TT this is
%trial and error but it's a row vector ?? what

end

