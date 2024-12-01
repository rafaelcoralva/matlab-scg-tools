function [max_diameter_over_volume] = anisotropy_calculator_avg_cycle_max_diameter_over_volume(scg_xyz)
% This function calculates a simple index of anisotropy of the input x/y/z acceleration trajectory (e.g., avg S1 complex):
% - The ratio of 3D acceleration trajetctory max diameter over its cubic root volume. 
% The greater this ratio, the greater the anisotropy.

% INPUTS:
% scg_xyz (3 column matrix) - The SCG acceleration trajectory for which we wish to compute the anisotropy (e.g., average S1 complex). 
%                             The three columns correspond to the 3 x, y, and z components.

% OUTPUT:
% max_diameter_over_volume (scalar) - Indicator of anisotropy.

% DESCRIPTION:
% Part 1. Calculation of the 3D hull surface and its volume using the MATLAB 'boundary' function.
% Part 2. Identification the maximum antipodal distance (i.e. diameter) of the surface.
% Part 3. Division by cubic root of volume.

% NOTE ON COORDINATE FRAMES:
% For the recorded SCG signals:
% - the SCG Cartesian coordinate system is: x = medial (left-) lateral axis.
%                                           y = long axis (cranial-positive).
%                                           z = posterior-anterior axis.
% - the SCG Spherical coordinate system is: azimuth (theta) = from the (SCG) x axis towards the (SCG) z axis.
%                                           elevation (phi) = from the (SCG) xz plane towards the (SCG) y axis.
%                                           radius = modulus of x,y, and z.
%
% NOTE: This SCG reference frame is NOT the same as that considered by various Cartesian and Spherical Matlab functions:
% - the Matlab Cartesian coordinate system is: x = medial-(left) lateral.
%                                              y = anterior-posterior.
%                                              z = long axis.
% - the Matlab Spherical coordinate system is: azimuth (theta) = from the (Matlab) x axis towards the (Matlab) y axis.
%                                              elevation (phi) = from the (Matlab) xy plane towards the (Matlab) z axis.
%                                              radius = modulus of x,y, and z.
%
% Hence, conversions are made in the code btween SCG and Matlab coordinate frames. These  conversions are:
% x_scg = x_mat
% y_scg = z_mat
% z_scg = -y_mat
% azimuth_scg = -1*azimuth_mat
% elevation_scg = elevation_mat
% radius_scg = radius_mat

diameter_sph_mat = zeros(181,180); % Pre-Initialization for speed.


%% 1. Calculation of the 3D hull surface using the MATLAB boundary function.

% In this entire part, the xyz coordinates are assumed to be in the SCG reference frame.
[triangles_idx_xzy, volume] = boundary(scg_xyz(:,1), scg_xyz(:,3), scg_xyz(:,2), 0); % A matrix of triangles that comprises the 3D hull surface. Each row represents a seperate triangle. Each element of a particular row represents the row index of one of the vertices of the triange. The '0' in the 'boundary' function call specifies that the boundary is a convex hull. If it is 1, it becomes the concave hull. Volume represents the volume enclosed by the surface.

tri_vert_1_xyz = zeros(size(triangles_idx_xzy,1),3); % Pre-Initialization for speed.
tri_vert_2_xyz = zeros(size(triangles_idx_xzy,1),3);
tri_vert_3_xyz = zeros(size(triangles_idx_xzy,1),3);

for i = 1:size(triangles_idx_xzy,1)
    tri_vert_1_xyz(i,:) = scg_xyz(triangles_idx_xzy(i,1),:); % The x,y,and z coordinates (in the SCG reference frame) of vertex 1 of the triangles comprising the 3D hull surface.
    tri_vert_2_xyz(i,:) = scg_xyz(triangles_idx_xzy(i,2),:); % The x,y,and z coordinates (in the SCG reference frame) of vertex 2 of the triangles comprising the 3D hull surface.
    tri_vert_3_xyz(i,:) = scg_xyz(triangles_idx_xzy(i,3),:); % The x,y,and z coordinates (in the SCG reference frame) of vertex 3 of the triangles comprising the 3D hull surface.
end


%% 2. Identification of the azimuthal and elevation angles that produce the maximum diameter.
% In this entire part, a distinction is made between the xyz coordinates and the polar and azimuthal angle between the SCG and Matlab reference frames.

% Calculating the "radius" of the 3D hull surface for all elevation angles (-90 to 90 degrees) and azimuthal angles (-180 to 180 degrees).
for el_mat = -90:1:90 % el_mat represents the elevation in the Matlab coordinate frame.
    for az_mat = -179:1:180 % az_mat represents the azimuthal angle in the Matlab coordinate frame.
        
        [xhat_mat, yhat_mat, zhat_mat] = sph2cart(az_mat*pi/180,el_mat*pi/180,1); % The unit-vectors (according to the Matlab reference frame) pointing the direction of the current elevation and azimuth angle combination.
        
        xhat_scg = xhat_mat;
        yhat_scg = zhat_mat;
        zhat_scg = -1*yhat_mat;
        
        % line([0 100*xhat_scg],[0 100*zhat_scg],[0 100*yhat_scg]) % Troubleshooting - Plotting a line primitive in the current axis.
        
        [idx_intersection, intersection_radius] = TriangleRayIntersection([0, 0, 0], [xhat_scg, yhat_scg, zhat_scg], tri_vert_1_xyz, tri_vert_2_xyz, tri_vert_3_xyz);
        
        % Calculating the "radius" of the 3D Hull Surface and storing it in radius_sph_mat - a matrix whose row and column indices represent elevation and azimuth, respectively (in the Matlab spherical reference frame), and whose elements represent the radii.
        largest_radius = max(intersection_radius(idx_intersection)); % There should only be 1 radius if the boundary is convex but for some reason there sometimes is two.
        radius_sph_mat(el_mat+91, az_mat+180) = largest_radius;
        
        clear largest_radius
        
        % Cartesian (x,y,z) coordinates of the (index-parametrised) radius_sph_mat, in the Matlab reference frame.
        [radius_x_mat(el_mat+91, az_mat+180), radius_y_mat(el_mat+91, az_mat+180), radius_z_mat(el_mat+91, az_mat+180)] = sph2cart(az_mat*pi/180, el_mat*pi/180, radius_sph_mat(el_mat+91, az_mat+180));
        
        clear xhat_mat yhat_mat zhat_mat xhat_scg yhat_scg zhat_scg idx_intersection intersection_radius
        
    end
end

clear az_mat el_mat radius_x_mat radius_y_mat radius_z_mat

% Calculating the 3D Hull Surface "diameter" (i.e. the antipodal distance) for each combination of elevation and azimuthal angle.
for el_mat = -90:1:90 % Considering the anterior hemisphere (only a single hemisphere needs to be considered as the "diameter" measured from the opposing hemisphere is the same).
    for az_mat = -179:1:0
        
        [el_anitpode_mat, az_antipode_mat] = antipodal_calculator(el_mat, az_mat); % Calculating the antipodal elevation and azimuth for the current elevation and azimuth angle combination.
        
        % Identifying the indices in radius_sph_mat for the radii corresponding to the current elevation and azimuth angle combination and its antipode.
        idx_el_mat = el_mat + 91;
        idx_az_mat = az_mat + 180;
        
        idx_el_antipode_mat = el_anitpode_mat + 91;
        idx_az_antipode_mat = az_antipode_mat + 180;
        
        % Calculating the "diameter" = sum of antipodal radii, for that particular elevation and azimuthal angle combination.
        diameter_sph_mat(el_mat + 91, az_mat + 180) = radius_sph_mat(idx_el_mat,idx_az_mat) + radius_sph_mat(idx_el_antipode_mat,idx_az_antipode_mat);
        
    end
end

clear el_mat idx_el_mat az_mat idx_az_mat el_anitpode_mat az_antipode_mat idx_el_antipode_mat idx_az_antipode_mat

% Finding maximum diameter (i.e. maximum radius antipodal sum) and then calculating the ratio to the enclosed surface volume.
max_diameter = max(diameter_sph_mat,[],'all');


%% 3. Dividing by cubic root of volume.
max_diameter_over_volume = max_diameter/(volume^1/3);


end


%% Sub-functions.
function [el_antipode, az_antipode] = antipodal_calculator(el, az)

el_antipode = -1*el;

if az < 0
    az_antipode = az + 180;
elseif az > 0
    az_antipode = az - 180;
elseif az == 0
    az_antipode = 180;
end
end