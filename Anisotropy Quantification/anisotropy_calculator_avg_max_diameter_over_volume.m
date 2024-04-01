function [avg_diameter_over_volume] = anisotropy_calculator_avg_max_diameter_over_volume(fs,r_pos,scg_xyz)
% This algorithm calculates a simple index of average S1 anisotropy - the ratio of S1 3D trajetctoy max diameter over its volume. 
% The greater this ratio, the greater the anisotropy.

% Part 1. Calculation of the 3D hull surface and its volume using the MATLAB 'boundary' function.
% Part 2. Identification the maximum antipodal distance (i.e. diameter) of the surface.

% INPUTS:
% fs - sampling frequency
% r_pos - ecg R wave sample positions
% scg_xyz - scg signals - the three columns correspond to the 3 x, y, and z components.

% OUTPUT:
% avg_diameter_over_volume = indicator of anisotropy.

% NOTE ON COORDINATE FRAMES:
% For the recorded SCG signals:
% - the SCG Cartesian coordinate system is: x = medial-(left) lateral axis.
%                                           y = long axis.
%                                           z = posterior-anterior axis.
% - the SCG Spherical coordinate system is: azimuth (theta) = from the (SCG) x axis towards the (SCG) z axis.
%                                           elevation (phi) = from the (SCG) xz plane towards the (SCG) y axis.
%                                           radius = modulus of x,y, and z.
%
% This SCG reference frame is NOT the same as that considered by various Cartesian and Spherical Matlab functions:
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

s1_xyz = zeros(1,0.2*fs); % Pre-Initialization for speed.
diameter_over_volume = zeros(1,length(r_pos));
diameter_sph_mat = zeros(181,180);

for i = 1:length(r_pos)
    
    s1_xyz = scg_xyz(r_pos(i):r_pos(i) + 0.2*fs,:); % The s1 complex is defined as the 200ms following each R-peak detection.
    
    %% Part 1.0 Calculation of the 3D hull surface using the MATLAB boundary function.
    
    % In this entire part, the xyz coordinates are assumed to be in the SCG reference frame.
    
    [triangles_idx_xzy, volume] = boundary(s1_xyz(:,1), s1_xyz(:,3), s1_xyz(:,2), 0); % A matrix of triangles that comprises the 3D hull surface. Each row represents a seperate triangle. Each element of a particular row represents the row index of one of the vertices of the triange. The '0' in the 'boundary' function call specifies that the boundary is a convex hull. If it is 1, it becomes the concave hull. Volume represents the volume enclosed by the surface.

    tri_vert_1_xyz = zeros(size(triangles_idx_xzy,1),3); % Pre-Initialization for speed.
    tri_vert_2_xyz = zeros(size(triangles_idx_xzy,1),3);
    tri_vert_3_xyz = zeros(size(triangles_idx_xzy,1),3);
    
    for j = 1:size(triangles_idx_xzy,1)
        tri_vert_1_xyz(j,:) = s1_xyz(triangles_idx_xzy(j,1),:); % The x,y,and z coordinates (in the SCG reference frame) of vertex 1 of the triangles comprising the 3D hull surface.
        tri_vert_2_xyz(j,:) = s1_xyz(triangles_idx_xzy(j,2),:); % The x,y,and z coordinates (in the SCG reference frame) of vertex 2 of the triangles comprising the 3D hull surface.
        tri_vert_3_xyz(j,:) = s1_xyz(triangles_idx_xzy(j,3),:); % The x,y,and z coordinates (in the SCG reference frame) of vertex 3 of the triangles comprising the 3D hull surface.
    end
    
    %% 2. Identification of the azimuthal and elevation angles that produce the maximum antipode.
    
    % In this entire part, a distinction is made between the xyz coordinates and the polar and azimuthal angle between the SCG and Matlab reference frames.
    
    %% Part 2.1 - Calculating the "radius" of the 3D hull surface for all elevation angles (-90 to 90 degrees) and azimuthal angles (-180 to 180 degrees).
    for el_mat = -90:1:90 % el_mat represents the elevation in the Matlab coordinate frame.
        for az_mat = -179:1:180 % az_mat represents the azimuthal angle in the Matlab coordinate frame.
            
            [xhat_mat, yhat_mat, zhat_mat] = sph2cart(az_mat*pi/180,el_mat*pi/180,1); % The unit-vectors (according to the Matlab reference frame) pointing the direction of the current elevation and azimuth angle combination.
            
            xhat_scg = xhat_mat;
            yhat_scg = zhat_mat;
            zhat_scg = -1*yhat_mat;
            
            % line([0 100*xhat_scg],[0 100*zhat_scg],[0 100*yhat_scg]) % Troubleshooting - Plotting a line primitive in the current axis.
            
            [intersection_idx, intersection_radius] = TriangleRayIntersection([0 0 0], [xhat_scg, yhat_scg, zhat_scg], tri_vert_1_xyz, tri_vert_2_xyz, tri_vert_3_xyz);
            
            % Calculating the "radius" of the 3D Hull Surface and storing it in radius_sph_mat - a matrix whose row and column indices represent elevation and azimuth, respectively (in the Matlab spherical reference frame), and whose elements represent the radii.
            largest_radius = max(intersection_radius(find(intersection_idx))); % There should only be 1 radius if the boundary is convex but for some reason there sometimes is two.
            radius_sph_mat(el_mat+91, az_mat+180) = largest_radius;
            
            clear largest_radius
            
            % Cartesian (x,y,z) coordinates of the (index-parametrised) radius_sph_mat, in the Matlab reference frame.
            [radius_x_mat(el_mat+91,az_mat+180),radius_y_mat(el_mat+91,az_mat+180),radius_z_mat(el_mat+91,az_mat+180)] = sph2cart(az_mat*pi/180, el_mat*pi/180, radius_sph_mat(el_mat+91, az_mat+180));
            
            clear xhat_mat yhat_mat zhat_mat xhat_scg yhat_scg zhat_scg intersection_idx intersection_radius
            
        end
    end
    
    clear az_mat el_mat s1_p_x_mat s1_p_y_mat s1_p_z_mat radius_x_mat radius_y_mat radius_z_mat
    
    %% Part 2.2 - Calculating the 3D Hull Surface "diameter" (i.e. the antipodal distance) for each combination of elevation and azimuthal angle.
    for el_mat = -90:1:90 % Considering the anterior hemisphere (only a single hemisphere needs to be considered as the "diameter" measured from the opposing hemisphere is the same).
        for az_mat = -179:1:0
            
            [antipodal_el_mat, antipodal_az_mat] = antipodal_calculator(el_mat, az_mat); % Calculating the antipodal elevation and azimuth for the current elevation and azimuth angle combination.
            
            % Identifying the indices in radius_sph_mat for the radii corresponding to the current elevation and azimuth angle combination and its antipode.
            el_mat_idx = el_mat + 91;
            az_mat_idx = az_mat + 180;
            
            antipodal_el_mat_idx = antipodal_el_mat + 91;
            antipodal_az_mat_idx = antipodal_az_mat + 180;
            
            % Calculating the "diameter" = sum of antipodal radii, for that particular elevation and azimuthal angle combination.
            diameter_sph_mat(el_mat + 91, az_mat + 180) = radius_sph_mat(el_mat_idx,az_mat_idx) + radius_sph_mat(antipodal_el_mat_idx,antipodal_az_mat_idx);
            
        end
    end
    
    clear el_mat el_mat_idx az_mat az_mat_idx antipodal_el_mat antipodal_az_mat antipodal_el_mat_idx antipodal_az_mat_idx
    
    %% Part 2.3 Finding maximum diameter.
    
    % Finding maximum diameter (i.e. maximum radius antipodal sum) and then calculating the ratio to the enclosed surface volume.
    max_diameter = max(max(diameter_sph_mat));
    
    diameter_over_volume(i) = max_diameter/volume;
       
end

avg_diameter_over_volume = mean(diameter_over_volume);

%% Sub-functions.
    function [antipodal_el, antipodal_az] = antipodal_calculator(el, az)
        
        antipodal_el = -1*el;
        
        if az < 0
            antipodal_az = az + 180;
        elseif az > 0
            antipodal_az = az - 180;
        elseif az == 0
            antipodal_az = 180;
        end
    end

end