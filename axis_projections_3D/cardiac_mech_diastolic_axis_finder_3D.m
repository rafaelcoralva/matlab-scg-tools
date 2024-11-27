function [axes_angles] = cardiac_mech_diastolic_axis_finder_3D(avg_s2_xyz)
% This algorithm identifies the cardiac mechanical diastolic (i.e. S2) axis in three dimensions using the 'boundary' function-derived hull surface.
% It consists of two parts:

% Part 1. Calculation of the 3D hull surface using the MATLAB 'boundary' function.
% Part 2. Identification of the azimuthal and elevation angles that produce the maximum antipodal distance of the 3D hull.

% INPUTS:
% avg_s2_xyz (Matrix: [length(avg_s2), 3]) - 3D average S2 acceleration trajectory matrix of size,
%                                            where the three columns correspond to the 3 scg xyz components

% OUTPUTS:
% axis_angles = [el_mech_axis, az_mech_axis, el_antimech_axis, az_antimech_axis], where:
%               el - refers to elevation angle from the SCG xz plane [i.e. anatomical transverse plane] TOWARDS the SCG y axis [i.e. anatomical long axis].
%               az - refers to azimuth angle from SCG x axis [i.e. anatomical left medial lateral axis] TOWARDS the SCG z axis [i.e. posterior-anterior axis].
%               mech_axis - mechanical axis maximising the 3D hull surface diameter.
%               antimech_axis - antimechanical axis minimising the 3D hull surface diameter - it is usually perpendicular to mech_axis

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

flag_plot = 1; % Plot flag. If 0, no plots are produced. If 1, only final plot is produced. If 2, all plots are produced.


%% 1.1 Calculation of the 3D hull surface using the MATLAB boundary function.
% In this entire part, the xyz coordinates are assumed to be in the SCG reference frame.
idx_triangles_xzy = boundary(avg_s2_xyz(:,1), avg_s2_xyz(:,3), avg_s2_xyz(:,2), 0); % A matrix of triangles that comprises the 3D hull surface. 
                                                                                    % Each row represents a seperate triangle. 
                                                                                    % Each element of a particular row represents the row index of one of the vertices of the triange. 
                                                                                    % The '0' in the 'boundary' function call is the shrink factor specifying that the boundary is a convex hull. If it is 1, it becomes the concave hull.


%% 1.2 Tri-force Interpolation of the 3D hull surface (interpolating each triangle into 4 smaller triangles)
% NOTE - I suspect this part isn't actually necessary because TriangleRayIntersection function used later doesn't care about the size of the triangles -> Check.
% In this entire part, the xyz coordinates are assumed to be in the SCG reference frame.
unos = [1, 1, 1]; % Required later computation of area of a triangle using coordinate geometry equation.

area_thr = 2.5; % The area limit (in mg^2) of the triangles that comprise the 3D hull surface. If the area of one of these triangles is greater than this threshold, it must be broken down (i.e. interpolated) into smaller triangles.
length_thr = 5; % The length limit (in mg) of the sides of the triangles that comprise the 3D hull surface. If the area of one of these triangles is greater than this threshold, it must be broken down (i.e. interpolated) into smaller triangles.

avg_s2_xyz_interp = avg_s2_xyz; % A copy of the average S2 component that will be extended with extra points should interpolation be required.
idx_triangles_xzy_interp = idx_triangles_xzy; % A copy of the matrix of triangles comprising the 3D hull surface that will be extended with extra points should interpolation be required.

counter_interps = 0; % Counter of times the interpolation is performed.

% Initialization of the three vertices of the four new triangles created upon each Tri-Force interpolation.
small_tri_1_vert_1_xyz = [nan, nan, nan]; small_tri_1_vert_2_xyz = [nan, nan, nan]; small_tri_1_vert_3_xyz = [nan, nan, nan];
small_tri_2_vert_1_xyz = [nan, nan, nan]; small_tri_2_vert_2_xyz = [nan, nan, nan]; small_tri_2_vert_3_xyz = [nan, nan, nan];
small_tri_3_vert_1_xyz = [nan, nan, nan]; small_tri_3_vert_2_xyz = [nan, nan, nan]; small_tri_3_vert_3_xyz = [nan, nan, nan];
small_tri_4_vert_1_xyz = [nan, nan, nan]; small_tri_4_vert_2_xyz = [nan, nan, nan]; small_tri_4_vert_3_xyz = [nan, nan, nan];

i = 1;
while (i <= size(idx_triangles_xzy_interp,1)) % Iterating through each large triangle (note, interpolating a large triangle may result in the creating of 1-4 smaller but still "large" triangles).
       
    large_tri_vert_1_xyz(i,:) = avg_s2_xyz_interp(idx_triangles_xzy_interp(i,1),:); % The x,y,and z coordinates (in the SCG reference frame) of vertex 1 of the triangles comprising the 3D hull surface.
    large_tri_vert_2_xyz(i,:) = avg_s2_xyz_interp(idx_triangles_xzy_interp(i,2),:); % The x,y,and z coordinates (in the SCG reference frame) of vertex 2 of the triangles comprising the 3D hull surface.
    large_tri_vert_3_xyz(i,:) = avg_s2_xyz_interp(idx_triangles_xzy_interp(i,3),:); % The x,y,and z coordinates (in the SCG reference frame) of vertex 3 of the triangles comprising the 3D hull surface.
    
    large_tri_vert_x_123(i,:) = [large_tri_vert_1_xyz(i,1) large_tri_vert_2_xyz(i,1) large_tri_vert_3_xyz(i,1)]; % The x coordinates (in the SCG reference frame) of vertex 1,2, and 3 of the triangles comprising the 3D hull surface.
    large_tri_vert_y_123(i,:) = [large_tri_vert_1_xyz(i,2) large_tri_vert_2_xyz(i,2) large_tri_vert_3_xyz(i,2)]; % The y coordinates (in the SCG reference frame) of vertex 1,2, and 3 of the triangles comprising the 3D hull surface.
    large_tri_vert_z_123(i,:) = [large_tri_vert_1_xyz(i,3) large_tri_vert_2_xyz(i,3) large_tri_vert_3_xyz(i,3)]; % The z coordinates (in the SCG reference frame) of vertex 1,2, and 3 of the triangles comprising the 3D hull surface.
    
    large_tri_max_length(i) = max([norm(large_tri_vert_2_xyz(i,:) - large_tri_vert_1_xyz(i,:)), norm(large_tri_vert_3_xyz(i,:) - large_tri_vert_2_xyz(i,:)), norm(large_tri_vert_3_xyz(i,:) - large_tri_vert_1_xyz(i,:))]); % The length of the longest side of that triangle.
    large_tri_area(i) = 0.5*sqrt(det([large_tri_vert_x_123(i,:);large_tri_vert_y_123(i,:);unos])^2 + det([large_tri_vert_y_123(i,:);large_tri_vert_z_123(i,:);unos])^2 + det([large_tri_vert_z_123(i,:);large_tri_vert_x_123(i,:);unos])^2); % The area of that triangle.
    
    if (large_tri_area(i) > area_thr) || (large_tri_max_length(i) > length_thr) % Check if large triangle requires interpolation (i.e., if it's actually "large).
        % Tri-force interpolation - creating 4 smaller triangles from the original triangle.
        counter_interps = counter_interps + 1; % Increasing the number of interpolations counter.
        
        % Finding the midpoints of the sides of the large triangle.
        large_tri_midpoint_1_2_xyz = (large_tri_vert_1_xyz(i,:) + large_tri_vert_2_xyz(i,:))/2;
        large_tri_midpoint_2_3_xyz = (large_tri_vert_2_xyz(i,:) + large_tri_vert_3_xyz(i,:))/2;
        large_tri_midpoint_3_1_xyz = (large_tri_vert_3_xyz(i,:) + large_tri_vert_1_xyz(i,:))/2;
        
        % Defining 4 new small triangles.        
        % Small Triangle 1.
        small_tri_1_vert_1_xyz(end+1,:) =  large_tri_vert_1_xyz(i,:);
        small_tri_1_vert_2_xyz(end+1,:) =  large_tri_midpoint_1_2_xyz;
        small_tri_1_vert_3_xyz(end+1,:) =  large_tri_midpoint_3_1_xyz;
        
        % Small Triangle 2.
        small_tri_2_vert_1_xyz(end+1,:) =  large_tri_midpoint_1_2_xyz;
        small_tri_2_vert_2_xyz(end+1,:) =  large_tri_vert_2_xyz(i,:);
        small_tri_2_vert_3_xyz(end+1,:) =  large_tri_midpoint_2_3_xyz;
        
        % Small Triangle 3.
        small_tri_3_vert_1_xyz(end+1,:) =  large_tri_midpoint_3_1_xyz;
        small_tri_3_vert_2_xyz(end+1,:) =  large_tri_midpoint_2_3_xyz;
        small_tri_3_vert_3_xyz(end+1,:) =  large_tri_vert_3_xyz(i,:);
        
        % Small Triangle 4 - the inverted middle triangle of the tri-force.
        small_tri_4_vert_1_xyz(end+1,:) =  large_tri_midpoint_1_2_xyz;
        small_tri_4_vert_2_xyz(end+1,:) =  large_tri_midpoint_2_3_xyz;
        small_tri_4_vert_3_xyz(end+1,:) =  large_tri_midpoint_3_1_xyz;
        
        % Adding the vertices of triangle 4 (the only "new" data points) to the interpolated average S2 xyz component.
        % This is necessary because the trisurf function will plot based on the indices of the average S2 xyz component (hence the new data points must first exist).
        % Although intiuitively it might not make sense to add new interpolated datapoints to the end of the avg s2 
        % (i.e., chronologically they may not make sense), it doesn't matter because trisurf just sees a cloud of points and doesn't care about their order.
        avg_s2_xyz_interp(end+1,:) = [small_tri_4_vert_1_xyz(end,1) small_tri_4_vert_1_xyz(end,2) small_tri_4_vert_1_xyz(end,3)]; % Smal Triangle 4 Vertex 1.
        avg_s2_xyz_interp(end+1,:) = [small_tri_4_vert_2_xyz(end,1) small_tri_4_vert_2_xyz(end,2) small_tri_4_vert_2_xyz(end,3)]; % Smal Triangle 4 Vertex 2.
        avg_s2_xyz_interp(end+1,:) = [small_tri_4_vert_3_xyz(end,1) small_tri_4_vert_3_xyz(end,2) small_tri_4_vert_3_xyz(end,3)]; % Smal Triangle 4 Vertex 3.
        
        % Deleting the triangle indices of the original large triangle.
        idx_triangles_xzy_interp(i,:) = [];
        
        % Adding the triangle indices of the new four small triangles.        
        % Small Triangle 1.
        [~,idx_small_tri_1_vert_1] = ismember([small_tri_1_vert_1_xyz(end,1), small_tri_1_vert_1_xyz(end,2), small_tri_1_vert_1_xyz(end,3)], avg_s2_xyz_interp,'rows'); % Index of small triangle 1's vertex 1.
        [~,idx_small_tri_1_vert_2] = ismember([small_tri_1_vert_2_xyz(end,1), small_tri_1_vert_2_xyz(end,2), small_tri_1_vert_2_xyz(end,3)], avg_s2_xyz_interp,'rows'); % Index of small triangle 1's vertex 2.
        [~,idx_small_tri_1_vert_3] = ismember([small_tri_1_vert_3_xyz(end,1), small_tri_1_vert_3_xyz(end,2), small_tri_1_vert_3_xyz(end,3)], avg_s2_xyz_interp,'rows'); % Index of small triangle 1's vertex 3.
       
        idx_small_tri_1_allvertices = [idx_small_tri_1_vert_1, idx_small_tri_1_vert_2, idx_small_tri_1_vert_3];
        
        % Small Triangle 2.
        [~,idx_small_tri_2_vert_1] = ismember([small_tri_2_vert_1_xyz(end,1), small_tri_2_vert_1_xyz(end,2), small_tri_2_vert_1_xyz(end,3)], avg_s2_xyz_interp,'rows'); % Index of small triangle 2's vertex 1.
        [~,idx_small_tri_2_vert_2] = ismember([small_tri_2_vert_2_xyz(end,1), small_tri_2_vert_2_xyz(end,2), small_tri_2_vert_2_xyz(end,3)], avg_s2_xyz_interp,'rows'); % Index of small triangle 2's vertex 2.
        [~,idx_small_tri_2_vert_3] = ismember([small_tri_2_vert_3_xyz(end,1), small_tri_2_vert_3_xyz(end,2), small_tri_2_vert_3_xyz(end,3)], avg_s2_xyz_interp,'rows'); % Index of small triangle 2's vertex 3.
   
        idx_small_tri_2_allvertices = [idx_small_tri_2_vert_1, idx_small_tri_2_vert_2, idx_small_tri_2_vert_3];
        
        % Small Triangle 3.
        [~,idx_small_tri_3_vert_1] = ismember([small_tri_3_vert_1_xyz(end,1) small_tri_3_vert_1_xyz(end,2) small_tri_3_vert_1_xyz(end,3)], avg_s2_xyz_interp,'rows'); % Index of small triangle 3's vertex 1.
        [~,idx_small_tri_3_vert_2] = ismember([small_tri_3_vert_2_xyz(end,1) small_tri_3_vert_2_xyz(end,2) small_tri_3_vert_2_xyz(end,3)], avg_s2_xyz_interp,'rows'); % Index of small triangle 3's vertex 2.
        [~,idx_small_tri_3_vert_3] = ismember([small_tri_3_vert_3_xyz(end,1) small_tri_3_vert_3_xyz(end,2) small_tri_3_vert_3_xyz(end,3)], avg_s2_xyz_interp,'rows'); % Index of small triangle 3's vertex 3.
      
        idx_small_tri_3_allvertices = [idx_small_tri_3_vert_1, idx_small_tri_3_vert_2, idx_small_tri_3_vert_3];
        
        % Small Triangle 4.
        [~,idx_small_tri_4_vert_1] = ismember([small_tri_4_vert_1_xyz(end,1) small_tri_4_vert_1_xyz(end,2) small_tri_4_vert_1_xyz(end,3)],avg_s2_xyz_interp,'rows'); % Index of small triangle 4's vertex 1.
        [~,idx_small_tri_4_vert_2] = ismember([small_tri_4_vert_2_xyz(end,1) small_tri_4_vert_2_xyz(end,2) small_tri_4_vert_2_xyz(end,3)],avg_s2_xyz_interp,'rows'); % Index of small triangle 4's vertex 2.
        [~,idx_small_tri_4_vert_3] = ismember([small_tri_4_vert_3_xyz(end,1) small_tri_4_vert_3_xyz(end,2) small_tri_4_vert_3_xyz(end,3)],avg_s2_xyz_interp,'rows'); % Index of small triangle 4's vertex 3.
       
        small_tri_4_allvertices_idx = [idx_small_tri_4_vert_1, idx_small_tri_4_vert_2, idx_small_tri_4_vert_3];
       
        idx_triangles_xzy_interp = insertrows(idx_triangles_xzy_interp,[idx_small_tri_1_allvertices; idx_small_tri_2_allvertices; idx_small_tri_3_allvertices; small_tri_4_allvertices_idx], i-1);
        i = i - 1; % This way, the four newly-created triangles will also be analysed for possible interpolation.
        
        % The interpolated large triangle has now become four small triangles. 
        % These four new small triangles are added to the surface.
        % On the next iteration the four new small triangles will be considered as "large" triangles which will themselves be analysed for further interpolation if necessary.
        
        clear large_tri_midpoint_1_2_xyz large_tri_midpoint_2_3_xyz large_tri_midpoint_3_1_xyz idx_small_tri_1_vert_1 idx_small_tri_1_vert_2 idx_small_tri_1_vert_3 idx_small_tri_1_allvertices idx_small_tri_2_vert_1 idx_small_tri_2_vert_2 idx_small_tri_2_vert_3 idx_small_tri_2_allvertices idx_small_tri_3_vert_1 idx_small_tri_3_vert_2 idx_small_tri_3_vert_3 idx_small_tri_3_allvertices idx_small_tri_4_vert_1 idx_small_tri_4_vert_2 idx_small_tri_4_vert_3 small_tri_4_allvertices_idx 
    end
        
    i = i + 1;
end

clear unos i large_tri_area large_tri_max_length small_tri_1_vert_1_xyz small_tri_1_vert_2_xyz small_tri_1_vert_3_xyz small_tri_1_vert_4_xyz small_tri_2_vert_1_xyz small_tri_2_vert_2_xyz small_tri_2_vert_3_xyz small_tri_3_vert_1_xyz small_tri_3_vert_2_xyz small_tri_3_vert_3_xyz small_tri_4_vert_1_xyz small_tri_4_vert_2_xyz small_tri_4_vert_3_xyz


%% 1.3 Plotting.
if flag_plot
    HullSurface = figure;
    title('S2 3D Hull Surface'); xlabel('x (mg)'); ylabel('z (mg)'); zlabel('y (mg)');

    % Plotting the original S2 acceleration trajectory.
    plot3(avg_s2_xyz(:,1), avg_s2_xyz(:,3), avg_s2_xyz(:,2),'k','LineWidth',2); hold on;

    % Plotting the vertices of the triangles defining the 3D Hull Surface.
    % hold on; plot3(large_tri_vert_1_xyz(:,1),large_tri_vert_1_xyz(:,3),large_tri_vert_1_xyz(:,2),'r.','MarkerSize',10); plot3(large_tri_vert_2_xyz(:,1),large_tri_vert_2_xyz(:,3),large_tri_vert_2_xyz(:,2),'r.','MarkerSize',10); plot3(large_tri_vert_3_xyz(:,1),large_tri_vert_3_xyz(:,3),large_tri_vert_3_xyz(:,2),'r.','MarkerSize',10);

    % Plotting the (possibly interpolated) 3D Hull Surface.
    trisurf(idx_triangles_xzy_interp, avg_s2_xyz_interp(:,1), avg_s2_xyz_interp(:,3), avg_s2_xyz_interp(:,2), sqrt(avg_s2_xyz_interp(:,1).^2 + avg_s2_xyz_interp(:,3).^2 + avg_s2_xyz_interp(:,2).^2) ,'FaceAlpha',0.5,'EdgeAlpha',0); colormap spring; xlabel('x (mg)'); ylabel('z (mg)'); zlabel('y (mg)'); hold on;

    min_lim = min(avg_s2_xyz_interp);
    max_lim = max(avg_s2_xyz_interp); 
    xlim([min_lim(1), max_lim(1)]);
    ylim([min_lim(1), max_lim(1)]); 
    zlim([min_lim(1), max_lim(1)]); % Setting all axes to equal span.

    legend('S2 Acceleration Trajectory','S2 Hull');
    hcb=colorbar;
    title(hcb,' Acceleration (mg)');
    title('S2')

    clear min_lim max_lim hcb
end


%% 2. Identification of the azimuthal and elevation angles that produce the maximum antipode.
% In this entire part, a distinction is made between the xyz coordinates and the polar and azimuthal angle between the SCG and Matlab reference frames.


%% 2.1 - Calculating the "radius" of the 3D hull surface for all elevation angles (-90 to 90 degrees) and azimuthal angles (-180 to 180 degrees).
for el_mat = -90:1:90 % el_mat represents the elevation in the Matlab coordinate frame.
    for az_mat = -179:1:180 % az_mat represents the azimuthal angle in the Matlab coordinate frame.
       
        [xhat_mat, yhat_mat, zhat_mat] = sph2cart(az_mat*pi/180,el_mat*pi/180,1); % The unit-vectors (according to the Matlab reference frame) pointing the direction of the current elevation and azimuth angle combination.

        xhat_scg = xhat_mat;
        yhat_scg = zhat_mat;
        zhat_scg = -1*yhat_mat;
        
        % line([0 100*xhat_scg],[0 100*zhat_scg],[0 100*yhat_scg]) % Troubleshooting - Plotting a line primitive in the current axis.
        
        [idx_intersection, intersection_radius] = TriangleRayIntersection([0, 0, 0], [xhat_scg, yhat_scg, zhat_scg], large_tri_vert_1_xyz, large_tri_vert_2_xyz, large_tri_vert_3_xyz); % Identifying the triangle with which the ray pointing in the current elevation+azimuth intersects with.

        % Calculating the "radius" of the 3D Hull Surface and storing it in radius_sph_mat - a matrix whose row and column indices represent elevation and azimuth, respectively (in the Matlab spherical reference frame), and whose elements represent the radii.
        largest_radius = max(intersection_radius(find(idx_intersection))); % There should only be 1 radius if the boundary is convex but for some reason there sometimes is two. 
        radius_sph_mat(el_mat+91, az_mat+180) = largest_radius; % The +91 and +180 are for correct Matlab indexing.
            
        clear largest_radius       
        
        % Cartesian (x,y,z) coordinates of the (index-parametrised) radius_sph_mat, in the Matlab reference frame. 
        [radius_x_mat(el_mat+91,az_mat+180), ...
         radius_y_mat(el_mat+91,az_mat+180), ...
         radius_z_mat(el_mat+91,az_mat+180)] = sph2cart(az_mat*pi/180, el_mat*pi/180, radius_sph_mat(el_mat+91, az_mat+180));     
        
        clear xhat_mat yhat_mat zhat_mat xhat_scg yhat_scg zhat_scg idx_intersection intersection_radius
        
    end
end

clear az_mat el_mat

% Converting from Cartesian Matlab reference frame to Cartesian SCG reference frame - this is needed later for plotting.
radius_x_scg = radius_x_mat;
radius_y_scg = radius_z_mat;
radius_z_scg = -1*radius_y_mat;

clear s2_p_x_mat s2_p_y_mat s2_p_z_mat radius_x_mat radius_y_mat radius_z_mat


%% 2.2 - Calculating the 3D Hull Surface "diameter" (i.e. the antipodal distance) for each combination of elevation and azimuthal angle.
diameter_sph_mat = zeros([181, 180]); % Initialization for efficiency
for el_mat = -90:1:90 % Considering the anterior hemisphere (only a single hemisphere needs to be considered as the "diameter" measured from the opposing hemisphere is the same).
    for az_mat = -179:1:0
        
        % Calculating the antipodal elevation and azimuth for the current elevation and azimuth angle combination.
        [el_antipode_mat, az_antipode_mat] = antipodal_calculator(el_mat, az_mat); % Calculating the antipodal elevation and azimuth for the current elevation and azimuth angle combination.

        % Identifying the indices in radius_sph_mat for the radii corresponding to the current elevation and azimuth angle combination and its antipode.
        idx_el_mat = el_mat + 91;
        idx_az_mat = az_mat + 180;
        
        idx_el_antipode_mat = el_antipode_mat + 91;
        idx_az_antipode_mat = az_antipode_mat + 180;
        
        % Calculating the "diameter" = sum of antipodal radii, for that particular elevation and azimuthal angle combination.
        diameter_sph_mat(el_mat + 91, az_mat + 180) = radius_sph_mat(idx_el_mat, idx_az_mat) + radius_sph_mat(idx_el_antipode_mat,idx_az_antipode_mat);
  
    end
end

clear el_mat idx_el_mat az_mat idx_az_mat el_antipode_mat az_antipode_mat idx_el_antipode_mat idx_az_antipode_mat


%% 2.3 Finding the angle combination maximising the 3D Hull Surface "diameter".

% Finding maximum diameter (i.e. maximum radius antipodal sum) and the elevation and azimuthal angles which produce it.
max_diameter = max(diameter_sph_mat, [], 'all');

[idx_max_diameter_el_sph, idx_max_diameter_az_sph] = find(diameter_sph_mat == max_diameter); % The row (elevation) and column (azimuth) indices of 'diameter_sph_mat' corresponding to the max diameter (in the anterior hemisphere). 
el_max_diameter_sph_mat = idx_max_diameter_el_sph - 91; % The actual elevation angle (in the Matlab reference frame) maximising diameter. 
az_max_diameter_sph_mat = idx_max_diameter_az_sph - 180; % The actual azimuthal angle (in the Matlab reference frame) maximising diameter.

% Finding the elevation and azimuth angles of the antipode (in the posterior hemisphere) of the point maximising diameter (referenced in the anterior hemisphere).
[el_max_diameter_antipode_sph_mat, az_max_diameter_antipode_sph_mat] = antipodal_calculator(el_max_diameter_sph_mat, az_max_diameter_sph_mat);


%% 2.4 Finding the angle combination minimising S2 p2p amplitude.

% Finding the minimum S2 p2p amplitude and the elevation and azimuthal angles which produce it (i.e. minimum S2 p antipodal sum)
min_diameter = min(min(diameter_sph_mat));

[idx_min_diameter_el_sph, idx_min_diameter_az_sph] = find(diameter_sph_mat == min_diameter); % The row (elevation) and column (azimuth) indices of 'diameter_sph_mat' corresponding to the min diameter (in the anterior hemisphere). 
min_diameter_el_sph_mat = idx_min_diameter_el_sph - 91; % The actual elevation angle minimising diameter. 
min_diameter_az_sph_mat = idx_min_diameter_az_sph - 180; % The actual azimuthal angle minimising diameter.

% Finding the antipode (in the posterior hemisphere) of the point minimising diameter (in the anterior hemisphere).
[min_diameter_el_antipode_sph_mat, min_diameter_az_antipode_sph_mat] = antipodal_calculator(min_diameter_el_sph_mat, min_diameter_az_sph_mat);


%% 3 Plotting
if flag_plot == 2
    
    % Plotting the 3D Hull Surface radius spherical projection map in 2D.
    figure; surf(radius_sph_mat); title('S2 3D Hull Surface Radius Spherical Projection Map'); xlabel('Azimuthal Angle (Degrees) {Matlab Reference Frame}'); xticks([0 40 80 120 160 200 240 280 320 360]); xticklabels({-180, -140, -100, -60, -20, 20, 60, 100, 140, 180}); ylabel('Elevation Angle (Degrees) {Matlab Reference Frame}'); yticks([1 21 41 61 81 101 121 141 161 181]); yticklabels({-90,-70,-50,-30,-10,10,30,50,70,90}); zlabel('S2 Peak Amplitude (mg)'); shading interp; colormap spring
    
    % Plotting the 3D Hull Surface hemispherical projection map (i.e.the diameter) in 2D.
    figure; surf(diameter_sph_mat); title('S2 3D Hull Surface Diameter (Anterior) Hemispherical Projection Map (= S2 Radius Antipodal Summation)');  xlabel('Azimuthal Angle (Degrees) {SCG Reference Frame}'); xticks([0 20 40 60 80 100 120 140 160 180]); xticklabels({180, 160, 140, 120, 100, 80, 60, 40, 20, 0}); ylabel('Elevation (Degrees) {SCG Reference Frame}'); yticks([0 11 21 31 41 51 61 71 81 91 101 111 121 131 141 151 161 171 181]); yticklabels({-90,-80,-70,-60,-50,-40,-30,-20,-10,0,10,20,30,40,50,60,70,80,90}); zlabel('S2 Peak-to-Peak Amplitude (mg)'); shading interp; colormap spring
    hold on; plot3(idx_max_diameter_az_sph, idx_max_diameter_el_sph, max_diameter,'r.','MarkerSize',30);
    hold on; plot3(idx_min_diameter_az_sph, idx_min_diameter_el_sph, min_diameter,'g.','MarkerSize',30); legend('S2 3D Hull Surface Diameter', 'Max S2 3D Hull Surface Diameter' , 'Min S2 3D Hull Surface Diameter')
       
end

% Adjusting the indices of max and min diameter and their antipodes for 'radius_sph' and 'radius_x/y/z'.
idx_max_diameter_el_antipode_sph = el_max_diameter_antipode_sph_mat + 91;
idx_max_diameter_az_antipode_sph = az_max_diameter_antipode_sph_mat + 180; 

idx_min_diameter_el_antipode_sph = min_diameter_el_antipode_sph_mat + 91; 
idx_min_diameter_az_antipode_sph = min_diameter_az_antipode_sph_mat + 180;

if flag_plot > 0
    
    % Highlighting the axes maximising and minimising diameter of the 3D Hull Surface on the plot displaying this surface and the S2 acceleration trajectory.
    figure(HullSurface)
    
    hold on; plot3([radius_x_scg(idx_max_diameter_el_sph, idx_max_diameter_az_sph); radius_x_scg(idx_max_diameter_el_antipode_sph, idx_max_diameter_az_antipode_sph)],[radius_z_scg(idx_max_diameter_el_sph, idx_max_diameter_az_sph); radius_z_scg(idx_max_diameter_el_antipode_sph, idx_max_diameter_az_antipode_sph)],[radius_y_scg(idx_max_diameter_el_sph, idx_max_diameter_az_sph); radius_y_scg(idx_max_diameter_el_antipode_sph, idx_max_diameter_az_antipode_sph)],'r--','LineWidth',3);
    hold on; plot3([radius_x_scg(idx_min_diameter_el_sph, idx_min_diameter_az_sph); radius_x_scg(idx_min_diameter_el_antipode_sph, idx_min_diameter_az_antipode_sph)],[radius_z_scg(idx_min_diameter_el_sph, idx_min_diameter_az_sph); radius_z_scg(idx_min_diameter_el_antipode_sph, idx_min_diameter_az_antipode_sph)],[radius_y_scg(idx_min_diameter_el_sph, idx_min_diameter_az_sph); radius_y_scg(idx_min_diameter_el_antipode_sph, idx_min_diameter_az_antipode_sph)],'b--','LineWidth',3); legend('S2 Acceleration Trajectory','S2 3D Hull Surface','Max 3D Hull Surface Diameter','Min 3D Hull Surface Diameter');
    
    clear max_lim min_lim el_max_diameter_antipode_sph_mat az_max_diameter_antipode_sph_mat min_diameter_el_antipode_sph_mat min_diameter_az_antipode_sph_mat

end
    

%% 4 Troubleshooting - Calculating max and min distances (should equal max S2 p2p and min S2 p2p, respectively).
max_distance = sqrt( (radius_x_scg(idx_max_diameter_el_sph, idx_max_diameter_az_sph) - radius_x_scg(idx_max_diameter_el_antipode_sph, idx_max_diameter_az_antipode_sph)).^2 + (radius_z_scg(idx_max_diameter_el_sph, idx_max_diameter_az_sph) - radius_z_scg(idx_max_diameter_el_antipode_sph, idx_max_diameter_az_antipode_sph)).^2 + (radius_y_scg(idx_max_diameter_el_sph, idx_max_diameter_az_sph) - radius_y_scg(idx_max_diameter_el_antipode_sph, idx_max_diameter_az_antipode_sph)).^2 );
min_distance = sqrt( (radius_x_scg(idx_min_diameter_el_sph, idx_min_diameter_az_sph) - radius_x_scg(idx_min_diameter_el_antipode_sph, idx_min_diameter_az_antipode_sph)).^2 + (radius_z_scg(idx_min_diameter_el_sph, idx_min_diameter_az_sph) - radius_z_scg(idx_min_diameter_el_antipode_sph, idx_min_diameter_az_antipode_sph)).^2 + (radius_y_scg(idx_min_diameter_el_sph, idx_min_diameter_az_sph) - radius_y_scg(idx_min_diameter_el_antipode_sph, idx_min_diameter_az_antipode_sph)).^2 );

clear idx_max_diameter_el_sph idx_max_diameter_az_sph idx_max_diameter_el_antipode_sph idx_max_diameter_az_antipode_sph idx_min_diameter_el_sph idx_min_diameter_az_sph idx_min_diameter_el_antipode_sph idx_min_diameter_az_antipode_sph radius_x_scg radius_y_scg radius_z_scg


%% 5 Converting the function output angles into the SCG reference frame.
el_mech_axis = el_max_diameter_sph_mat;
az_mech_axis = -1*az_max_diameter_sph_mat;
el_antimech_axis = min_diameter_el_sph_mat;
az_antimech_axis = -1*min_diameter_az_sph_mat;

clear el_max_diameter_sph_mat az_max_diameter_sph_mat min_diameter_el_sph_mat min_diameter_az_sph_mat

axes_angles = [el_mech_axis, az_mech_axis, el_antimech_axis, az_antimech_axis];
end


%% Service-functions.
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