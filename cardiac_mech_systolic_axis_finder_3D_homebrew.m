function [axis_angles] = cardiac_mech_systolic_axis_finder_3D_homebrew(avg_s1_xyz)
% This algorithm identifies the cardiac mechanical systolic axis in three dimensions. It consists of two parts.

% 1. Calculation of the 3D hull using a home-brew method.
% 2. Identification of the azimuthal and elevation angles that produce the maximum 3D hull diameter.

% INPUTS:
% avg_s1_xyz - 3D average S1 acceleration trajectory matrix of size [length(avg_s1) 3],where the three columns correspond to the 3 scgxyz components.
%
% OUTPUTS:
% axis_angles = [elevation azimuth], where:
%               elevation - elevation angle of max 3D hull diameter (angle from the SCG xz plane TOWARDS the SCG y axis).
%               azimuth - azimuthal angle of max 3D hull diameter (angle from SCG x axis TOWARDS the SCG z axis).

% NOTE ON COORDINATE FRAMES: 
% For the recorded SCG signals:
% - the SCG Cartesian coordinate system is: x = medial-(left) lateral.
%                                           y = long axis.
%                                           z = posterior-anterior.
% - the SCG Spherical coordinate system is: azimuth (theta) = from x towards z.
%                                           elevation (phi) = from the xz plane towards y.
%                                           radius = modulus of x,y, and z.
%
% This SCG reference frame is NOT the same as that considered by various Matlab functions:
% - the Matlab Cartesian coordinate system is: x = medial-(left) lateral.
%                                              y = anterior-posterior.
%                                              z = long axis.
% - the Matlab Spherical coordinate system is: azimuth (theta) = from x towards y.
%                                              elevation (phi) = from the xy plane towards y.
%                                              radius = modulus of x,y, and z.
%
% Hence, conversions must be made in the code btween SCG and Matlab coordinate frames. These  conversions are:
% - x_scg = x_mat
% - y_scg = z_mat
% - z_scg = -y_mat
% - azimuth_scg = { 180 - azimuth_mat   if azimuth_mat > 0
%                   abs(azimuth)mat     if azimuth)mat < 0
% - elevation_scg = elevation_mat
% - radius_scg = radius_mat
%% 1. Calculation of the hull.

s1_AzimuthElevationR_mat = zeros(size(avg_s1_xyz,1),3); % Initialization - matrix of azimuthal angle (column 1), elevation angle (column 2), and radius (column 3), in the Matlab reference frame.

[s1_AzimuthElevationR_mat(:,1), s1_AzimuthElevationR_mat(:,2), s1_AzimuthElevationR_mat(:,3)] = cart2sph(avg_s1_xyz(:,1), -1*avg_s1_xyz(:,3), avg_s1_xyz(:,2)); % Conversion from SCG cartesian (x,y,z) to Matlab spherical (azimuth, elevation, radius).

% Converting azimuth and elevation angles to degrees
s1_AzimuthElevationR_mat(:,1) = (180/pi)*s1_AzimuthElevationR_mat(:,1);
s1_AzimuthElevationR_mat(:,2) = (180/pi)*s1_AzimuthElevationR_mat(:,2);

angular_resolution = 2; % In degrees. This is the angular precision with which we may identify the cardiac mechanical systolic axis. The smaller this is, the greater the computational cost.

% Rounding according to the angular resolution
s1_AzimuthElevationR_mat(:,1) = round(s1_AzimuthElevationR_mat(:,1)/angular_resolution)*angular_resolution;
s1_AzimuthElevationR_mat(:,2) = round(s1_AzimuthElevationR_mat(:,2)/angular_resolution)*angular_resolution;

% Sorting in order of increasing elevation angle and then increasing azimuthal angle.
s1_AzimuthElevationR_mat = sortrows(s1_AzimuthElevationR_mat,[2,1]);

s1_hull_AzimuthElevationR_mat = s1_AzimuthElevationR_mat; % Initialization of the 3D hull.

% Discarding "interior" elements of the hull (i.e. rows with identical elevation and azimuthal angles but non-maximal radii)
for i=size(s1_AzimuthElevationR_mat,1)-1:-1:1 % For all the rows of the 3D S1 acceleration trajectory...
    if s1_AzimuthElevationR_mat(i,2) ==  s1_AzimuthElevationR_mat(i+1,2) && s1_AzimuthElevationR_mat(i,1) ==  s1_AzimuthElevationR_mat(i+1,1) % If elevation and azimuthal angles of current and previous rows are identical...
        
        % Delete hull row with smallest radius.
        if s1_AzimuthElevationR_mat(i,3) > s1_AzimuthElevationR_mat(i+1,3)
            s1_hull_AzimuthElevationR_mat(i+1,:)=[];
        elseif s1_AzimuthElevationR_mat(i,3) <= s1_AzimuthElevationR_mat(i+1,3) 
            s1_hull_AzimuthElevationR_mat(i,:)=[];
        end
        
    end  
end

%% Plotting 
[hull_x_mat, hull_y_mat, hull_z_mat] = sph2cart((pi/180)*s1_hull_AzimuthElevationR_mat(:,1),(pi/180)*s1_hull_AzimuthElevationR_mat(:,2),s1_hull_AzimuthElevationR_mat(:,3));
hull_mat = [hull_x_mat, hull_y_mat, hull_z_mat];

% Converting to the SCG Cartesian reference frame
hull_x_scg = hull_x_mat;
hull_y_scg = hull_z_mat;
hull_z_scg = -1*hull_y_mat;
hull_scg = [hull_x_scg, hull_y_scg, hull_z_scg];

figure;
title('S1 Hull (home-brew)'); xlabel('x (mg)'); ylabel('z (mg)'); zlabel('y (mg)');

% Plotting the original S1 acceleration trajectory
plot3(avg_s1_xyz(:,1),avg_s1_xyz(:,3),avg_s1_xyz(:,2),'k','LineWidth',2); hold on;

% Plotting the isolated elements (i.e. vertices) of the hull
plot3(hull_scg(:,1), hull_scg(:,3), hull_scg(:,2),'b.'); 

% Plotting the hull surface
hull_idx = boundary(hull_x_scg,hull_z_scg,hull_y_scg,0.99);
trisurf(hull_idx,hull_x_scg,hull_z_scg,hull_y_scg,'Facecolor','blue','FaceAlpha',0.2,'EdgeAlpha',0); 

min_lim = min(min(avg_s1_xyz)); max_lim = max(max(avg_s1_xyz));
xlim([min_lim(1) max_lim(1)]); ylim([min_lim(1) max_lim(1)]); zlim([min_lim(1) max_lim(1)]);
legend('S1 Acceleration Trajectory','S1 Hull Vertices','S1 Hull Surface')

%% 2. Identification of the azimuthal and elevation angles that produce the maximum antipode.

azimuth=180;
elevation=45;

axis_angles = [azimuth elevation];

end