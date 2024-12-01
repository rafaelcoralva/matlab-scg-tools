function [axes_angles] = cardiac_max_s2p2p_axis_finder_3D(avg_s2_xyz)
% This algorithm projects the average S2 complex along every combination of azimuthal (-180 to 180) and elevation (-90 to 90) angles.
% It then outputs the combination of these angles that maximises the average s2 complex's peak-to-peak (p2p) amplitude.

% INPUTS:
% avg_s2_xyz (Matrix: [length(s1), 3]) - Average x,y, and z SCG S2 components. The three columns correspond to the 3 scg x/y/z components.

% OUTPUTS:
% axes_angles = [max_s2_p2p_angles, min_s2_p2p_angles], where:
%                max_s2_p2p_angles = [max_s2_p2p_el, max_s2_p2p_az];
%                min_s2_p2p_angles = [min_s2_p2p_el, min_s2_p2p_az]; where:
%                max_s2_p2p - refers to the axis maximising s2 p2p amplitude.
%                min_s2_p2p - refers to the axis minimising s2 p2p amplitude.
%                el - refers to the elevation angle in the SCG reference frame (angle from the SCG xz plane [i.e. anatomical transverse plane] TOWARDS the SCG y axis [i.e. anatomical long axis]).
%                ax - referst o the azimuth angle in the SCG reference frame (angle from SCG x axis [i.e. anatomical left medial lateral axis] TOWARDS the SCG z axis [i.e. posterior-anterior axis]).

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

% NOTE 2: The axis maximising average S2 p2p is defined by convention as pointing towards the anterior hemipshere.

flag_plot = 1; % Plot flag. If 0, no plots are produced. If 1, only final plot is produced. If 2, all plots are produced.


%% 1. S2 Peak Amplitude Spherical Projection

% Considering s2 peak amplitude (not p2p amplitude) to then create a spherical mapping of s2 peak amplitude across elevation and azimuth.
% s2 p2p amplitude is then the sum of the antipodal s2 peak amplitudes.
for el_mat = -90:1:90 % el_mat represents the elevation in the Matlab coordinate frame.
    for az_mat = -179:1:180 % az_mat represents the azimuthal angle in the Matlab coordinate frame.

        [xhat_mat,yhat_mat,zhat_mat] = sph2cart(az_mat*pi/180,el_mat*pi/180,1); % The unit-vectors (according to the Matlab reference frame) pointing the direction of the current elevation and azimuth angle combination.

        % Converting to the SCG convention reference frame (this depends on the accelerometer and its placement orientation! - they should be such that they match the conventional SCG reference frame).
        xhat_scg = xhat_mat;
        yhat_scg = zhat_mat;
        zhat_scg = -1*yhat_mat;

        projection_axis_scg = [xhat_scg, yhat_scg, zhat_scg]; % The unit-vector axis (according to the SCG reference frame) pointing in the direction of the current elevation and azimuth angle combination.

        for i = 1:size(avg_s2_xyz,1)
            projection(i) = dot(avg_s2_xyz(i,:), projection_axis_scg); % Using the dot product to perform the projection of the average SCG s2 component along the current azimuthal and elevation angles combination.
        end

        clear i xhat_scg yhat_scg zhat_scg xhat_mat yhat_mat zhat_mat projection_axis_scg

        % Calculating the s2 peak amplitude of the projection and storing it in s2_p_sph_mat - a matrix whose row and column indices represent elevation and azimuth, respectively (in the Matlab spherical reference frame), and whose elements represent the s2 peak ampltudes (i.e. the radii).
        s2_p_sph_mat(el_mat+91, az_mat+180) = max(projection);

        % Cartesian (x,y,z) coordinates of the (index-parametrised) s2_p_sph_mat, in the Matlab reference frame.
        [s2_p_x_mat(el_mat+91,az_mat+180), s2_p_y_mat(el_mat+91,az_mat+180), s2_p_z_mat(el_mat+91,az_mat+180)] = sph2cart(az_mat*pi/180, el_mat*pi/180, max(projection));
        clear projection

    end
end

clear az_mat el_mat

% Converting from Cartesian Matlab reference frame to Cartesian SCG reference frame - this is needed later for plotting.
s2_p_x_scg = s2_p_x_mat;
s2_p_y_scg = s2_p_z_mat;
s2_p_z_scg = -1*s2_p_y_mat;

clear s2_p_x_mat s2_p_y_mat s2_p_z_mat


%% 2. Calculating the s2 peak-to-peak Amplitude (i.e. the antipodal distance) for each combination of elevation and azimuthal angle.
for el_mat = -90:1:90 % Considering the anterior hemisphere (only a single hemisphere needs to be considered as the p2p amplitude captures the opposing hemisphere).
    for az_mat = -179:1:0

        [antipodal_el_mat, antipodal_az_mat] = antipodal_calculator(el_mat, az_mat); % Calculating the antipodal elevation and azimuth for the current elevation and azimuth angle combination.

        % Identifying the indices in s2_p_sph_mat for the s2 peak amplitudes corresponding to the current elevation and azimuth angle combination and its antipode.
        idx_el_mat = el_mat + 91;
        idx_az_mat = az_mat + 180;

        idx_el_antipode_mat = antipodal_el_mat + 91;
        idx_az_antipode_mat = antipodal_az_mat + 180;

        % Calculating the s2 p2p = sum of antipodal s2 peak amplitudes, for that particular elevation and azimuthal angle combination.
        s2_p_antipodal_sum_sph_mat(el_mat + 91, az_mat + 180) = s2_p_sph_mat(idx_el_mat, idx_az_mat) + s2_p_sph_mat(idx_el_antipode_mat, idx_az_antipode_mat);

    end
end

clear el_mat idx_el_mat az_mat idx_az_mat antipodal_el_mat antipodal_az_mat idx_el_antipode_mat idx_az_antipode_mat


%% 3. Finding the angle combination maximising s2 p2p amplitude.

% Finding maximum s2 p2p amplitude and the elevation and azimuthal angles which produce it (i.e.maximum s2 p antipodal sum).
max_s2_p2p = max(s2_p_antipodal_sum_sph_mat,[],'all');

[idx_el_max_s2_p2p_sph, idx_az_max_s2_p2p_sph] = find(s2_p_antipodal_sum_sph_mat == max_s2_p2p); % The row (elevation) and column (azimuth) indices of 's2_p_antipodal_sum_sph_mat' corresponding to the max s2 p2p (in the anterior hemisphere).
el_max_s2_p2p_sph_mat = idx_el_max_s2_p2p_sph - 91; % The actual elevation angle (in the Matlab reference frame) maximising s2 p2p.
az_max_s2_p2p_sph_mat = idx_az_max_s2_p2p_sph - 180; % The actual azimuthal angle (in the Matlab reference frame) maximising s2 p2p.

% Finding the elevation and azimuth angles of the antipode (in the posterior hemisphere) of the point maximising s2 p2p amplitude (in the anterior hemisphere).
[el_max_s2_p2p_antipode_sph_mat, az_max_s2_p2p_antipode_sph_mat] = antipodal_calculator(el_max_s2_p2p_sph_mat,az_max_s2_p2p_sph_mat);


%% 4. Finding the angle combination minimising s2 p2p amplitude.

% Finding the minimum s2 p2p amplitude and the elevation and azimuthal angles which produce it (i.e. minimum s2 p antipodal sum)
min_s2_p2p = min(s2_p_antipodal_sum_sph_mat, [], 'all');

[idx_el_min_s2_p2p_sph, idx_az_min_s2_p2p_sph] = find(s2_p_antipodal_sum_sph_mat == min_s2_p2p); % The row (elevation) and column (azimuth) indices of 's2_p_antipodal_sum_sph_mat' corresponding to the min s2 p2p (in the anterior hemisphere).
el_min_s2_p2p_sph_mat = idx_el_min_s2_p2p_sph - 91; % The actual elevation angle minimising s2 p2p. The '- 1' is to compensate for 1-based indexing.
az_min_s2_p2p_sph_mat = idx_az_min_s2_p2p_sph - 180; % The actual azimuthal angle minimising s2 p2p.

% Finding the antipode (in the posterior hemisphere) of the point minimising s2 p2p amplitude (in the anterior hemisphere).
[el_min_s2_p2p_antipode_sph_mat, az_min_s2_p2p_antipode_sph_mat] = antipodal_calculator(el_min_s2_p2p_sph_mat,az_min_s2_p2p_sph_mat);

% Adjusting the indices of max and min s2 peak amplitudes and their antipodes for 's2_p_sph' and 's2_p_x/y/z'.
idx_el_max_s2_p2p_antipode_sph = el_max_s2_p2p_antipode_sph_mat + 91;
idx_az_max_s2_p2p_antipode_sph = az_max_s2_p2p_antipode_sph_mat + 180;

idx_el_min_s2_p2p_antipode_sph = el_min_s2_p2p_antipode_sph_mat + 91;
idx_az_min_s2_p2p_antipode_sph = az_min_s2_p2p_antipode_sph_mat + 180;


%% 5. Plotting
if flag_plot == 2

    % Plotting the s2 peak amplitude spherical projection map (i.e. in 2D)
    figure;
    surf(s2_p_sph_mat); title('S2 Peak Amplitude Spherical Projection Map'); xlabel('Azimuthal Angle (Degrees) {Matlab Reference Frame}'); xticks([0 40 80 120 160 200 240 280 320 360]); xticklabels({-180, -140, -100, -60, -20, 20, 60, 100, 140, 180}); ylabel('Elevation Angle (Degrees) {Matlab Reference Frame}'); yticks([1 21 41 61 81 101 121 141 161 181]); yticklabels({-90,-70,-50,-30,-10,10,30,50,70,90}); zlabel('s2 Peak Amplitude (mg)'); shading interp; colormap spring;

    % Plotting the s2 peak-to-peak amplitude hemispherical projection map (i.e.the s2 peak amplitude antipodal summation).
    figure;
    surf(s2_p_antipodal_sum_sph_mat); title('S2 P2P Amplitude (Anterior) Hemispherical Projection Map (= S2 Peak Antipodal Summation)');  xlabel('Azimuthal Angle (Degrees) {SCG Reference Frame}'); xticks([0 20 40 60 80 100 120 140 160 180]); xticklabels({180, 160, 140, 120, 100, 80, 60, 40, 20, 0}); ylabel('Elevation (Degrees) {SCG Reference Frame}'); yticks([0 11 21 31 41 51 61 71 81 91 101 111 121 131 141 151 161 171 181]); yticklabels({-90,-80,-70,-60,-50,-40,-30,-20,-10,0,10,20,30,40,50,60,70,80,90}); zlabel('S2 Peak-to-Peak Amplitude (mg)'); shading interp; colormap spring;
    hold on; plot3(idx_az_max_s2_p2p_sph, idx_el_max_s2_p2p_sph, max_s2_p2p,'r.','MarkerSize',30); legend('S2 P2P Amplitude', 'Max S2 P2P Amplitude')
    hold on; plot3(idx_az_min_s2_p2p_sph, idx_el_min_s2_p2p_sph, min_s2_p2p,'g.','MarkerSize',30); legend('S2 P2P Amplitude', 'Max S2 P2P Amplitude' , 'Min S2 P2P Amplitude')

end

if flag_plot > 0 % Plotting the s2 acceleration trajectory, the s2 peak amplitude spherical projection (in 3D), and highlighting the axes maximising and minimising s2 p2p.

    figure; plot3(avg_s2_xyz(:,1),avg_s2_xyz(:,3),avg_s2_xyz(:,2),'k','LineWidth',2); % s2 acceleration trajectory.
    hold on; surface(s2_p_x_scg, s2_p_z_scg, s2_p_y_scg, sqrt(s2_p_x_scg.^2 + s2_p_z_scg.^2 + s2_p_y_scg.^2),'FaceAlpha',0.5,'EdgeAlpha',0); title('S2 Peak Amplitude Projection'); xlabel('x (mg)'); ylabel('z (mg)'); zlabel('y (mg)'); hcb = colorbar; hcb.Title.String = "S2 Peak Amplitude (mg)"; colormap spring; hold on;  % s2 Peak Amplitude Projection.

    hold on; plot3([s2_p_x_scg(idx_el_max_s2_p2p_sph, idx_az_max_s2_p2p_sph); s2_p_x_scg(idx_el_max_s2_p2p_antipode_sph, idx_az_max_s2_p2p_antipode_sph)],[s2_p_z_scg(idx_el_max_s2_p2p_sph, idx_az_max_s2_p2p_sph); s2_p_z_scg(idx_el_max_s2_p2p_antipode_sph, idx_az_max_s2_p2p_antipode_sph)],[s2_p_y_scg(idx_el_max_s2_p2p_sph, idx_az_max_s2_p2p_sph); s2_p_y_scg(idx_el_max_s2_p2p_antipode_sph, idx_az_max_s2_p2p_antipode_sph)],'r--','LineWidth',3);
    hold on; plot3([s2_p_x_scg(idx_el_min_s2_p2p_sph, idx_az_min_s2_p2p_sph); s2_p_x_scg(idx_el_min_s2_p2p_antipode_sph, idx_az_min_s2_p2p_antipode_sph)],[s2_p_z_scg(idx_el_min_s2_p2p_sph, idx_az_min_s2_p2p_sph); s2_p_z_scg(idx_el_min_s2_p2p_antipode_sph, idx_az_min_s2_p2p_antipode_sph)],[s2_p_y_scg(idx_el_min_s2_p2p_sph, idx_az_min_s2_p2p_sph); s2_p_y_scg(idx_el_min_s2_p2p_antipode_sph, idx_az_min_s2_p2p_antipode_sph)],'b--','LineWidth',3); legend('S2 Acceleration Trajectory','S2 Peak Amplitude Projection','Max S2 p2p Axis','Min S2 p2p Axis');
    clear hcb

    min_lim = min(min([s2_p_x_scg; s2_p_y_scg; s2_p_z_scg])); max_lim = max(max(([s2_p_x_scg; s2_p_y_scg; s2_p_z_scg])));
    xlim([min_lim(1) max_lim(1)]); ylim([min_lim(1) max_lim(1)]); zlim([min_lim(1) max_lim(1)]);

    clear max_lim min_lim

end


%% 6. Troubleshooting - Calclating max and min distances (should equal max s2 p2p and min s2 p2p, respectively).
max_distance = sqrt( (s2_p_x_scg(idx_el_max_s2_p2p_sph, idx_az_max_s2_p2p_sph) - s2_p_x_scg(idx_el_max_s2_p2p_antipode_sph, idx_az_max_s2_p2p_antipode_sph)).^2 + (s2_p_z_scg(idx_el_max_s2_p2p_sph, idx_az_max_s2_p2p_sph) - s2_p_z_scg(idx_el_max_s2_p2p_antipode_sph, idx_az_max_s2_p2p_antipode_sph)).^2 + (s2_p_y_scg(idx_el_max_s2_p2p_sph, idx_az_max_s2_p2p_sph) - s2_p_y_scg(idx_el_max_s2_p2p_antipode_sph, idx_az_max_s2_p2p_antipode_sph)).^2 );
min_distance = sqrt( (s2_p_x_scg(idx_el_min_s2_p2p_sph, idx_az_min_s2_p2p_sph) - s2_p_x_scg(idx_el_min_s2_p2p_antipode_sph, idx_az_min_s2_p2p_antipode_sph)).^2 + (s2_p_z_scg(idx_el_min_s2_p2p_sph, idx_az_min_s2_p2p_sph) - s2_p_z_scg(idx_el_min_s2_p2p_antipode_sph, idx_az_min_s2_p2p_antipode_sph)).^2 + (s2_p_y_scg(idx_el_min_s2_p2p_sph, idx_az_min_s2_p2p_sph) - s2_p_y_scg(idx_el_min_s2_p2p_antipode_sph, idx_az_min_s2_p2p_antipode_sph)).^2 );


%% 7. Converting the function output angles into the SCG reference frame.
el_max_s2_p2p_sph_scg = el_max_s2_p2p_sph_mat;
az_max_s2_p2p_sph_scg = -1*az_max_s2_p2p_sph_mat;
el_min_s2_p2p_sph_scg = el_min_s2_p2p_sph_mat;
az_min_s2_p2p_sph_scg = -1*az_min_s2_p2p_sph_mat;

clear el_max_s2_p2p_sph_mat az_max_s2_p2p_sph_mat el_min_s2_p2p_sph_mat az_min_s2_p2p_sph_mat

max_s2_p2p_angles = [el_max_s2_p2p_sph_scg, az_max_s2_p2p_sph_scg];
min_s2_p2p_angles = [el_min_s2_p2p_sph_scg, az_min_s2_p2p_sph_scg];

axes_angles = [max_s2_p2p_angles, min_s2_p2p_angles]; % The output of the function is the max s2 elevation and azimuth angles and the min s2 elevation and azimuth angles in the SCG reference frame.

end


%% Service-functions.
function [el_antipode, az_antipode] = antipodal_calculator(el, az) % Function to return the antipode of an input angle given by its elevation and azimuth.

el_antipode = -1*el;

if az < 0
    az_antipode = az + 180;

elseif az > 0
    az_antipode = az - 180;
    
elseif az == 0
    az_antipode = 180;
end

end
