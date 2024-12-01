function [axes_angles] = cardiac_max_s1p2p_axis_finder(avg_s1_xyz)
% This algorithm projects the average S1 complex along every combination of azimuthal (-180 to 180) and elevation (-90 to 90) angles.
% It then outputs the angles that maximise the average S1 complex's peak-to-peak (p2p) amplitude, as well as the angles that minimise it.

% INPUTS:
% avg_s1_xyz (Matrix: [length(s1), 3]) - Average x,y, and z SCG S1 components. The three columns correspond to the 3 x/y/z components.

% OUTPUTS:
% axes_angles = [max_s1_p2p_angles, min_s1_p2p_angles], where:
%                max_s1_p2p_angles = [max_s1_p2p_el, max_s1_p2p_az];
%                min_s1_p2p_angles = [min_s1_p2p_el, min_s1_p2p_az]; where:
%                max_s1_p2p - refers to the axis maximising S1 p2p amplitude.
%                min_s1_p2p - refers to the axis minimising S1 p2p amplitude.
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

% NOTE 2: The axis maximising average S1 p2p is defined by convention as pointing towards the anterior hemipshere.

flag_plot = 1; % Plot flag. If 0, no plots are produced. If 1, only final plot is produced. If 2, all plots are produced.


%% 1. S1 Peak Amplitude Spherical Projection
    
% Considering S1 peak amplitude (not p2p amplitude) to then create a spherical mapping of S1 peak amplitude across elevation and azimuth. 
% S1 p2p amplitude is then the sum of the antipodal S1 peak amplitudes.

% Initializations for speed
s1_p_sph_mat = zeros([181,360]); % Contains the S1 peak amplitude for each elevation + azimuth combination.
s1_p_x_mat = zeros([181,360]); % Contains the cartesian coordinates in the Matlab reference frame of the S1 peak amplitude point.
s1_p_y_mat = zeros([181,360]);
s1_p_z_mat = zeros([181,360]);

for el_mat = -90:1:90 % el_mat represents the elevation in the Matlab coordinate frame.
    for az_mat = -179:1:180 % az_mat represents the azimuthal angle in the Matlab coordinate frame.

        [xhat_mat, yhat_mat, zhat_mat] = sph2cart(az_mat*pi/180,el_mat*pi/180,1); % The unit-vectors (according to the Matlab reference frame) pointing the direction of the current elevation and azimuth angle combination.

        % Converting to the SCG convention reference frame (this depends on the accelerometer and its placement orientation! - they should be such that they match the conventional SCG reference frame).
        xhat_scg = xhat_mat;
        yhat_scg = zhat_mat;
        zhat_scg = -1*yhat_mat;

        projection_axis_scg = [xhat_scg, yhat_scg, zhat_scg]; % The unit-vector axis (according to the SCG reference frame) pointing in the direction of the current elevation and azimuth angle combination.

        projection = zeros(length(avg_s1_xyz),1); % Initializing for speed.
        for i = 1:size(avg_s1_xyz,1)
            projection(i) = dot(avg_s1_xyz(i,:), projection_axis_scg); % Using the dot product to perform the projection of the average SCG S1 component along the current azimuthal and elevation angles combination.
        end

        clear i xhat_scg yhat_scg zhat_scg xhat_mat yhat_mat zhat_mat projection_axis_scg

        % Calculating the S1 peak amplitude of the projection and storing it in s1_p_sph_mat 
        % - a matrix whose row and column indices represent elevation and azimuth, respectively (in the Matlab spherical reference frame), 
        % and whose elements represent the S1 peak ampltudes (i.e. the radii).
        s1_p_sph_mat(el_mat+91, az_mat+180) = max(projection); % The +91 and +180 are to allow for correct indexing.

        % Cartesian (x,y,z) coordinates of the s1_p_sph_mat, in the Matlab reference frame.
        [s1_p_x_mat(el_mat+91, az_mat+180), s1_p_y_mat(el_mat+91, az_mat+180), s1_p_z_mat(el_mat+91,az_mat+180)] = sph2cart(az_mat*pi/180, el_mat*pi/180, max(projection));
        clear projection

    end
end

clear az_mat el_mat

% Converting from Cartesian Matlab reference frame to Cartesian SCG reference frame - this is needed later for plotting.
s1_p_x_scg = s1_p_x_mat;
s1_p_y_scg = s1_p_z_mat;
s1_p_z_scg = -1*s1_p_y_mat;

clear s1_p_x_mat s1_p_y_mat s1_p_z_mat


%% 2. Calculating the S1 peak-to-peak Amplitude (i.e. the antipodal distance) for each combination of elevation and azimuthal angle.
s1_p_antipodal_sum_sph_mat = zeros(181, 180); % Initialization for speed

for el_mat = -90:1:90 % Considering the anterior hemisphere (only a single hemisphere needs to be considered as the p2p amplitude captures also the opposing hemisphere).
    for az_mat = -179:1:0

        [el_antipode_mat, az_antipode_mat] = antipodal_calculator(el_mat, az_mat); % Calculating the antipodal elevation and azimuth for the current elevation and azimuth angle combination.

        % Identifying the indices in s1_p_sph_mat for the S1 peak amplitudes corresponding to the current elevation and azimuth angle combination and its antipode.
        idx_el_mat = el_mat + 91;
        idx_az_mat = az_mat + 180;

        idx_el_antipode_mat = el_antipode_mat + 91;
        idx_az_antipode_mat = az_antipode_mat + 180;

        % Calculating the S1 p2p = sum of antipodal S1 peak amplitudes, for that particular elevation and azimuthal angle combination.
        s1_p_antipodal_sum_sph_mat(el_mat + 91, az_mat + 180) = s1_p_sph_mat(idx_el_mat, idx_az_mat) + s1_p_sph_mat(idx_el_antipode_mat, idx_az_antipode_mat);

    end
end

clear el_mat idx_el_mat az_mat idx_az_mat el_antipode_mat az_antipode_mat idx_el_antipode_mat idx_az_antipode_mat


%% 3. Finding the angle combination maximising S1 p2p amplitude.

% Finding maximum S1 p2p amplitude and the elevation and azimuthal angles which produce it (i.e.maximum S1 p antipodal sum).
max_s1_p2p = max(max(s1_p_antipodal_sum_sph_mat));

[idx_el_max_s1_p2p_sph, idx_az_max_s1_p2p_sph] = find(s1_p_antipodal_sum_sph_mat == max_s1_p2p); % The row (elevation) and column (azimuth) indices of 's1_p_antipodal_sum_sph_mat' corresponding to the max s1 p2p (in the anterior hemisphere).
el_max_s1_p2p_sph_mat = idx_el_max_s1_p2p_sph - 91; % The actual elevation angle (in the Matlab reference frame) maximising S1 p2p.
az_max_s1_p2p_sph_mat = idx_az_max_s1_p2p_sph - 180; % The actual azimuthal angle (in the Matlab reference frame) maximising S1 p2p.

% Finding the elevation and azimuth angles of the antipode (in the posterior hemisphere) of the point maximising S1 p2p amplitude (in the anterior hemisphere).
[el_max_s1_p2p_antipode_sph_mat, az_max_s1_p2p_antipode_sph_mat] = antipodal_calculator(el_max_s1_p2p_sph_mat, az_max_s1_p2p_sph_mat);


%% 4. Finding the angle combination minimising S1 p2p amplitude.

% Finding the minimum S1 p2p amplitude and the elevation and azimuthal angles which produce it (i.e. minimum S1 p antipodal sum)
min_s1_p2p = min(s1_p_antipodal_sum_sph_mat, [], 'all');

[idx_min_s1_p2p_el_sph, idx_min_s1_p2p_az_sph] = find(s1_p_antipodal_sum_sph_mat == min_s1_p2p); % The row (elevation) and column (azimuth) indices of 's1_p_antipodal_sum_sph_mat' corresponding to the min s1 p2p (in the anterior hemisphere).
el_min_s1_p2p_sph_mat = idx_min_s1_p2p_el_sph - 91; % The actual elevation angle minimising S1 p2p. 
az_min_s1_p2p_sph_mat = idx_min_s1_p2p_az_sph - 180; % The actual azimuthal angle minimising S1 p2p.

% Finding the antipode (in the posterior hemisphere) of the point minimising S1 p2p amplitude (in the anterior hemisphere).
[el_min_s1_p2p_antipode_sph_mat, az_min_s1_p2p_antipode_sph_mat] = antipodal_calculator(el_min_s1_p2p_sph_mat, az_min_s1_p2p_sph_mat);

% Adjusting the indices of max and min S1 peak amplitudes and their antipodes for 's1_p_sph' and 's1_p_x/y/z'.
idx_el_max_s1_p2p_antipode_sph = el_max_s1_p2p_antipode_sph_mat + 91;
idx_az_max_s1_p2p_antipode_sph = az_max_s1_p2p_antipode_sph_mat + 180;

idx_el_min_s1_p2p_antipode_sph = el_min_s1_p2p_antipode_sph_mat + 91;
idx_az_min_s1_p2p_antipode_sph = az_min_s1_p2p_antipode_sph_mat + 180;


%% 5. Plotting
if flag_plot == 2

    % Plotting the S1 peak amplitude spherical projection map (i.e. in 2D)
    figure; 
    surf(s1_p_sph_mat); title('S1 Peak Amplitude Spherical Projection Map'); xlabel('Azimuthal Angle (Degrees) {Matlab Reference Frame}'); xticks([0 40 80 120 160 200 240 280 320 360]); xticklabels({-180, -140, -100, -60, -20, 20, 60, 100, 140, 180}); ylabel('Elevation Angle (Degrees) {Matlab Reference Frame}'); yticks([1 21 41 61 81 101 121 141 161 181]); yticklabels({-90,-70,-50,-30,-10,10,30,50,70,90}); zlabel('S1 Peak Amplitude (mg)'); shading interp

    % Plotting the S1 peak-to-peak amplitude hemispherical projection map (i.e.the S1 peak amplitude antipodal summation).
    figure; 
    surf(s1_p_antipodal_sum_sph_mat); title('S1 P2P Amplitude (Anterior) Hemispherical Projection Map (= S1 Peak Antipodal Summation)');  xlabel('Azimuthal Angle (Degrees) {SCG Reference Frame}'); xticks([0 20 40 60 80 100 120 140 160 180]); xticklabels({180, 160, 140, 120, 100, 80, 60, 40, 20, 0}); ylabel('Elevation (Degrees) {SCG Reference Frame}'); yticks([0 11 21 31 41 51 61 71 81 91 101 111 121 131 141 151 161 171 181]); yticklabels({-90,-80,-70,-60,-50,-40,-30,-20,-10,0,10,20,30,40,50,60,70,80,90}); zlabel('S1 Peak-to-Peak Amplitude (mg)'); shading interp;
    hold on; plot3(idx_az_max_s1_p2p_sph, idx_el_max_s1_p2p_sph, max_s1_p2p,'r.','MarkerSize',30); legend('S1 P2P Amplitude', 'Max S1 P2P Amplitude')
    hold on; plot3(idx_min_s1_p2p_az_sph, idx_min_s1_p2p_el_sph, min_s1_p2p,'g.','MarkerSize',30); legend('S1 P2P Amplitude', 'Max S1 P2P Amplitude' , 'Min S1 P2P Amplitude')

end

if flag_plot == 1 % Plotting the S1 acceleration trajectory, the S1 peak amplitude spherical projection (in 3D), and highlighting the axes maximising and minimising S1 p2p.
 
    figure; 
    plot3(avg_s1_xyz(:,1),avg_s1_xyz(:,3),avg_s1_xyz(:,2),'k','LineWidth',2); % S1 acceleration trajectory.
    hold on; surface(s1_p_x_scg, s1_p_z_scg, s1_p_y_scg, sqrt(s1_p_x_scg.^2 + s1_p_z_scg.^2 + s1_p_y_scg.^2),'FaceAlpha',0.5,'EdgeAlpha',0); title('S1 Peak Amplitude Projection'); xlabel('x (mg)'); ylabel('z (mg)'); zlabel('y (mg)'); hcb = colorbar; hcb.Title.String = "S1 Peak Amplitude (mg)"; hold on; % S1 Peak Amplitude Projection.
    hold on; plot3([s1_p_x_scg(idx_el_max_s1_p2p_sph, idx_az_max_s1_p2p_sph); s1_p_x_scg(idx_el_max_s1_p2p_antipode_sph, idx_az_max_s1_p2p_antipode_sph)],[s1_p_z_scg(idx_el_max_s1_p2p_sph, idx_az_max_s1_p2p_sph); s1_p_z_scg(idx_el_max_s1_p2p_antipode_sph, idx_az_max_s1_p2p_antipode_sph)],[s1_p_y_scg(idx_el_max_s1_p2p_sph, idx_az_max_s1_p2p_sph); s1_p_y_scg(idx_el_max_s1_p2p_antipode_sph, idx_az_max_s1_p2p_antipode_sph)],'r--','LineWidth',3);
    hold on; plot3([s1_p_x_scg(idx_min_s1_p2p_el_sph, idx_min_s1_p2p_az_sph); s1_p_x_scg(idx_el_min_s1_p2p_antipode_sph, idx_az_min_s1_p2p_antipode_sph)],[s1_p_z_scg(idx_min_s1_p2p_el_sph, idx_min_s1_p2p_az_sph); s1_p_z_scg(idx_el_min_s1_p2p_antipode_sph, idx_az_min_s1_p2p_antipode_sph)],[s1_p_y_scg(idx_min_s1_p2p_el_sph, idx_min_s1_p2p_az_sph); s1_p_y_scg(idx_el_min_s1_p2p_antipode_sph, idx_az_min_s1_p2p_antipode_sph)],'b--','LineWidth',3); legend('S1 Acceleration Trajectory','S1 Peak Amplitude Projection','Max S1 p2p Axis','Min S1 p2p Axis');
    clear hcb

    min_lim = min(min([s1_p_x_scg; s1_p_y_scg; s1_p_z_scg])); max_lim = max(max(([s1_p_x_scg; s1_p_y_scg; s1_p_z_scg])));
    xlim([min_lim(1) max_lim(1)]); ylim([min_lim(1) max_lim(1)]); zlim([min_lim(1) max_lim(1)]);

    clear max_lim min_lim
end


%% 6. Troubleshooting - Calclating max and min distances (should equal max S1 p2p and min S1 p2p, respectively).
max_distance = sqrt( (s1_p_x_scg(idx_el_max_s1_p2p_sph, idx_az_max_s1_p2p_sph) - s1_p_x_scg(idx_el_max_s1_p2p_antipode_sph, idx_az_max_s1_p2p_antipode_sph)).^2 + (s1_p_z_scg(idx_el_max_s1_p2p_sph, idx_az_max_s1_p2p_sph) - s1_p_z_scg(idx_el_max_s1_p2p_antipode_sph, idx_az_max_s1_p2p_antipode_sph)).^2 + (s1_p_y_scg(idx_el_max_s1_p2p_sph, idx_az_max_s1_p2p_sph) - s1_p_y_scg(idx_el_max_s1_p2p_antipode_sph, idx_az_max_s1_p2p_antipode_sph)).^2 );
min_distance = sqrt( (s1_p_x_scg(idx_min_s1_p2p_el_sph, idx_min_s1_p2p_az_sph) - s1_p_x_scg(idx_el_min_s1_p2p_antipode_sph, idx_az_min_s1_p2p_antipode_sph)).^2 + (s1_p_z_scg(idx_min_s1_p2p_el_sph, idx_min_s1_p2p_az_sph) - s1_p_z_scg(idx_el_min_s1_p2p_antipode_sph, idx_az_min_s1_p2p_antipode_sph)).^2 + (s1_p_y_scg(idx_min_s1_p2p_el_sph, idx_min_s1_p2p_az_sph) - s1_p_y_scg(idx_el_min_s1_p2p_antipode_sph, idx_az_min_s1_p2p_antipode_sph)).^2 );


%% 7. Converting the function output angles into the SCG reference frame.
el_max_s1_p2p = el_max_s1_p2p_sph_mat;
az_max_s1_p2p = -1*az_max_s1_p2p_sph_mat;
el_min_s1_p2p = el_min_s1_p2p_sph_mat;
az_min_s1_p2p = -1*az_min_s1_p2p_sph_mat;

clear el_max_s1_p2p_sph_mat az_max_s1_p2p_sph_mat el_min_s1_p2p_sph_mat az_min_s1_p2p_sph_mat

max_s1_p2p_angles = [el_max_s1_p2p, az_max_s1_p2p];
min_s1_p2p_angles = [el_min_s1_p2p, az_min_s1_p2p];

axes_angles = [max_s1_p2p_angles, min_s1_p2p_angles]; % The output of the function is the max S1 elevation and azimuth angles and the min S1 elevation and azimuth angles in the SCG reference frame.
end


%% Service functions.
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

