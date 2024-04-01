function plotting_preprocessed_scg_xyz(preprocessed_scg_xyz, raw_scg_xyz, fs)
% Inside of Wrapped SCG Preprocessing Plotting Function to plot the raw and
% preprocessed SCGs.

t = 0:1/fs:((size(preprocessed_scg_xyz,1)-1))/fs;

fig_scg_preprocessing = figure; 
% SCG x
ax(1) = subplot(6,1,1); plot(t,raw_scg_xyz(:,1)); ylabel('Acceleration (mg)'); box off; title('Raw SCG x'); 
ax(2) = subplot(6,1,2); plot(t,preprocessed_scg_xyz(:,1),'r'); ylabel('Acceleration (mg)'); box off; title('Preprocessed SCG x'); 

% SCG y
ax(3) = subplot(6,1,3); plot(t,raw_scg_xyz(:,2)); ylabel('Acceleration (mg)'); box off; title('Raw SCG y'); 
ax(4) = subplot(6,1,4); plot(t,preprocessed_scg_xyz(:,2),'r'); ylabel('Acceleration (mg)'); box off; title('Preprocessed SCG y'); 

% SCG z
ax(5) = subplot(6,1,5); plot(t,raw_scg_xyz(:,3)); ylabel('Acceleration (mg)'); box off; title('Raw SCG z'); 
ax(6) = subplot(6,1,6); plot(t,preprocessed_scg_xyz(:,3),'r'); ylabel('Acceleration (mg)'); box off; title('Preprocessed SCG z'); 

xlabel('Time (sec)');  linkaxes(ax,'x')

end

