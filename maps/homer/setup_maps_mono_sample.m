function setup_maps_mono_sample
% Set options for monochromatic sample homering

% Remove existing options global variable:
% ----------------------------------------
ixf_global_var('homer_mono_sample','remove');

% Set options
% -----------
% Default normalisation is with monitor 1
ixf_global_var('homer_mono_sample','set','normalisation',1);
ixf_global_var('homer_mono_sample','set','range',[1000,2000]);

% Normalisation scales (mon_scale is instrument dependent)
ixf_global_var('homer_mono_sample','set','mon_scale',1.7016e8);
ixf_global_var('homer_mono_sample','set','peak_scale',1000000);
ixf_global_var('homer_mono_sample','set','uamp_scale',1000);

% Particular to mono_sample operation
ixf_global_var('homer_mono_sample','set','det_units','$w');
ixf_global_var('homer_mono_sample','set','background',[12000, 18000]);
