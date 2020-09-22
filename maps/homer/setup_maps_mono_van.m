function setup_maps_mono_van
% Set options for monochromatic vanadium homering

% Remove existing options global variable:
% ----------------------------------------
ixf_global_var('homer_mono_van','remove');

% Set options
% -----------
% Default normalisation is with monitor 1
ixf_global_var('homer_mono_van','set','normalisation',1);
ixf_global_var('homer_mono_van','set','range',[1000,2000]);

% Normalisation scales (mon_scale is instrument dependent)
ixf_global_var('homer_mono_van','set','mon_scale',1.7016e8);
ixf_global_var('homer_mono_van','set','peak_scale',1000000);
ixf_global_var('homer_mono_van','set','uamp_scale',1000);

% Particular to mono_van operation
ixf_global_var('homer_mono_van','set','det_units','$w');
ixf_global_var('homer_mono_van','set','background',[12000, 18000]);
ixf_global_var('homer_mono_van','set','d_int',[-1 1]);
