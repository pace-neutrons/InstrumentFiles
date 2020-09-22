function setup_maps_white_beam
% Set options for white beam vanadium homering

% Remove existing options global variable:
% ----------------------------------------
ixf_global_var('homer_white_beam','remove');

% Set options
% -----------
% Default normalisation is with monitor 1
ixf_global_var('homer_white_beam','set','normalisation',1);
ixf_global_var('homer_white_beam','set','range',[1000,2000]);

% Normalisation scales (mon_scale is instrument dependent)
ixf_global_var('homer_white_beam','set','mon_scale',1000);  % *** why different to mono_van and mono_sample?
ixf_global_var('homer_white_beam','set','peak_scale',1000000);
ixf_global_var('homer_white_beam','set','uamp_scale',1000);

% Particular to white_van operation
ixf_global_var('homer_white_beam','set','ei','white');
ixf_global_var('homer_white_beam','set','det_units','$e');
ixf_global_var('homer_white_beam','set','d_int',[20,300]);
