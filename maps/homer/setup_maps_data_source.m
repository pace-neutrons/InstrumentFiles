function DSO=setup_maps_data_source
% Creates default data source object

DSO=IXTdata_source;

% Required by higher level homer commands:
% ----------------------------------------
% (instrument population requires them; moderator for its location for Ei estimates, source for frequency?)
DSO=add_item(DSO,'inst_nxs:::source.nxs','source');
DSO=add_item(DSO,'inst_nxs:::moderator.nxs','moderator');

% Optional instrument population:
% -------------------------------
%  DSO=add_item(DSO,'inst_nxs:::attenuator.nxs','attenuator')
%  DSO=add_item(DSO,'inst_nxs:::aperture.nxs','aperture')

