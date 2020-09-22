function [ eres_tot_ave]=LET_resolution2(ei,ef,freq,shape,specres)
% Monte Carlo for determining  pathlength variance in detectors for LET
% ei is incident energy
% ef is final energy of neutron , ef=ei for elastic
% freq is frequency of resolution choppers
% shape is sample shape =flat, for flatplate or cyl for cyclindrical
% specres.......if it is a number then it gives the resolution for that
% spectra number (assigns run below)
% specres.......if it is a file then it opens a map file 

Lmch=15.67; %mod to chop 5 distance in m. mod is taken at chop 1
Lchs=1.5; % distance from chop 5 to sample
neutron=1000; % number of neutrons for mc simulation of sample shape 
pixel=.0156; % length of a pixel in m assuming a 4 to 1 mapping
tube_res=0.025; % position resolution of the 4m tubes in m
thi=-20; % if a flat plate sample thi is the angle between the normal to the sample and ki
mod_sam=25.0; % mod to sample distance
x_size=0.03;  %sample width in m or diameter if a cylinder
z_size=0.03; % sample height in m
x_range=x_size*cos((thi/180)*pi);
chop5_rad=280; % disk radius in mm
tchop=1e6*10/(2*2*pi*chop5_rad*freq); % chop5 opening time in us, 10 is the chopper opening in mm
tmod=1e6*40/(2*2*pi*chop5_rad*freq); % chop1 or effective moderator width in us
Lami=sqrt(81.81/ei);  % wavelength of incident neutron
Lamf=sqrt(81.81/ef);  % wavelength of neutron after interaction with sample
tmch=252.82.*Lmch*Lami; %  time from chop1 (mod) to chopper
pressure=10; %10 atmospheres in LETs tubes
k=(pressure/10)*0.74*Lamf ; % k is N *sigma for a 10 atm He tube
diam=25.4;   %diameter of deteco of detector

vi=3955.4/Lami; %incident velocity
vf=3955.4/Lamf; % final velocity

to=1e6*(23.5/vi); %opening time of chopper in us
vmax_chop= 1e6*23.5/(to-0.5*tchop); %max velocity of neutrons through chopper
vmin_chop= 1e6*23.5/(to+0.5*tchop); % min velocity of neutrons through chopper
et=ei-ef; % transfered energy to sample
emax_in=5.227e-6*vmax_chop^2;
emin_in=5.227e-6*vmin_chop^2;
emax_out=emax_in-et;
emin_out=emin_in-et;
vmax_f=sqrt(emax_out)*437.3949;
vmin_f=sqrt(emin_out)*437.3949;

tmin_det=1e6*(((23.5+1.4)/vmax_chop)+(3.5/vmax_f));
tmax_det=1e6*(((23.5+1.4)/vmin_chop)+(3.5/vmin_f));
dif_t=(tmax_det-tmin_det);

tom=1e6*(7.83/vi); %opening time of mod chopper in us
vmin_mod= 1e6*15.67/(to-(tom-(0.5*tmod))); %max velocity of neutrons through chopper
vmax_mod= 1e6*15.67/(to-(tom+(0.5*tmod))); % min velocity of neutrons through chopper
et=ei-ef; % transfered energy to sample
emax_in=5.227e-6*vmax_mod^2;
emin_in=5.227e-6*vmin_mod^2;
emax_out=emax_in-et;
emin_out=emin_in-et;
vmax_f=sqrt(emax_out)*437.3949;
vmin_f=sqrt(emin_out)*437.3949;
tmin_det_m=to+1e6*(((1.4)/vmax_mod)+(3.5/vmax_f));
tmax_det_m=to+1e6*(((1.4)/vmin_mod)+(3.5/vmin_f));
dif_tm=(tmax_det_m-tmin_det_m);


ass 1084
sp = gget('spec');
[sp,ind]=sort(sp);   % sort the spectra into ascending number
% get first index where spectrum number=1
ilo=lower_index(sp,1);
if ilo>length(sp) % all spectra are zero length
    delta = [];
    twotheta = [];
    azimuth = [];
    x2 = [];
else

   [delta, twotheta, azimuth, x2] = get_secondary;

   
% ******************** get the average angles and lengths to each spectrum
    nspec = sp(end);  % number of spectra
    ndet  = accumarray(sp(ilo:end)',ones(size(ind(ilo:end))),[nspec,1])';
    index=find(ndet==0);
    ndet(index)=1;
    twotheta= accumarray(sp(ilo:end)',twotheta(ind(ilo:end)),[nspec,1])'./ndet;
    check= accumarray(sp(ilo:end)',sign(azimuth(ind(ilo:end))),[nspec,1])'./ndet;
    azimuth = accumarray(sp(ilo:end)',azimuth(ind(ilo:end)),[nspec,1])'./ndet;
    x2      = accumarray(sp(ilo:end)',x2(ind(ilo:end)),[nspec,1])'./ndet;
    spectra   = accumarray(sp(ilo:end)',sp(ilo:end),[nspec,1])'./ndet;
end

%ind=find(x2<3.505);
%sps=spectra(ind);
%ind=find(sps<24576);
%sps=sps(ind);
%specres=sps;







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isnumeric(specres)
    groups=1;
     map_out = cell(1,1); % create cell array
     map_out{1} = specres;

else
    map_out = map_read (specres);
    groups=length(map_out);
end

for group=1:groups
    specgroup=map_out{group};
    numspec=length(specgroup);
    eres_tot=0;
    fwhm_sam_tot=0;
    for specindex=1:numspec
        indexspec=find(spectra==specgroup(specindex)); %get spectrum 


% *************** convert coordinates from spherical to cartesian

[x,z,y]=sph2cart((azimuth(indexspec(1))./57.3),((twotheta(indexspec(1))-90)./57.3),x2(indexspec(1)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% this section calculates the position error due to the pixel size
zd=z+((rand(neutron,1).*pixel)-(pixel/2))+((tube_res/2.35).*(randn(neutron,1)));  % this gives the random position along the pixel the neutron arrives at
% the second term above is due to error in position from the electronics
% approx 25mm fwhm

%%%%% this section calculates the position error due to the detector shape
  radius=diam/2;
    pos=zeros(neutron*2,2);
    pos(:,1)=(rand(neutron*2,1)*diam)-diam/2;   %determines random x position neutron arrives at
    pos(:,2)=sqrt((radius*radius)-(pos(:,1).*pos(:,1))); %halve thickness of tube at position pos
    depth=10*(-log(1.0-rand(neutron*2,1))/k); % gives depth in mm neutron travels in tube
    a=find(depth<2.0*pos(:,2)); % these are indexs where neutron has been absorbed within the tube thickness
     pos=pos(a,:);
     depth=depth(a);
    pos(:,2)=(-pos(:,2))+depth;
    pos=pos(1:neutron,:); %just take 1st 1000 values
    thi=-atan(x/y); % angle of tube in horizontal plane relative to ki 
    rot=[cos(thi) -sin(thi);sin(thi) cos(thi)]; % rotation matrix around z by thi
    pos=rot*pos';  % the new x,y position errors when rotated by thi
    pos=pos'/1000.0; % put it back into meters
   xd=x+pos(:,1);
   yd=y+pos(:,2);


%%%% this section will calculate time error due to sample shape


    xs=((rand(neutron,1).*x_range)-(x_range/2));   %determines random position neutron arrives at in x
    ys=tan((thi/180)*pi).*xs;
    zs=((rand(neutron,1).*z_size)-(z_size/2)); 
    sam_det=sqrt((xd-xs).^2 +(yd-ys).^2 +(zd-zs).^2);  % this is the range of sample to detecot distances
   

sigma_sam=sqrt(var(sam_det)) ;  %variance of distances to detector
FWHM_sam=2.35*sigma_sam;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% now work out the time variance due to the choppers

Lsd=x2(indexspec(1)); % sample to detector distance
FWHM_chop=tchop*(Lsd+Lmch+Lchs)/Lmch; % FWHM time at detector due to chopper at elastic line
FWHM_mod=tmod*(Lsd+Lchs)/Lmch; % FWHM time at detector due to moderator 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% energy resolutions due to various components

eres_chop=(2*tchop/tmch)*((Lsd+Lmch+Lchs)/Lsd)*((ef/ei)^1.5);
eres_mod=(2*tmod/tmch)*(1+((Lchs/Lsd)*((ef/ei)^1.5)));
eres_sam=(2*FWHM_sam/Lsd)*((ef/ei)^1.5);

eres_tot=eres_tot+sqrt(eres_chop^2 + eres_mod^2 +eres_sam^2);  % % energy resolution due to everything
fwhm_sam_tot=fwhm_sam_tot+FWHM_sam;


    end  % end spectra
    eres_tot_ave(group)=(1000*ei*eres_tot)/numspec;  % gives average resolution in uev for the mapping group
    group;
    fwhm_sam_ave=fwhm_sam_tot/numspec;
end  % end groups

%fprintf(1,'FWHM resolution is %6.2f ueV \n',eres_tot_ave)



