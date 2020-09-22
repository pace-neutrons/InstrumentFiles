function Multi_rep(efocus,efocus2, freq1)

%efocus is the energy all choppers are focused on
%freq1 is frequency of main resolution choppers
% efocus2 is what the frame overlap chopper focuses on

freqcont=freq1/2;  % this is freq of contaminant removal chopper
%if (int16(freqcont/7)~=freqcont/7);
%    fprintf(1,'resolution frequency /2 must be a multiple of 5Hz');
%    return
%end


dist=[9.3 10.1]; % distance to each chopper (m)
freq=[50 freq1]; % frequency of each chopper
nslot=[1 2]; % number of slots in each chopper assumed equally spaced
width=[950 10]; % width of windows (mm) in each disk (assumed the same)
g_wid=[64 10]; % width of guide opening at each chopper (mm)
radius=[250 290]; % radius in mm of each disk at centre of window
counter=[1 1]; % 2 means a counter rotating chopper and 1 is a single disk



% wavelength to focus choppers on
lamfoc=sqrt(81.82/efocus);
lamfoc2=sqrt(81.82/efocus2);
lf=[lamfoc2 lamfoc]


N=length(dist);  %number of choppers
samp_det=2.5  ;% sample to detector distance
chop_samp=1.8 ;%final chopper to sample distance


tm=5 ;% maximum moderator time us
% temp=10K standard
temp=5;% sample temperature in K
%temp=500;% sample temperature in K
frac_energy=0.80;%fraction of Ei to plot energy loss lines, i.e this is 80% energy loss
source_rep=50; %rep rate of source in Hz
num_frame=1; %number of frames to look at
tf=num_frame*1e6*1/source_rep ;% time frame to look at us
close all; %closes all previous figure windows


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this section calculates all possible frequencies of the last chopper for
% each frequency of the pulse removal chopper. This is not for the program
% but for information only

pulse_rem=1;% this shows which chopper is the pulse removal chopper
chop1_max=150 ;%max frequency of pulse removal chopper
chop2_max=300 ;%max freq of final chopper
freq_max1=nslot(pulse_rem)*chop1_max;
freq_max2=nslot(N)*chop2_max;
freq_min1=source_rep/nslot(pulse_rem);
freq_min2=source_rep/nslot(N);
figure
axis([0 chop1_max 0 chop2_max]);
xlabel('Pulse removal chopper freq (Hz)');
ylabel('Final chopper freq (HZ)');
title('Allowed chopper frequencies');
hold on
count=1;

for fr1=freq_min1:freq_min1:chop1_max;  %actual chopper frequencies
   % for fr2=freq_min2:freq_min2:chop2_max; 
    for fr2=10:10:chop2_max; 
    n=(dist(N)*fr2*nslot(N))/(dist(pulse_rem)*fr1*nslot(pulse_rem)); %if n is an integer then puse is transmitted
 %n=(dist(pulse_rem)*fr1)/(dist(N)*fr2);
        if (round(n)~=0)   
       df1=(fr1)*(abs(n-round(n))/n); %error in frequency of chopper
       dd1=2*3.1415*radius(pulse_rem)*df1*dist(pulse_rem)*252.6e-6 ;%error in chopper window position(mm) per angstrom
        if (dd1<5) 
            if (fr1)>15
            plot(fr1,fr2,'rs')
            text(fr1,fr2,num2str(fr2),'Fontsize',6)
            diff=(abs(n-round(n))/n)*100;
           % text((fr1/nslot(pulse_rem)),(fr2/nslot(N))-6,num2str(dd1,3),'Fontsize',6)
            text(fr1,fr2-5,num2str(fr1),'Fontsize',6)
            count=count+1;
            end
        end    
        
        end
    end
end









%for fr1=30:freq_min1:freq_max1;
 %   for fr2=30:freq_min2:freq_max2; 
  %  n=(dist(N)*fr2)/(dist(pulse_rem)*fr1); %if n is an integer then puse is transmitted
 %n=(dist(pulse_rem)*fr1)/(dist(N)*fr2);
   %     if (round(n)~=0)   
    %   df1=(fr1/nslot(pulse_rem))*(abs(n-round(n))/n); %error in frequency of chopper
     %  dd1=2*3.1415*radius(pulse_rem)*df1*dist(pulse_rem)*252.6e-6 ;%error in chopper window position(mm) per angstrom
      %  if (dd1<5) 
       %     if (fr1/nslot(pulse_rem))>15
        %    plot((fr1/nslot(pulse_rem)),(fr2/nslot(N)),'rs')
%            text((fr1/nslot(pulse_rem)),(fr2/nslot(N)),num2str(fr2/nslot(N)),'Fontsize',6)
 %           diff=(abs(n-round(n))/n)*100;
           % text((fr1/nslot(pulse_rem)),(fr2/nslot(N))-6,num2str(dd1,3),'Fontsize',6)
  %          text((fr1/nslot(pulse_rem))+10,(fr2/nslot(N)),num2str(fr1/nslot(pulse_rem)),'Fontsize',6)
   %         count=count+1;
    %        end
     %   end    
        
      %  end
%    end
%end
count;
hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


freq_tot=freq.*nslot;
freq_max=max(freq_tot);
chop_times=zeros(round(freq_max*tf*1e-6)+1,3,N); %This will be used to store the chopper opening and closing times for all N choppers


tot_dist=dist(N)+samp_det+chop_samp; %moderator to detector distance
 figure
hold on
axis([0 tf 0 tot_dist]);
 xlabel('Time (us)');
 ylabel('Distance from moderator');
 title('Multiple Chopper system');

for loop=1:N
    x=[0 tf];
    y=[dist(loop) dist(loop)]; 
    plot(x,y,'-k')  % plots the positions of choppers
end

for loop=1:N
    lam=lf(loop)
t_open=252.6*lam*dist(loop); %opening time of chopper such that it is open for wavelength lam
t_wid=1e6*(width(loop)+g_wid(loop))/(2*pi*radius(loop)*counter(loop)*freq(loop));   % full opening time of choper
t_wid2=1e6*(width(loop)-g_wid(loop))/(2*pi*radius(loop)*counter(loop)*freq(loop));   % opening time of choper for 100% transmition
x=[ (t_open-t_wid/2)  (t_open+t_wid/2)];
x2=[ (t_open-t_wid2/2)  (t_open+t_wid2/2)];
 y=[dist(loop) dist(loop)];
t_shift=1e6/(nslot(loop)*freq(loop)); %time difference from one window to the next in the chopper
num=fix(t_open/t_shift);
x=x-num*t_shift;
x2=x2-num*t_shift;
if (x(1)<0) 
    if loop~=1
         x=x+t_shift;
         x2=x2+t_shift;
     else
        x(1)=100;
        x2(1)=100;
    end
end
count=0;
    while (x > 0) & (x < tf)
        count=count+1;
        plot(x,y,'-y')
        plot(x2,y,'-w')
        chop_times(count,1:2,loop)=x;  
        x=x+t_shift; 
        x2=x2+t_shift; 
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the next section scans to see what neutrons can pass through the whole
% chopper system. It starts with the final chopper nearest the sample and
% scans every opening

count=1;
while chop_times(count,1,N)~=0  ; %loops around every opening in last chopper

    %count=16;
    lam_max=chop_times(count,1,N)/(252.6*dist(N)); %min wavelength to pass through fiinal chopper
    lam_min=(chop_times(count,2,N)-tm)/(252.6*dist(N)); % max wavelength to pass through final chopper
    chop_times(:,3,:)=0;
    %determine number of openings each chopper has within the final chopper
    %chopper opening time/lambda space
 
    for chop=1:N    %loop around all choppers up to last chopper
         tmax=(252.6*lam_min*dist(chop))+tm;
         tmin=252.6*lam_max*dist(chop);
        count2=1;
            while chop_times(count2,1,chop)~=0  %loop around all openings in chopper
                if     ((chop_times(count2,1,chop)>=tmin) & (chop_times(count2,1,chop)<=tmax)) | ((chop_times(count2,2,chop)>=tmin) & (chop_times(count2,2,chop)<=tmax)) | ((chop_times(count2,1,chop)<=tmin) & (chop_times(count2,2,chop)>=tmax))               
                chop_times(count2,3,chop)=1;  
                end 
                count2=count2+1;
            end  
    end
    
    
    
    a=sum(chop_times(:,3,1:N-1));
    if (sum(find(a==0))==0) %if true this means all choppers have at least one opening within final chopper time frame
    
            index=find(chop_times(:,3,N)); %get index of final chopper opening
            %next line gets 2 lines of gradient m and intercept c defining chopper opening in wavelength time space 
            [m,c]=chop_line(chop_times(index(1),1,N),chop_times(index(1),2,N),dist(N));
            % the 4 points defining the vertices of the slanted rectangle in lam,t space are given by
            x(1)=0; x(2)=0.1 ;x(3)=tm ;x(4)=tm+0.1;
            y(1)=c(1); y(2)=c(2) ;y(3)=(m*tm)+c(2); y(4)=(m*tm)+c(1); 
            % plot(x,y,'rs')
            % this line converts the points x,y into lines joining the points
            lines_last=point_line(x,y);
            
            % find out the number of openings each chopper has
            for chop=1:N
                num_open(chop)=length(find(chop_times(:,3,chop)));
            end
            % this returns a matrix of all possible combinations of chopper openings, i.e 1st opening of chop1 with 3rd of chop2 and 2 of chop3 ect
            [all_permutations,tot]=permut(num_open,N); % tot is the total number of permutations
            flag=0;
            
            
        for loop=1:tot  % loop around all permutations of openings up to last chopper opening
            loop;
            flag=0;
            lines_old=lines_last;  % resets lam,t space to that of last chopper
            for chop=1:N %N-1
                 open=find(chop_times(:,3,chop))  ;
                 index=all_permutations(loop,chop);
                 
            [m,c]=chop_line(chop_times(open(index),1,chop),chop_times(open(index),2,chop),dist(chop));
            % determine intersecting points of these choppers lam,t space with last chopper (lines_old)
            [xi,yi]=intersect(lines_old,m,c);
            %plot(xi,yi,'rs')
                if sum(xi)==0;
                     flag=1;
                     break % i.e no neutrons can ge through these choppers
          
                end
                 lines_old=point_line(xi,yi); % new area of lam, time space is given by common area intersection  
             end
             if (flag~=1)
                 all_permutations(loop,:);  %this is the successful permutation
                 draw_ray(xi,yi,samp_det,chop_samp,dist(N),temp,frac_energy,freq1,source_rep,num_frame);
                 flag=0;
             end    
        end
        
        
    end
   count=count+1; 
end   
    hold off
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [pe,tot]=permut(num_open,N)
len=length(num_open);
tot=1;
for loop=1:len
    tot=tot*num_open(loop) ; %this is the total number of permutations
end
for loop=1:N
    count=1;
    while count < tot+1
        for loop2=1:num_open(loop)
        pe(count,loop)=loop2;
        count=count+1;
        end
    end
end    
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xi,yi]=intersect(lines_old,m,c);
xi=[];yi=[];
len1=size(lines_old,1);
count=1;
for loop1=1:len1
    m1=lines_old(loop1,1);
    c1=lines_old(loop1,2);
    x1_old=lines_old(loop1,3);
    x2_old=lines_old(loop1,4);
    for line_new=1:2
        m2=m;
        c2=c(line_new);
        
        xint=(c2-c1)/(m1-m2);
        yint=(m1*xint)+c1;
        if ((xint > x1_old) & (xint < x2_old)) | ((xint < x1_old) & (xint > x2_old)) ;%i.e intercept is between x points of old line
            xi(count)=xint;
            yi(count)=yint;
            count=count+1;
        end  
    end
  
end
for line=1:len1
    for point=1:2
     x=lines_old(line,(2+point));
     y=lines_old(line,(4+point));
     ylower=(m*x)+c(1);
     yupper=(m*x)+c(2);
        if (y <= yupper) & (y >= ylower) & ((sum(find(xi==x))+sum(find(yi==y)))==0)
         xi(count)=x;yi(count)=y;
         count=count+1;
       end
   end   
end


  
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [line]=point_line(x,y);%determines the lines outling the points x,y in time lambda space
len=length(x);
line=zeros(len,6);
num=1;
for point1=1:2:len %pick first point
   
    for point2=1:len %pick 2nd point to make a line with
        if point2~=point1 
            ax=x(point1)-x(point2);ay=y(point1)-y(point2);
            count=1;
            flag=0;
          for point3=1:len % scroll through rest of points
                if ((point3~=point2) & (point3~=point1))
          
                    bx=x(point1)-x(point3);by=y(point1)-y(point3);
                    up_down=sign((ax*by)-(ay*bx));
                    
                    if count==1;side=up_down;end
                    if up_down~=side
                        flag=1;
                        break
                    end
                    count=count+1;
                end
          end
          if flag==0; %means all other points are on one side of the line created by point1 and point2
              if ax==0
                  ax=1e-6;
              end
              grad=ay/ax;
              inter=y(point1)-grad.*x(point1);
              line(num,:)=[grad inter x(point1) x(point2) y(point1) y(point2)];
              num=num+1;
          end
       end
   end
end
%plot(x,y,'rs')

%len=size(line,1);
%for loop=1:len
 %   if (line(loop,3)<line(loop,4))
  %      step=abs(line(loop,3)-line(loop,4))/100;
  % else
  %      step=-abs(line(loop,3)-line(loop,4))/100;
  % end
  %  count=1;
  %for x_co=line(loop,3):step:line(loop,4)
   %   xp(count)=x_co;
    %  yp(count)=(line(loop,1)*x_co)+line(loop,2);
     % count=count+1;
     %end   
  
  % plot(xp,yp);
  %end
  return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [m,c]=chop_line(t1,t2,r);%function determines the two lines outlining the shape in time lambda space by the chopper 
cmin=t1/(252.6*r);
c(1)=cmin;
cmax=t2/(252.6*r);
c(2)=cmax;
m=-1/(252.6*r);
count=1;
%for x_co=1:20:4000
 %     xp(count)=x_co;
  %    y1(count)=(m*x_co)+c(1);
   %   y2(count)=(m*x_co)+c(2);
    %  count=count+1;
    %end   
  
 % plot(xp,y1); 
  %plot(xp,y2);
   %  hold on
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function draw_ray(x,y,samp_det,chop_samp,dist_N,temp,energy_frac,freq1,source_rep,num_frame)
area=polyarea(x,y);
energy_gain=0.86*temp;
len=length(x);
x_tot=0;y_tot=0;
for loop=1:len
    x_tot=x_tot+x(loop); y_tot=y_tot+y(loop);
end
x_ave=x_tot/len;y_ave=y_tot/len;
string=num2str(y_ave,3);
energy=81.81/(y_ave*y_ave);
string2=num2str(energy,3);
if (x_ave<700 && energy <250)
    fprintf(1,'Energy transmitted by choppers is %s mev \n',string2);
end
energy_max=energy+energy_gain;
y_min=sqrt(81.81/energy_max);
new_energy=(1-energy_frac)*energy;
y_max=sqrt(81.81/(new_energy));
t(1)=x_ave;dist(1)=0;

if (energy <250)
tshift=1e6/source_rep;
for loop=1:num_frame+1
    
    if (t(1)-(loop*tshift))>1500
        t(2)=(252.6*y_ave*(dist_N+chop_samp+samp_det))+t(1);dist(2)=dist_N+chop_samp+samp_det;
        if (loop >1)
        t=t+tshift;
        end
        plot(t,dist,'-r')
    else
        t(2)=(252.6*y_ave*(dist_N+chop_samp))+t(1);dist(2)=dist_N+chop_samp;
        t(3)=t(2)+(252.6*y_min*(samp_det));dist(3)=dist_N+chop_samp+samp_det;
        if (loop > 1)
        t=t+tshift;
        end
        plot(t,dist,'-b')
        t2(1)=t(2)+(252.6*y_max*(samp_det));dist2(1)=dist(2)+samp_det;
        t2(2)=t(2);dist2(2)=dist(2);
        plot(t2,dist2,'-b')
       % text(t(2)+2000,dist_N+chop_samp,string,'Fontsize',6);
        text(t(2),dist_N+chop_samp-1,string2,'Fontsize',6);
    end
    
    
    
end

end

return


