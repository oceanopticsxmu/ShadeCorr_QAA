function [Rrs_corr, epsilon,a,bb]=shade_corr(wl,Rrs_shade,radius,theta)
%############################################################# 
% Created by Xiaolong Yu on June 22, 2020.

%#############################################################
% This script is designed for shading error correction of SBA measured Rrs
% not applicable for shallow water.
% Reference: Yu, X., Lee, Z., Shang, Z., Lin, H., & Lin, G. (2021). A simple and robust 
           % shade correction scheme for remote sensing reflectance obtained by the 
           % skylight-blocked approach. Optics Express, 29, 470-486.
%#################################
%% input parameters 
% wl          : wavelength (nm)
% Rrs_shade   : obtained Rrs from SBA (i.e., Rrs with shade error, match wl)
% radius      : radius of the cone attached radiance radiometer, in meters
% theta       : solar zenith above sea surface 

%% output parameters
% epsilon     : shadow error
% Rrs_corr    : Shadow-corrected Rrs ('true' Rrs)

%% call functions
% aw_lee2015;    % get aw, pure water absorption, adapted from Lee et al., 2015
% h2o_iops;      % get bbw, pure water backscattering, adapted from Morel et al., 2001,  used in Hydrolight
% get_bbp_ref    % get bbp at 750 nm using optimization, assume a750=aw750;

%% start 

wl_ref0=750; % start wavlenegth to initial the shade correction scheme.

 % ##### index of reference wavelength at 750 nm ##### 
[~,id440] = min(abs(wl-440)); 
[~,id550] = min(abs(wl-550)); 
[~,id_ref0] = min(abs(wl-wl_ref0));

if wl(end) < 750
    disp('750 nm is required');
    return;
end

if wl(end) > 750
    wl=wl(1:id_ref0);
    Rrs_shade=Rrs_shade(1:id_ref0);
end




  n=1.34;
  theta_w=asin(sin(theta*pi/180)/n); % solar zenith just beneath sea surface 
  
% s1 s2 s3 s4 : parameters for calculating the attenuation coefficient and epsilon
  s1=1.15+3.15*sin(theta_w);
  s2=-1.57;
  s3=5.62*sin(theta_w)-0.23;
  s4=-0.5;  
  
aw=aw_lee2015(wl);  % pure seawater absorption coefficients, Lee et al., 2015
bbw=0.5*h2o_iops(wl,'b');  % pure seawater backscattering coefficients

aw_ref0=aw(id_ref0);
bbw_ref0=bbw(id_ref0);

%% get bbp at ref 750 nm
  Rrs_shade_ref=Rrs_shade(id_ref0);   
  a_ref0=aw_ref0;      
  [bbp_ref0,eps,~] = get_bbp_ref(Rrs_shade_ref,radius,theta_w,s1,s2,s3,s4,a_ref0,bbw_ref0);     
  epsilon(id_ref0)=eps; 
    
%% get a at given wavelength
  rrs_shade = Rrs_shade./(0.52+1.7*Rrs_shade); 
  
  for j=1:length(wl)-1
        wl_ref=wl(j);
        [~,id_ref]=min(abs(wl-wl_ref));
        bbw_ref=bbw(id_ref);
        Rrs_shade_ref2=Rrs_shade(id_ref);           
        Y= 2.0*(1-1.2.*exp(-0.9*rrs_shade(id440)./rrs_shade(id550)));  
        [a_ref,bb_ref,eps,~]=get_a_ref(Rrs_shade_ref2,wl_ref,radius,theta_w,s1,s2,s3,s4,Y,bbp_ref0,bbw_ref);          
        epsilon(j)=eps; 
        a_temp(j)=a_ref;
        bb_temp(j)=bb_ref;
  end
    a=[a_temp,aw_ref0];
    bb=[bb_temp,bbp_ref0+bbw_ref0];
    Rrs_corr=Rrs_shade./(1-epsilon');
end








