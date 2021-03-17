function [a_ref,bb_ref,epsilon,fval]=get_a_ref(Rrs_shade_ref,wl_ref,radius,theta_w,s1,s2,s3,s4,Y,bbp750,bbw_ref)
%############################################################# 
% Xiaolong Yu on June 22, 2020, calculate the bbp at ref(750) assuming
% a(ref) = aw(ref)

%#################################
%% input parameters
% aw          : pure water absorption, adapted from Lee et al., 2015
% bbw         : pure water backscattering, adapted from Morel et al., 2001,
%               used in Hydrolight
% theta       : solar zenith above sea surface 
% radius      : radius of the cone attached radiance radiometer, in meters
% Rrs_shade   : raw Rrs data from SBA (Rrs with shadow errors)
% bbp_ref     : bbp at reference wavelength (e.g., 750 nm)
%% other paramters
% theta_w     : solar zenith just below sea surface
% s1 s2 s3 s4 : parameters for calculating the attenuation coefficient and epsilon

%% output parameters
% a_ref       : a at the ref wavelength (e.g., 400 or 440 nm)
% Rrs_corr    : Shadow-corrected Rrs ('true' Rrs) 

%% start  
rrs = Rrs_shade_ref./(0.52+1.7*Rrs_shade_ref);   
g0 = 0.089;
g1 = 0.125;
u_ref  = (-g0 + (g0^2 + 4*g1*rrs).^0.5)/(2*g1); 

bbp_ref = bbp750*(750./wl_ref)^Y;
% bbw_ref=0.5*h2o_iops(wl_ref,'b');
% bbw_ref=bbw(wl_ref);
bb_ref=bbp_ref+bbw_ref;
a_ref_temp = bb_ref.*(1-u_ref)./u_ref;

x0=[a_ref_temp];
options=optimset('largescale','on','display','iter','tolx',1e-8,'tolfun',1e-12,'Algorithm','active-set','MaxFunEvals',10000,'MaxIter',10000);

LB=[a_ref_temp*0.1];            %lower boundary
UB=[a_ref_temp*10];             %upper boundary  

[xx,fval]=fmincon(@(x0) optimization2(Rrs_shade_ref,x0,radius,theta_w,s1,s2,s3,s4,bb_ref,bbw_ref),x0,[],[],[],[],LB,UB,[],options);

a_ref=xx; 

epsilon=get_epsilon(a_ref,bb_ref,radius,theta_w,s1,s2,s3,s4); 

end








