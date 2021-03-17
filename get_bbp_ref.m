function [bbp_ref,epsilon,fval]=get_bbp_ref(Rrs_shade_ref,radius,theta_w,s1,s2,s3,s4,a_ref,bbw_ref)
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
%% other paramters
% theta_w     : solar zenith just below sea surface
% s1 s2 s3 s4 : parameters for calculating the attenuation coefficient and epsilon

%% output parameters
% bbp_ref     : bbp at the ref wavelength (ref=750 nm)
% Rrs_corr    : Shadow-corrected Rrs ('true' Rrs) 

%% start  
rrs = Rrs_shade_ref./(0.52+1.7*Rrs_shade_ref);  

g0 = 0.089;
g1 = 0.125;
u_ref  = (-g0 + (g0^2 + 4*g1*rrs).^0.5)/(2*g1);
bbp_ref_temp = u_ref.*a_ref/(1-u_ref) - bbw_ref;

x0=[bbp_ref_temp];
options=optimset('largescale','on','display','iter','tolx',1e-8,'tolfun',1e-12,'Algorithm','active-set','MaxFunEvals',10000,'MaxIter',10000);

LB=[bbp_ref_temp*0.1];            %lower boundary
UB=[bbp_ref_temp*10];             %upper boundary  

[xx,fval]=fmincon(@(x0) optimization1(Rrs_shade_ref,x0,radius,theta_w,s1,s2,s3,s4,a_ref,bbw_ref),x0,[],[],[],[],LB,UB,[],options);

bbp_ref=xx;
bb_ref=bbp_ref+bbw_ref;
epsilon=get_epsilon(a_ref,bb_ref,radius,theta_w,s1,s2,s3,s4); 

end








