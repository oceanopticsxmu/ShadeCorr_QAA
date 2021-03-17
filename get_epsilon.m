function [epsilon]=get_epsilon(a,bb,radius,theta_w,s1,s2,s3,s4)
% calculate shading error, epsilon 
Kappa=s1.*a.*exp(s2.*bb)+s3.*bb.*exp(s4.*a);     
epsilon=1-exp(-Kappa.*radius./tan(theta_w)); 
     
end