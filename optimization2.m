function f=optimization2(Rrs_shade_ref0,a_ref,radius,theta_w,s1,s2,s3,s4,bb_ref,bbw_ref)

%% optimization process to get bb at reference band
bbp_ref=bb_ref-bbw_ref;

rrsw=0.113.*bbw_ref./(a_ref+bb_ref); 

u_temp=bbp_ref/(a_ref+bb_ref); 

rrsdp=0.2*(1-0.63*exp(-2.448*u_temp))*u_temp; 
rrs=rrsdp+rrsw;
Rrs_ref=0.52.*rrs./(1-1.7*rrs);  

epsilon_ref=get_epsilon(a_ref,bb_ref,radius,theta_w,s1,s2,s3,s4); 
Rrs_shade_ref1=Rrs_ref.*(1-epsilon_ref); 
diff=abs(Rrs_shade_ref1-Rrs_shade_ref0)./Rrs_shade_ref0; 

f=diff; % calculated error between simulate Rrs_shade and measured Rrs_shade at ref

end
