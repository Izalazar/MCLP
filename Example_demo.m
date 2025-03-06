clear all;
close all;
clc;

data=load('Male_A_110Hz__D_0.mat');%50ms voice segmet , extracted from Repository II of the OPENGLOT database.  
s=data.data.s;
vg=data.data.vg;
vg=vg/max(abs(vg(30:end-30)));

alpha=0.99;% Pre-emphasis filter
P=11;% Vocal tract filter order
epsilon1=0.0001;% Thresholds
epsilon2=0.00001;
[a_vt1,h_e]=MCLP(s',P,epsilon1,epsilon2,alpha);

a_vt= remove_spurious_poles(a_vt1);
vg_est=filter(a_vt,1,s); %Inverse filtering
vg_est=vg_est/max(abs(vg_est(30:end-30)));


plot(vg);hold on;plot((vg_est))
plot(h_e)
legend('$v_g$','$\hat{v}_g$','$h_e$','Interpreter','latex')
