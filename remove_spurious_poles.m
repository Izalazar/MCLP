function [a_vt] = remove_spurious_poles(a_vt)
%Remove spurious poles of the vocal tract filter. 

a_vt=a_vt./sum(a_vt);%Normalization of the filter

[z,p,gain] = tf2zpk(1,a_vt);

dist=abs(p);
p(dist>1)=1./p(dist>1);
aux=real(p)>0 & abs(p)<0.8;
p=p(aux==0); % non-spurious poles selection. 

[~,a_vt]= zp2tf(z,p,gain);
a_vt=a_vt./sum(a_vt);

end