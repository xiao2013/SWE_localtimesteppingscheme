function [F, amax]=solver(hl,hr,ul,ur,vl,vr,sn,cn)

global grav

%Compute Roe averages
duml=hl^0.5;
dumr=hr^0.5;
cl=(grav*hl)^0.5;
cr=(grav*hr)^0.5;
hhat=duml*dumr;
uhat=(duml*ul + dumr*ur)/(duml+dumr);
vhat=(duml*vl + dumr*vr)/(duml+dumr);
chat=(0.5*grav*(hl+hr))^0.5;
uperp=uhat*cn+vhat*sn;
dh=hr-hl;
du=ur-ul;
dv=vr-vl;
dupar=-du*sn+dv*cn;
duperp=du*cn+dv*sn;
dW=[0.5*(dh-hhat*duperp/chat); hhat*dupar; 0.5*(dh+hhat*duperp/chat)];

uperpl=ul*cn+vl*sn;
uperpr=ur*cn+vr*sn;
al1=uperpl-cl;
al3=uperpl+cl;
ar1=uperpr-cr;
ar3=uperpr+cr;
R=[1 0 1;
    uhat-chat*cn -sn uhat+chat*cn;
    vhat-chat*sn cn vhat+chat*sn];
da1=max([0 2*(ar1-al1)]);
da3=max([0 2*(ar3-al3)]);
a1=abs(uperp-chat);
a2=abs(uperp);
a3=abs(uperp+chat);

%Critical flow fix
if a1 < da1,
     a1=0.5*(a1*a1/da1+da1);
end
if a3 < da3,
     a3=0.5*(a3*a3/da3+da3);
end

%Compute interface flux
A=diag([a1 a2 a3]);
FL=[uperpl*hl; ul*uperpl*hl + 0.5*grav*hl*hl*cn; vl*uperpl*hl + 0.5*grav*hl*hl*sn];
FR=[uperpr*hr; ur*uperpr*hr + 0.5*grav*hr*hr*cn; vr*uperpr*hr + 0.5*grav*hr*hr*sn];
%ur
%hr
%uperpr
%ur*uperpr*hr
%cn
%FL
%FR
%R*A*dW
%R
%A
%dW
F=0.5*(FL + FR - R*A*dW);
amax=chat+abs(uperp);




