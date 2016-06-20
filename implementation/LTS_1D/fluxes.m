function fluxes(U,n,level)

global F levelf lambdaf
%Compute fluxes in x-direction
%levelf
sn=0; cn=1;
vl=0;
vr=0;
for j=2:n,
    if levelf(1,j) <= level,
        hl=U(1,j-1,1);
        hr=U(1,j,1);
        ul=U(1,j-1,2)/hl;
        ur=U(1,j,2)/hr;
        [F(1,j,:), lambdaf(1,j)]=solver(hl,hr,ul,ur,vl,vr,sn,cn);
    end
end
%Enforce wall boundary on left side of wall
if levelf(1,1) <= level,
hr=U(1,1,1);
ur=U(1,1,2)/hr;
[F(1,1,:), lambdaf(1,1)]=solver(hr,hr,-ur,ur,vr,vr,sn,cn);
end

%Enforce wall boundary on right side of wall
if levelf(1,n+1) <= level,
hl=U(1,n,1);
ul=U(1,n,2)/hl;
[F(1,n+1,:), lambdaf(1,n+1)]=solver(hl,hl,ul,-ul,vl,vl,sn,cn);
end
