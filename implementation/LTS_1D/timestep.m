function [dt dt0 dt_real crmax]=timestep(LTS_levels,cellsize,n,dt0)

global levelc levelf lambdaf

crmax=0;
dt=dt0*ones(1,n);
dt_real=dt0*ones(1,n);
cr=zeros(1,n);
levelc=LTS_levels*ones(1,n);
levelf=LTS_levels*ones(1,n+1);
if (LTS_levels > 1),
    %Compute Courant number and assign level id to cells    
    for j=1:n,
        crdum=dt0*max([lambdaf(1,j) lambdaf(1,j+1)])/cellsize;
        dt_real(1,j)=dt0*0.8/crdum;

        %crdum
        for i=1:LTS_levels-1,
            crtest=0.8*2^(-i);
            %crtest
            if crdum >= crtest,
                levelc(1,j)=i;
                dt(1,j)=dt0*2^(i-1);
                cr(1,j)=crdum*dt(1,j)/dt0;
                break
            end
            levelc(1,j)=LTS_levels;
            dt(1,j)=dt0*2^(LTS_levels-1);
            cr(1,j)=crdum*dt(1,j)/dt0;
        end
        if cr(1,j) > crmax,
            crmax=cr(1,j);
        end
    end
    
    %Assign a level id to edges (faces) of mesh

    for j=2:n,
        levelf(1,j)=min([levelc(1,j-1) levelc(1,j)]);
    end
    levelf(1,1)=levelc(1,1);
    levelf(1,n+1)=levelc(1,n);
    
    %Create a layer of buffer cells
    for j=1:n,
        levelc(1,j)=min([levelf(1,j) levelf(1,j+1)]);
        dt(1,j)=dt0*2^(levelc(1,j)-1);
    end
else %Compute maximum courant number if using global time stepping
    for j=1:n,
        cr(1,j)=dt0*(max([lambdaf(1,j) lambdaf(1,j+1)]))/cellsize;
        if cr(1,j) > crmax,
            crmax=cr(1,j);
            %umax=crmax/dt0;
        end        
    end
%    dt0=0.8/umax;
end

           