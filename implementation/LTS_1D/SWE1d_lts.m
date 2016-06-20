clc
clear
close all
global grav F lambdaf

grav=9.806;
% aviobj = avifile('swe1D_LTS.avi','compression','None');
cellsize=1; %Cell size (m)
n=50; %number of cells across box
L=n*cellsize; %box size (m)

x=cellsize/2:cellsize:L-cellsize/2; %x coordinate of cell center
y=cellsize/2:cellsize:L-cellsize/2; %y coordinate of cell center

%LTS_levels=1;
% for k = 1:10
k = 1;
% for LTS_levels=1:2:3
LTS_levels = 3
nt(k)=200; %number of time steps
% nt(k) = nt(k) + (k-1)*50;
% nt = 50;
ntplot=1; %plotting interval
dt0=0.1; %time step (s)

%Initial condition
t=0; %start time
h=0.1*ones(1,n); %depth
h(1,1:10)=1; %higher water level in corner

u=zeros(1,n);

F=zeros(1,n+1,3); %fluxes in the $x$ direction
lambdaf=zeros(1,n+1); %maximum wave speeds in $x$ direction
U(1,:,1)=h;
U(1,:,2)=u.*h;

level_B=LTS_levels;
for tstep=1:nt(k),
    
    level_A=level_B;
    %Assign level_B
    if(mod(tstep-1,2) == 0)
        level_B=1;
    else
        for i=1:LTS_levels
            idum=2^(i-1);
            if (mod(tstep,idum) == 0),
                level_B=i;
            end
        end
    end

    %Compute fluxes
    fluxes(U,n,level_A);
    %Assign LTS levels
    if (level_A == LTS_levels),
        [dt dt0 dt_real crmax]=timestep(LTS_levels,cellsize,n,dt0);
    end
    %Update solution to next time level
    U = corrector(U,n,dt,cellsize,level_B);

    fprintf('Time Step = %d, Courant Number = %g, level_B= %d \n',tstep,crmax,level_B)

    if (level_B == LTS_levels),
        subplot(2,1,1)
        U_d=U(1,:,1);        
        bar(x,U_d,1)
        xlabel('x (m)'); ylabel('h (m)')
        title('1D shallow water equation - dam break numerical simulation')
        axis([0 L 0 1])
        subplot(2,1,2)
        bar(x,dt(1,:,1),1)
        xlabel('x - domain x'); ylabel('LTS timestepping')
        axis([0 L 0 0.8])

%         subplot(3,1,3)
%         bar(x,dt_real(1,:,1),1)
%         xlabel('x - domain x'); ylabel('necessary timestepping')
%         axis([0 L 0 0.8])
%         aviobj = addframe(aviobj,gcf);
        pause(0.2)
    end
end
% viobj = close(aviobj)









