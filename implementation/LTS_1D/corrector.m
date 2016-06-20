function U=corrector(U,n,dt,cellsize,level)

global F levelc

for j=1:n,
    if levelc(1,j) <= level,
        for l=1:2,
            U(1,j,l)=U(1,j,l)-dt(1,j)*(F(1,j+1,l)-F(1,j,l))/cellsize;
        end
    end
end
