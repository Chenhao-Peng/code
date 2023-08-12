function [A,sWeight,B,Efinal]=BH(setU,settx,setty,setmu,setny,setN,setlocd,setphi,chi)

%%%%% Set simulation options
% chi = 10; % maximum bond dimension
Nsites = setN; % number of lattice sites
OPTS.numsweeps = 15; % number of DMRG sweeps
OPTS.display = 2; % level of output display   ??????
OPTS.updateon = 1; % update MPS tensors       
OPTS.maxit = 2; % iterations of Lanczos method
OPTS.krydim = 4; % dimension of Krylov subspace  

phi=setphi;
mu=setmu;
ny=setny;
U=setU;
tx=settx;
ty = setty;
chid=setlocd;
b = zeros(chid,chid);
n = zeros(chid,chid);
for i=1:chid-1
    b(i,i+1)=sqrt(i);
    n(i+1,i+1)=i;
end

omg=U/2*n*n-(U/2+mu)*n;
M = zeros(4,4,chid,chid);
M(1,1,:,:) = eye(chid); M(4,2,:,:) = -tx*b;
M(2,1,:,:) = b'; M(4,3,:,:) = -tx*b';
M(3,1,:,:) = b; M(4,4,:,:) = eye(chid);
M(4,1,:,:) = omg;

%define M1 operator, to get the cos term.
M1 = zeros(4,4,chid,chid);
M1(4,1,:,:) = -2*ty*n;


%For Periodic boundary
 ML = zeros(4,4,1,1); %left MPO boundary
 ML(:,:,1,1)=eye(4);
 MR = zeros(4,4,1,1); %right MPO boundary
 MR(:,:,1,1)=eye(4);
 T=zeros(4,4,chid,chid);
T(1,2,:,:) = -tx*b;
T(1,3,:,:) = -tx*b';
T(1,4,:,:) = eye(chid);
MT=ncon({M+M1*cos(2*pi*ny/Nsites-2*pi*phi*(Nsites)),T},{[-1,1,2,-4],[1,-2,-3,2]});

%For open boundary 
%open boundary condition
% ML = reshape([0,0,0,1],[4,1,1]);
% MR = reshape([1,0,0,0],[4,1,1]);
% MT=M;

% Initialize MPS tensors
A = {};
A{1} = rand(1,chid,min(chi,chid));
for k = 2:Nsites
    A{k} = rand(size(A{k-1},3),chid,min(min(chi,size(A{k-1},3)*chid),chid^(Nsites-k)));
end
[A,sWeight,B,Ekeep2] = doDMRG_MPO(A,ML,M,MT,MR,M1,phi,ny,chi,OPTS);
Efinal=Ekeep2(end);
end
