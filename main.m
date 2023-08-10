% %% test of orth
% locald=10;
% Nsites=10;
% res={};
% op={};
% %def local n oper
% n = zeros(locald,locald);
% for i=1:locald-1
%     n(i+1,i+1)=i;
% end
% %
% res{1}=ones(1,1);
% for u = 1:Nsites
%     op{u}=eye(locald);
% end
% for u = 1:Nsites  
%     res{u+1} = ncon({res{u},op{u},A{u},conj(A{u})},{[1,2],[3,4],[1,3,-1],[2,4,-2]});
% end
% result=res{Nsites+1};

%%Test if the GS energy is correct;
clear all;
locald=10;
Nsites=10;
chi = 10;
ty=0;
ny=1;
phi=1;
U=1;
tx=0.25;
mu=0.25;

x=1.0;
ty = sin(pi/2*x);
U=cos(pi/2*x);
[A,sWeight,B,E] = BH(U,tx,ty,mu,ny,Nsites,locald,phi,chi);
%% test sample: get entropy,totenergy,sumnum
clear all;
locald=10;
Nsites=10;
simnum=100;
scaleForSimArea=0.5;
chi = 10;
ty=0;
ny=1;
phi=1;
U=1;
nbetween=9;
TotEnerAll={};
sumnumAll={};
entAll={};

for ntry = 1:1:nbetween
    theta = pi/(2*(nbetween-1))*(ntry-1);
    ty=sin(theta);
    U=cos(theta);
    
    n = zeros(locald,locald);
    % These three things are what we wanna to calculate
    TotEner=zeros(simnum,simnum);
    sumnum=zeros(simnum,simnum);
    ent=zeros(simnum,simnum);
    for i=1:locald-1
        n(i+1,i+1)=i;
    end
    
    for y=1:simnum
      mu=0+y/simnum*scaleForSimArea;
      for x=1:simnum
        tx=0+x/simnum*scaleForSimArea;
        [A,sWeight,B,E] = BH(U,tx,ty,mu,ny,Nsites,locald,phi,chi);
        %totenergy
        TotEner(x,y)=E;
        %sumnum of particle
        for i=1:Nsites
          res={};
          op={};
          res{1}=ones(1,1);
          for u = 1:Nsites
            op{u}=eye(locald);
          end
          op{i}=n;
          for u = 1:Nsites  
            res{u+1} = ncon({res{u},op{u},A{u},conj(A{u})},{[1,2],[3,4],[1,3,-1],[2,4,-2]});
          end
          sumnum(x,y)=sumnum(x,y)+res{Nsites+1};
        end
        %entanentropy
        entanS=sWeight{Nsites/2};
        for i=1:length(entanS)
          if entanS(i,i)>10e-4
            ent(x,y)=ent(x,y)-entanS(i,i)*log(entanS(i,i));
          end
        end
      end
    end
    TotEnerAll{ntry}=TotEner;
    sumnumAll{ntry}=sumnum;
    entAll{ntry}=ent;
end

%% Plot phase diagram
tgrid=1/simnum*scaleForSimArea:1/simnum*scaleForSimArea:1*scaleForSimArea;
mugrid=1/simnum*scaleForSimArea:1/simnum*scaleForSimArea:1*scaleForSimArea;
[tt,mumu]=meshgrid(tgrid,mugrid);
surf(mumu,tt,TotEnerAll{4});
xlabel('Constant t');
ylabel('Constant mu');
zlabel('total gs energy')
title('total gs energy change with t,mu')
%% Plot2
tgrid=1/simnum*scaleForSimArea:1/simnum*scaleForSimArea:1*scaleForSimArea;
mugrid=1/simnum*scaleForSimArea:1/simnum*scaleForSimArea:1*scaleForSimArea;
[tt,mumu]=meshgrid(tgrid,mugrid);
surf(mumu,tt,entAll{4});
xlabel('Constant t');
ylabel('Constant mu');
zlabel('entanglement entropy')
title('entanglement entropy change with t,mu')
%% Plot3
tgrid=1/simnum*scaleForSimArea:1/simnum*scaleForSimArea:1*scaleForSimArea;
mugrid=1/simnum*scaleForSimArea:1/simnum*scaleForSimArea:1*scaleForSimArea;
[tt,mumu]=meshgrid(tgrid,mugrid);
surf(mumu,tt,sumnumAll{4});
xlabel('Constant t');
ylabel('Constant mu');
zlabel('Total particle number in this chains')
title('Total particle number change with t,mu')

% %% test of superfluid density
% clear all;
% locald=10;
% Nsites=10;
% simnum=50;
% chi=10;
% sfd=zeros(simnum,simnum);
% for y=1:simnum
%     mu=0.2+0.5*y/simnum;
% for x=1:simnum
%     t=0+0.3*x/simnum;
% [A,sWeight,B,E] = BH(1,t,mu,Nsites,locald,0,chi);
% [A2,sWeight2,B2,E2] = BH(1,t,mu,Nsites,locald,0.1,chi);
% %totenergy
% sfd(x,y)=(E2-E)/(0.1*0.1*t);
% end
% end

%% Plot4
tgrid=0.3/simnum:0.3/simnum:0.3;
mugrid=0.2+9*0.5/simnum:0.5/simnum:0.7;
[tt,mumu]=meshgrid(tgrid,mugrid);
sfdp=sfd(9:1:simnum,1:1:simnum);
surf(tt,mumu,sfdp);
xlabel('Constant t');
ylabel('Constant mu');
zlabel('superfluid density for this chains')
title('superfluid density changed with t,mu')

%% Energy Plot1
clear all
Elist=zeros(10,1);
for n=1:10
[A,sWeight,B,E] = BH(1,0.1*n,0.4,10,10,0,10);
Elist(n)=E;
end
%% Energy Plot2
x=0.1:0.1:1;
hold on
plot(x,Elist)
plot(x,Elist1)
hold off
%% test of locnum
% locald=10;
% Nsites=10;
% %def local n oper
% n = zeros(locald,locald);
% for i=1:locald-1
%     n(i+1,i+1)=i;
% end
% %sum all
% result=0;
% i=9;
% res={};
% op={};
% res{1}=ones(1,1);
% for u = 1:Nsites
%     op{u}=eye(locald);
% end
% op{i}=n;
% for u = 1:Nsites  
%     res{u+1} = ncon({res{u},op{u},A{u},conj(A{u})},{[1,2],[3,4],[1,3,-1],[2,4,-2]});
% end
% result=result+res{Nsites+1};
% %% Compute 2-site reduced density matrices, local energy density
% rhotwo = {};
% hamloc = reshape(kron(sP,sM) + kron(sM,sP),[2,2,2,2]);
% for k = 1:Nsites-1
%     rhotwo{k} = ncon({A{k},conj(A{k}),A{k+1},conj(A{k+1}),sWeight{k+2},sWeight{k+2}},...
%         {[1,-3,2],[1,-1,3],[2,-4,4],[3,-2,5],[4,6],[5,6]});
%     Enloc(k) = ncon({hamloc,rhotwo{k}},{[1:4],[1:4]});
% end