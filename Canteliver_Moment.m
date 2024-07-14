clc
clear
close all
%% MinCompliance Volume
% Phi is constant Phi = 0.1;
% ngp = Number of guss points
% dx = Length of each element at x direction
% dy = Length of each element at y direction
% E = Modulus of elasticity
% v = Possion ratio
% nor = Number of rectangle
% lx = Length of each element at x direction
% ly = Length of each element at y direction
% X0 = X-coordinate of left-up for each rectangle
% Y0 = Y-coordinate of left-up for each rectangle
% EBc,NBc: Essensual and Necessary boundary conditions
%% Canteliver example
MP = struct("Emax",1,"Emin",1e-3,"v",0.3);
GP = struct("ngp",4,"nor",1,"dx",1,"dy",1,"X0",0,"Y0",30,"lx",60,"ly",30);
OPP = struct("tho",0.5,"Noi",150,"dt",0.05,"SP",1);
%%
Emax = MP.Emax; Emin = MP.Emin; v = MP.v;
ngp = GP.ngp; nor = GP.nor; dx=GP.dx; dy=GP.dy;
X0 = GP.X0; Y0 = GP.Y0;lx = GP.lx; ly = GP.ly;
%%
X = cell(nor,1); Y = cell(nor,1); % Nodes coordinates
NodeCoord = []; nelx = zeros(nor,1); nely = zeros(nor,1);
coord = cell(nor,1);
for i=1:nor
    x = X0(i):dx:X0(i)+lx(i); y = Y0(i):-dy:Y0(i)-ly(i);
    [X{i},Y{i}] = meshgrid(x,y);
    nelx(i) = numel(x)-1; nely(i) = numel(y)-1;
    xNode = repmat(x',nely(i)+1,1); 
    yNode = reshape(repmat(y,nelx(i)+1,1),(nelx(i)+1)*(nely(i)+1),1);
    coord{i} = round([xNode,yNode],5);
    NodeCoord = [NodeCoord;coord{i}];
end
NodeCoord = unique(NodeCoord,'stable','rows');
Xmin = min(NodeCoord(:,1));Xmax = max(NodeCoord(:,1));
Ymin = min(NodeCoord(:,2));Ymax = max(NodeCoord(:,2));
%%
ElNodes = [];IC=cell(nor,1);
for i=1:nor
    nodes = (1:size(coord{i},1))';
    [~,ic,iC] = intersect(coord{i},NodeCoord,'stable','rows');
    nodes(ic) = iC;IC{i}=iC;
    eleN1 = (repmat((1:nelx(i)),nely(i),1)+nelx(i)*(0:nely(i)-1)')'; % Index of elements
    eleNode = repmat(eleN1(:),1,4)+repmat([0,1,nelx(i)+...
        [1,2]],nelx(i)*nely(i),1)+kron((0:nely(i)-1)',ones(nelx(i),1));
    eleNode = nodes(eleNode);
    ElNodes = [ElNodes;eleNode]; %Nodes of elements
end
%% Display Initial design
x=NodeCoord(:,1);y=NodeCoord(:,2);
%% Full solid
Phi=0.1*ones(size(NodeCoord,1),1); %Const phi
Phimax=1;
filename = 'Fullsolid';
%% Case1
% Phi(y<10 | y>20)=-0.1;
% filename = 'Case1';
%% Case2
% Phi((y>10 & y<20) & x<50)=-0.1;
% filename = 'Case2';
%% Case3
% indy = (x>12.5 & x<17.5) | (x>27.5 & x<32.5) | (x>42.5 & x<47.5);
% indx = y>12.5 & y<17.5; Phi(indx & indy)=-0.1;
% filename = 'Case3';
%% Case4
% indy = (x>5 & x<10) | (x>20 & x<25) | (x>35 & x<40) | (x>50 & x<55); 
% indx = (y>5 & y<10) | (y>20 & y<25);
% Phi(indx & indy)=-0.1;
% filename = 'Case4';
%% Case5
% xc = 30;yc = 15;
% Phi=max(min(min(sqrt((x(:)-(xc(:))').^2+(y(:)-(yc(:))').^2)-8,[],2),Phimax),-1);
% filename = 'Case5';
%% Case6
% xc=5:10:55;yc=5:10:25;
% [xc,yc]=meshgrid(xc,yc);
% Phi=max(min(min(sqrt((x(:)-(xc(:))').^2+(y(:)-(yc(:))').^2)-4,[],2),Phimax),-Phimax);
% filename = 'Case6';
%% Case7
% n=3;
% for i=1:nor
%     x=NodeCoord(IC{i},1);y=NodeCoord(IC{i},2);
%     wx=2*pi/(n*dx);wy=2*pi/(n*dy);
%     Phi(IC{i}) = round(sin((x-X0(i))*wx+pi/2).*sin((Y0(i)-y)*wy+pi/2),4);
% end
% filename = 'Case7';
%%
Phi0=Phi;
%%
figure;
view(3)
% Define the 8 corners of the cube
vertices = [Xmin Ymin -1;
            Xmax Ymin -1;
            Xmax Ymax -1;
            Xmin Ymax -1;
            Xmin Ymin 1;
            Xmax Ymin 1;
            Xmax Ymax 1;
            Xmin Ymax 1];
% Define the faces by connecting the vertices
faces = [1 2 6 5;  % bottom face
         2 3 7 6;  % right face
         3 4 8 7;  % top face
         4 1 5 8;  % left face
         1 2 3 4;  % front face
         5 6 7 8]; % back face
patch('Vertices',vertices,'Faces',faces,'FaceColor','none',...
    'EdgeColor','black',"linewidth",1.2)
hold on
for i=1:nor
    surf(X{i},Y{i},reshape(Phi(IC{i}),(nelx(i)+1),(nely(i)+1))');
end
axis tight
box on
% axis equal
hold off
figure
for i=1:nor
    contourf(X{i},Y{i},reshape(Phi(IC{i}),(nelx(i)+1),(nely(i)+1))',[0,0])
    hold on
end
axis equal
%%
noe = sum(nelx.*nely); ndof = 2*max(ElNodes,[],'all');
eDofMat = kron(ElNodes,[2,2])+repmat([-1,0],1,4);
iK = reshape(kron(eDofMat,ones(8,1))',64*noe,1);
jK = reshape(kron(eDofMat,ones(1,8))',64*noe,1);
iT = reshape(kron(ElNodes,ones(4,1))',16*noe,1);
jT = reshape(kron(ElNodes,ones(1,4))',16*noe,1);
%% Stifness calculation for each element
[r,W]=makegussianpoint(ngp); % guss points and their wheights
ke = zeros(8); ks = ke;
Emat = 1/(1-v^2)*[1,v,0;...
    v,1,0;...
    0,0,(1-v)/2]; % Elasticity matrix
eXcor = [0,dx,0,dx]';eYcor = [0,0,-dy,-dy]'; % Element X-Y coordinates
Bmatk = zeros(3,8,noe,ngp^2); Ar = zeros(ngp^2,1);
t1e = zeros(4);
t2e = zeros(4);
for i=1:ngp
    n = r(ngp-i+1);
    for j=1:ngp
        e=r(j);
        N=0.25*[(1-e)*(1+n),(1+e)*(1+n),(1-e)*(1-n),(1+e)*(1-n)]; % 4-node shape functions
        N_e = 0.25*[-(1+n),(1+n),(n-1),(1-n)];
        N_n = 0.25*[(1-e),(1+e),(e-1),-(1+e)];
        X_e = N_e*eXcor; X_n = N_n*eXcor;
        Y_e = N_e*eYcor;Y_n = N_n*eYcor;
        J_matrix = [X_e Y_e;X_n Y_n];
        detJ = det(J_matrix);
        N_XY=J_matrix\[N_e;N_n];
        b = [kron(N_XY(1,:),[1,0]);kron(N_XY(2,:),[0,1]);...
            kron(N_XY(1,:),[0,1])+kron(N_XY(2,:),[1,0])];
        Bmatk(:,:,:,ngp*(i-1)+j) = repmat(b,1,1,noe,1);
        b=[kron(N_XY,[1,0]);kron(N_XY,[0,1])];
        Ar(ngp*(i-1)+j) = detJ*W(ngp-i+1)*W(j);
        ke = ke+Bmatk(:,:,1,ngp*(i-1)+j)'*Emat*Bmatk(:,:,1,ngp*(i-1)+j)*Ar(ngp*(i-1)+j);
        t1e = t1e + N*N'*Ar(ngp*(i-1)+j);
        t2e = t2e + N_XY'*N_XY*Ar(ngp*(i-1)+j);
    end
end
Emat = repmat(Emat,1,1,noe,ngp^2);
st1 = reshape(t1e(:)*ones(1,noe),16*noe,1);
st2 = reshape(t2e(:)*ones(1,noe),16*noe,1);
T1 = sparse(iT,jT,st1);
T2 = sparse(iT,jT,st2);
Ar = repmat(reshape(Ar,1,1,1,ngp^2),1,1,noe,1);
%% NBC and EBC Condition
x = NodeCoord(:,1); y = NodeCoord(:,2);
indEBC = find(NodeCoord(:,1)==0);
EBC = [indEBC,zeros(numel(indEBC),2)];
indNBC = find(abs(NodeCoord(:,2)-15)<=0.1 & NodeCoord(:,1)==60);
NBC = [indNBC,zeros(numel(indNBC),1),-1/numel(indNBC)*ones(numel(indNBC),1)];
%% Boundary condition implementation
U = nan(ndof,1); F = zeros(ndof,1); AdjV = zeros(ndof,1); 
U(EBC(:,1)*2-1) = EBC(:,2);U(EBC(:,1)*2) = EBC(:,3);
F(NBC(:,1)*2-1) = NBC(:,2);F(NBC(:,1)*2) = NBC(:,3);
ukdis = isnan(U); kdis = ~isnan(U); % known and unknown displacement
FixPhi = NBC(:,1);
% Phi(FixPhi)=1;
%% initialize some parameters
tho = OPP.tho;
Noi = OPP.Noi;
dt = OPP.dt;
[s,t]=meshgrid(linspace(-1,1,20));
%% PrepheralPhi
do = find(y==0);
up = find(y==30);
ri = find(x==60);
le = find(x==0);
PrPhi = setdiff([do;up;ri;le],FixPhi);
%% initialize parameter of constraints
Lambda=0;
Mu = 5;
G0=0.4;
Gamma = 0.5;
dGamma = 0.05;
maxGamma = 5;
nRelax = 30;
%% Optimization Loop
J = 0;
G = 0;
k=1;
% vsave = VideoWriter(filename);
% open(vsave);
while k<=OPP.Noi
    % Volume fraction calculation for each element
    tmPhi=(1-s(:)).*(1-t(:))*Phi(ElNodes(:,3))'/4+(1+s(:)).*(1-t(:))*...
        Phi(ElNodes(:,4))'/4+(1-s(:)).*(1+t(:))*...
        Phi(ElNodes(:,1))'/4+(1+s(:)).*(1+t(:))*Phi(ElNodes(:,2))'/4;
    vfe = (sum(tmPhi>=0)/numel(s))';
    G(k)= sum(vfe)/noe;
% Global stiffness calculation
    Ee = vfe'*(Emax-Emin)+Emin;
    sK = reshape(ke(:)*(Ee.*ones(1,noe)),64*noe,1);
    K = sparse(iK,jK,sK); 
    K = (K+K')/2;
% Solving system of equlibrium equations
    U(ukdis)=K(ukdis,ukdis)\(F(ukdis)-K(ukdis,kdis)*U(kdis));
% Objective function sensitivity wrt phi
    ElemComp=sum(0.5*(ke*U(eDofMat)').*(U(eDofMat)'.*Ee));
    J(k) = sum(ElemComp);
    %% CONVERGENCE CHECK
    if k>nRelax && abs(G(k)-G0)/G0<1e-4 && all(abs(J(k)-J(k-9:k-1))/J(k)<1e-4)
        break;
    end
    %%
    J_phi = sparse(ElNodes,ones(noe,4),0.25*ElemComp'.*ones(noe,4));
    if k<nRelax
        Lambda=Mu*(G(k)-G(1)+(G(1)-G0)*k/nRelax);
    else
        Lambda =Lambda+Gamma*(G(k)-G0);
        Gamma = min(Gamma+dGamma,maxGamma);
    end
    V = -J_phi/mean(abs(J_phi))+Lambda;
%     V = normalize(-J_phi,'zscore')+Lambda;
 % Update scheme
    T = (T1/dt+tho*T2);
    Yy = (T1*(Phi/dt-V));
    Phi=T\Yy;
    Phi = min(max(Phi,-Phimax),Phimax);
    Phi(FixPhi)=1;
    ind = find(Phi(PrPhi)>-1e-3);
    Phi(PrPhi(ind))=-1e-3;
    if OPP.SP==1
        Xdef = NodeCoord(:,1);
        Ydef = NodeCoord(:,2);        
        for j=1:nor
            x = reshape(Xdef(IC{j}),(nelx(j)+1),(nely(j)+1))';
            y = reshape(Ydef(IC{j}),(nelx(j)+1),(nely(j)+1))';
            z = reshape(Phi(IC{j}),(nelx(j)+1),(nely(j)+1))';
            contourf(x,y,z,[0,0])
            hold on
            axis('equal')
        end
        title('Optimum Shape');
        hold off
%         writeVideo(vsave,getframe)
        drawnow
    end
    
    [k,J(k),G(k)]
    k=k+1;
end
% close(vsave)
