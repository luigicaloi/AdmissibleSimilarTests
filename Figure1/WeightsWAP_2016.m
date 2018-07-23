%This script replicates Figure 3 in p. 19 of Montiel-Olea (2016)

%More generally, this file approximates the density of the weights in 
%(3.9) and (3.10) for different values of the reduced-form correlation
%and different number of instruments

%This version: March 19th, 2016
%Last revision: May 24th, 2016.


%% 1) Inputs
rho=.5; %Correlation parameter for the matrix Psi
k=4;    %Number of instruments
beta0=0; b0=[1;-beta0]; a0=[beta0;1];
Psimatrix=[1,rho;rho,1]; %Matrix Psi (assumed to have 1 in the diagonals)
C0=[((b0'*Psimatrix*b0)^(-1/2))*b0';((a0'*(Psimatrix^(-1))*a0)^(-1/2))*a0'*Psimatrix];

%% 2) Draws from rho and phi

%Since phi is uniform on the sphere, phi can be represented as
%phi=[cos(\theta),sin(\theta)]', with theta \in [0,2pi] (See Stroock(1999) p.87)
anginf=0; 
angsup=+2*pi;
 
grid=(rand(50000,1).*(angsup-anginf))+anginf; %set the angle theta
phi=[cos(grid),sin(grid)]'; %uniform on the circle
z=sum(randn(50000,k).^2,2)'; %chi^2_k
%Grid for the values of phi and the values of the chi-square part of lambda


%% 3) Simulation of the weights for \lambda and \lambda^.5 (\beta-\beta0)

weightl=bsxfun(@times,([0,1]*Psimatrix*C0'*phi).^2,z); %Formula 3.10 in the paper
weightlp= ((b0'*Psimatrix*b0)^(1/2))*((z.^.5).*phi(1,:)); %Formula 3.9 in the paper



%% 4) Compute the joint density

%To match the simulations in AMS(06), I report joint density over
%\lambda^.5 (\beta-beta0) and \lambda. This requires the program gkde2


p=gkde2([weightlp',weightl']);

%Joint Density
figure
surf(p.x,p.y,p.pdf);
view(202,34);
axis([floor(min(weightlp))-1 ceil(max(weightlp))+1 0 20 0 .1]);
xlabel('\lambda^{1/2} (\beta-\beta_0)');
ylabel('\lambda');
colormap(jet);
nam=strcat('k=',num2str(k),'_','rho=',num2str(rho),'_','FigWeightsFinal_2016.eps');
name=strcat(nam);
print(gcf,'-depsc2',name);







