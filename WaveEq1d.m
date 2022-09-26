%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%       A conservative FD solver for a System of Hyperbolic PDE's
%                  by Manuel Diaz, ENSMA, 2021.02.26
%
%   A Numerical solver based on weno scheme for solving the second order
%   wave equation q_tt + c^2 q_xx = 0, with c > 0, as a 1st order system
%
%                    [p]     [ 0 , r c^2 ][p]   
%                    [u]_t + [ 1 / r , 0 ][u]_x = 0,
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Refs:
% [1] Brady, Peter T., and Daniel Livescu. "High-order, stable, and
%     conservative boundary schemes for central and compact finite
%     differences." Computers & Fluids 183 (2019): 84-101.  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc; close all;
global c r Dx periodic_x idL idR

% Domain size and CFL 
	  N = [ 21, 41, 81,161,321]; % number of cells
    CFL = [0.9,0.9,0.9,0.9,0.9]; % CFL condition

% Propagation Media properties
	c = 1.0;	% speed of sound
	r = 1.0;	% media density

% Parameters
      t0 = 0.;           % Initial time
    tEnd = 6.1;          % Output time
FDscheme ='lele643';     % lele643, pade43 (only!)
RKscheme = 04;           % 3:ERK3, 4:ERK4
plot_fig = true;%false   % visualize evolution of the equation.
    PATH = './figures/'; % path to folder of figures

% Build figures folder if it nos exists
if not(isfolder(PATH)), mkdir(PATH); end

% Error container
L2norm = []; OOA = true;

for n=1:numel(N), nx=N(n); cfl=CFL(n);
% Discretize spatial domain
periodic_x = false;
switch periodic_x
    case 0, ax=0; bx=1; dx=(bx-ax)/(nx-1); % non-periodic
    case 1, ax=0; bx=1; dx=(bx-ax)/( nx ); % periodic
end, x=linspace(ax,bx,nx)'; xc=ax+0.75*(bx-ax);

% Build initial condition (IC)
p0 = -3*pi/2*sin(3*pi*x/2);
u0 = zeros(size(x));
q0 = [p0,u0];

% Bluid exact solution (ES)
pE = @(t) -3*pi/4*(sin(3*pi*(x-t)/2)+sin(3*pi*(x+t)/2));
uE = @(t) -3*pi/4*(sin(3*pi*(x-t)/2)-sin(3*pi*(x+t)/2));

% Initial time step
dt0 = cfl*dx/c;

% print to terminal dx and dt
fprintf('using dx: %1.4f and dt: %1.4f\n',dx,dt0);

%% Build numerical scheme

% Build scheme
FD = compactSchemes(FDscheme,nx,periodic_x);

% Diff-operators
Dx = FD.Dx/dx;

% Boundary masks
idL = FD.index_L;
idR = FD.index_R;

% Define RHS
RHS = @WaveEq1d_RHS;

%% Solver Loop

% Load IC
q=q0; t=t0; it=0; dt=dt0; 

% Prepare visualization figures
if plot_fig~=false
    figure(1); region=[ax,bx,-5.2,5.2];
    title_string = sprintf('%s, time %1.2f',FDscheme,t);
    subplot(2,1,1); hu=plot(x,p0,'.k',x,pE(t),':r'); 
    axis(region); grid on; grid minor;
    ylabel('$\wp(x,t)$','interpreter','latex','fontsize',20);
    hT=title(title_string,'interpreter','latex','fontsize',20);
    subplot(2,1,2); hp=plot(x,u0,'.k',x,uE(t),':r'); 
    axis(region); grid on; grid minor;
    xlabel('$x$','interpreter','latex','fontsize',20);
    ylabel('$u(x,t)$','interpreter','latex','fontsize',20);
    legend(hu,{FDscheme,'Exact'},'location','northwest','interpreter','latex'); 
    legend(hp,{FDscheme,'Exact'},'location','northwest','interpreter','latex');
end

% for max norm plot
time=[]; max_norm=[];

tic
while t < tEnd
    % Set iteration time
    if (t+dt)>tEnd; dt=tEnd-t; end
    
    switch RKscheme
        case 3 % SSP-RK33
            qo= q;
            q = qo+dt*dF(q,t);
            q = 0.75*qo+0.25*(q+dt*dF(q,t+dt/2));
            q = (qo+2*(q+dt*dF(q,t+dt)))/3;
        case 4 % ERK4
            qo= q;          L1 = RHS(q,t);
            q = qo+dt/2*L1; L2 = RHS(q,t+0.5*dt);
            q = qo+dt/2*L2; L3 = RHS(q,t+0.5*dt);
            q = qo+dt*L3;   L4 = RHS(q,t+dt);
            q = qo+dt*(L1+2*(L2+L3)+L4)/6;
        otherwise, error('ERROR: RK scheme not set :P');
    end
    
    % Update time
    t=t+dt;
        
    % compute conserved properties
    p=q(:,1); u=q(:,2);
    
    % Update dt and iteration counter
    dt=cfl*dx/c; it=it+1;
    
    % Compute max-norm
    time=[time,t]; max_norm=[max_norm,norm(pE(t)-p,Inf)]; %#ok<AGROW>
    
    % Plot figure
    if plot_fig && (rem(it,10)==0)
        hu(1).YData = u;   hu(2).YData = uE(t);
        hp(1).YData = p;   hp(2).YData = pE(t);
        hT.String = sprintf('%s, time %1.2f',FDscheme,t);
        drawnow
    end
end
cputime = toc; fprintf('CPUtime: %1.2f\n',cputime);

%% Post-process 

% get flow properties
p=q(:,1); u=q(:,2);

if OOA
    % Compute error L2-norm (and display it!)
    err = (p-pE(t)).^2;
    err = sqrt(sum(err))/ N(n); 
    disp(['L2norm = ',num2str(err)]);
    L2norm = [L2norm,err]; %#ok<AGROW>
end

if plot_fig 
    hu(1).YData = u;   hu(2).YData = uE(t);
    hp(1).YData = p;   hp(2).YData = pE(t);
    hT.String = sprintf('%s, time %1.2f',FDscheme,t);
    drawnow;
end

% Stability graph
fig=figure(2);
semilogy(time,max_norm); hold on; xlim([0,tEnd]); ylim([1E-8,2]);
end, hold off;
title(FDscheme,'interpreter','latex','fontsize',20);
xlabel('$t$','interpreter','latex','fontsize',20);
ylabel('$L_\infty$','interpreter','latex','fontsize',20);
legend({['N=',num2str(N(1))],...
        ['N=',num2str(N(2))],...
        ['N=',num2str(N(3))],...
        ['N=',num2str(N(4))],...
        ['N=',num2str(N(5))]},'location','southeast','interpreter','latex'); 
switch periodic_x
    case 0, print(fig,[PATH,'WaveEq1d_',FDscheme,'_RK',num2str(RKscheme),'_stability'],'-dpng');
    case 1, print(fig,[PATH,'WaveEq1d_',FDscheme,'_RK',num2str(RKscheme),'_stability_periodic'],'-dpng');
end

% Save OOA data
if OOA, save([PATH,'WaveEq1d_',FDscheme,'_RK',num2str(RKscheme),'_OOA.mat'],'N','L2norm'); end