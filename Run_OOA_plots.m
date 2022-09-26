% Visualization parameters
set(groot,'defaultTextInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'DefaultTextFontName','Times',...
'DefaultTextFontSize',20,...
'DefaultAxesFontName','Times',...
'DefaultAxesFontSize',20,...
'DefaultLineLineWidth',1.5,...
'DefaultAxesBox','on',...
'defaultAxesLineWidth',1.5,...
'DefaultFigureColor','w',...
'DefaultLineMarkerSize',7.75)

% Saving Path
PATH = './figures/';

% List *.mat files
file = dir([PATH,'*.mat']);

for i = 1:size(file,1)
    % data file into memory
    load([PATH,file(i).name]);
    
    % get data file name
    fullfilename = file(i).name(1:end-4);
    [path,name,ext]=fileparts(fullfilename);
    
    % get data range
    N = N(1:end);
    L2norm = L2norm(1:end);

    % Plot OOA figure
    h = figure(999); colormap winter;
    loglog(N,L2norm,'o-b'); 
    title('Error - log plot');
    hold on
    loglog(N,L2norm(1)*(N(1)./N).^2,'--c');
    loglog(N,L2norm(1)*(N(1)./N).^4,'--m');
    loglog(N,L2norm(1)*(N(1)./N).^6,'--r');
    hold off
    legend({'observed','order 2','order 4','order 6'},'location','best');
    legend boxoff
    xlabel('N'); ylabel('error(N)');
    drawnow;
    %fig = gcf; fig.PaperUnits = 'inches'; fig.PaperPosition = [0 0 10 8];
    %print(h,[PATH,name],'-depsc');
    print(h,[PATH,name,'.png'],'-dpng');
end
    