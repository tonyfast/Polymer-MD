%% Principal Components Analysis of Polymer Molecular Dynamics Simulations
%
%% Workflow
%
% # Converted Raw Data into Structure data -
% <https://github.com/tonyfast/Polymer-MD/blob/master/CODES/Analytics/Converted_Data.m
% Convert_Data.m>
% # Create Tree structures from point cloud data -
% <https://github.com/tonyfast/Polymer-MD/blob/master/CODES/Analytics/Tree_Build.m
% Tree_Build.m>
% # Partition the Volume using the trees -
% <https://github.com/tonyfast/Polymer-MD/blob/master/CODES/Analytics/Tree2Partition_Stats.m
% Tree2Partition_Stats.m>
% # Compute the Spatial Statistics and PCA Embedding -
% <https://github.com/tonyfast/Polymer-MD/blob/master/CODES/Analytics/Analyze_Stats.m
% Analyze_Stats.m >

%% Reporting

% Load PCA scores and metadata for analysis
load PCA_embed.mat

% Visualize local Results
gplot3( G(vals, vals), U2(:,1:3) )
hold on
scatter3(  U2(:,1),U2(:,2), U2(:,3), 50,info2(:,1), 'filled' )
text(U2(:,1),U2(:,2), U2(:,3), num2str(info2(:,2)))
hold off
close all
% dcm = datacursormode;
% set(dcm,'UpdateFcn',@(x,y)view_molecules(x,y,U2, AA))

%% Create Plotly
%
% * Download plot.ly
% <https://plot.ly/api/MATLAB/getting-started/installation Matlab Api>
% * Follow the <https://plot.ly/api/MATLAB/getting-started Getting Started>
% instructions.

plotdata{1} = struct;
co = cbrewer( 'qual','Paired',20);
ct = 0;
for ii = 1 : max(info2(:,1))
    
    b = info2(:,1) == ii;
    if any(b)
        ct = ct + 1;
        plotdata{ct}.x = U2(b,1);
        plotdata{ct}.y = U2(b,2);
        plotdata{ct}.name = num2str(unique(info2(b,1)));
        plotdata{ct}.line = struct('width',0,'color','');
        plotdata{ct}.marker = struct('opacity',.9, ...
            'symbol','circle', ...
            'size',16, ...
            'color', sprintf('rgb(%i,%i,%i)',round(co(ct,:)*255)), ...
            'line', struct('width',3,'color','black') );
    end
end
layout = struct( 'showlegend' , true, ...
    'xaxis' , struct( 'title', 'First Principal Components',...
    'autorange' , true ), ...
    'yaxis' , struct( 'title', 'Second Principal Components',...
    'autorange' , true ) ...
    );

% response = plotly( plotdata, struct('layout',layout) );

%% 
% <html>
% <iframe height="600" id="igraph" scrolling="no" seamless="seamless" src="https://plot.ly/~TonyFast/8" width="100%"></iframe>
% </html>
% This plot shows the PCA embedding of some of the Polymer Molecular
% Dynamics Simulations.

