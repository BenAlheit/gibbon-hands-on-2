clear; close all; clc;

%% Plot settings
markerSize = 25;

%% Control parameter

sampleSize = 2;
numSplitSteps = 2;

%% Mesh
v = [0 0 0;
     1 0 0;
     1 1 0;
     0 1 0;
     0 0 1;
     1 0 1;
     1 1 1;
     0 1 1] * sampleSize;

E = [1 2 3 4 5 6 7 8]; % Hex element

% [E,V,C,CV]=subHex(E,V,numSplitSteps,1);


[E,v]=subHex(E,v,numSplitSteps,1);

% F = [1 2 3 4; 5 6 7 8];% Face / element array

[F] = element2patch(E);% Face / element array

indBoundary = tesBoundary(F, v);

Fb = F(indBoundary, :);

%% Anfgular threshold

a = (45/180)*pi;


Cb = patchFeatureDetect(Fb, v, a);
C = [1:1:size(F,1)]';

%% Plot
cFigure;

% plotV(v, 'k.', 'MarkerSize', markerSize)

gpatch(Fb, v, Cb, 'k', 1, 3);
colormap viridis; icolorbar;
axisGeom;
drawnow;
%% Logic top nodes

L = Cb == 6;

dsfsd

%%

L = C==2;

F_top = F(L, :);
indTopNodes = unique(F_top);

%%
Z = v(:, 3);

logicTopNodes = Z>=(max(Z)-eps(0));
indTopNodes = find(logicTopNodes);

logicTopFaces = all(ismember(F, indTopNodes), 2);
F_top = F(logicTopFaces, :);
%% Plot
cFigure;

% plotV(v, 'k.', 'MarkerSize', markerSize)


% plotV(v(indTopNodes, :), 'r', 'Markersize', 25);
plotV(v(logicTopNodes, :), 'r.', 'Markersize', 25);
gpatch(F, v, 'w', 'k', 1, 3);
patchNormPlot(F, v)
gpatch(F_top, v, 'b', 'k', 1, 3);
colormap viridis; icolorbar;
axisGeom;
drawnow;