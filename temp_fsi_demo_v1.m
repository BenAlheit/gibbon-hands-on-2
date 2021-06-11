clear; close all; clc; 

%%
% Plot settings
fontSize=15;
faceAlpha1=0.3;
faceAlpha2=1;
cMap=gjet(4); 
patchColor=cMap(1,:);
markerSize=25; 

%% Control parameters

pointSpacing=3; 
cylRadius=10; 
cylHeight=cylRadius*7;
nt=ceil(cylHeight/pointSpacing);
numHeight=nt+iseven(nt); 

% Creating input structure
cylInputStruct.cylRadius=cylRadius;
cylInputStruct.numRadial=ceil((2*pi*cylRadius)/pointSpacing);
cylInputStruct.cylHeight=cylHeight;
cylInputStruct.numHeight=numHeight;
cylInputStruct.meshType='tri';
cylInputStruct.closeOpt=1;

layerThickness=cylRadius/5; 
numSteps=2; 

%%
% Derive patch data for a cylinder
[F,V,C]=patchcylinder(cylInputStruct); 
R=euler2DCM([0 0.5*pi 0]);
V=V*R; 

V_regions=getInnerPoint(F,V); %Define region points
V_holes=[]; %Define hole points
[regionTetVolumes]=tetVolMeanEst(F,V); %Volume estimate for regular tets
stringOpt='-pq1.2AaY'; %Options for tetgen

%%
% Mesh using TetGen

%Create tetgen input structure
inputStruct.stringOpt=stringOpt; %Tetgen options
inputStruct.Faces=F; %Boundary faces
inputStruct.Nodes=V; %Nodes of boundary
inputStruct.faceBoundaryMarker=C; 
inputStruct.regionPoints=V_regions; %Interior points for regions
inputStruct.holePoints=V_holes; %Interior points for holes
inputStruct.regionA=regionTetVolumes; %Desired tetrahedral volume for each region

% Mesh model using tetrahedral elements using tetGen 
[meshOutput]=runTetGen(inputStruct); %Run tetGen 

%% 
% Access mesh output structure

E=meshOutput.elements; %The elements
V=meshOutput.nodes; %The vertices or nodes
CE=meshOutput.elementMaterialID; %Element material or region id
Fb=meshOutput.facesBoundary; %The boundary faces
Cb=meshOutput.boundaryMarker; %The boundary markers

%%
% Visualization

hf=cFigure; 
subplot(1,2,1); hold on;
title('Input boundaries','FontSize',fontSize);
hp(1)=gpatch(Fb,V,Cb,'k',faceAlpha1);
hp(2)=plotV(V_regions,'r.','MarkerSize',markerSize);
legend(hp,{'Input mesh','Interior point(s)'},'Location','NorthWestOutside');
axisGeom(gca,fontSize); camlight headlight;
colormap(cMap); icolorbar;

hs=subplot(1,2,2); hold on;
title('Tetrahedral mesh','FontSize',fontSize);

% Visualizing using |meshView|
optionStruct.hFig=[hf,hs];
meshView(meshOutput,optionStruct);

axisGeom(gca,fontSize); 
gdrawnow;

%%
% Visualizing meshed regions

cFigure; 
gpatch(Fb,V,Cb);
colormap(gjet(3)); icolorbar; 
axisGeom; 
camlight headlight; 
drawnow; 

%% Create pentahedral skin
[Ep,Vp,Fq1,Fq2]=patchThick(Fb(Cb==1,:),V,-1,layerThickness,numSteps);

%Use element2patch to get patch data 
Fp=element2patch(Ep,[],'penta6');

%%
% Visualizing meshed regions

cFigure; 
gpatch(Fb,V,'w','none',0.5);
gpatch(Fp,Vp,'rw','k',1);
axisGeom; 
camlight headlight; 
drawnow; 

%%

F_inlet=Fb(Cb==2,:);
Eb_inlet=patchBoundary(F_inlet,V); 
indBoundaryNodes_inlet=edgeListToCurve(Eb_inlet);
F_outlet=Fb(Cb==3,:);
Eb_outlet=patchBoundary(F_outlet,V); 
indBoundaryNodes_outlet=edgeListToCurve(Eb_outlet);

%%
% Visualizing meshed regions

cFigure; hold on; 
gpatch(Fb,V,'w','none',0.5);
hp1=gpatch(F_inlet,V,'rw','k',1);
hp2=gpatch(F_outlet,V,'bw','k',1);
hp3=plotV(V(indBoundaryNodes_inlet,:),'r-','LineWidth',3);
hp4=plotV(V(indBoundaryNodes_outlet,:),'b-','LineWidth',3);
legend([hp1 hp2 hp3 hp4],{'Inlet faces','Outlet faces','Inlet boundary nodes','Outlet boundary nodes'})
axisGeom; 
camlight headlight; 
drawnow; 

