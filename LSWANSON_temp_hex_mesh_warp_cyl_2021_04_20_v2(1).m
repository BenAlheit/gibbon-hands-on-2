clear; close all; clc; 

%% PLOT SETTINGS

fontSize=15;
edgeWidth=0.5; 

%%
% Path names
defaultFolder = fileparts(mfilename('fullpath'));
saveFolder = fileparts(mfilename('fullpath'));
saveName=fullfile(saveFolder,'hexMeshBaseLine_heart_basic.mat');
saveOn=0; 

interpMethod='linear';
extrapMethod='linear';

%% Building a quadrilateral circular mesh

ne=6; %Elements in radius
f=0.5; %Fraction (with respect to outer radius) where central square appears

%Create the mesh
[F_template,V_template]=discQuadMesh(ne,1,f);
V_template(:,3)=0;
F_template=fliplr(F_template);
Eb=patchBoundary(F_template,V_template);
indB=edgeListToCurve(Eb);
indB=indB(1:end-1);

cPar.n=50;
cPar.Method='LAP';
cPar.RigidConstraints=indB;
[V_template]=patchSmooth(F_template,V_template,[],cPar);

%% Create guide curves
pointSpacing = 1; %spacing on cross-section
axisSpacing = 2; %spacing on axis
radii = readmatrix(fullfile(defaultFolder, 'radii.txt')); %outerRadius array 1xn
radiiOuter = radii + 1;
radiiInner = radii; %inner radius array 1xn

nRadii = size(radiiOuter, 2); %number axis points
vesselLength = axisSpacing * nRadii; %length of vessel

axisPoints = (0 : nRadii-1) * axisSpacing; %array of axis coordinates
nRad = ne*4+1; %Number of radial steps based on a mean radius between max an min radii

t = linspace(0, 2 * pi, nRad)'; %Angles
t = t(1:end-1); %take away last which equals start

% Inner Wall boundary locations
x_inner = axisPoints.*ones(size(t));
y_inner = radiiInner.*sin(t); %array nRad-1 x nRadii where column i is the x coordinate of each point on the cross section
z_inner = radiiInner.*cos(t);%array nRad-1 x nRadii where column j is the y coordinate of each point on the cross section

% Building array of polygons in space
v_Inner = zeros(nRad-1, 3, nRadii);
V_control=cell(1,nRadii);
V_control_template=V_template(indB,:);

for radPoint=1:1:nRadii %Loop over each radiuss
    v_Inner(:, :, radPoint) = [x_inner(:, radPoint) y_inner(:, radPoint) z_inner(:, radPoint)];
    V_control{radPoint} = [x_inner(:, radPoint) y_inner(:, radPoint) z_inner(:, radPoint)];
end       

%%

% Cp = gjet(nRadii);
% 
% cFigure; 
% subplot(1,2,1); hold on;
% title('template');
% hp1=gpatch(F_template,V_template,'bw','k',1,edgeWidth);
% hp2=plotV(V_control_template,'r-','LineWidth',3,'MarkerSize',25);
% legend([hp1 hp2],{'Template mesh','Template control curve'});
% axisGeom(gca,fontSize);
% camlight headlight; 
% 
% subplot(1,2,2); hold on;
% title('Model boundary polygons','FontSize',fontSize);
% 
% for radPoint=1:1:nRadii
%     hp = plotV(v_Inner(: ,:, radPoint));
%     hp.Color = Cp(radPoint, :);
% end
% 
% axisGeom(gca,fontSize);
% camlight headlight; 
% drawnow;

%% Morph template quad mesh 

FT=cell(1,nRadii);
VT=cell(1,nRadii);

for q=1:1:nRadii
    FT{q}=F_template;
    VT{q}=interpMorph(V_template,V_control_template,V_control{q},interpMethod,extrapMethod);            
end

[FT,VT,CT]=joinElementSets(FT,VT);

%%
% %Visualizing mesh
% 
% cFigure; hold on;
% title('Morphed template meshes');
% gpatch(FT,VT,CT,'k',1,edgeWidth);
% 
% for q=1:1:nRadii
%     plotV(V_control{q},'g-','LineWidth',3,'MarkerSize',25);
% end
% 
% axisGeom(gca,fontSize);
% camlight headlight; 
% colormap gjet;  icolorbar;
% drawnow;

%% Build hexahedral elements

E1=[];
C1=[];
for q=1:1:nRadii-1
    e=[FT(CT==q,:) FT(CT==q+1,:)];
    E1=[E1; e];
    C1=[C1; q*ones(size(e,1),1)];
end

[F1,CF1]=element2patch(E1,C1);
V1=VT;

%%
% cFigure; hold on;
% title('Hex mesh fluid');
% 
% gpatch(F1,V1,CF1,'k',1,edgeWidth);
% patchNormPlot(F1,V1);
% 
% axisGeom(gca,fontSize);
% camlight headlight; 
% colormap gjet;  icolorbar;
% drawnow;

%% Zoning Cells on the wall for thickening
indBoundaryFaces=tesBoundary(F1,V1);
Fb1=F1(indBoundaryFaces,:);
Cb1=CF1(indBoundaryFaces,:);
Cb1_V=faceToVertexMeasure(Fb1,V1,Cb1);

logicSides1=all(Cb1_V(Fb1)==1,2);
logicSides2=all(Cb1_V(Fb1)==nRadii-1,2);
logicSides_wall = zeros(size(logicSides1));
logicSides_wall = all(logicSides_wall==logicSides1, 2);
logicSides_wall = ~all(logicSides_wall==logicSides2, 2);

Fb1_sides=Fb1(logicSides_wall,:);
Cb1_sides=Cb1(logicSides_wall,:);

Fb_inlet_fluid = Fb1(logicSides1,:);
Cb1_inlet_fluid=Cb1(logicSides1,:);
Fb_outlet_fluid = Fb1(logicSides2,:);
Cb1_outlet_fluid=Cb1(logicSides2,:);

%% Thickening the wall
dirSet=-1;
layerThickness=1;
numSteps=4;

[Ft,Vt]=patchCleanUnused(Fb1_sides,V1);
[E2,V2,F2_1,F2_2]=patchThick(Ft,Vt,dirSet,layerThickness,numSteps);
C2=repmat(Cb1_sides,[numSteps,1]);
[F2,CF2]=element2patch(E2,C2);

%% Plot the mesh
% cFigure; hold on;
% title('Hex mesh wall');
% % gpatch(Fb1_sides,V1,'kw','k',1,edgeWidth);
% gpatch(F2,V2,'b','k',1,edgeWidth);
% % patchNormPlot(F1,V1);
% 
% axisGeom(gca,fontSize);
% camlight headlight; 
% colormap gjet;  icolorbar;
% drawnow;

%% Construct Mesh

E=[E1;E2+size(V1,1)];
M=[ones(size(E1,1),1); 2*ones(size(E2,1),1)]; %Domain label
V=[V1;V2];
C=[C1;C2];
[F,CF,CF_type]=element2patch(E,C);
[~,MF]=element2patch(E,M);

indBoundaryFaces=tesBoundary(F,V);

%% Construct Mesh Struct

meshStruct.elements=E;
meshStruct.nodes=V;
meshStruct.faces=F;
meshStruct.elementMaterialID=M;

hFig=cFigure; 

subplot(1,2,1); hold on; 
title('Hex mesh');
gpatch(F,V,'bw','k',edgeWidth);
axisGeom(gca,fontSize);
camlight headlight; 
colormap gjet;

hs=subplot(1,2,2); hold on; 
title('Cut view of solid mesh','FontSize',fontSize);
optionStruct.hFig=[hFig hs];
meshView(meshStruct,optionStruct);
axisGeom(gca,fontSize);
drawnow; 

%% Get fluid domain inlet/outlet/outer faces

logicFluidFaces= MF==1; 
F_fluid=F(logicFluidFaces,:);
CF_fluid=CF(logicFluidFaces,:);
CF_type_fluid=CF_type(logicFluidFaces,:);
indBoundaryFacesFluid=tesBoundary(F_fluid,V);

F_fluid_boundary=F_fluid(indBoundaryFacesFluid,:);
CF_type_fluid_boundary=CF_type_fluid(indBoundaryFacesFluid,:);

F_inlet=F_fluid_boundary(CF_type_fluid_boundary==1,:);
F_outlet=F_fluid_boundary(CF_type_fluid_boundary==2,:);
F_fsi=F_fluid_boundary(CF_type_fluid_boundary==5,:);

%% Solid front/back node faces

logicSolidFaces= MF==2; 
F_solid=F(logicSolidFaces,:);
CF_solid=CF(logicSolidFaces,:);
CF_type_solid=CF_type(logicSolidFaces,:);

indBoundaryFacesSolid=tesBoundary(F_solid,V);

F_solid_boundary=F_solid(indBoundaryFacesSolid,:);
CF_type_solid_boundary=CF_type_solid(indBoundaryFacesSolid,:);

F_inlet_solid=F_solid_boundary(CF_type_solid_boundary==3,:);
F_outlet_solid=F_solid_boundary(CF_type_solid_boundary==4,:);

%%
% cFigure; hold on; 
% 
% gpatch(F(indBoundaryFaces,:),V,'w','none',0.25);
% gpatch(F_solid_boundary,V,CF_type_solid_boundary,'k',1);
% axisGeom(gca,fontSize);
% icolorbar;
% camlight headlight; 
% drawnow;

%%
cFigure; hold on; 

gpatch(F(indBoundaryFaces,:),V,'w','none',0.25);

hp1=gpatch(F_inlet,V,'rw','k',1);
hp2=gpatch(F_outlet,V,'gw','k',1);
hp3=gpatch(F_fsi,V,'bw','k',1);
hp4=gpatch(F_inlet_solid,V,'yellow','k',1);
hp5=gpatch(F_outlet_solid,V,'magenta','k',1);
% gpatch(F_fluid_boundary,V,CF_type_fluid_boundary,'k',1);

axisGeom(gca,fontSize);
legend([hp1 hp2 hp3 hp4 hp5],{'Inlet Flow','Outlet Flow','FSI Boundary','Inlet Solid', 'Outlet Solid'})
camlight headlight; 
drawnow; 

%%
% Eb_inlet=patchBoundary(F_inlet,V); 
% indBoundaryNodes_inlet=edgeListToCurve(Eb_inlet);
% 
% Eb_outlet=patchBoundary(F_outlet,V); 
% indBoundaryNodes_outlet=edgeListToCurve(Eb_outlet);
% 
% %%
% % Visualizing meshed regions
% 
% cFigure; hold on; 
% gpatch(F(indBoundaryFaces,:),V,'w','none',0.5);
% hp1=gpatch(F_inlet,V,'rw','k',1);
% hp2=gpatch(F_outlet,V,'bw','k',1);
% hp3=plotV(V(indBoundaryNodes_inlet,:),'r-','LineWidth',3);
% hp4=plotV(V(indBoundaryNodes_outlet,:),'b-','LineWidth',3);
% legend([hp1 hp2 hp3 hp4],{'Inlet faces','Outlet faces','Inlet boundary nodes','Outlet boundary nodes'})
% axisGeom; 
% camlight headlight; 
% drawnow; 

%%

if saveOn==1
    save(saveName,'meshStruct');
end

%%

function [V_target]=interpMorph(V_template,V_control_template,V_control_target,interpMethod,extrapMethod)

U_control=V_control_target-V_control_template;

interpFunction_Ux=scatteredInterpolant(V_control_template(:,[1 2]),U_control(:,1),interpMethod,extrapMethod);
interpFunction_Uy=scatteredInterpolant(V_control_template(:,[1 2]),U_control(:,2),interpMethod,extrapMethod);
interpFunction_Uz=scatteredInterpolant(V_control_template(:,[1 2]),U_control(:,3),interpMethod,extrapMethod);

Ux_template=interpFunction_Ux(V_template(:,[1 2])); %X-displacement interpolated
Uy_template=interpFunction_Uy(V_template(:,[1 2])); %Y-displacement interpolated
Uz_template=interpFunction_Uz(V_template(:,[1 2])); %Z-displacement interpolated

V_target=V_template+[Ux_template Uy_template Uz_template];

end