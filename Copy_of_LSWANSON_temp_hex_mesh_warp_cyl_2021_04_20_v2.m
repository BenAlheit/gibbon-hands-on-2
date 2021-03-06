
clear; close all; clc

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

%% Join and Merging Mesh

E2 = E2+size(V1,1);
E=[E1;E2];
M=[ones(size(E1,1),1); 2*ones(size(E2,1),1)]; %Domain label
V=[V1;V2];
C=[C1;C2];
[F,CF,CF_type]=element2patch(E,C);
[~,MF]=element2patch(E,M);

[F, V, ~, indFix] = mergeVertices(F, V);

E1 = indFix(E1);
E2 = indFix(E2);
E = indFix(E);


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

%% Check elements

F2 = element2patch(E2);
F1 = element2patch(E1);

%% Plot faces

cFigure; hold on; 

hp1=gpatch(F1,V,'rw','k',1);
hp2=gpatch(F2,V,'bw','k',1);
axisGeom; 
camlight headlight; 
drawnow; 

%%
Eb_inlet=patchBoundary(F_inlet,V); 
indBoundaryNodes_inlet=edgeListToCurve(Eb_inlet);

Eb_outlet=patchBoundary(F_outlet,V); 
indBoundaryNodes_outlet=edgeListToCurve(Eb_outlet);
% 
% %%
% Visualizing meshed regions

cFigure; hold on; 

gpatch(F(indBoundaryFaces,:),V,'w','none',0.5);
hp1=gpatch(F_inlet,V,'rw','k',1);
hp2=gpatch(F_outlet,V,'bw','k',1);
hp3=plotV(V(indBoundaryNodes_inlet,:),'r-','LineWidth',3);
hp4=plotV(V(indBoundaryNodes_outlet,:),'b-','LineWidth',3);
legend([hp1 hp2 hp3 hp4],{'Inlet faces','Outlet faces','Inlet boundary nodes','Outlet boundary nodes'})
axisGeom; 
camlight headlight; 
drawnow; 

%%

if saveOn==1
    save(saveName,'meshStruct');
end


%% Control parameters

t_end = 5;
bpm = 60;
t = linspace(0, t_end, 1000);


% V = rand([10,3]);

E_solid = E1;
E_fluid = E2;

E_solid=E_solid(:,[5 6 7 8 1 2 3 4]);
E_fluid=E_fluid(:,[5 6 7 8 1 2 3 4]);

F_i_solid = fliplr(F_inlet_solid);
F_o_solid = fliplr(F_outlet_solid);

F_i_fluid = fliplr(F_inlet);
F_o_fluid = fliplr(F_outlet);

% F_FSI = fliplr(F_fsi);

F_FSI = F_fsi;

inletInd = unique(F_i_solid);
outletInd = unique(F_o_solid);

defaultFolder = 'C:\Users\alhei\dev\gibbon-hands-on';
savePath=fullfile(defaultFolder,'data','temp');

% Defining file names
febioFebFileNamePart='aorta';
febioFebFileName=fullfile(savePath,[febioFebFileNamePart,'.feb']); %FEB file name
febioLogFileName=[febioFebFileNamePart,'.txt']; %FEBio log file name
febioLogFileName_disp=[febioFebFileNamePart,'_disp_out.txt']; %Log file name for exporting displacement
febioLogFileName_strain=[febioFebFileNamePart,'_strain_out.txt']; %Log file name for exporting strain

%% Region names

inlet_solid_surface = 'inlet-solid';
inlet_fluid_surface = 'inlet-fluid';
outlet_solid_surface = 'outlet-solid';
outlet_fluid_surface = 'outlet-fluid';
FSI_surface = 'FSIInterface';

inlet_node_set = 'inlet-node-set';
outlet_node_set = 'outlet-node-set';

%% Visualize mesh


cFigure; hold on; 

gpatch(F(indBoundaryFaces,:),V,'w','none',0.5);
hp1=gpatch(F_i_solid,V,'rw','k',1);
hp2=gpatch(F_o_solid,V,'bw','k',1);
hp3=gpatch(F_i_fluid,V,'rw','k',1);
hp4=gpatch(F_o_fluid,V,'bw','k',1);
hp5=gpatch(F_FSI,V,'bw','k',1);
hp6=patchNormPlot(F_FSI,V);

% hp2=gpatch(F_outlet,V,'bw','k',1);
% hp3=plotV(V(indBoundaryNodes_inlet,:),'r-','LineWidth',3);
% hp4=plotV(V(indBoundaryNodes_outlet,:),'b-','LineWidth',3);
% legend([hp1 hp2 hp3 hp4],{'Inlet faces','Outlet faces','Inlet boundary nodes','Outlet boundary nodes'})
axisGeom; 
camlight headlight; 
drawnow; 

%% FSI Solver

%Get a template with default settings 
[febio_spec]=febioStructTemplate;

%febio_spec version 
febio_spec.ATTR.version='3.0'; 

%Module section
febio_spec.Module.ATTR.type='solid'; 

% febio_spec.Control.analysis='STATIC';

%Control section

% febio_spec.Module.ATTR.type='fluid-FSI'; 

% febio_spec.Module.ATTR.type='fluid-FSI'; 
% 
febio_spec.Control.analysis='DYNAMIC'; 
febio_spec.Control.time_steps=1550;
febio_spec.Control.step_size=0.001; 

febio_spec.Control.solver.max_refs = 5;
febio_spec.Control.solver.max_ups = 50;
febio_spec.Control.solver.diverge_reform = 0;
febio_spec.Control.solver.reform_each_time_step = 0;
febio_spec.Control.solver.dtol = 0.001;
% febio_spec.Control.solver.vtol = 0.001;
% febio_spec.Control.solver.ftol = 0.001;
febio_spec.Control.solver.etol = 0.01;
febio_spec.Control.solver.rtol = 0.001;
febio_spec.Control.solver.lstol = 0.9;
febio_spec.Control.solver.min_residual = 1.e-16;
febio_spec.Control.solver.max_residual = 1.e+10;
% febio_spec.Control.solver.rhoi = 0;
febio_spec.Control.solver.qnmethod = 'BROYDEN';
febio_spec.Control.solver.symmetric_stiffness = 0;


febio_spec.Control.time_stepper.dtmin = 0.0001;
febio_spec.Control.time_stepper.dtmax = 0.001;
febio_spec.Control.time_stepper.max_retries = 9;
febio_spec.Control.time_stepper.opt_iter = 53;

febio_spec.Globals.Constants.T = 0;
febio_spec.Globals.Constants.R = 0;
febio_spec.Globals.Constants.Fc = 0;

%% Material specification

% Solid
febio_spec.Material.material{1}.ATTR.name='artery-wall';
febio_spec.Material.material{1}.ATTR.type='neo-Hookean';
febio_spec.Material.material{1}.ATTR.id=1;
febio_spec.Material.material{1}.E=11700;
febio_spec.Material.material{1}.v=0.3;
febio_spec.Material.material{1}.density=1000;
% TODO change fluid to solid

% TODO check density and units in general
febio_spec.Material.material{2}.ATTR.name='blood';
febio_spec.Material.material{2}.ATTR.type = 'neo-Hookean';
% febio_spec.Material.material{2}.solid.ATTR.type = 'neo-Hookean';
febio_spec.Material.material{2}.ATTR.id=2;
febio_spec.Material.material{2}.E=11700;
febio_spec.Material.material{2}.v=0.3;
febio_spec.Material.material{2}.density=1000;

% TODO Paramter out material name
% Fluid
% febio_spec.Material.material{2}.ATTR.name='blood';
% febio_spec.Material.material{2}.ATTR.type='fluid-FSI';
% febio_spec.Material.material{2}.ATTR.id=2;
% febio_spec.Material.material{2}.fluid.ATTR.type = 'fluid';
% febio_spec.Material.material{2}.fluid.density = 1060;
% febio_spec.Material.material{2}.fluid.k = 2.2e9;
% febio_spec.Material.material{2}.fluid.viscous.ATTR.type = 'Carreau';
% febio_spec.Material.material{2}.fluid.viscous.mu0 = 0.056;
% febio_spec.Material.material{2}.fluid.viscous.mui = 0.00345;
% febio_spec.Material.material{2}.fluid.viscous.lambda = 3.313;
% febio_spec.Material.material{2}.fluid.viscous.n = 0.3568;
% febio_spec.Material.material{2}.solid.ATTR.type = 'neo-Hookean';
% febio_spec.Material.material{2}.solid.density = 0.;
% febio_spec.Material.material{2}.solid.E = 1.e-9;
% febio_spec.Material.material{2}.solid.v = 0;

%% Mesh

febio_spec.Mesh.Nodes{1}.ATTR.name = 'nodes';
febio_spec.Mesh.Nodes{1}.node.ATTR.id=(1:size(V,1))'; %The node id's
febio_spec.Mesh.Nodes{1}.node.VAL=V; %The nodel coordinates


febio_spec.Mesh.NodeSet{1}.ATTR.name=inlet_node_set;
febio_spec.Mesh.NodeSet{1}.node.ATTR.id=inletInd(:);


febio_spec.Mesh.NodeSet{2}.ATTR.name=outlet_node_set;
febio_spec.Mesh.NodeSet{2}.node.ATTR.id=outletInd(:);
% TODO add nodes

% Elemenets for artery blood flow
febio_spec.Mesh.Elements{1}.ATTR.type = 'hex8';
febio_spec.Mesh.Elements{1}.ATTR.name = 'fluid';
febio_spec.Mesh.Elements{1}.elem.ATTR.id=(1:size(E_fluid,1))'; 
febio_spec.Mesh.Elements{1}.elem.VAL=E_fluid; 

% TODO add elemenets

% Elemenets for artery wall
febio_spec.Mesh.Elements{2}.ATTR.type = 'hex8';
febio_spec.Mesh.Elements{2}.ATTR.name = 'solid';
febio_spec.Mesh.Elements{2}.elem.ATTR.id=(1:size(E_solid,1))'+size(E_fluid,1)'; 
febio_spec.Mesh.Elements{2}.elem.VAL=E_solid; 
% TODO add elemenets

% Surfaces
febio_spec.Mesh.Surface{1}.ATTR.name = inlet_solid_surface;
febio_spec.Mesh.Surface{1}.quad4.ATTR.id = (1:size(F_i_solid,1))';
febio_spec.Mesh.Surface{1}.quad4.VAL = F_i_solid;

febio_spec.Mesh.Surface{2}.ATTR.name = inlet_fluid_surface;
febio_spec.Mesh.Surface{2}.quad4.ATTR.id = (1:size(F_i_fluid,1))';
febio_spec.Mesh.Surface{2}.quad4.VAL = F_i_fluid;

febio_spec.Mesh.Surface{3}.ATTR.name = outlet_solid_surface;
febio_spec.Mesh.Surface{3}.quad4.ATTR.id = (1:size(F_o_solid,1))';
febio_spec.Mesh.Surface{3}.quad4.VAL = F_o_solid;

febio_spec.Mesh.Surface{4}.ATTR.name = outlet_fluid_surface;
febio_spec.Mesh.Surface{4}.quad4.ATTR.id = (1:size(F_o_fluid,1))';
febio_spec.Mesh.Surface{4}.quad4.VAL = F_o_fluid;


febio_spec.Mesh.Surface{5}.ATTR.name = FSI_surface;
febio_spec.Mesh.Surface{5}.quad4.ATTR.id = (1:size(F_FSI,1))';
febio_spec.Mesh.Surface{5}.quad4.VAL = F_FSI;

% Mesh domains
% TODO CHeck that name matches up 
febio_spec.MeshDomains.SolidDomain{1}.ATTR.name = 'solid';
febio_spec.MeshDomains.SolidDomain{1}.ATTR.mat = 'artery-wall';

febio_spec.MeshDomains.SolidDomain{2}.ATTR.name = 'fluid';
febio_spec.MeshDomains.SolidDomain{2}.ATTR.mat = 'blood';
% TODO Add mesh domains

%% Boundary conditions

% febio_spec.Boundary.bc{1}.ATTR.name = 'fix-inlet-solid';
febio_spec.Boundary.bc{1}.ATTR.type = 'fix';
% TODO Refer to node set instead of @surface
febio_spec.Boundary.bc{1}.ATTR.node_set = inlet_node_set;
febio_spec.Boundary.bc{1}.dofs = 'x,y,z';


febio_spec.Boundary.bc{2}.ATTR.type='prescribe';
febio_spec.Boundary.bc{2}.ATTR.node_set=outlet_node_set;
febio_spec.Boundary.bc{2}.dof='x';
febio_spec.Boundary.bc{2}.scale.ATTR.lc=1;
febio_spec.Boundary.bc{2}.scale.VAL=0.0002;
febio_spec.Boundary.bc{2}.relative=0;

% 
% febio_spec.Boundary.bc{2}.ATTR.name = 'fix-outlet-solid';
% febio_spec.Boundary.bc{2}.ATTR.type = 'fix';
% febio_spec.Boundary.bc{2}.ATTR.node_set = outlet_node_set;
% febio_spec.Boundary.bc{2}.dofs = 'x,y,z';


%% Loads
% febio_spec.Loads.surface_load{1}.ATTR.name = 'FSIInterfaceTraction';
% febio_spec.Loads.surface_load{1}.ATTR.type = 'fluid-FSI traction';
% febio_spec.Loads.surface_load{1}.ATTR.surface = FSI_surface;
% 
% febio_spec.Loads.surface_load{2}.ATTR.name = 'inlet-pressure';
% febio_spec.Loads.surface_load{2}.ATTR.type = 'fluid pressure';
% febio_spec.Loads.surface_load{2}.ATTR.surface = inlet_fluid_surface;
% febio_spec.Loads.surface_load{2}.pressure.ATTR.lc = 1;
% febio_spec.Loads.surface_load{2}.pressure.VAL = 500;
% 
% febio_spec.Loads.surface_load{3}.ATTR.name = 'inlet-backflow-stabilization';
% febio_spec.Loads.surface_load{3}.ATTR.type = 'fluid backflow stabilization';
% febio_spec.Loads.surface_load{3}.ATTR.surface = inlet_fluid_surface;
% febio_spec.Loads.surface_load{3}.beta.ATTR.lc = 2;
% febio_spec.Loads.surface_load{3}.beta.VAL  = 1;
% 
% febio_spec.Loads.surface_load{4}.ATTR.name = 'inlet-tangential-stabilization';
% febio_spec.Loads.surface_load{4}.ATTR.type = 'fluid tangential stabilization';
% febio_spec.Loads.surface_load{4}.ATTR.surface = inlet_fluid_surface;
% febio_spec.Loads.surface_load{4}.beta.ATTR.lc = 3;
% febio_spec.Loads.surface_load{4}.beta.VAL  = 1;
% 
% febio_spec.Loads.surface_load{5}.ATTR.name = 'outlet-pressure';
% febio_spec.Loads.surface_load{5}.ATTR.type = 'fluid pressure';
% febio_spec.Loads.surface_load{5}.ATTR.surface = outlet_fluid_surface;
% febio_spec.Loads.surface_load{5}.pressure.ATTR.lc = 4;
% febio_spec.Loads.surface_load{5}.pressure.VAL  = 0;
% 
% febio_spec.Loads.surface_load{6}.ATTR.name = 'outlet-backflow-stabilization';
% febio_spec.Loads.surface_load{6}.ATTR.type = 'fluid backflow stabilization';
% febio_spec.Loads.surface_load{6}.ATTR.surface = outlet_fluid_surface;
% febio_spec.Loads.surface_load{6}.beta.ATTR.lc = 5;
% febio_spec.Loads.surface_load{6}.beta.VAL  = 1;
% 
% febio_spec.Loads.surface_load{7}.ATTR.name = 'outlet-tangential-stabilization';
% febio_spec.Loads.surface_load{7}.ATTR.type = 'fluid tangential stabilization';
% febio_spec.Loads.surface_load{7}.ATTR.surface = outlet_fluid_surface;
% febio_spec.Loads.surface_load{7}.beta.ATTR.lc = 6;
% febio_spec.Loads.surface_load{7}.beta.VAL= 1;


%% Load data
febio_spec.LoadData.load_controller{1}.ATTR.id=1;
febio_spec.LoadData.load_controller{1}.ATTR.type='loadcurve';
febio_spec.LoadData.load_controller{1}.interpolate='SMOOTH';
febio_spec.LoadData.load_controller{1}.extend='REPEAT';
febio_spec.LoadData.load_controller{1}.points.point.VAL=[t; (sin(2*pi*t/(bpm/60) - pi/2) + 1)/2]';

febio_spec.LoadData.load_controller{2}.ATTR.id=2;
febio_spec.LoadData.load_controller{2}.ATTR.type='loadcurve';
febio_spec.LoadData.load_controller{2}.interpolate='STEP';
febio_spec.LoadData.load_controller{2}.points.point.VAL=[0 0; 1 1];

febio_spec.LoadData.load_controller{3}.ATTR.id=3;
febio_spec.LoadData.load_controller{3}.ATTR.type='loadcurve';
febio_spec.LoadData.load_controller{3}.interpolate='STEP';
febio_spec.LoadData.load_controller{3}.points.point.VAL=[0 0; 1 1];

febio_spec.LoadData.load_controller{4}.ATTR.id=4;
febio_spec.LoadData.load_controller{4}.ATTR.type='loadcurve';
febio_spec.LoadData.load_controller{4}.interpolate='SMOOTH';
febio_spec.LoadData.load_controller{4}.points.point.VAL=[0 0; 1 1];

febio_spec.LoadData.load_controller{5}.ATTR.id=5;
febio_spec.LoadData.load_controller{5}.ATTR.type='loadcurve';
febio_spec.LoadData.load_controller{5}.interpolate='STEP';
febio_spec.LoadData.load_controller{5}.points.point.VAL=[0 0; 1 1];

febio_spec.LoadData.load_controller{6}.ATTR.id=6;
febio_spec.LoadData.load_controller{6}.ATTR.type='loadcurve';
febio_spec.LoadData.load_controller{6}.interpolate='STEP';
febio_spec.LoadData.load_controller{6}.points.point.VAL=[0 0; 1 1];
% TODO CHeck capital letters


febioStruct2xml(febio_spec,febioFebFileName); %Exporting to file and domNode

%% Running the FEBio analysis
% To run the analysis defined by the created FEBio input file the
% |runMonitorFEBio| function is used. The input for this function is a
% structure defining job settings e.g. the FEBio input file name. The
% optional output runFlag informs the user if the analysis was run
% succesfully. 

febioAnalysis.run_filename=febioFebFileName; %The input file name
febioAnalysis.run_logname=febioLogFileName; %The name for the log file
febioAnalysis.disp_on=1; %Display information on the command window
febioAnalysis.runMode='internal';%'external' or 'internal';
febioAnalysis.maxLogCheckTime=10; %Max log file checking time

[runFlag]=runMonitorFEBio(febioAnalysis);%START FEBio NOW!!!!!!!!


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
