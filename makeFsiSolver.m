%% Control parameters
clear; close all; clc;

t_end = 5;
bpm = 60;
t = linspace(0, t_end, 1000);

n_rand_els = 10;
n_rand_nodes = 10;

V = rand([10,3]);

E_solid = randi([0 n_rand_nodes],n_rand_els,4);
E_fluid = randi([0 n_rand_nodes],n_rand_els,4);

F_i_solid = randi([0 n_rand_nodes],n_rand_els,3);
F_o_solid = randi([0 n_rand_nodes],n_rand_els,3);

F_i_fluid = randi([0 n_rand_nodes],n_rand_els,3);
F_o_fluid = randi([0 n_rand_nodes],n_rand_els,3);

F_FSI = randi([0 n_rand_nodes],n_rand_els, 4);

% Path names
% defaultFolder = fileparts(fileparts(mfilename('C:\Users\alhei\dev\gibbon-hands-on')));

defaultFolder = fileparts(fileparts(mfilename('fullpath'))); 
savePath=fullfile(defaultFolder,'data','temp');

% Defining file names
febioFebFileNamePart='tempModel';
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


%% FSI Solver

%Get a template with default settings 
[febio_spec]=febioStructTemplate;

%febio_spec version 
febio_spec.ATTR.version='3.0'; 

%Module section
febio_spec.Module.ATTR.type='fluid-FSI'; 

febio_spec.Control.analysis='Dynamic'; 
febio_spec.Control.time_steps=1550;
febio_spec.Control.step_size=0.001; 

febio_spec.Control.solver.max_refs = 5;
febio_spec.Control.solver.max_ups = 50;
febio_spec.Control.solver.diverge_reform = 0;
febio_spec.Control.solver.reform_each_time_step = 0;
febio_spec.Control.solver.dtol = 0.001;
febio_spec.Control.solver.vtol = 0.001;
febio_spec.Control.solver.ftol = 0.001;
febio_spec.Control.solver.etol = 0.01;
febio_spec.Control.solver.rtol = 0.001;
febio_spec.Control.solver.lstol = 0.9;
febio_spec.Control.solver.min_residual = 1.e-16;
febio_spec.Control.solver.max_residual = 1.e+10;
febio_spec.Control.solver.rhoi = 0;
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

% Fluid
febio_spec.Material.material{2}.ATTR.name='blood';
febio_spec.Material.material{2}.ATTR.type='fluid-FSI';
febio_spec.Material.material{2}.ATTR.id=2;
febio_spec.Material.material{2}.fluid.ATTR.type = 'fluid';
febio_spec.Material.material{2}.fluid.density = 1060;
febio_spec.Material.material{2}.fluid.k = 2.2e9;
febio_spec.Material.material{2}.fluid.viscous.ATTR.type = 'Carreau';
febio_spec.Material.material{2}.fluid.viscous.mu0 = 0.056;
febio_spec.Material.material{2}.fluid.viscous.mui = 0.00345;
febio_spec.Material.material{2}.fluid.viscous.lambda = 3.313;
febio_spec.Material.material{2}.fluid.viscous.n = 0.3568;

%% Mesh

febio_spec.Mesh.Nodes{1}.ATTR.name = 'nodes';
febio_spec.Mesh.Nodes{1}.node.ATTR.id=(1:size(V,1))'; %The node id's
febio_spec.Mesh.Nodes{1}.node.VAL=V; %The nodel coordinates

% TODO add nodes

% Elemenets for artery blood flow
febio_spec.Mesh.Elements{1}.ATTR.type = 'tet4';
febio_spec.Mesh.Elements{1}.ATTR.name = 'fluid';
febio_spec.Mesh.Elements{1}.elem.ATTR.id=(1:size(E_fluid,1))'; 
febio_spec.Mesh.Elements{1}.elem.VAL=E_fluid; 

% TODO add elemenets

% Elemenets for artery wall
febio_spec.Mesh.Elements{2}.ATTR.type = 'tet4';
febio_spec.Mesh.Elements{2}.ATTR.name = 'solid';
febio_spec.Mesh.Elements{2}.elem.ATTR.id=(1:size(E_solid,1))'; 
febio_spec.Mesh.Elements{2}.elem.VAL=E_solid; 
% TODO add elemenets

% Surfaces
febio_spec.Mesh.Surface{1}.ATTR.name = inlet_solid_surface;
febio_spec.Mesh.Surface{1}.tri3.ATTR.id = (1:size(F_i_solid,1))';
febio_spec.Mesh.Surface{1}.tri3.VAL = F_i_solid;

febio_spec.Mesh.Surface{2}.ATTR.name = inlet_fluid_surface;
febio_spec.Mesh.Surface{2}.tri3.ATTR.id = (1:size(F_i_fluid,1))';
febio_spec.Mesh.Surface{2}.tri3.VAL = F_i_fluid;

febio_spec.Mesh.Surface{3}.ATTR.name = outlet_solid_surface;
febio_spec.Mesh.Surface{3}.tri3.ATTR.id = (1:size(F_o_solid,1))';
febio_spec.Mesh.Surface{3}.tri3.VAL = F_o_solid;

febio_spec.Mesh.Surface{4}.ATTR.name = outlet_fluid_surface;
febio_spec.Mesh.Surface{4}.tri3.ATTR.id = (1:size(F_o_fluid,1))';
febio_spec.Mesh.Surface{4}.tri3.VAL = F_o_fluid;


febio_spec.Mesh.Surface{5}.ATTR.name = FSI_surface;
febio_spec.Mesh.Surface{5}.tri3.ATTR.id = (1:size(F_FSI,1))';
febio_spec.Mesh.Surface{5}.tri3.VAL = F_FSI;

% Mesh domains
febio_spec.MeshDomains.SolidDomain{1}.ATTR.name = 'solid';
febio_spec.MeshDomains.SolidDomain{1}.ATTR.mat = 'artery-wall';

febio_spec.MeshDomains.SolidDomain{2}.ATTR.name = 'fluid';
febio_spec.MeshDomains.SolidDomain{2}.ATTR.mat = 'blood';
% TODO Add mesh domains

%% Boundary conditions

febio_spec.Boundary.bc{1}.ATTR.name = 'fix-inlet-solid';
febio_spec.Boundary.bc{1}.ATTR.type = 'fix';
febio_spec.Boundary.bc{1}.ATTR.node_set = strcat('@surface:',inlet_solid_surface);
febio_spec.Boundary.bc{1}.dofs = 'x, y, z';


febio_spec.Boundary.bc{2}.ATTR.name = 'fix-outlet-solid';
febio_spec.Boundary.bc{2}.ATTR.type = 'fix';
febio_spec.Boundary.bc{2}.ATTR.node_set = strcat('@surface:',outlet_solid_surface);
febio_spec.Boundary.bc{2}.dofs = 'x, y, z';


%% Loads
febio_spec.Loads.surface_load{1}.ATTR.name = 'FSIInterfaceTraction';
febio_spec.Loads.surface_load{1}.ATTR.type = 'fluid-FSI traction';
febio_spec.Loads.surface_load{1}.ATTR.surface = FSI_surface;

febio_spec.Loads.surface_load{2}.ATTR.name = 'inlet-pressure';
febio_spec.Loads.surface_load{2}.ATTR.type = 'fluid pressure';
febio_spec.Loads.surface_load{2}.ATTR.name = inlet_fluid_surface;
febio_spec.Loads.surface_load{2}.pressure.ATTR.lc = '1';
febio_spec.Loads.surface_load{2}.pressure.VAL = 500;

febio_spec.Loads.surface_load{3}.ATTR.name = 'inlet-backflow-stabilization';
febio_spec.Loads.surface_load{3}.ATTR.type = 'fluid backflow stabilization';
febio_spec.Loads.surface_load{3}.ATTR.name = inlet_fluid_surface;
febio_spec.Loads.surface_load{3}.beta.ATTR.lc = '2';
febio_spec.Loads.surface_load{3}.beta.VAL  = 1;

febio_spec.Loads.surface_load{4}.ATTR.name = 'inlet-tangential-stabilization';
febio_spec.Loads.surface_load{4}.ATTR.type = 'fluid tangential stabilization';
febio_spec.Loads.surface_load{4}.ATTR.name = inlet_fluid_surface;
febio_spec.Loads.surface_load{4}.beta.ATTR.lc = '3';
febio_spec.Loads.surface_load{4}.beta.VAL  = 1;

febio_spec.Loads.surface_load{5}.ATTR.name = 'outlet-pressure';
febio_spec.Loads.surface_load{5}.ATTR.type = 'fluid pressure';
febio_spec.Loads.surface_load{5}.ATTR.name = outlet_fluid_surface;
febio_spec.Loads.surface_load{5}.pressure.ATTR.lc = '4';
febio_spec.Loads.surface_load{5}.pressure.VAL  = 0;

febio_spec.Loads.surface_load{6}.ATTR.name = 'outlet-backflow-stabilization';
febio_spec.Loads.surface_load{6}.ATTR.type = 'fluid backflow stabilization';
febio_spec.Loads.surface_load{6}.ATTR.name = outlet_fluid_surface;
febio_spec.Loads.surface_load{6}.beta.ATTR.lc = '5';
febio_spec.Loads.surface_load{6}.beta.VAL  = 1;

febio_spec.Loads.surface_load{7}.ATTR.name = 'outlet-tangential-stabilization';
febio_spec.Loads.surface_load{7}.ATTR.type = 'fluid tangential stabilization';
febio_spec.Loads.surface_load{7}.ATTR.name = outlet_fluid_surface;
febio_spec.Loads.surface_load{7}.beta.ATTR.lc = '6';
febio_spec.Loads.surface_load{7}.beta.VAL= 1;


%% Load data
febio_spec.LoadData.load_controller{1}.ATTR.id='1';
febio_spec.LoadData.load_controller{1}.ATTR.type='loadcurve';
febio_spec.LoadData.load_controller{1}.interpolate='SMOOTH';
febio_spec.LoadData.load_controller{1}.extend='REPEAT';
febio_spec.LoadData.load_controller{1}.points.point.VAL=[t; (sin(2*pi*t/(bpm/60) - pi/2) + 1)/2]';

febio_spec.LoadData.load_controller{2}.ATTR.id='2';
febio_spec.LoadData.load_controller{2}.ATTR.type='loadcurve';
febio_spec.LoadData.load_controller{2}.interpolate='STEP';
febio_spec.LoadData.load_controller{2}.points.point.VAL=[0 0; 1 1];

febio_spec.LoadData.load_controller{3}.ATTR.id='3';
febio_spec.LoadData.load_controller{3}.ATTR.type='loadcurve';
febio_spec.LoadData.load_controller{3}.interpolate='STEP';
febio_spec.LoadData.load_controller{3}.points.point.VAL=[0 0; 1 1];

febio_spec.LoadData.load_controller{4}.ATTR.id='4';
febio_spec.LoadData.load_controller{4}.ATTR.type='loadcurve';
febio_spec.LoadData.load_controller{4}.interpolate='SMOOTH';
febio_spec.LoadData.load_controller{4}.points.point.VAL=[0 0; 1 1];

febio_spec.LoadData.load_controller{5}.ATTR.id='5';
febio_spec.LoadData.load_controller{5}.ATTR.type='loadcurve';
febio_spec.LoadData.load_controller{5}.interpolate='STEP';
febio_spec.LoadData.load_controller{5}.points.point.VAL=[0 0; 1 1];

febio_spec.LoadData.load_controller{6}.ATTR.id='6';
febio_spec.LoadData.load_controller{6}.ATTR.type='loadcurve';
febio_spec.LoadData.load_controller{6}.interpolate='STEP';
febio_spec.LoadData.load_controller{6}.points.point.VAL=[0 0; 1 1];


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
febioAnalysis.runMode='external';%'internal';
febioAnalysis.maxLogCheckTime=10; %Max log file checking time

[runFlag]=runMonitorFEBio(febioAnalysis);%START FEBio NOW!!!!!!!!

