clear; close all; clc;

%%
% Plot settings
fontSize=10;
faceAlpha1=1;
faceAlpha2=0.3;
cmap=gjet(250);

%% Control parameters

if ispc
    dataDir='D:\Experimental_Data\NUIG\INSIST\2019_03_26';
else
    dataDir='/mnt/data/Experimental_Data/NUIG/INSIST/2019_03_26';
end
dirContentStruct = dir(fullfile(dataDir));

%Remove hidden folders containing dot in path
logicKeep=false(1,numel(dirContentStruct));
for q=1:1:numel(dirContentStruct)
    if isempty(strfind(dirContentStruct(q).name,'.'))
        logicKeep(q)=1;
    end
end
dirContentStruct=dirContentStruct(logicKeep);

defaultFolder = fileparts(mfilename('fullpath'));
savePath=fullfile(defaultFolder,'data','temp');

%%
% pressure_offsets=[8000 12000];
% velocityVal=[0.4 0.6];
%%

anisoType=1;
includeClot=0;

for qSim=1%:1:2
    close all; 

    % Defining file names
    switch includeClot
        case 0
            febioFebFileNamePart=['vesselModel_aniso_circum_noclot_aniso',num2str(qSim)];
        case 1
            febioFebFileNamePart=['vesselModel_aniso_circum_clot_aniso',num2str(qSim)];
    end
    
febioFebFileName=fullfile(savePath,[febioFebFileNamePart,'.feb']); %FEB file name
febioLogFileName=fullfile(savePath,[febioFebFileNamePart,'.txt']); %FEBio log file name
febioLogFileName_disp=[febioFebFileNamePart,'_disp_out.txt']; %Log file name for exporting displacement

%Path and file names for load curve data
loadNameLoadCurve='/mnt/data/MATLAB/Projects/vessel_segment/data/Cornelissen_2018_data1.mat';



%%

pressure_offset=10000;%pressure_offsets(qSim);

velocityData.min=0.1;
velocityData.max=0.5;%velocityVal(qSim);
numWaves=1;

%%

%Material parameter set

%Bulk modulus factors
k_factor_vessel=1e2; 
k_factor_clot=1e2;

%--> Vessel
c1_vessel=200000;%1000000; %Shear-modulus-like parameter
m1_vessel=2; %Material parameter setting degree of non-linearity
k_vessel=c1_vessel*k_factor_vessel; %Bulk modulus
ksi=c1_vessel/8;
alphaPar=2;
beta=2;

numMaterialsEnd=50;
enhancementFactor=100;
c1_vessel_range=linspace(c1_vessel,c1_vessel*enhancementFactor,numMaterialsEnd);
k_vessel_range=c1_vessel_range*k_factor_vessel;
materialIndexFluid=numMaterialsEnd+1;

%--> Fluid, fluid part
density_fluid  = 1060;
k_fluid        = 2.2e9;
mu0_Carreau    = 0.056;
mui_Carreau    = 0.00345;
lambda_Carreau = 3.313;
n_Carreau      = 0.3568;

%--> Fluid, solid part
density_neo=0;
E_neo=1e-9;
v_neo=0;

% Clot
clotType=1;
switch clotType
    case 1
        c1_clot=7.2689377051225247e2;
        m1_clot=1.9925054551032964e+00;
    case 2
        c1_clot=8.0355113780266512e2;
        m1_clot=2.1201008264322891e+00;
    case 3
        c1_clot=2.8099638619828640e2;
        m1_clot=2.3852821240150659e+00;
end
k_clot=c1_clot*k_factor_clot; %Bulk modulus


constraintType=1;

% FEA control settings
step_size_desired=0.004; %0.003;
numOutput=200;
max_refs=5;
max_ups=50;
diverge_reform=0;
reform_each_time_step=0;
dtol=0.001;
vtol=0.001;
ftol=0.001;
etol=0.01;
rtol=0.001;
lstol=0.9;
min_residual=1e-16;
max_residual=1e+10;
rhoi=0;
qnmethod=1;
max_retries=5;
opt_iter=53;
analysis_type='dynamic';
symmetric_stiffness=0;

%%
%Parse each valid folder
% [3 5 6 7 8 18 28];
for qf=18%[3 5 6 7 8 18 28]; %1:1:numel(dirContentStruct)
    
    close all;
    
    %%
    
    dirNameNow=dirContentStruct(qf).name;
    
    %% Import surface model data
    if includeClot==1
        loadNameSurf=fullfile(dataDir,dirNameNow,['modelGeometryLocal',num2str(includeClot),'.mat']);
        dataStructOut=load(loadNameSurf);
    else
        try
            loadNameSurf=fullfile(dataDir,dirNameNow,['modelGeometryLocal',num2str(includeClot),'.mat']);
            dataStructOut=load(loadNameSurf);
        catch
            loadNameSurf=fullfile(dataDir,dirNameNow,'modelGeometryLocal.mat');
            dataStructOut=load(loadNameSurf);
        end
    end
    
    
    
    %% Access data
    
    distExtend=dataStructOut.distExtend/1000;
    E_fluid=dataStructOut.E_fluid;
    V=dataStructOut.V;
    Fb_fluid=dataStructOut.Fb_fluid;
    Cb_fluid=dataStructOut.Cb_fluid;
    E_vessel=dataStructOut.E_vessel;
    F_vessel_all=dataStructOut.F_vessel_all;
    F_vessel_inner=dataStructOut.F_vessel_inner;
    F_vessel_outer=dataStructOut.F_vessel_outer;
    F_vessel_sides=dataStructOut.F_vessel_sides;
    E_levelset_cut=dataStructOut.E_levelset_cut;
    V_levelset=dataStructOut.V_levelset;
    if includeClot==1
        E_clot=dataStructOut.E_clot;
        Fb_clot=dataStructOut.Fb_clot;
        Cb_clot=dataStructOut.Cb_clot;
        F_vessel_inner_fluid=dataStructOut.F_vessel_inner_fluid;
    end
    
    %% Convert units to meters
    V=V/1000;
    V_levelset=V_levelset/1000;
    
    %%
    % Visualize model geometry
    
    cFigure; hold on;
    gpatch(Fb_fluid,V,Cb_fluid,'k',1);
    gpatch(F_vessel_all,V,'bw','k',1);
    if includeClot==1
        gpatch(Fb_clot,V,Cb_clot,'k',1);
    end
    axisGeom;
    camlight headlight;
    colormap(cmap); icolorbar;
    drawnow;

    %%
    
    F_inlet=fliplr(Fb_fluid(Cb_fluid==2,:));
    F_outlet1=fliplr(Fb_fluid(Cb_fluid==3,:));
    if includeClot==1        
        F_fsi=[fliplr(F_vessel_inner_fluid); fliplr(Fb_clot(Cb_clot==3,:))];
    else
        F_outlet2=fliplr(Fb_fluid(Cb_fluid==4,:));
        F_fsi=fliplr(F_vessel_inner);
    end
    
    %%
    cFigure; hold on;
    gpatch(F_fsi,V,'bw','k',1);
    patchNormPlot(F_fsi,V);
    
    axisGeom;
    camlight headlight;    
    drawnow;

    %%
    % Visualize surface
    
    cFigure; hold on;
    %     gpatch(F_vessel_all,V,'kw','none',0.1);
    gpatch(Fb_fluid,V,'kw','r',0.5);
    gpatch(F_inlet,V,'rw','k',1);
    patchNormPlot(F_inlet,V);
    gpatch(F_outlet1,V,'gw','k',1);
    patchNormPlot(F_outlet1,V);
    if includeClot~=1
        gpatch(F_outlet2,V,'bw','k',1);
        patchNormPlot(F_outlet2,V);
    end
    gpatch(F_vessel_inner,V,'kw','g',0.5);
    patchNormPlot(F_vessel_inner,V);
    axisGeom;
    camlight headlight;
    drawnow;
        
    %% Rotate model so inlet is origin
    
    indInlet=unique(F_inlet(:));
    mean_V=mean(V(indInlet,:),1);
    V=V-mean_V;
    N1=-vecnormalize(mean(patchNormal(F_inlet,V),1));
    N2p=[0 1 0];
    N3=vecnormalize(cross(N1,N2p));
    N2=vecnormalize(cross(N3,N1));
    
    R=[N1(:) N2(:) N3(:)];
    V=V*R;
    
    V_levelset=V_levelset-mean_V;
    V_levelset=V_levelset*R;
    
    %%
    cFigure; hold on;
    gpatch(Fb_fluid,V,'kw','r',0.5);
    gpatch(F_inlet,V,'rw','k',1);
    patchNormPlot(F_inlet,V);
    gpatch(F_outlet1,V,'gw','k',1);
    patchNormPlot(F_outlet1,V);
    if includeClot~=1
        gpatch(F_outlet2,V,'bw','k',1);
        patchNormPlot(F_outlet2,V);
    end
    gpatch(F_vessel_inner,V,'kw','g',0.5);
    patchNormPlot(F_vessel_inner,V);
    quiverVec(zeros(1,3),N1,0.01,'r');
    quiverVec(zeros(1,3),N2,0.01,'g');
    quiverVec(zeros(1,3),N3,0.01,'b');
    axisGeom;
    camlight headlight;
    drawnow;
    
    %%
    
    clear optionStruct
    optionStruct.outputType='label';
    G=tesgroup(F_vessel_sides,optionStruct);
    Vc=patchCentre(F_vessel_sides,V);
    
    d=nan(1,max(G(:)));
    for q=1:1:max(G(:))
        d(q)=mean(sqrt(sum(Vc(G==q,:).^2,2)));
    end
    [~,inletLabel]=min(d);
    
    %%
    
    cFigure; hold on;
    gpatch(F_vessel_sides,V,G,'k',1);
    gpatch(F_vessel_inner,V,'kw','none',0.5);
    axisGeom;
    camlight headlight;
    colormap(gjet(250)); icolorbar;
    drawnow;

    
    %%
    switch constraintType
        case 1
            
            if includeClot==1
%                 f=F_vessel_sides(G~=2,:); %Exclude clot side
                f=F_vessel_sides(:); %All sides
                bcSupportList=unique([F_inlet(:); F_outlet1(:); f(:)]);
            else
                bcSupportList=unique([F_inlet(:); F_outlet1(:); F_outlet2(:); F_vessel_sides(:)]);
            end
        case 2
            if includeClot==1
                f=F_vessel_sides(G~=inletLabel,:);
                bcSupportList=unique([F_outlet1(:); f(:)]);
                
                f=F_vessel_sides(G==inletLabel,:);
                bcSupport_X=unique([F_inlet(:);f(:)]);
            else
                f=F_vessel_sides(G~=inletLabel,:);
                bcSupportList=unique([F_outlet1(:); F_outlet2(:); f(:)]);
                
                f=F_vessel_sides(G==inletLabel,:);
                bcSupport_X=unique([F_inlet(:);f(:)]);
            end
    end
    
    bcFluidFix=unique(F_fsi(:));
    bcDilatationPrescribe=unique(patchBoundary(F_inlet,V));
    
    %%
    % Visualize surface
    
    cFigure; hold on;
    gpatch(F_vessel_all,V,'kw','none',0.1);
    gpatch(Fb_fluid,V,'kw','none',0.1);
    hp(1)=plotV(V(bcSupportList,:),'k.','MarkerSize',35);
    hp(2)=plotV(V(bcFluidFix,:),'g.','MarkerSize',15);
    hp(3)=plotV(V(bcDilatationPrescribe,:),'r.','MarkerSize',40);
    switch constraintType
        case 1
            legend(hp,{'Supported nodes','Fixed fluid velocity','ef'});
        case 2
            hp(4)=plotV(V(bcSupport_X,:),'c.','MarkerSize',35);
            legend(hp,{'Supported nodes','Fixed fluid velocity','ef','X support'});
    end
    axisGeom;
    camlight headlight;
    drawnow;
    
    %% Propagate distance from inlet to inform spatially varying stiffness
    
    indStart=bcDilatationPrescribe; %Index of the start point
    optionStruct.toleranceLevel=0; %Tolerance for convergence
    optionStruct.waitBarOn=0; %Turn on/off waitbar
    
    %Compute distances on mesh description
    distMarch=meshDistMarch(F_vessel_inner,V,indStart,optionStruct);
    distMarch(isnan(distMarch))=max(distMarch);
    distMarch=vertexToFaceMeasure(F_vessel_inner,distMarch);
    
    materialIndex=distMarch;
    materialIndex=materialIndex-min(materialIndex(:));
    materialIndex=materialIndex./distExtend;
    materialIndex(materialIndex>1)=1;
    materialIndex=abs(materialIndex-1);
    materialIndex=round((materialIndex.*(numMaterialsEnd-1))+1);
    
    %%
    cFigure; hold on;
    gpatch(F_vessel_inner,V,materialIndex,'none',1);
%     plotV(V(bcDilatationPrescribe,:),'r.','MarkerSize',40);
    axisGeom(gca,75);
    camlight headlight;
    colormap(gjet(50)); h=colorbar; 
    drawnow;
    
    %%
    VM=patchCentre(F_vessel_inner,V);
    VE=patchCentre(E_vessel,V);
    [~,indMin]=minDist(VE,VM);
    materialIndexElements=materialIndex(indMin);
    
    %     %% Visualize
    %     [Fp,Cp,CFp]=element2patch(E_vessel,materialIndexElements,'penta6');
    %
    %     cFigure; hold on;
    %     for q=1:1:numel(Fp)
    %         gpatch(Fp{q},V,Cp{q},'none',1);
    %     end
    %     axisGeom;
    %     camlight headlight;
    %     colormap(gjet(250)); colorbar;
    %     drawnow;
    
    %%
        
    VP=patchCentre(E_levelset_cut,V_levelset);    
    N_fib_axial_levelset=V_levelset(E_levelset_cut(:,2),:)-V_levelset(E_levelset_cut(:,1),:);

    [~,indMinLevelSet]=minDist(VE,VP);
        
    N_fib_axial=vecnormalize(N_fib_axial_levelset(indMinLevelSet,:));
    D_fib=vecnormalize(VP(indMinLevelSet,:)-VE);
    N_fib_circum=vecnormalize(cross(N_fib_axial,D_fib));
    D_fib=vecnormalize(cross(N_fib_axial,N_fib_circum));
    
    cFigure; hold on;    
    gpatch(F_vessel_all,V,'kw','none',0.1);
    gpatch(E_levelset_cut,V_levelset,'none','k',1,3);
%     quiverVec(VE,N_fib_axial,a);
%     quiverVec(VP,N_fib_axial_levelset,0.001,'r');
%     quiverVec(VE,N_fib_axial,0.001,'r');
%     quiverVec(VE,D,0.001,'g');
    quiverVec(VE,N_fib_circum,0.0005,'b');
    axisGeom;
    camlight headlight;    
    drawnow;
    
    %% set-up load curves
    
    % Load load-curve data
    loadCurveData=load(loadNameLoadCurve);
    t=loadCurveData.x; %Time
    a=loadCurveData.y; %Amplitude
    
    tDevelopStart=(max(distMarch)./velocityData.min);
    tWaitStart=0.5*tDevelopStart;
    tWaitIntermediate=(max(t(:))-min(t(:)))/6;
    
    sigmoidRampDurationVelocity=tDevelopStart/2;
    sigmoidRampDuration_P0=sigmoidRampDurationVelocity+(tDevelopStart-sigmoidRampDurationVelocity)/2;
    
    w=velocityData.max-velocityData.min;
    
    %in-vivo part
    
    a=a-min(a(:));
    a=a./max(a(:));
    a=a.*w;
    a=a+velocityData.min;
    
    %Sigmoid part
    s=10;
    T=linspace(-sigmoidRampDurationVelocity/2,sigmoidRampDurationVelocity/2,100);
    
    S=tanh(T.*s);
    S=S-min(S(:));
    S=S./max(S(:));
    S=S.*velocityData.min;
    S(end+1)=S(end);
    T=T-T(1);
    T(end+1)=T(end)+tWaitStart;
    
    t=t+T(end);
    t(end+1)=t(end)+tWaitIntermediate;
    a(end+1)=a(end);
    
    TT=[T(:);];
    AA=[S(:);];
    
    for q=1:1:numWaves
        TT=[TT;t(:)-min(t(:))+max(TT(:))];
        AA=[AA;a(:)];
    end
    
    [TT,indUni]=unique(TT);
    AA=AA(indUni);
    loadCurve_velocity=[TT(:) AA(:)];
    
    dt=mean(diff(loadCurve_velocity(:,1)));
    n=ceil(max(loadCurve_velocity(:,1))/dt);
    loadCurve_velocity=evenlySampleCurve(loadCurve_velocity,n,'pchip',0);
    
    %Sigmoid part
    s=15;
    T=linspace(-sigmoidRampDuration_P0/2,sigmoidRampDuration_P0/2,100);
    
    S=tanh(T.*s);
    S=S-min(S(:));
    S=S./max(S(:));
    %     S=S.*velocityData.min;
    %     S(end+1)=S(end);
    T=T-T(1);
    %     T(end+1)=T(end)+tWaitStart;
    
    dt=mean(diff(loadCurve_velocity(:,1)));
    n=ceil(max(T(:))/dt);
    loadCurve_other=[T(:) S(:)];
    loadCurve_other=evenlySampleCurve(loadCurve_other,n,'pchip',0);
    
    %%
    
    cFigure;
    subplot(1,3,1);hold on;
    plot(loadCurve_velocity(:,1),loadCurve_velocity(:,2),'r.-','MarkerSize',15,'LineWidth',2);
    set(gca,'FontSize',35);
    grid on; box on; axis square;
    set(gca,'XTick',0:0.5:max(loadCurve_velocity(:,1)));
    
    subplot(1,3,2);hold on;
    plot(loadCurve_other(:,1),loadCurve_other(:,2),'g.-','MarkerSize',15,'LineWidth',2);
    set(gca,'FontSize',35);
    grid on; box on; axis square;
    set(gca,'XTick',0:0.1:max(loadCurve_other(:,1)));
    
    subplot(1,3,3);hold on;
    plot(loadCurve_velocity(:,1),loadCurve_velocity(:,2),'r.-','MarkerSize',15,'LineWidth',2);
    plot(loadCurve_other(:,1),loadCurve_other(:,2).*velocityData.min,'g.-','MarkerSize',15,'LineWidth',2);
    set(gca,'FontSize',35);
    grid on; box on; axis square;
    set(gca,'XTick',0:0.5:max(loadCurve_velocity(:,1)));
    
    drawnow;

    %%
    
    cFigure; hold on;
    title('Inlet fluid velocity profile');
    plot(loadCurve_velocity(:,1),loadCurve_velocity(:,2),'r-','LineWidth',5);
    set(gca,'FontSize',35);
    grid on; box on; %axis square;
    axis tight; 
    set(gca,'XTick',0:0.25:max(loadCurve_velocity(:,1)));
    set(gca,'YTick',round(100*linspace(0,max(loadCurve_velocity(:,2)),11))/100);    
    ylim([0 0.5]);
    xlabel('Time [s]');
    ylabel('Velocity [m/s]');    
    ha=gca;
    ha.LineWidth=2;
    drawnow; 
    
    %%
    
    cFigure;
    hold on;
    title('Inlet fluid velocity profile');
    xlabel('Time [s]'); ylabel('Velocity [m/s]')
    plot(loadCurve_velocity(:,1),loadCurve_velocity(:,2),'r-','LineWidth',4);
    set(gca,'FontSize',35);
    grid on; box on; axis tight;
    ylim([0 0.5])
    set(gca,'YTick',0:0.05:0.5);
    set(gca,'XTick',0:0.5:max(loadCurve_velocity(:,1)));
    drawnow; 
    
    %% Define time stepping
    
    tTotal=max(loadCurve_velocity(:,1));
    time_steps=round(tTotal/step_size_desired);
    step_size=tTotal/time_steps;
    dtmin=step_size/100;
    dtmax=step_size;
    plot_stride=round(time_steps/numOutput);
    
    %%
    cFigure; hold on;
    gpatch(F_vessel_all,V,'kw','none',0.1);
%     gpatch(F_inlet,V,'rw','r',0.5);
%     gpatch(F_outlet1,V,'gw','g',0.5);
    gpatch(F_fsi,V,'bw','b',0.5);
    
    axisGeom;
    camlight headlight;    
    drawnow;

    %% Defining the FEBio input structure
    % See also |febioStructTemplate| and |febioStruct2xml| and the FEBio user
    % manual.
    
    %Get a template with default settings
    [febio_spec]=febioStructTemplate_v2p5;
    
    %febio_spec version
    febio_spec.ATTR.version='2.5';
    
    %Module section
    febio_spec.Module.ATTR.type='fluid-FSI';
    
    %Control section
    febio_spec.Control.time_steps=time_steps;
    febio_spec.Control.step_size=step_size;
    febio_spec.Control.max_refs=max_refs;
    febio_spec.Control.max_ups=max_ups;
    febio_spec.Control.diverge_reform=diverge_reform;
    febio_spec.Control.reform_each_time_step=reform_each_time_step;
    febio_spec.Control.dtol=dtol;
    febio_spec.Control.vtol=vtol;
    febio_spec.Control.ftol=ftol;
    febio_spec.Control.etol=etol;
    febio_spec.Control.rtol=rtol;
    febio_spec.Control.lstol=lstol;
    febio_spec.Control.min_residual=min_residual;
    febio_spec.Control.max_residual=max_residual;
    febio_spec.Control.plot_stride=plot_stride;
    febio_spec.Control.rhoi=rhoi;
    febio_spec.Control.qnmethod=qnmethod;
    febio_spec.Control.time_stepper.dtmin=dtmin;
    febio_spec.Control.time_stepper.dtmax=dtmax;
    febio_spec.Control.time_stepper.max_retries=max_retries;
    febio_spec.Control.time_stepper.opt_iter=opt_iter;
    febio_spec.Control.analysis.ATTR.type=analysis_type;
    febio_spec.Control.symmetric_stiffness=symmetric_stiffness;
    
    %Material section
    % -> Material 1: The vessel
    switch anisoType
        case 0
            for q=1:1:numMaterialsEnd
                febio_spec.Material.material{q}.ATTR.name=['Vessel_mat',num2str(q)];
                febio_spec.Material.material{q}.ATTR.type='Ogden unconstrained';  
                febio_spec.Material.material{q}.ATTR.id=q;
                febio_spec.Material.material{q}.c1=c1_vessel_range(q);
                febio_spec.Material.material{q}.m1=m1_vessel;
                febio_spec.Material.material{q}.c2=c1_vessel_range(q);
                febio_spec.Material.material{q}.m2=-m1_vessel;
                febio_spec.Material.material{q}.cp=k_vessel_range(q);
                febio_spec.Material.material{q}.density=1000;
            end
        otherwise
            for q=1:1:numMaterialsEnd
                febio_spec.Material.material{q}.ATTR.name=['Vessel_mat',num2str(q)];
                febio_spec.Material.material{q}.ATTR.type='solid mixture';
                febio_spec.Material.material{q}.ATTR.id=q;                
                febio_spec.Material.material{q}.mat_axis.ATTR.type='user';
                
                febio_spec.Material.material{q}.solid{1}.ATTR.type='Ogden unconstrained';                
                febio_spec.Material.material{q}.solid{1}.c1=c1_vessel_range(q);
                febio_spec.Material.material{q}.solid{1}.m1=m1_vessel;
                febio_spec.Material.material{q}.solid{1}.c2=c1_vessel_range(q);
                febio_spec.Material.material{q}.solid{1}.m2=-m1_vessel;
                febio_spec.Material.material{q}.solid{1}.cp=k_vessel_range(q);
                febio_spec.Material.material{q}.solid{1}.density=1000;
                
                febio_spec.Material.material{q}.solid{2}.ATTR.type='fiber-exp-pow';
                febio_spec.Material.material{q}.solid{2}.ksi=ksi;
                febio_spec.Material.material{q}.solid{2}.alpha=alphaPar;
                febio_spec.Material.material{q}.solid{2}.beta=beta;
                febio_spec.Material.material{q}.solid{2}.theta=0;
                febio_spec.Material.material{q}.solid{2}.phi=0;
            end
    end
    
    % ----> The fluid part
    febio_spec.Material.material{materialIndexFluid}.ATTR.name='Blood';
    febio_spec.Material.material{materialIndexFluid}.ATTR.type='fluid-FSI';
    febio_spec.Material.material{materialIndexFluid}.ATTR.id=materialIndexFluid;
    febio_spec.Material.material{materialIndexFluid}.fluid{1}.ATTR.type='fluid';
    febio_spec.Material.material{materialIndexFluid}.fluid{1}.density=density_fluid;
    febio_spec.Material.material{materialIndexFluid}.fluid{1}.k=k_fluid;
    
    febio_spec.Material.material{materialIndexFluid}.fluid{1}.viscous.ATTR.type='Carreau';
    febio_spec.Material.material{materialIndexFluid}.fluid{1}.viscous.mu0=mu0_Carreau;
    febio_spec.Material.material{materialIndexFluid}.fluid{1}.viscous.mui=mui_Carreau;
    febio_spec.Material.material{materialIndexFluid}.fluid{1}.viscous.lambda=lambda_Carreau;
    febio_spec.Material.material{materialIndexFluid}.fluid{1}.viscous.n=n_Carreau;
    
    % ----> The solid part
    febio_spec.Material.material{materialIndexFluid}.solid{1}.ATTR.type='neo-Hookean';
    febio_spec.Material.material{materialIndexFluid}.solid{1}.E=E_neo;
    febio_spec.Material.material{materialIndexFluid}.solid{1}.v=v_neo;
    febio_spec.Material.material{materialIndexFluid}.solid{1}.density=density_neo;
    
    if includeClot==1
        materialIndexClot=numel(febio_spec.Material.material)+1;       
        febio_spec.Material.material{materialIndexClot}.ATTR.name='Clot';
        febio_spec.Material.material{materialIndexClot}.ATTR.type='Ogden';
        febio_spec.Material.material{materialIndexClot}.ATTR.id=materialIndexClot;
        febio_spec.Material.material{materialIndexClot}.c1=c1_clot;
        febio_spec.Material.material{materialIndexClot}.m1=m1_clot;
        febio_spec.Material.material{materialIndexClot}.c2=c1_clot;
        febio_spec.Material.material{materialIndexClot}.m2=-m1_clot;
        febio_spec.Material.material{materialIndexClot}.k=k_clot;
    end
        
    %Geometry section
    % -> Nodes
    febio_spec.Geometry.Nodes{1}.ATTR.name='nodeSet_all'; %The node set name
    febio_spec.Geometry.Nodes{1}.node.ATTR.id=(1:size(V,1))'; %The node id's
    febio_spec.Geometry.Nodes{1}.node.VAL=V; %The nodel coordinates
    
    % -> Elements
    startInd=1;
    endInd=0;
    for q=1:1:numMaterialsEnd
        
        E_now=E_vessel(materialIndexElements==q,:);
        endInd=endInd+size(E_now,1);
        
        febio_spec.Geometry.Elements{q}.ATTR.type='penta6'; %Element type of this set
        febio_spec.Geometry.Elements{q}.ATTR.mat=q; %material index for this set
        febio_spec.Geometry.Elements{q}.ATTR.name='Vessel'; %Name of the element set
        
        febio_spec.Geometry.Elements{q}.elem.ATTR.id=(startInd:1:endInd)'; %Element id's
        febio_spec.Geometry.Elements{q}.elem.VAL=E_now;
        
        startInd=endInd+1;
    end
    
    febio_spec.Geometry.Elements{materialIndexFluid}.ATTR.type='tet4'; %Element type of this set
    febio_spec.Geometry.Elements{materialIndexFluid}.ATTR.mat=materialIndexFluid; %material index for this set
    febio_spec.Geometry.Elements{materialIndexFluid}.ATTR.name='Blood'; %Name of the element set
    febio_spec.Geometry.Elements{materialIndexFluid}.elem.ATTR.id=(size(E_vessel,1)+1:1:size(E_vessel,1)+size(E_fluid,1))'; %Element id's
    febio_spec.Geometry.Elements{materialIndexFluid}.elem.VAL=E_fluid;
    
    if includeClot==1
        materialIndexClot=materialIndexFluid+1;
        febio_spec.Geometry.Elements{materialIndexClot}.ATTR.type='tet4'; %Element type of this set
        febio_spec.Geometry.Elements{materialIndexClot}.ATTR.mat=materialIndexClot; %material index for this set
        febio_spec.Geometry.Elements{materialIndexClot}.ATTR.name='Clot'; %Name of the element set
        febio_spec.Geometry.Elements{materialIndexClot}.elem.ATTR.id=(size(E_vessel,1)+size(E_fluid,1)+1:1:size(E_vessel,1)+size(E_fluid,1)+size(E_clot,1))'; %Element id's
        febio_spec.Geometry.Elements{materialIndexClot}.elem.VAL=E_clot;
    end
        
    % -> NodeSets
    febio_spec.Geometry.NodeSet{1}.ATTR.name='bcSupportList';
    febio_spec.Geometry.NodeSet{1}.node.ATTR.id=bcSupportList(:);
    
    febio_spec.Geometry.NodeSet{2}.ATTR.name='FixedFluidVelocity';
    febio_spec.Geometry.NodeSet{2}.node.ATTR.id=bcFluidFix(:);
    
    febio_spec.Geometry.NodeSet{3}.ATTR.name='PrescribedFluidDilatation';
    febio_spec.Geometry.NodeSet{3}.node.ATTR.id=bcDilatationPrescribe(:);
    
    switch constraintType
        case 1
        case 2
            febio_spec.Geometry.NodeSet{4}.ATTR.name='bcSupportList_X';
            febio_spec.Geometry.NodeSet{4}.node.ATTR.id=bcSupport_X(:);
    end
    
    switch anisoType
        case 0
        otherwise
            % -> ElementSets
            febio_spec.Geometry.ElementSet{1}.ATTR.name='elementSetTransiso';
            febio_spec.Geometry.ElementSet{1}.elem.ATTR.id=(1:size(E_vessel,1))';
            
            %MeshData section
            % -> ElementData
            febio_spec.MeshData.ElementData{1}.ATTR.elem_set=febio_spec.Geometry.ElementSet{1}.ATTR.name;
            febio_spec.MeshData.ElementData{1}.ATTR.var='mat_axis';
            
            for q=1:1:size(E_vessel,1)
                febio_spec.MeshData.ElementData{1}.elem{q}.ATTR.lid=q;
                switch anisoType
                    case 1
                        febio_spec.MeshData.ElementData{1}.elem{q}.a=N_fib_axial(q,:);
                        febio_spec.MeshData.ElementData{1}.elem{q}.d=D_fib(q,:);
                    case 2
                        febio_spec.MeshData.ElementData{1}.elem{q}.a=N_fib_circum(q,:);
                        febio_spec.MeshData.ElementData{1}.elem{q}.d=D_fib(q,:);
                end
            end
    end
    
    if includeClot==1
        % -> Surfaces
        febio_spec.Geometry.Surface{1}.ATTR.name='Fluid_inlet';
        febio_spec.Geometry.Surface{1}.tri3.ATTR.lid=(1:1:size(F_inlet,1))';
        febio_spec.Geometry.Surface{1}.tri3.VAL=F_inlet;
        
        febio_spec.Geometry.Surface{2}.ATTR.name='Fluid_outlet1';
        febio_spec.Geometry.Surface{2}.tri3.ATTR.lid=(1:1:size(F_outlet1,1))';
        febio_spec.Geometry.Surface{2}.tri3.VAL=F_outlet1;
        
        febio_spec.Geometry.Surface{3}.ATTR.name='FSI';
        febio_spec.Geometry.Surface{3}.tri3.ATTR.lid=(1:1:size(F_fsi,1))';
        febio_spec.Geometry.Surface{3}.tri3.VAL=F_fsi;
        
        %Boundary condition section
        % -> Fix boundary conditions
        febio_spec.Boundary.fix{1}.ATTR.bc='x,y,z';
        febio_spec.Boundary.fix{1}.ATTR.node_set=febio_spec.Geometry.NodeSet{1}.ATTR.name;
        febio_spec.Boundary.fix{2}.ATTR.bc='wx,wy,wz';
        febio_spec.Boundary.fix{2}.ATTR.node_set=febio_spec.Geometry.NodeSet{2}.ATTR.name;
        
        switch constraintType
            case 1
            case 2
                febio_spec.Boundary.fix{3}.ATTR.bc='x';
                febio_spec.Boundary.fix{3}.ATTR.node_set=febio_spec.Geometry.NodeSet{4}.ATTR.name;
        end
        
        febio_spec.Boundary.prescribe{1}.ATTR.bc='ef';
        febio_spec.Boundary.prescribe{1}.ATTR.node_set=febio_spec.Geometry.NodeSet{3}.ATTR.name;
        febio_spec.Boundary.prescribe{1}.scale.ATTR.lc=1;
        febio_spec.Boundary.prescribe{1}.scale.VAL=-4.54545e-06;
        febio_spec.Boundary.prescribe{1}.relative=0;
        
        %Loads
        febio_spec.Loads.surface_load{1}.ATTR.type='fluid normal velocity';
        febio_spec.Loads.surface_load{1}.ATTR.surface=febio_spec.Geometry.Surface{1}.ATTR.name;
        febio_spec.Loads.surface_load{1}.velocity.ATTR.lc=2;
        febio_spec.Loads.surface_load{1}.velocity.VAL=-1;
        febio_spec.Loads.surface_load{1}.prescribe_nodal_velocities=1;
        febio_spec.Loads.surface_load{1}.parabolic=1;
        
        febio_spec.Loads.surface_load{2}.ATTR.type='fluid resistance';
        febio_spec.Loads.surface_load{2}.ATTR.surface=febio_spec.Geometry.Surface{2}.ATTR.name;
        febio_spec.Loads.surface_load{2}.R.ATTR.lc=3;
        febio_spec.Loads.surface_load{2}.R.VAL=4e+08;
        febio_spec.Loads.surface_load{2}.pressure_offset.ATTR.lc=1;
        febio_spec.Loads.surface_load{2}.pressure_offset.VAL=pressure_offset;
        
        febio_spec.Loads.surface_load{3}.ATTR.type='fluid backflow stabilization';
        febio_spec.Loads.surface_load{3}.ATTR.surface=febio_spec.Geometry.Surface{2}.ATTR.name;
        febio_spec.Loads.surface_load{3}.beta.ATTR.lc=3;
        febio_spec.Loads.surface_load{3}.beta.VAL=1;
        
        febio_spec.Loads.surface_load{4}.ATTR.type='fluid tangential stabilization';
        febio_spec.Loads.surface_load{4}.ATTR.surface=febio_spec.Geometry.Surface{2}.ATTR.name;
        febio_spec.Loads.surface_load{4}.beta.ATTR.lc=3;
        febio_spec.Loads.surface_load{4}.beta.VAL=1;
        
        febio_spec.Loads.surface_load{5}.ATTR.type='fluid-FSI traction';
        febio_spec.Loads.surface_load{5}.ATTR.surface=febio_spec.Geometry.Surface{3}.ATTR.name;

    else
        
        % -> Surfaces
        febio_spec.Geometry.Surface{1}.ATTR.name='Fluid_inlet';
        febio_spec.Geometry.Surface{1}.tri3.ATTR.lid=(1:1:size(F_inlet,1))';
        febio_spec.Geometry.Surface{1}.tri3.VAL=F_inlet;
        
        febio_spec.Geometry.Surface{2}.ATTR.name='Fluid_outlet1';
        febio_spec.Geometry.Surface{2}.tri3.ATTR.lid=(1:1:size(F_outlet1,1))';
        febio_spec.Geometry.Surface{2}.tri3.VAL=F_outlet1;
        
        febio_spec.Geometry.Surface{3}.ATTR.name='Fluid_outlet2';
        febio_spec.Geometry.Surface{3}.tri3.ATTR.lid=(1:1:size(F_outlet2,1))';
        febio_spec.Geometry.Surface{3}.tri3.VAL=F_outlet2;
        
        febio_spec.Geometry.Surface{4}.ATTR.name='FSI';
        febio_spec.Geometry.Surface{4}.tri3.ATTR.lid=(1:1:size(F_fsi,1))';
        febio_spec.Geometry.Surface{4}.tri3.VAL=F_fsi;
        
        %Boundary condition section
        % -> Fix boundary conditions
        febio_spec.Boundary.fix{1}.ATTR.bc='x,y,z';
        febio_spec.Boundary.fix{1}.ATTR.node_set=febio_spec.Geometry.NodeSet{1}.ATTR.name;
        febio_spec.Boundary.fix{2}.ATTR.bc='wx,wy,wz';
        febio_spec.Boundary.fix{2}.ATTR.node_set=febio_spec.Geometry.NodeSet{2}.ATTR.name;
        
        switch constraintType
            case 1
            case 2
                febio_spec.Boundary.fix{3}.ATTR.bc='x';
                febio_spec.Boundary.fix{3}.ATTR.node_set=febio_spec.Geometry.NodeSet{4}.ATTR.name;
        end
        
        febio_spec.Boundary.prescribe{1}.ATTR.bc='ef';
        febio_spec.Boundary.prescribe{1}.ATTR.node_set=febio_spec.Geometry.NodeSet{3}.ATTR.name;
        febio_spec.Boundary.prescribe{1}.scale.ATTR.lc=1;
        febio_spec.Boundary.prescribe{1}.scale.VAL=-4.54545e-06;
        febio_spec.Boundary.prescribe{1}.relative=0;
        
        %Loads
        febio_spec.Loads.surface_load{1}.ATTR.type='fluid normal velocity';
        febio_spec.Loads.surface_load{1}.ATTR.surface=febio_spec.Geometry.Surface{1}.ATTR.name;
        febio_spec.Loads.surface_load{1}.velocity.ATTR.lc=2;
        febio_spec.Loads.surface_load{1}.velocity.VAL=-1;
        febio_spec.Loads.surface_load{1}.prescribe_nodal_velocities=1;
        febio_spec.Loads.surface_load{1}.parabolic=1;
        
        febio_spec.Loads.surface_load{2}.ATTR.type='fluid resistance';
        febio_spec.Loads.surface_load{2}.ATTR.surface=febio_spec.Geometry.Surface{2}.ATTR.name;
        febio_spec.Loads.surface_load{2}.R.ATTR.lc=3;
        febio_spec.Loads.surface_load{2}.R.VAL=4e+08;
        febio_spec.Loads.surface_load{2}.pressure_offset.ATTR.lc=1;
        febio_spec.Loads.surface_load{2}.pressure_offset.VAL=pressure_offset;
        
        febio_spec.Loads.surface_load{3}.ATTR.type='fluid resistance';
        febio_spec.Loads.surface_load{3}.ATTR.surface=febio_spec.Geometry.Surface{3}.ATTR.name;
        febio_spec.Loads.surface_load{3}.R.ATTR.lc=3;
        febio_spec.Loads.surface_load{3}.R.VAL=4e+08;
        febio_spec.Loads.surface_load{3}.pressure_offset.ATTR.lc=1;
        febio_spec.Loads.surface_load{3}.pressure_offset.VAL=pressure_offset;
        
        febio_spec.Loads.surface_load{4}.ATTR.type='fluid backflow stabilization';
        febio_spec.Loads.surface_load{4}.ATTR.surface=febio_spec.Geometry.Surface{2}.ATTR.name;
        febio_spec.Loads.surface_load{4}.beta.ATTR.lc=3;
        febio_spec.Loads.surface_load{4}.beta.VAL=1;
        
        febio_spec.Loads.surface_load{5}.ATTR.type='fluid backflow stabilization';
        febio_spec.Loads.surface_load{5}.ATTR.surface=febio_spec.Geometry.Surface{3}.ATTR.name;
        febio_spec.Loads.surface_load{5}.beta.ATTR.lc=3;
        febio_spec.Loads.surface_load{5}.beta.VAL=1;
        
        febio_spec.Loads.surface_load{6}.ATTR.type='fluid tangential stabilization';
        febio_spec.Loads.surface_load{6}.ATTR.surface=febio_spec.Geometry.Surface{2}.ATTR.name;
        febio_spec.Loads.surface_load{6}.beta.ATTR.lc=3;
        febio_spec.Loads.surface_load{6}.beta.VAL=1;
        
        febio_spec.Loads.surface_load{7}.ATTR.type='fluid tangential stabilization';
        febio_spec.Loads.surface_load{7}.ATTR.surface=febio_spec.Geometry.Surface{3}.ATTR.name;
        febio_spec.Loads.surface_load{7}.beta.ATTR.lc=3;
        febio_spec.Loads.surface_load{7}.beta.VAL=1;
        
        febio_spec.Loads.surface_load{8}.ATTR.type='fluid-FSI traction';
        febio_spec.Loads.surface_load{8}.ATTR.surface=febio_spec.Geometry.Surface{4}.ATTR.name;

    end
    
    %LoadData section
    febio_spec.LoadData.loadcurve{1}.ATTR.id=1;
    febio_spec.LoadData.loadcurve{1}.ATTR.type='linear';
    febio_spec.LoadData.loadcurve{1}.point.VAL=loadCurve_other;
    
    febio_spec.LoadData.loadcurve{2}.ATTR.id=2;
    febio_spec.LoadData.loadcurve{2}.ATTR.type='linear';
    febio_spec.LoadData.loadcurve{2}.point.VAL=loadCurve_velocity;
    
    febio_spec.LoadData.loadcurve{3}.ATTR.id=3;
    febio_spec.LoadData.loadcurve{3}.ATTR.type='step';
    febio_spec.LoadData.loadcurve{3}.point.VAL=[0 0; 1 1];
    
    %Output section
    % -> log file
    % febio_spec.Output.logfile.ATTR.file=febioLogFileName;
    % febio_spec.Output.logfile.node_data{1}.ATTR.file=febioLogFileName_disp;
    % febio_spec.Output.logfile.node_data{1}.ATTR.data='ux;uy;uz';
    % febio_spec.Output.logfile.node_data{1}.ATTR.delim=',';
    % febio_spec.Output.logfile.node_data{1}.VAL=1:size(V,1);
    
    febio_spec.Output.plotfile.var{end+1}.ATTR.type='velocity';
    febio_spec.Output.plotfile.var{end+1}.ATTR.type='fluid pressure';
    febio_spec.Output.plotfile.var{end+1}.ATTR.type='fluid density';
    febio_spec.Output.plotfile.var{end+1}.ATTR.type='fluid dilatation';
    febio_spec.Output.plotfile.var{end+1}.ATTR.type='fluid stress';
    febio_spec.Output.plotfile.var{end+1}.ATTR.type='fluid velocity';
    febio_spec.Output.plotfile.var{end+1}.ATTR.type='fluid volume ratio';
    febio_spec.Output.plotfile.var{end+1}.ATTR.type='fluid shear viscosity';
    febio_spec.Output.plotfile.var{end+1}.ATTR.type='fluid vorticity';
    febio_spec.Output.plotfile.var{end+1}.ATTR.type='nodal fluid velocity';
    febio_spec.Output.plotfile.var{end+1}.ATTR.type='nodal relative fluid velocity';
    febio_spec.Output.plotfile.var{end+1}.ATTR.type='relative fluid velocity';
    
    %%
    
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
    febioAnalysis.disp_log_on=1; %Display convergence information in the command window
    febioAnalysis.runMode='internal';%'internal';
    febioAnalysis.t_check=0.25; %Time for checking log file (dont set too small)
    febioAnalysis.maxtpi=1e99; %Max analysis time
    febioAnalysis.maxLogCheckTime=1000; %Max log file checking time
    
    [runFlag]=runMonitorFEBio(febioAnalysis);%START FEBio NOW!!!!!!!!
end
end