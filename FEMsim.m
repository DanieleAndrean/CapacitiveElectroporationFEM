%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AUTHOR: Daniele Andrean, PhD candidate Department of Information Engineering, University of Padua
% Date: 12/11/2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clear all
clc
%%
scal=100;
SaveFile="Results.mat";
%% -------------------- Parameters --------------------
eps0 = 8.854187817e-12;
% Material properties (electrolyte bath)
epsr_bath = 102; 
epsilon_bath=epsr_bath*eps0;
sigma_bath = 1.79;    % electrolyte

% Cell
epsr_mem = 3.0*scal; 
epsilon_mem=epsr_mem*eps0;
sigma_mem = 1e-6*scal;

epsr_cyt=80;
epsilon_cyt=epsr_cyt*eps0;
sigma_cyt=0.28;

% Electrode metal (if included in geometry, else electrode is external face)
epsr_elec = 11.4; 
epsilon_elec=epsr_elec*eps0;
sigma_elec = 10;
% Thin dielectric layer (characterized by surface capacitance)
epsr_diel = 3.9*scal; 
epsilon_diel=epsr_diel*eps0;
sigma_diel=1e-16*scal;

%% ----------------------- Time 
t_rise = 10e-9;   % 10 ns voltage rise taken from 33250A 80 MHz, Agilent Technologies datasheet
t_fine = 10e-6;
t_on = 0.1e-3;
t_pre  = 5e-9;       % start time
t_end  = t_pre+0.1e-3;    % final time 1 ms

Npre=3;
Nrise  = 10;
Nfall = 10;
Non = 10; % fine points during rise
Nend = 3;    % coarse points after rise

t_prestim= linspace(0,t_pre,Npre);
t_rise1   = linspace(t_pre, t_pre+2*t_rise, Nrise);
t_fall1 = linspace(t_rise1(end),t_rise1(end)+t_fine, Nfall);
t_stim = linspace(t_fall1(end),t_fall1(end)+t_on, Non);
t_rise2 = linspace(t_stim(end),t_stim(end)+t_rise, Nrise);
t_fall2 = linspace(t_rise2(end),t_rise2(end)+t_fine, Nfall);
t_post = linspace(t_rise2(end),t_end, Nend);
tvec = unique([t_prestim,t_rise1,t_fall1,t_stim,t_post]);
dt=diff(tvec);
num_steps=length(tvec);


%% ------------------------ SIGNAL
% % Applied voltage waveform (function of time)
Amp=3;

Vapplied = @(t) Amp*((( (t > t_pre) .* (t <= t_rise+t_pre) ) .* ( (t-t_pre)/ t_rise )) + ... % 1. Rise (0 to A)
    (((t > t_rise+t_pre) .* (t <= t_rise+t_pre+t_on)).* (1 - (2 *(t - t_rise-t_pre)/ t_on))) + ... % 2. Fall (A to -A)
    (((t > t_rise+t_pre+t_on) & (t <= 2*t_rise+t_pre+t_on) ) .* (-1 + (t - t_rise-t_pre-t_on) / t_rise))) + ...% 3. Return (-A to 0)
    ( (t < 0) + (t > 2*t_rise+t_pre+t_on) ) .* 0; % 4. Zero (Outside the pulse duration)

%Check signal
figure
plot(tvec,Vapplied(tvec),'o-')
Vchem=0.07;
%% --------------------------Geometries
Lelec=30e-6;
Helec=1e-6;
Hdiel=15e-9*scal;
Rcell=10e-6;
Hmem=5e-9*scal;
c_e_dist=50e-9;

%Electrolyte
bkg=multicuboid(2e-4,2e-4,2e-4,Zoffset=-1e-4);
%Membrane
cub1=fegeometry(multicuboid(2*Rcell,2*Rcell,Rcell,ZOffset=-Rcell));
cub2=multicuboid(2*Rcell,2*Rcell,Rcell+Hmem,ZOffset=-Rcell);
sph1=fegeometry(multisphere(Rcell));
sph2=multisphere(Rcell-Hmem);
sph1=subtract(sph1,cub1);
sph2=subtract(sph2,cub2);

sph1_m=generateMesh(sph1);
sph2_m=generateMesh(sph2);
modelcast=createpde();
modelcast2=createpde();
geometryFromMesh(modelcast,sph1_m.Mesh.Nodes,sph1_m.Mesh.Elements);
sph1=modelcast.Geometry;
geometryFromMesh(modelcast2,sph2_m.Mesh.Nodes,sph2_m.Mesh.Elements);
sph2=modelcast2.Geometry;

mem=addCell(sph1,sph2);
mem=translate(mem,[0,0,Helec+Hdiel+c_e_dist]);

%electrode
electrode=multicuboid(Lelec,Lelec,[Helec,Hdiel],ZOffset=[0,Helec]);


g=addCell(bkg,mem);
g=addCell(g,electrode);

%Geometry check
pdegplot(g,"CellLabels","on","FaceAlpha",0.5)


%C1: back
%C2: out sph
%C3: in sph
%C4: big cub
%C5: small cub

%% === Create PDE Model ===
model = createpde();
model.Geometry=g;
generateMesh(model, Hmin=1e-7);
figure()
pdegplot(model.Geometry,FaceAlpha=0.5, ...
                        FaceLabels="on")
% -- Solution Storage --
P = model.Mesh.Nodes;
num_nodes = size(P, 2);
V_results = zeros(num_nodes, num_steps);

% -- Initial Condition (V=0 everywhere at t=0) --
load("InitConds.mat");
V_prev_solution = zeros(num_nodes, 1);
V_prev_solution=inconds;
V_results(:, 1) = V_prev_solution;

% -- Electrode & Wall IDs (Example) --
signal_edge_ID = 12;
insulating_wall_IDs = [11, 13:21]; % All other outer walls
faceBath=1:6;

cellIDs=1:5;% bath, mem, cyt, elec, dielec
sigmas=[sigma_bath,sigma_mem,sigma_cyt,sigma_elec,sigma_diel];
epsilons=[epsilon_bath,epsilon_mem, epsilon_cyt,epsilon_elec,epsilon_diel];
thick=[1,Hmem,1,1,Hdiel];
disp('Starting Backward Euler simulation...');

for i = 2:num_steps
    disp("Step "+(i-1)+" / "+(num_steps-1));
    t = tvec(i); % Current time
    
    % --- A. Create Interpolant for Previous Solution V^{N-1} ---
    % This is how we pass the 'f' term to the solver
    V_prev_interpolant = scatteredInterpolant(P', V_prev_solution);

    % --- B. Define PDE Coefficients FOR THIS STEP ---
    
    % 'c' coefficient 
    c_func = @(location, state) coeff_c_selector(location, ...
        sigmas,cellIDs);
    
    % 'a' coefficient 
    a_func = @(location, state) coeff_a_selector(location, dt(i-1),thick, ...
        epsilons,cellIDs);

    % 'f' coefficient 
    f_func = @(location, state) coeff_f_selector(location, dt(i-1),thick, ...
        [0,epsilon_mem,0,0,epsilon_diel], V_prev_interpolant,[0,sigma_mem,0,0,0],Vchem,cellIDs);

    specifyCoefficients(model, 'm', 0, 'd', 0, 'c', c_func, 'a', a_func, 'f', f_func);

    % --- C. Apply Boundary Conditions FOR THIS STEP ---
    % Clear all old BCs
    model.BoundaryConditions = [];
    %Potential zero at bath boudaries
    applyBoundaryCondition(model, 'dirichlet', 'Face', faceBath, 'u',0);
    % Apply signal
    applyBoundaryCondition(model, 'dirichlet', 'Face', signal_edge_ID, 'u', Vapplied(tvec(i)));
    
    % Apply insulating walls
    applyBoundaryCondition(model, 'neumann', 'Face', insulating_wall_IDs, 'g', 0, 'q', 0);

    % --- D. Solve the Static Problem ---
   
    result = solvepde(model);
    
    % --- E. Store and Update ---
    V_current_solution = result.NodalSolution;
    V_results(:, i) = V_current_solution;
    V_prev_solution = V_current_solution; % Setup for the next loop
    
    if mod(i, 100) == 0
        fprintf('Progress: %.1f%%\n', (t/t_end)*100);
    end
    disp("Done")
    disp("")
end

save(SaveFile,"V_results","model","tvec","scal","g");
disp('Simulation complete.');

elems = findElements(model.Mesh,"Region",Cell=[2:g.NumCells]);

% for i=1:numel(tvec)
%     figure;
%     pdeplot3D(model.Mesh.Nodes,model.Mesh.Elements(:,elems),ColorMapData=V_results(:,i), FaceAlpha=0.7)
%     clim([-3,3])
%     axis([-10e-6,10e-6,-10e-6,10e-6,0,10e-6])
%     title("Phi t: "+tvec(i));
% end

%% ---------------------------------------- 4. Post-Processing ---

%% AVERAGE VM

% upper compart points
theta_h=linspace(0,2*pi,20);
theta_v=linspace(0,pi/2,10);

[th_h,th_v]=meshgrid(theta_h,theta_v);

xup_in=(Rcell-(Hmem+1e-6)).*sin(th_v).*cos(th_h);
yup_in=(Rcell-(Hmem+1e-6)).*sin(th_v).*sin(th_h);
zup_in=(Helec+Hdiel+c_e_dist+(Hmem+1e-6))+(Rcell-2*(Hmem+1e-6)).*cos(th_v);

xup_out=(Rcell).*sin(th_v).*cos(th_h);
yup_out=(Rcell).*sin(th_v).*sin(th_h);
zup_out=(Helec+Hdiel+c_e_dist+Hmem)+(Rcell).*cos(th_v);

pup_in=[xup_in(:),yup_in(:),zup_in(:)];
pup_out=[xup_out(:),yup_out(:),zup_out(:)];

node_in_up = findNodes(model.Mesh, 'nearest', pup_in');
node_out_up = findNodes(model.Mesh, 'nearest', pup_out');

% down compart points

R_in=linspace(0,Rcell-(Hmem+1e-6),20);
[R,th]=meshgrid(R_in,theta_h);

xdw=R.*cos(th);
ydw=R.*sin(th);
zdw_in=(Helec+Hdiel+c_e_dist+Hmem+1e-6)*ones(size(xdw));
zdw_out=(Helec+Hdiel+c_e_dist)*ones(size(xdw));

pdw_in=[xdw(:),ydw(:),zdw_in(:)];
pdw_out=[xdw(:),ydw(:),zdw_out(:)];

node_in_dw = findNodes(model.Mesh, 'nearest', pdw_in');
node_out_dw = findNodes(model.Mesh, 'nearest', pdw_out');


figure
plot3(xup_in,yup_in,zup_in,'og')
hold on
plot3(xup_out,yup_out,zup_out,'or')
plot3(xdw,ydw,zdw_in,'og')
plot3(xdw,ydw,zdw_out,'or')
pdegplot(g,"FaceAlpha",0.5)

% Get the potential time-series at these two nodes
V_in_up = mean(V_results(node_in_dw, :));
V_out_up = mean(V_results(node_out_dw, :));
V_in_dw = mean(V_results(node_in_up, :));
V_out_dw = mean(V_results(node_out_up, :));

exclude=[]; %exclude nan or numerical errors
idxplt=setdiff([1:length(tvec)],exclude);
tplot=(tvec(idxplt))*1e6;
figure
plot(tplot,V_in_up(idxplt)-V_out_up(idxplt),'b',"LineWidth",3)
hold on
plot(tplot,V_in_dw(idxplt)-V_out_dw(idxplt),'r',"LineWidth",3)
theme("light")
xlim([-1,50])
ylim([-1,1])
xlabel("time [\mus]")
legend("Free Membrane","Attached Membrane")
title("Membrane Potential")
ylabel("Potential [V]")
fontsize(30,"points")

figure
plot(tplot,V_in_up(idxplt),'b',"LineWidth",3)
hold on
plot(tplot,V_in_dw(idxplt),'r',"LineWidth",3)
plot(tplot,V_out_up(idxplt),'b--',"LineWidth",3)
plot(tplot,V_out_dw(idxplt),'r--',"LineWidth",3)
theme("light")
title("Extracellular and Intracellular Potential")
xlabel("Time [\mus]")
ylabel("\Phi [V]")
legend("Free Intracellular","Attached Intracellular","Free Extracellular","Attached Extracellular")
xlim([-1,50])
ylim([-0.1,1.5])
fontsize(30,"points")



%3D plot at interesting points
tscal=(tvec)*1e6;
tidx=[find(tscal==0),find(tscal>=0.004,1,'first'),find(tscal>=0.009,1,'first'),...
    find(tscal>=0.66,1,'first'),find(tscal>1.1,1,'first'),find(tscal>10,1,'first'),find(tscal>30,1,'first')];
tv=tscal(tidx);
elems = findElements(model.Mesh,"Region",Cell=[2:g.NumCells]);
for i=1:numel(tv)
    figure;
    pdeplot3D(model.Mesh.Nodes,model.Mesh.Elements(:,elems),ColorMapData=V_results(:,tidx(i)), FaceAlpha=0.7)
    clim([0.7,1.1])
    colormap("parula")
    axis([-10e-6,10e-6,-10e-6,10e-6,0,10e-6])
    title("\Phi at t = "+round(tv(i),1,"significant")+" \mus");
    theme("light")
    fontsize(30,"points")
end

%% --- 3. Helper Functions (Define these in separate .m files or at end) ---

function c = coeff_c_selector(location, sigmas,cellIDs)
    num_pts = length(location.x);
    c = zeros(1, num_pts);
    for i=1:length(cellIDs)
        c(location.subdomain==cellIDs(i))=sigmas(cellIDs(i));
    end
   
   
end

function a = coeff_a_selector(location, dt, thick,epsilons,cellIDs)
    num_pts = length(location.x);
    a = zeros(1, num_pts);
     for i=1:length(cellIDs)
        a(location.subdomain==cellIDs(i))=epsilons(cellIDs(i))/(thick(cellIDs(i))^2*dt);
     end   
end

function f = coeff_f_selector(location, dt, thick, epsilons, V_prev_interp,sigmas,Vchem,cellIDs)
    num_pts = length(location.x);
    f = zeros(1, num_pts);
    
    % Evaluate previous potential at all requested locations
    V_prev_vals = V_prev_interp(location.x, location.y,location.z);
    
     for i=1:length(cellIDs)
        f(location.subdomain==cellIDs(i))=epsilons(cellIDs(i))/(thick(cellIDs(i))^2*dt) * V_prev_vals(location.subdomain==cellIDs(i))...
                                            -sigmas(cellIDs(i))/(thick(cellIDs(i))^2)*Vchem;
    end
   
end