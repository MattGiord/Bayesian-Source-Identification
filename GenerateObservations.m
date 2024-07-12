Copyright (c) 2024 Matteo Giordano
%
% Codes accompanying the article "A Bayesian approach with Gaussian priors to 
% the inverse problem of source identification in elliptic PDEs" 
% by Matteo Giordano

%%
% Discrete noisy observations of elliptic PDE solution
%
% Let O be a smooth domain in R^2. Consider the 2D elliptic PDE in
% divergence form with homogeneous Dirichlet boundary conditions:
%
%   div(c grad u)=f, in O
%   u=0, on boundary of O
%
% where f:O -> R is the unknown (sufficiently smooth) source function and 
% c:O -> [0,+infty) is a given smooth diffusivity.
%
% There exists a unique classical solution G(f)=u_f, that depends linearly on
% f. We observe data
%
%   Y_i = u_f(x_i) + sigma W_i,   i = 1,...,n
%
% where x_i are design points in O, sigma>0 and W_i are i.i.d. 
% N(0,1) random variables.
%
% The following code generates n observations Y_1,...,Y_n.

%%
% Create rotate elliptically-shaped domain and triangular mesh

% Display more digits
format long

ax_h = 1; 
    % length of horizontal semiaxis
ax_v = .75; 
    % length of vertical semiaxis
rot = pi/6;
    % angle of rotation
t = linspace(0,2*pi,1000);
pgon = polyshape({ax_h*cos(t)*cos(rot) - ax_v*sin(t)*sin(rot)},...
    {ax_v*sin(t)*cos(rot) + ax_h*cos(t)*sin(rot)});
vol = pi*ax_h*ax_v;

% Create a triangulation representation of pgon
tr = triangulation(pgon);

% Create a PDE model
model = createpde;

% With the triangulation data as a mesh, use the geometryFromMesh function
% to create a geometry
tnodes = tr.Points';
telements = tr.ConnectivityList';
geometryFromMesh(model,tnodes,telements);
%pdegplot(model)

% Generate and plot triangular mesh
generateMesh(model,'Hmax',.05);
figure()
axes('FontSize', 20, 'NextPlot','add')
pdemesh(model)
xticks([-1,-.5,0,.5,1])
yticks([-1,-.5,0,.5,1])
xlabel('x', 'FontSize', 20);
ylabel('y', 'FontSize', 20);
mesh_nodes = model.Mesh.Nodes; 
    % 2 x mesh_size matrix whose columns contain the (x,y) coordinates 
    % of the nodes in the mesh
mesh_nodes_num = size(mesh_nodes); 
mesh_nodes_num=mesh_nodes_num(2); 
    % number of nodes in the mesh
mesh_elements = model.Mesh.Elements; 
    % 6 x mesh_elements_num whose columns contain the 6 node indices 
    % identifying each triangle. The first 3 elements of each column contain 
    % the indices of the 3 vertices of the triangle 
mesh_elements_num = size(mesh_elements); 
mesh_elements_num = mesh_elements_num(2); 
    % number of triangles in the mesh
[~,mesh_elements_area] = area(model.Mesh); 

% Compute barycenters of triangular mesh elements
barycenters = zeros(2,mesh_elements_num);
for i=1:mesh_elements_num
    barycenters(:,i) = mean(mesh_nodes(:,mesh_elements(1:3,i)),2);
end


%%
% Specify diffusivity c and true source f_0

% Specify the diffusivity c as a function of (x,y)
c_min = 1;
c = @(x,y) c_min +1 +5*exp(-(5*x-2).^2-(5*y-2).^2) ...
   + 5*exp(-(5*x+2).^2-(5*y+2).^2);
c_mesh = c_min + 1 + 5*exp(-(5*mesh_nodes(1,:)-2).^2 ...
    - (5*mesh_nodes(2,:)-2).^2) ...
    + 5*exp(-(5*mesh_nodes(1,:)+2).^2-(5*mesh_nodes(2,:)+2).^2);

% Plot c
figure()
axes('FontSize', 20, 'NextPlot','add')
pdeplot(model,'XYData',c_mesh)
colorbar('Fontsize',20)
% pdeplot(model,'XYData',f0_num,'ZData',f0_num,'ColorMap',jet) 
    % 3D plot
%title('Known diffusivity c','FontSize',20);
xticks([-1,-.5,0,.5,1])
yticks([-1,-.5,0,.5,1])
xlabel('x', 'FontSize', 20);
ylabel('y', 'FontSize', 20);
crameri vik

% Specify the unknown source f_0 as a function of (x,y)
f0 = @(x,y) exp(-(5*x-2.5).^2-(5*y).^2)+exp(-(7.5*x).^2-(2.5*y).^2)...
    +exp(-(5*x+2.5).^2-(5*y).^2);
f0_mesh=exp(-(5*mesh_nodes(1,:)-2.5).^2-(5*mesh_nodes(2,:)).^2)...
    +exp(-(7.5*mesh_nodes(1,:)).^2-(2.5*mesh_nodes(2,:)).^2)...
    +exp(-(5*mesh_nodes(1,:)+2.5).^2-(5*mesh_nodes(2,:)).^2);

% Plot f0
figure()
axes('FontSize', 20, 'NextPlot','add')
pdeplot(model,'XYData',f0_mesh,'ColorMap',jet)
colorbar('Fontsize',20)
%title('True source f_0','FontSize',20);
xticks([-1,-.5,0,.5,1])
yticks([-1,-.5,0,.5,1])
xlabel('x', 'FontSize', 20);
ylabel('y', 'FontSize', 20);
crameri vik

%%
% Elliptic PDE solution corresponding to diffusivity c and true source f_0

% Specify diffusivity and source as a functions of (location,state) 
% to pass to elliptic PDE solver
c_fun=@(location,state) c(location.x,location.y);
f0_fun=@(location,state) f0(location.x,location.y);

% Specify zero Dirichlet boundary conditions on all edges
applyBoundaryCondition(model,'dirichlet','Edge',1:model.Geometry.NumEdges,'u',0);

% Specify the PDE coefficients
specifyCoefficients(model,'m',0,'d',0,'c',c_fun,'a',0,'f',f0_fun);

% Solve the PDE and plot solution
results = solvepde(model);
u0 = results.NodalSolution; 
figure()
axes('FontSize', 20, 'NextPlot','add')
pdeplot(model,'XYData',u0)
%title('PDE solution G(f_0)\equiv u_{f_0}','FontSize', 20)
xlabel('x','FontSize', 20)
ylabel('y', 'FontSize', 20)
xticks([-1,-.5,0,.5,1])
yticks([-1,-.5,0,.5,1])
crameri batlowW

u0_interp=scatteredInterpolant(mesh_nodes(1,:)',mesh_nodes(2,:)',u0);
u0_bary=u0_interp(barycenters(1,:),barycenters(2,:));
% Approximate L^2 distance between f_0 and posterior mean
u0_norm = sqrt(sum((u0_bary).^2.*mesh_elements_area));
disp(['Norm of u0 = ', num2str(u0_norm)])

%%
% Noisy observations of PDE solution

rng(1)

% Sample design points and noise variales
sample_size=4500; 
    % number of observations
sigma=0.0005;
disp(['Signal to noise ratio = ', num2str(u0_norm/sigma)])
    % noise standard deviation
rand_index=sort(randsample(mesh_nodes_num,sample_size)); 
    % random indices in the mesh
rand_mesh=mesh_nodes(:,rand_index); 
    % random sample of mesh points drawn uniformly at random
figure()
axes('FontSize', 20, 'NextPlot','add')
scatter(rand_mesh(1,:),rand_mesh(2,:),'filled') 
    % plot of sampled locations
%title('Design points X_1,...,X_n','FontSize', 20)
xticks([-1,-.5,0,.5,1])
yticks([-1,-.5,0,.5,1])
xlabel('x','FontSize', 20)
ylabel('y', 'FontSize', 20)

% Sample observations
observations=u0(rand_index)+(mvnrnd(zeros(sample_size,1),...
    sigma^2*eye(sample_size)))'; 
    % add i.i.d N(0,sigma^2) noise to the observation
u0_noisy=u0;
u0_noisy(rand_index)=observations;

% Plot corrupted PDE solution
figure()
axes('FontSize', 20, 'NextPlot','add')
pdeplot(model,'XYData',u0_noisy)
%title('Observations Y_i=u_{f_0}(X_i)+\sigma W_i','FontSize', 20);
xlabel('x','FontSize', 20)
ylabel('y', 'FontSize', 20)
xticks([-1,-.5,0,.5,1])
yticks([-1,-.5,0,.5,1])
crameri batlowW
