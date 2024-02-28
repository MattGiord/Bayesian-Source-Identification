% Copyright (c) 2024 Matteo Giordano
%
% Codes accompanying the article "A Bayesian approach with Gaussian priors to 
% the inverse problem of source identification in elliptic PDEs" 
% by Matteo Giordano

%%
% Bayesian nonparametric inference for the source f:O -> R 
% with the Matérn process priors via conjugate formulae
%
% Requires output of GenerateObservations.m (including f_0, observations 
% and geometry) and K_mat.m (Matérn covariance kernel).

%%
% Mesh for discretisation of parameter space via piecewise linear functions

model_prior = createpde(); 
geometryFromMesh(model_prior,tnodes,telements);
generateMesh(model_prior,'Hmax',0.1);
mesh_nodes_prior=model_prior.Mesh.Nodes;
mesh_nodes_num_prior=size(mesh_nodes_prior); 
mesh_nodes_num_prior=mesh_nodes_num_prior(2);
mesh_elements_prior = model_prior.Mesh.Elements;
mesh_elements_num_prior = size(mesh_elements_prior); 
mesh_elements_num_prior = mesh_elements_num_prior(2); 
    % number of triangular mesh elements
[~,mesh_elements_area_prior] = area(model_prior.Mesh); 
    % area of triangular mesh elements

% Compute barycenters of triangular mesh elements
barycenters_prior = zeros(2,mesh_elements_num_prior);
for i=1:mesh_elements_num_prior
    barycenters_prior(:,i) = mean(mesh_nodes_prior(:,mesh_elements_prior(1:3,i)),2);
end

f0_mesh_prior = f0(mesh_nodes_prior(1,:),mesh_nodes_prior(2,:));
f0_bary_prior = f0(barycenters_prior(1,:),barycenters_prior(2,:));
f0_norm=sqrt(sum((f0_bary_prior).^2.*mesh_elements_area_prior));

%%
% Piecewise linear basis functions

f_lin=zeros(mesh_nodes_num_prior,1);
f_lin(115)=1;

figure()
pdeplot(model_prior,'XYData',f_lin,'ZData',f_lin,'ColorMap','parula');
legend('e_{115}','Fontsize',25)

%%
% Discretisation of solution operator

% Discretisation of operator for the elements of basis over full mesh
fwd_operator=zeros(sample_size,mesh_nodes_num_prior); 
    % sample_size x mesh_elements_num_prior, whose (i,j)-th element is 
    % L^{-1}e_j((x_i,y_i)) where L^{-1} is the solution operator, e_j is
    % the j-th basis function and (x_i,y_i) is the i-th design point.

tic

for j=1:mesh_nodes_num_prior
    %ej_mesh_prior=exp(-(mesh_nodes_prior(1,:)-barycenters_prior(1,j)).^2 ...
    %    -(mesh_nodes_prior(2,:)-barycenters_prior(2,j)).^2);
    ej_mesh_prior=zeros(1,mesh_nodes_num_prior);
    ej_mesh_prior(j)=1;
    ej_fun=@(location,state) griddata(mesh_nodes_prior(1,:), ...
        mesh_nodes_prior(2,:),ej_mesh_prior,location.x,location.y);
    specifyCoefficients(model,'m',0,'d',0,'c',c_fun,'a',0,'f',ej_fun);
    results = solvepde(model);
    u=results.NodalSolution;
    fwd_operator(:,j) = u(rand_index);
end

toc

%%
% Prior covariance matrix for Matern process prior

% Compute covariance matrix for Matérn process on mesh
prior_regularity = 10; prior_len = .25;
prior_cov = zeros(mesh_nodes_num_prior,mesh_nodes_num_prior);
for i=1:mesh_nodes_num_prior
    for j=i:mesh_nodes_num_prior
        prior_cov(i,j) = K_mat(mesh_nodes_prior(:,i),...
            mesh_nodes_prior(:,j),prior_regularity,prior_len);
        prior_cov(j,i) = prior_cov(i,j);
    end
end

%%
% Sample and plot a prior draw

rng(1)

f_rand = mvnrnd(zeros(mesh_nodes_num_prior,1),prior_cov,1)';

figure()
pdeplot(model_prior,'XYData',f_rand,'ZData',f_rand,'ColorMap',parula);
title('f\sim\Pi(\cdot)','FontSize',15)
legend('Matérn, \nu=10, \ell=.25','FontSize',25)

%%
% Conjugate formulae for posterior variance and mean

prior_cov_inv=inv(prior_cov); 
    % inverse of prior covariance matrix

% Posterior covariance matrix
post_cov=inv(fwd_operator'*fwd_operator/(sigma^2)+prior_cov_inv);
% Posterior mean
f_mean_mesh_prior = post_cov*fwd_operator'*observations/(sigma^2);

% Plot posterior mean and estimation erorr
figure()
subplot(1,2,1)
pdeplot(model_prior,'XYData',f0_mesh_prior,'ColorMap',jet)
title('True f_0','FontSize',20)
clim([min(f0_mesh_prior),max(f0_mesh_prior)])
subplot(1,2,2)
pdeplot(model_prior,'XYData',f_mean_mesh_prior,'ColorMap',jet)
title('Posterior mean estimate','FontSize',20)
clim([min(f0_mesh_prior),max(f0_mesh_prior)])

figure()
pdeplot(model_prior,'XYData',f0_mesh_prior-f_mean_mesh_prior','ColorMap',hot)
title('Estimation error','FontSize',20)

% Compute piecewise constant approximations of posterior mean at
% the triangle baricenters
f_mean_bary_prior = griddata(mesh_nodes_prior(1,:),mesh_nodes_prior(2,:),...
    f_mean_mesh_prior,barycenters_prior(1,:),barycenters_prior(2,:));

% Approximate L^2 distance between f_0 and posterior mean
estim_error = sqrt(sum((f0_bary_prior-f_mean_bary_prior).^2.*mesh_elements_area_prior));
disp(['L^2 estim error = ', num2str(estim_error)])
disp(['L^2 norm of F_0 = ', num2str(f0_norm)])
