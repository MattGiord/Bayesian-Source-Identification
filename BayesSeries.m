% Copyright (c) 2024 Matteo Giordano
%
% Codes accompanying the article "A Bayesian approach with Gaussian priors to 
% the inverse problem of source identification in elliptic PDEs" 
% by Matteo Giordano

%%
% Bayesian nonparametric inference for the source f:O -> R 
% with truncated Gaussian series priors on the Dirichlet Laplacian eigebasis
% via conjugate formulae
%
% Requires output of GenerateObservations.m (including f_0, observations 
% and geometry).

%%
% Mesh for computation of the Dirichlet-Laplacian eigenpairs

model_prior = createpde(); 
geometryFromMesh(model_prior,tnodes,telements);
generateMesh(model_prior,'Hmax',0.075);
mesh_nodes_prior=model_prior.Mesh.Nodes;
mesh_nodes_num_prior=size(mesh_nodes_prior); 
mesh_nodes_num_prior=mesh_nodes_num_prior(2);
mesh_elements_prior = model_prior.Mesh.Elements;
mesh_elements_num_prior = size(mesh_elements_prior); 
mesh_elements_num_prior = mesh_elements_num_prior(2); 
[~,mesh_elements_area_prior] = area(model_prior.Mesh); 

% Compute barycenters of triangular mesh elements
barycenters_prior = zeros(2,mesh_elements_num_prior);
for i=1:mesh_elements_num_prior
    barycenters_prior(:,i) = mean(mesh_nodes_prior(:,mesh_elements_prior(1:3,i)),2);
end

%%
% Solve elliptic eigenvalue problem for the Dirichlet-Laplacian

tic

% Specity homogeneous Dirichlet boundary conditions
applyBoundaryCondition(model_prior,'dirichlet','Edge', ...
    1:model.Geometry.NumEdges,'u',0); 
% Specify coefficients for eigenvalue equation
specifyCoefficients(model_prior,'m',0,'d',1,'c',1,'a',0,'f',0);
range = [-1,500]; 
    % range of search for eigenvalues
results = solvepdeeig(model_prior,range); 
    % solve eigenvalue equation
lambdas_basis = results.Eigenvalues; 
    % extract eigenvalues
J_basis = length(lambdas_basis); 
    % number of eigenvalues (dimension of discretised parameter space)
e_basis = results.Eigenvectors; 
    % extract eigenfunctions

toc

figure() 
subplot(1,3,1)
pdeplot(model_prior,'XYData',e_basis(:,1)); 
    % plot first eigenfunction
title('e_0','FontSize',15)
subplot(1,3,2)
pdeplot(model_prior,'XYData',e_basis(:,2)); 
    % plot second eigenfunction
title('e_2','FontSize',15)
subplot(1,3,3)
pdeplot(model_prior,'XYData',e_basis(:,J_basis)); 
% plot eigenfunction corresponding to the largest found eigenvalue
title('e_J','FontSize',15)

% Plot the eigenvalues
figure()
axes('FontSize', 15, 'NextPlot','add')
plot(lambdas_basis,'*','Linewidth',3)
xlabel('j', 'FontSize', 15);
ylabel('\lambda_j', 'FontSize', 15);
legend('\lambda_j=O(j)','FontSize',25)

%%
%L2 normalisation of eigenfunction

% The eigenfunction evaluations over the mesh returned by the solver (in MATLAB 
% B2023a) are normalised in a way that the sum of the squared values are equal 
% to 1. Hence, the value of the eigenfunctions returned is connect to the mesh 
% size: the more points in the mesh, the smaller the values
%for i=1:J_basis
%    sum(e_basis(:,i).^2) % = 1 for each eigenfunction
%    mean(e_basis(:,i).^2) % = 1/mesh_size for each eigenfunction
%end

% L^2-normalisation of eigenfunctions
%for j=1:J
    % Put e_basis(:,j).^2 in polar coordinates. Need to add 1 to 
    % avoid numerical instability
    %ejpolar2 = @(thet,rad) (1+griddata(mesh_nodes(1,:),mesh_nodes(2,:),e_basis(:,j),...
    %rad.*cos(thet),rad.*sin(thet)).^2).*rad; 
    % Compute squared L2 norm of e_basis(:,j). Need to remove pi*(r*.999)^2 
    % to account for addition of 1 above
    %sqnorm = integral2(ejpolar2,0,2*pi,0,r*.9995)- pi*(r*.9995)^2 
    %e_basis(:,j) = e_basis(:,j)/sqrt(sqnorm); % normalise e_basis(:,j) to 
    % have unit L2 norm
%end

% For a faster normalisation, normalise the eigenfunctions so that
% the mean(e_basis(:,j).^2) = 1/sqrt(vol), thereby approximating the unit
% L^2-norm normalisation over the disk. Since originally 
% mean(e_basis(:,i).^2) = 1/mesh_prior_size, we need to multiply 
% e_basis(:,j) by sqrt(mesh_prior_size/area)
e_basis(:,1:J_basis) = e_basis(:,1:J_basis)*sqrt(mesh_nodes_num_prior/vol);

%for i=1:J_basis
%    mean(e_basis(:,i).^2) % = 1/sqrt(vol)
%end

figure() 
subplot(1,3,1)
pdeplot(model_prior,'XYData',e_basis(:,1));
title('Normalised e_1','FontSize',15)
%legend('e_1','FontSize',25)
subplot(1,3,2)
pdeplot(model_prior,'XYData',e_basis(:,2)); 
title('Normalised e_2','FontSize',15)
subplot(1,3,3)
pdeplot(model_prior,'XYData',e_basis(:,J_basis)); 
title('Normalised e_D','FontSize',15)

%%
% Projection of f_0 onto the Dirichlet-Laplacian eigenbasis

f0_mesh_prior = f0(mesh_nodes_prior(1,:),mesh_nodes_prior(2,:));
f0_interp=scatteredInterpolant(mesh_nodes_prior(1,:)',mesh_nodes_prior(2,:)',f0_mesh_prior');
f0_bary_prior=f0_interp(barycenters_prior(1,:),barycenters_prior(2,:));

f0_coeff=zeros(J_basis,1); 
    % initialises vector to store the Fourier coefficients of F0 in the 
    % Dirichlet-Laplacian eigenbasis

for j=1:J_basis
    ej_basis_interp=scatteredInterpolant(mesh_nodes_prior(1,:)',...
        mesh_nodes_prior(2,:)',e_basis(:,j));
    ej_basis_bary=ej_basis_interp(barycenters_prior(1,:),...
        barycenters_prior(2,:));
    f0_coeff(j)=sum(mesh_elements_area_prior.*f0_bary_prior.*ej_basis_bary);
end

f0_proj=zeros(1,mesh_nodes_num_prior);
for j=1:J_basis
    f0_proj = f0_proj+f0_coeff(j)*e_basis(:,j)';
end

figure()
subplot(1,3,1)
pdeplot(model_prior,'XYData',f0_mesh_prior,'ColorMap',jet)
title('True f_0','FontSize',20)
colorbar('Fontsize',15)

subplot(1,3,2)
pdeplot(model_prior,'XYData',f0_proj,'ColorMap',jet)
title('Projection of f_0','FontSize',20)
colorbar('Fontsize',15)

subplot(1,3,3)
pdeplot(model_prior,'XYData',f0_mesh_prior-f0_proj,'ColorMap',jet)
title('Approximation error','FontSize',20)
colorbar('Fontsize',15)

% Compute piecewise constant approximations of the projection of f_0 at
% the triangle baricenters
f0_proj_bary_prior = griddata(mesh_nodes_prior(1,:),mesh_nodes_prior(2,:),f0_proj,...
    barycenters_prior(1,:),barycenters_prior(2,:));

% Approximate L^2 distance between f_0 and posterior mean
approx_error = sqrt(sum((f0_bary_prior-f0_proj_bary_prior).^2.*mesh_elements_area_prior));
disp(['L^2 approximation error via projection = ', num2str(approx_error)])

% L^2 norm of f0
f0_norm=norm(f0_coeff);
disp(['L^2 norm of f_0 = ', num2str(f0_norm)])

%%
% Discretisation of solution operator

% Discretisation of operator for the elements of basis over full mesh
fwd_operator=zeros(sample_size,J_basis); 
    % mesh_size x coeff_num matrix, whose (i,j)-th element is 
    % L^{-1}e_j((x_i,y_i)) where L^{-1} is the solution operator, e_j is
    % the j-th basis function and (x_i,y_i) is the i-th design point.

tic

for j=1:J_basis
    ej_fun=@(location,state) griddata(mesh_nodes_prior(1,:), ...
        mesh_nodes_prior(2,:),e_basis(:,j),location.x,location.y);
    specifyCoefficients(model,'m',0,'d',0,'c',c_fun,'a',0,'f',ej_fun);
    results = solvepde(model);
    u=results.NodalSolution;
    fwd_operator(:,j) = u(rand_index);
end

toc

%%
% Prior covariance matrix for Gaussian series prior draws

prior_regularity=.75; 
prior_cov=diag(lambdas_basis.^(-2*prior_regularity)); 
    % diagonal prior covariance matrix

%%
% Sample and plot a prior draw

rng(1)

theta_rand=mvnrnd(zeros(J_basis,1),prior_cov,1)'; 
% sample Fourier coefficients from prior

f_rand=zeros(1,mesh_nodes_num_prior);

for j=1:J_basis
    f_rand = f_rand+theta_rand(j)*e_basis(:,j)';
end

figure()
axis equal
pdeplot(model_prior,'XYData',f_rand,'ColorMap',jet);
%title('f\sim\Pi(\cdot)','FontSize',15)
legend('s=1.5','FontSize',25)

%%
% Conjugate formulae for posterior variance and mean

prior_cov_inv=diag(lambdas_basis.^(2*prior_regularity)); 
    % inverse of prior covariance matrix

% Posterior covariance matrix
post_cov=inv(fwd_operator'*fwd_operator/(sigma^2)+prior_cov_inv);
% Posterior mean
post_mean = post_cov*fwd_operator'*observations/(sigma^2);

f_mean_mesh_prior=zeros(mesh_nodes_num_prior,1);
for j=1:J_basis
    f_mean_mesh_prior = f_mean_mesh_prior+post_mean(j)*e_basis(:,j);
end

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

% Approximate L^2 distance between f_0 and posterior mean
estim_error = norm(post_mean-f0_coeff);
disp(['L^2 estimation error = ', num2str(estim_error)])
disp(['L^2 norm of f_0 = ', num2str(f0_norm)])

%%
% Uncertainty quantification along x-axis

% x-values
xvals = linspace(-.9,.9,1000);
yvals= zeros(1,length(xvals));

% Computation of values of true f0 along x-axis
f0_xaxis=zeros(1,length(xvals));
for i=1:length(xvals)
    f0_xaxis(i) = f0(xvals(i),yvals(i));
end

% Computation of values of posterior mean along principal diagonal
f_mean_xaxis=griddata(mesh_nodes_prior(1,:),mesh_nodes_prior(2,:),...
    f_mean_mesh_prior,xvals,yvals);

% Plot f_0 and the posterior mean along the principal diagonal
figure()
hold on
plot(xvals,f0_xaxis,'Linewidth',2,'Color','k')
plot(xvals,f_mean_xaxis,'Linewidth',2,'Color','r')
legend('f_0(\cdot,0)','posterior mean')

rng(1)

% Draw and plot samples from posterior distribution along the x-axis
n_post_samples = 1000;
for k=1:n_post_samples
    theta_post=mvnrnd(post_mean,post_cov,1);
        % sample the vector of coefficients from posterior

    f_post_mesh_prior=zeros(mesh_nodes_num_prior,1);
    for j=1:J_basis
        f_post_mesh_prior = f_post_mesh_prior+theta_post(j)*e_basis(:,j);
    end

    f_post_xaxis=griddata(mesh_nodes_prior(1,:),mesh_nodes_prior(2,:),...
        f_post_mesh_prior,xvals,yvals);
    
    plot(xvals,f_post_xaxis,'Linewidth',.01,'Color','green')
end

% re-plot f0 and posterior mean
plot(xvals,f0_xaxis,'Linewidth',2,'Color','k')
plot(xvals,f_mean_xaxis,'Linewidth',2,'Color','r')
hold on
legend('f_0(\cdot,0)','posterior mean','Fontsize',15)
hold off

%%
% Estimation for increasing sample sizes

rng(1)

sigma=0.0005;
    % noise standard deviation

sample_sizes=[50;100;250;500;750;1000;1500;2000;3000;4500];
n_exp=length(sample_sizes);

post_mean=zeros(J_basis,n_exp);
estim_error=zeros(1,n_exp);

for s=1:n_exp
    rand_index=sort(randsample(mesh_nodes_num,sample_sizes(s))); 
        % random indices in the mesh
    rand_mesh=mesh_nodes(:,rand_index); 
        % random sample of mesh points drawn uniformly at random

    % Sample observations
    observations=u0(rand_index)+(mvnrnd(zeros(sample_sizes(s),1),...
        sigma^2*eye(sample_sizes(s))))'; 
    
    % Forward operator
    fwd_op_exp=fwd_operator(rand_index,:);

    % Posterior covariance matrix
    post_cov=inv(fwd_op_exp'*fwd_op_exp/(sigma^2)+prior_cov_inv);
    % Posterior mean
    post_mean(:,s) = post_cov*fwd_op_exp'*observations/(sigma^2);

    f_mean_mesh_prior=zeros(mesh_nodes_num_prior,1);
    for j=1:J_basis
        f_mean_mesh_prior = f_mean_mesh_prior+post_mean(j,s)*e_basis(:,j);
    end

    % Plot posterior mean and estimation erorr
    figure()
    %subplot(1,2,1)
    %pdeplot(model_prior,'XYData',f0_mesh_prior,'ColorMap',jet)
    %title('True f_0','FontSize',20)
    %clim([min(f0_mesh_prior),max(f0_mesh_prior)])
    %subplot(1,2,2)
    pdeplot(model_prior,'XYData',f_mean_mesh_prior,'ColorMap',jet)
    title(['Posterior mean, n=' num2str(sample_sizes(s))],'FontSize',20)
    clim([min(f0_mesh_prior),max(f0_mesh_prior)])

    % Approximate L^2 distance between f_0 and posterior mean
    estim_error(s) = norm(post_mean(:,s)-f0_coeff);
    disp(['L^2 estimation error = ', num2str(estim_error(s))])
    %disp(['L^2 norm of f_0 = ', num2str(f0_norm)])
end

figure()
scatter(sample_sizes,estim_error,'*','linewidth',5)
hold on
plot(sample_sizes,estim_error)
yline(approx_error,'r')
ylabel('L^2-estimation error','FontSize',15)
xlabel('n','FontSize',15)

