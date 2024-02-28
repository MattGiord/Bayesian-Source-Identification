% Copyright (c) 2024 Matteo Giordano
%
% Codes accompanying the article "A Bayesian approach with Gaussian priors to 
% the inverse problem of source identification in elliptic PDEs" 
% by Matteo Giordano

%%
% Semiparametric Bayesian inference for the source f:O -> R with truncated 
% Gaussian series priors defined on the Dirichlet Laplacian eigebasis
%
% Requires output of GenerateObservations.m (including f_0, observations 
% and geometry).


%%
% Find the Dirichlet-Laplacian eigenpairs

% Create mesh
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
barycenters_prior = zeros(2,mesh_elements_num_prior);
for i=1:mesh_elements_num_prior
    barycenters_prior(:,i) = mean(mesh_nodes_prior(:,mesh_elements_prior(1:3,i)),2);
end

% Solve elliptic eigenvalue problem for the Dirichlet-Laplacian on the
% created mesh
applyBoundaryCondition(model_prior,'dirichlet','Edge', ...
    1:model.Geometry.NumEdges,'u',0); 
specifyCoefficients(model_prior,'m',0,'d',1,'c',1,'a',0,'f',0);
range = [-1,250]; 
results = solvepdeeig(model_prior,range); 
lambdas_basis = results.Eigenvalues; 
J_basis = length(lambdas_basis); 
e_basis = results.Eigenvectors; 
e_basis(:,1:J_basis) = e_basis(:,1:J_basis)*sqrt(mesh_nodes_num_prior/vol);

%%
% Projection of f_0 onto the Dirichlet-Laplacian eigenbasis

f0_mesh_prior = f0(mesh_nodes_prior(1,:),mesh_nodes_prior(2,:));
f0_interp=scatteredInterpolant(mesh_nodes_prior(1,:)',mesh_nodes_prior(2,:)',f0_mesh_prior');
f0_bary_prior=f0_interp(barycenters_prior(1,:),barycenters_prior(2,:));
f0_coeff=zeros(J_basis,1); 
for j=1:J_basis
    psi_interp=scatteredInterpolant(mesh_nodes_prior(1,:)',...
        mesh_nodes_prior(2,:)',e_basis(:,j));
    psi_bary_prior=psi_interp(barycenters_prior(1,:),...
        barycenters_prior(2,:));
    f0_coeff(j)=sum(mesh_elements_area_prior.*f0_bary_prior.*psi_bary_prior);
end
f0_proj=zeros(1,mesh_nodes_num_prior);
for j=1:J_basis
    f0_proj = f0_proj+f0_coeff(j)*e_basis(:,j)';
end
f0_proj_bary_prior = griddata(mesh_nodes_prior(1,:),mesh_nodes_prior(2,:),f0_proj,...
    barycenters_prior(1,:),barycenters_prior(2,:));

% Approximation error
approx_error = sqrt(sum((f0_bary_prior-f0_proj_bary_prior).^2.*mesh_elements_area_prior));
disp(['L^2 approximation error via projection = ', num2str(approx_error)])
f0_norm=norm(f0_coeff);
disp(['L^2 norm of f_0 = ', num2str(f0_norm)])

%%
% Discretisation of solution operator

% Discretisation of operator for the elements of basis over full mesh
fwd_operator=zeros(mesh_nodes_num,J_basis); 
for j=1:J_basis
    ej_fun=@(location,state) griddata(mesh_nodes_prior(1,:), ...
        mesh_nodes_prior(2,:),e_basis(:,j),location.x,location.y);
    specifyCoefficients(model,'m',0,'d',0,'c',c_fun,'a',0,'f',ej_fun);
    results = solvepde(model);
    fwd_operator(:,j)=results.NodalSolution;
end

%%
% Prior covariance matrix for Gaussian series prior draws

prior_regularity=.75; 
prior_cov=diag(lambdas_basis.^(-2*prior_regularity)); 
    % diagonal prior covariance matrix
prior_cov_inv=diag(lambdas_basis.^(2*prior_regularity)); 
    % inverse of prior covariance matrix


%%
% Find SVD basis to diagonalise differential (and forward) operator

% Specify coefficients for eigenvalue equation associated to differential 
% operator 
specifyCoefficients(model_prior,'m',0,'d',1,'c',c_fun,'a',0,'f',0);
range = [-1,500]; 
results = solvepdeeig(model_prior,range); 
lambdas_SVD = results.Eigenvalues; 
J_SVD = length(lambdas_SVD); 
e_SVD = results.Eigenvectors; 
e_SVD(:,1:J_SVD) = e_SVD(:,1:J_SVD)*sqrt(mesh_nodes_num_prior/vol);

% Plot the eigenfunctions
figure() 
subplot(1,3,1)
pdeplot(model_prior,'XYData',e_SVD(:,1));
title('e_0','FontSize',15)
subplot(1,3,2)
pdeplot(model_prior,'XYData',e_SVD(:,2));
title('e_2','FontSize',15)
subplot(1,3,3)
pdeplot(model_prior,'XYData',e_SVD(:,J_SVD));
title('e_J','FontSize',15)

% Plot the eigenvalues
figure()
axes('FontSize', 15, 'NextPlot','add')
plot(lambdas_basis,'.','Linewidth',3)
xlabel('j', 'FontSize', 15);
ylabel('\lambda_j', 'FontSize', 15);

%%
% Comparison of finite element methods and spectral formula for forward 
% and differential operators

f0_SVD_coeff=zeros(J_SVD,1); 
    % initialises vector to store the Fourier coefficients of f0 in the 
    % SVD basis
for j=1:J_SVD
    ej_SVD_interp=scatteredInterpolant(mesh_nodes_prior(1,:)',...
        mesh_nodes_prior(2,:)',e_SVD(:,j));
    ej_SVD_bary_prior=ej_SVD_interp(barycenters_prior(1,:),...
        barycenters_prior(2,:));
    f0_SVD_coeff(j)=sum(mesh_elements_area_prior.*f0_bary_prior.*ej_SVD_bary_prior);
end

% Compute u_f0 via the spectral characterisation
u0_SVD=zeros(1,mesh_nodes_num_prior);
for j=1:J_SVD
    u0_SVD = u0_SVD+f0_SVD_coeff(j)*e_SVD(:,j)'/lambdas_SVD(j);
end

% Compare the two solutions
figure()
subplot(1,2,1)
pdeplot(model,'XYData',u0,'ColorMap',jet)
title('u_{f0} via FEM','FontSize',20)
colorbar('Fontsize',15)
subplot(1,2,2)
pdeplot(model_prior,'XYData',u0_SVD,'ColorMap',jet)
title('u_{f0} via SVD','FontSize',20)
colorbar('Fontsize',15)

% Apply differential operator to u_f0 (defined on original mesh)
u0_interp=scatteredInterpolant(mesh_nodes(1,:)',...
        mesh_nodes(2,:)',u0);
u0_bary_prior=u0_interp(barycenters_prior(1,:),...
        barycenters_prior(2,:));
u0_SVD_coeff=zeros(J_SVD,1); 
for j=1:J_SVD
    ej_SVD_interp=scatteredInterpolant(mesh_nodes_prior(1,:)',...
        mesh_nodes_prior(2,:)',e_SVD(:,j));
    ej_SVD_bary_prior=ej_SVD_interp(barycenters_prior(1,:),...
        barycenters_prior(2,:));
    u0_SVD_coeff(j)=sum(mesh_elements_area_prior.*u0_bary_prior.*ej_SVD_bary_prior);
end
L_cu0_SVD=zeros(1,mesh_nodes_num_prior);
for j=1:J_SVD
    L_cu0_SVD = L_cu0_SVD+u0_SVD_coeff(j)*e_SVD(:,j)'*lambdas_SVD(j);
end

% Compare f0 and L_c(u0)
figure()
subplot(1,2,1)
pdeplot(model,'XYData',f0_mesh,'ColorMap',jet)
title('f_0','FontSize',20)
colorbar('Fontsize',15)
subplot(1,2,2)
pdeplot(model_prior,'XYData',L_cu0_SVD,'ColorMap',jet)
title('\nabla\cdot(c\nabla u_{f0}) via SVD','FontSize',20)
colorbar('Fontsize',15)

%%
% Asymptotic variance for one-dimensional marginal posteriors

% Fix test function (in the Dirichlet-Laplacian eigenbasis)
index_test_function=2;
psi_interp=scatteredInterpolant(mesh_nodes_prior(1,:)',...
        mesh_nodes_prior(2,:)',e_basis(:,index_test_function));
psi_bary_prior=psi_interp(barycenters_prior(1,:),barycenters_prior(2,:));

% Projection of test function onto the SVD eigenbasis of the operator
psi_SVD_coeff=zeros(J_SVD,1); 
for k=1:J_SVD
    ej_SVD_interp=scatteredInterpolant(mesh_nodes_prior(1,:)',...
        mesh_nodes_prior(2,:)',e_SVD(:,k));
    ej_SVD_bary_prior=ej_SVD_interp(barycenters_prior(1,:),...
        barycenters_prior(2,:));
    psi_SVD_coeff(k)=sum(mesh_elements_area_prior.*psi_bary_prior.*ej_SVD_bary_prior);
end

% Plot true test function and SVD projection
%psi_SVD_proj=zeros(1,mesh_nodes_num_prior);
%for k=1:J_SVD
%    psi_SVD_proj = psi_SVD_proj+psi_SVD_coeff(k)*e_SVD(:,k)';
%end
%figure()
%subplot(1,2,1)
%pdeplot(model_prior,'XYData',e_basis(:,index_test_function),'ColorMap',jet)
%title('True test function \psi','FontSize',20)
%colorbar('Fontsize',15)
%subplot(1,2,2)
%pdeplot(model_prior,'XYData',psi_SVD_proj,'ColorMap',jet)
%title('SVD projection of \psi','FontSize',20)
%colorbar('Fontsize',15)

% Apply differential operator to test function
L_cpsi_coeff=zeros(J_SVD,1);
%L_cpsi=zeros(1,mesh_nodes_num_prior);
for k=1:J_SVD
    L_cpsi_coeff(k)=lambdas_SVD(k)*psi_SVD_coeff(k);
    %L_cpsi = L_cpsi+L_cpsi_coeff(k)*e_SVD(:,k)';
end

% Plot differential operator applied to test function
%figure()
%pdeplot(model_prior,'XYData',L_cpsi,'ColorMap',jet)
%title('\nabla\cdot(c\nabla\psi)','FontSize',20)
%colorbar('Fontsize',15)

asymptotic_var=norm(L_cpsi_coeff)^2*sigma^2/sample_size;
asymptotic_sd=sqrt(asymptotic_var);

disp(['Asymptotic standard deviation = ', num2str(asymptotic_sd)])

%%
% Posterior samples and comparison to asymptotic variance

% Posterior covariance matrix
post_cov=inv(fwd_operator'*fwd_operator/(sigma^2)+prior_cov_inv);
% Posterior mean
post_mean = post_cov*fwd_operator'*observations/(sigma^2);

post_means_j=post_mean(index_test_function);
post_var_j=post_cov(index_test_function,index_test_function);

% Histogram of posterior samples
rng(1)
thetaj_post_samples=normrnd(post_means_j,sqrt(post_var_j),250);
figure()
histogram(thetaj_post_samples)
xline(f0_coeff(index_test_function),'r','LineWidth',2)
xline(post_mean(index_test_function),'blue','LineWidth',2)
xline(post_mean(index_test_function)-1.96*asymptotic_sd,'green','LineWidth',2)
xline(post_mean(index_test_function)+1.96*asymptotic_sd,'green','LineWidth',2)
legend('\Pi(\theta_j|Y^{\epsilon})','\theta_{0,j}','Posterior mean','Fontsize',15)

%%
% Repeated experiments to evaluate coverage and variance of posterior mean
% estimator

index_test_function=2;
    %index of basis function

% Compute asymptotic variance
psi_interp=scatteredInterpolant(mesh_nodes_prior(1,:)',...
        mesh_nodes_prior(2,:)',e_basis(:,index_test_function));
psi_bary_prior=psi_interp(barycenters_prior(1,:),barycenters_prior(2,:));
psi_SVD_coeff=zeros(J_SVD,1); 
for k=1:J_SVD
    ej_SVD_interp=scatteredInterpolant(mesh_nodes_prior(1,:)',...
        mesh_nodes_prior(2,:)',e_SVD(:,k));
    ej_SVD_bary_prior=ej_SVD_interp(barycenters_prior(1,:),...
        barycenters_prior(2,:));
    psi_SVD_coeff(k)=sum(mesh_elements_area_prior.*psi_bary_prior.*ej_SVD_bary_prior);
end
% Apply differential operator to test function
L_cpsi_coeff=zeros(J_SVD,1);
for k=1:J_SVD
    L_cpsi_coeff(k)=lambdas_SVD(k)*psi_SVD_coeff(k);
end
asymptotic_var=norm(L_cpsi_coeff)^2*sigma^2/sample_size;
asymptotic_sd=sqrt(asymptotic_var);
disp(['Asymptotic standard deviation = ', num2str(asymptotic_sd)])

tic 

% Repeat experiments
rng(1)
sigma=0.0005;
sample_size=1000;
n_exp=1000;
    % number of repeated experiments
post_means=zeros(J_basis,n_exp);
post_vars=zeros(J_basis,n_exp);
coverage=0;

for s=1:n_exp
    % Sample design points
    rand_index=sort(randsample(mesh_nodes_num,sample_size)); 
    rand_mesh=mesh_nodes(:,rand_index); 

    % Sample observations
    observations=u0(rand_index)+(mvnrnd(zeros(sample_size,1),...
        sigma^2*eye(sample_size)))'; 
    
    % Discretise forward operator
    fwd_op_exp=fwd_operator(rand_index,:);

    % Posterior covariance matrix
    post_cov=inv(fwd_op_exp'*fwd_op_exp/(sigma^2)+prior_cov_inv);
    % Posterior mean
    post_means(:,s) = post_cov*fwd_op_exp'*observations/(sigma^2);
    post_vars(s) = post_cov(index_test_function,index_test_function);

    % Check if true coefficient is in credible interval
    if abs(f0_coeff(index_test_function) ...
            -post_means(index_test_function,s))<1.96*sqrt(post_vars(s))
        coverage=coverage+1;
    end

    % Estimation error (check uniformity across experiments)
    %estim_error = norm(post_mean(:,s)-f0_coeff);
    %disp(['L^2 estimation error = ', num2str(estim_error)])
end

toc

coverage=coverage/n_exp;
disp(['Obtained coverage = ', num2str(coverage)])

%%
% Histogram of realisation of posterior mean estimator with asymptotic
% variance

post_means_j=post_means(index_test_function,:);
post_means_mean=mean(post_means_j);
histogram(post_means_j)
xline(f0_coeff(index_test_function),'r','LineWidth',2)
xline(f0_coeff(index_test_function)-1.96*asymptotic_sd,'green','LineWidth',2)
xline(f0_coeff(index_test_function)+1.96*asymptotic_sd,'green','LineWidth',2)
%legend('Posterior mean realisations','\langle f_0,e_4\rangle_{L^2}','Fontsize',20)
