%% Inputing Initial Data 
clc
clear all
close all
%clear classes
%%  Making Structure of the data
data=Problem();

%Initializing parameter
l=1; % initial width of the domain
epsilon=2; % random number to initialize the loop
error=0.01; % desired RMAE (stopping criteria)
b_initial_temp=data.p_trial; % initializing the b
count=0; %for counting number of loops
bounds=[];
b_initial_vec=[]; %storing initial b's from all iterations
%b_opt_orig_vec=[]; %storing optimal b's from original model
b_opt_surr_vec=[]; %storing optimal b's from surrogate model
psi0_surr_vec1=[];
psi0_initial_surr_vec=[]; %storing initial psi0 for all iterations
%psi0_orig_vec=[]; %storing optimal psi0 from original model
psi0_surr_vec=[]; %storing optimal psi0 from surrogate model
psi0_surr_vec1=[];
psi1_initial_vec=[]; %storing initial psi1 for all iterations
%psi1_orig_vec=[]; %storing optimal psi1 from original model
psi1_surr_vec=[]; %storing optimal psi1 from surrogate model
epsilon_vec=[]; %RMAE of each iteration
%%
tic;
while epsilon>error
tic;
    %%
    % Placing the optimized values at 
    u=data.interpolate(1,b_initial_temp);
    for i=1:size(data.params_keys,2)
    if data.params_keys(i)~=0   
        data.p(i)=u(1);
        u(1)=[];
    end
    end
    %%  creating new bounds
    p_lower=data.p-l; p_upper=data.p+l;
    for i=1:size(data.p,2)
        data.params_values(i,1)=p_lower(i);
        data.params_values(i,2)=p_upper(i);
    end
%% Applying conditions to keep the bounds within the limit
%if the count is equals 0 (for the first iteration, it stores the values as
%initial value set
%if count is different(for all other iterations), it tellies the upper and
%lower bounds and checks if the bounds have exceeded the initial bounds; if
%yes, it restricts the bound to be within. 

     if count==0
         data.params_values(3,1)=0.1;
         data.params_values(3,2)=0.6;
        initial_value_set=data.params_values;
    else 
        for i=1:size(data.params_keys,2)
            if data.params_values(i,1)<initial_value_set(i,1)
               data.params_values(i,1)=initial_value_set(i,1);
            end
            if data.params_values(i,2)>initial_value_set(i,2)
               data.params_values(i,2)=initial_value_set(i,2);
            end
        end
     end

     %%
    data.p_trial=[];
    %Adding a value of p_trial
    for i=1:size(data.params_keys,2)
        if data.params_keys(i)~=0
            data.p_trial=[data.p_trial,repmat(data.p(i),1,data.n(i))];
        end
    end
%% generating sampling points using latin hypercube
p_train=create_p_samplingmethod(data,data.points);
data.p_train=p_train;
%% creating snapshot matrix
%initializing t vector
data.tvector=create_t_vector(data);
tic;
%Filling the solution of ODE in snapshot matrices
[snapshot_matrix,snapshot_matrix_y1,snapshot_matrix_y2]= ...
create_snapshots(data,p_train);
fprintf('Time Elapsed for creating snapshots is %d\n',toc)

%% SVD
[svd_u,svd_s,svd_v]=svd(snapshot_matrix);
% Here, svd_u and svd_v are orthogonal matrices, whereas svd_s is a
% diagonal matrix with singular values in leading diagonal. 

%% Calculation of Amplitudes

%After observing the singular value chart and commulative energy chart, the
%rank is selected.
rank=7;
data.rank=rank;

%Once we have the rank, the first k=rank columns of orthogonal matrix svd_u
%from the SVD decomposition make our basis phi
phi=svd_u(:,1:data.rank);
data.phi=phi;

%After having the basis, it is easy to calculate amplitudes A; since the
%basis is orthogonal, the inverse of matrix is same as transpose.
A= phi'*snapshot_matrix;


%% Radial Basis Functions
%Calculation of B matrix, which is matrix of basis of RBF interpolation
B=data.rbf(A,p_train); %p_train is the matrix of our sample points
data.B=B;

%% Calculation of G, the matrix of interpolation coefficients, 
% that depend on each sample point. For now, we are using our p_trial to 
%test the results.  

%calculation of g(p) %we need pX as input of podGvec and we create it
%seperately for each test point 
pX=create_pX(data,data.p_trial);
G=data.pod_G_vec(p_train,pX);
[p_point_approx_y1,p_point_real_y1,p_point_approx_y2,p_point_real_y2]=...
       create_real_approx(data,phi,B,G);

%% Solving optimization problem using ode for a single points

opt = {[true true], [false true], [true, false], [false false]};
% The above cell of options "opt" represent the combination of contraints
% available. 1: [true true] respresnet both optimization and bilateral 
% constraints are present. 2: [false true] represent that optimization
% constraint is absent, but bilateral constraint is present and so on. 

% To make is easier, we allow to select the value of desired combination of
% constraints by assigning its index to the variable "setno"
setno=1; 

%Our trial sample is set as b0 and b initial
% data.p_trial=[];
% %Adding a value of p_trial
% for i=1:size(data.params_keys,2)
%     if data.params_keys(i)~=0
%         data.p_trial=[data.p_trial,repmat(data.p(i),1,data.n(i))];
%     end
% end
if count==0
    tic;
b_initial= data.p_trial;
b=b_initial;
data.b = b;
data.b0 = b;
data.response_type= 2; %for direct method
psi0_initial_orig=data.criteria();
psi1_initial=data.constraint();
data.optimize(opt{setno}(1), opt{setno}(2));
psi0_orig = data.criteria();
psi1_orig = data.constraint();
b_opt_orig=data.b;

fprintf('Time Elapsed for calculating original optimization criteria is %d\n',toc)
end
%% Solving optimization problem using surrogate model for a single point
% This one is same as above, only the response type changes from 2=direct
% to 1=surrogate. Also, some values like phi, B, p_train, that are used in
% surrogate modeeling are redefined for the data structed q that we have
% made now. 
tic;
b=b_initial;
data.b = b;
data.b0 = b;
data.phi=phi;
data.p_train=p_train;
data.B=B;
data.response_type=1;
psi0_initial_surr=data.criteria();
data.optimize(opt{setno}(1), opt{setno}(2));
psi0_surr = data.criteria();
psi1_surr = data.constraint();
b_opt_surr = data.b;
% responses of surrogate point using initial model
data.response_type=2;
psi0_surr1=data.criteria;
psi1_surr1=data.constraint();


fprintf('Time Elapsed for calculating surrogate optimization criteria is %d\n',toc)
epsilon=max(data.relmaxabs(psi0_surr,psi0_surr1));
fprintf('The relative maximum absolute error for psi0 for test points is %d\n',epsilon);
fprintf('Time Elapsed for estimating error for this iteration is %d\n',toc)

%% Updating parameters
count=count+1;
b_initial_temp=data.b;
b_initial_vec=[b_initial_vec;b_initial]; 
%b_opt_orig_vec=[b_opt_orig_vec;b_opt_orig]; 
b_opt_surr_vec=[b_opt_surr_vec;b_opt_surr]; 
%psi0_initial_orig_vec=[psi0_initial_orig_vec;psi0_initial_orig]; 
psi0_initial_surr_vec=[psi0_initial_surr_vec;psi0_initial_surr];
%psi0_orig_vec=[psi0_orig_vec;psi0_orig]; 
psi0_surr_vec=[psi0_surr_vec;psi0_surr]; 
psi0_surr_vec1=[psi0_surr_vec1;psi0_surr1]; 
psi1_initial_vec=[psi1_initial_vec;psi1_initial]; 
%psi1_orig_vec=[psi1_orig_vec;psi1_orig]; 
psi1_surr_vec=[psi1_surr_vec;psi1_surr]; 
epsilon_vec=[epsilon_vec; epsilon];
bounds=[bounds;data.params_values];
l=0.5*l;
end
fprintf('Time elapsed for optimization through surrogate model is %d\n',toc);
fprintf('The total number of iterations done for achieving the tolerance level are %d\n',count);
%% Clearing intermediate parameters
clear u; clear setno; clear p_upper; clear p_lower; clear i; clear b; clear opt; clear b_initial_temp;
clear b_initial; clear b_opt_orig; clear b_opt_surr; clear psi0_initial; clear psi0_orig;
clear psi0_surr; clear psi1_initial; clear psi1_orig; clear psi1_surr; clear epsilon;