classdef Problem < handle
    properties
        p = [0.734,0.175,-0.500,-0.246,-0.425,0.635,-0.912] %model parameters
        p_train; %p_train represents the vector of sample points created after LHS design
        t0 = 0; %initial period of time
        Tend = 10; % end period of time
        y0 = [1 2]; %initial value of y
        points = 80; %total number of initial sampling points or ns
        t_points = 100; %total number of time instances at which system responses(snapshots) are recorded
        tvector; %vector with "t_points" equally distanced time periods 
        change_var = 2; %total number of model parameters that are allowed to change 
        params_keys=[0 0 0 0 1 0 1];
        params_values; %set of bounds for all the parameters
        p_trial; %if interested in checking model for one particular trial point eg. [-0.2565]%
        rank; %Chosen rank of approximation
        ng=10; %Total number of sample points for evaluation of model in online phase
        ode_options=odeset('RelTol',1e-7,'AbsTol',1e-7); %ODEoptions
        rbf_type='C';%'L' for linear spline rbf, 'C' for cubic spline rbf, 'T' for thin plate spline
        sampling_type='LHS'; %'LHS' for Latin Hypercube sampling, 'SLHS' for Symmetric LHS and 'Corner' for corner points
        B; %B is the matrix of interpolation basis obtained from rbf. It is calculated once for all.
        phi; % Phi is the matrix of basis obtained from SVD and it is also fixed and calculated once for all.
        response_type; % 1  for surrogate model, 2 for direct method - The response type determines which kind of model response is expected, original or surrogate
        
        % Parameters specific for Optimization Problem
        b;  % current opt. var values
        b0; % initial guess to be used as a constraint
        n=[0,0,0,0,2,0,2]; % number of points taken for interpolation of control function
        method = 'linear'; % method used for the interpolation of control function
        optY = 2; % index of y-vector to minimize for the current model
        yd=[5 5]; %desired value of y
        yMax=[6 6]; 
        cache={}; %save direct problem solution
    end
    methods
%% METHODS FOR SIMPLE SURROGATE MODEL   
 % The functions defined in this subsection correspond only to simple
 % simulation and solution of ODEs, construction of Surrogate models, with
 % no mention of optimization
        function p_vec=create_p_lhsdesign(data,ns)
            % This function is used to make a vector of sampling points
            % using LHS design sampling technique. 
            % Along with the data structure, it takes "ns" as input,
            % where ns is number of sampling points desired by the user. 
            xn = lhsdesign(ns,sum(data.n));% generate normalized design
            %Since the values of upper bounds and lower bounds are stored
            %in a cell, the are retrieved using a for loop for each of the
            %parameters that is varying. 
                       

            for i=1:size(data.n,2)
                lb1(i)=data.params_values(i,1); %lowebounds
                ub1(i)=data.params_values(i,2); %upperbounds
            end
            lb=[];ub=[];
            % This loops replicates the upper bounds and lower bound
            % according to the desired number of optimization parameters
            % given by data.n
            for i=1:size(data.n,2)
                if data.n(i)~=0
                    lb=[lb,repmat(lb1(i),1,data.n(i))];
                    ub=[ub,repmat(ub1(i),1,data.n(i))];
                end
            end
            %bsxfun is used to rescale the sampling points to fit in the
            %desired upper and lower bounds. 
            p_vec = bsxfun(@plus,lb,bsxfun(@times,xn,(ub-lb)));
            p_vec=transpose(p_vec);
        end
        function p_vec=create_p_SLHS(data,ns)
            delta = (1/ns)*ones(1,sum(data.n));
            
            X = zeros(ns,sum(data.n));
            for j = 1:sum(data.n)
                for i = 1:ns
                    X(i,j) = ((2*i-1)/2)*delta(j);
                end
            end
            
            P = zeros(ns,sum(data.n));
            P(:,1) = (1:ns)';
            if (mod(ns,2) == 0)
                k = ns/2;
            else
                k = (ns-1)/2;
                P(k+1,:) = (k+1)*ones(1,sum(data.n));
            end
            
            for j = 2:sum(data.n)
                P(1:k,j) = randperm(k)';
                for i = 1:k
                    if (rand(1) <= 0.5)
                        P(ns+1-i,j) = ns+1-P(i,j);
                    else
                        P(ns+1-i,j) = P(i,j);
                        P(i,j) = ns+1-P(i,j);
                    end
                end
            end
            
            xn = zeros(ns,sum(data.n));
            for j = 1:sum(data.n)
                for i = 1:ns
                    xn(i,j) = X(P(i,j),j);
                end
            end
            for i=1:size(data.n,2)
                lb1(i)=data.params_values(i,1); %lowebounds
                ub1(i)=data.params_values(i,2); %upperbounds
            end
            lb=[];ub=[];
            % This loops replicates the upper bounds and lower bound
            % according to the desired number of optimization parameters
            % given by data.n
            for i=1:size(data.n,2)
                if data.n(i)~=0
                    lb=[lb,repmat(lb1(i),1,data.n(i))];
                    ub=[ub,repmat(ub1(i),1,data.n(i))];
                end
            end
            p_vec = bsxfun(@plus,lb,bsxfun(@times,xn,(ub-lb)));
            p_vec=p_vec';
        end 
        function p_vec=create_vec_corner(data,ns)
            if ns > 2^sum(data.n)+1
                error('There are not enough corners to create the requested number of starting points')
            end
            
            for i=1:size(data.n,2)
                lb1(i)=data.params_values(i,1); %lowebounds
                ub1(i)=data.params_values(i,2); %upperbounds
            end
            lb=[];ub=[];
            % This loops replicates the upper bounds and lower bound
            % according to the desired number of optimization parameters
            % given by data.n
            for i=1:size(data.n,2)
                if data.n(i)~=0
                    lb=[lb,repmat(lb1(i),1,data.n(i))];
                    ub=[ub,repmat(ub1(i),1,data.n(i))];
                end
            end
            
            S=zeros(2^sum(data.n),sum(data.n)); %initialize sample site matrix
            for ii = 1:sum(data.n)
                S(:,ii)=repmat([repmat(lb(ii),2^(sum(data.n)-ii),1);...
                    repmat(ub(ii),2^(sum(data.n)-ii),1)],2^(ii-1),1);
            end
            
            
            %if there are more corner points than desired number of starting  points,
            %randomly select corner points
            if size(S,1) >ns-1
                r=randperm(size(S,1));
                use=S(r(1:ns-1),:); %select random corners
                samplesites=[use;lb+0.5*(ub-lb)]; %add center point
            else
                samplesites=[S;lb+0.5*(ub-lb)];%add center point
            end
            p_vec=samplesites';
        end
        function p_vec=create_p_samplingmethod(data,ns)
            %Setting up bounds
            for i=1:size(data.n,2)
                lb1(i)=data.params_values(i,1); %lowebounds
                ub1(i)=data.params_values(i,2); %upperbounds
            end
            lb=[];ub=[];
            % This loops replicates the upper bounds and lower bound
            % according to the desired number of optimization parameters
            % given by data.n
            for i=1:size(data.n,2)
                if data.n(i)~=0
                    lb=[lb,repmat(lb1(i),1,data.n(i))];
                    ub=[ub,repmat(ub1(i),1,data.n(i))];
                end
            end
            
            if strcmp(data.sampling_type,'LHS')
                xn = lhsdesign(ns,sum(data.n));
                p_vec = (bsxfun(@plus,lb,bsxfun(@times,xn,(ub-lb))))';
            elseif strcmp(data.sampling_type,'SLHS')
                delta = (1/ns)*ones(1,sum(data.n));
                
                X = zeros(ns,sum(data.n));
                for j = 1:sum(data.n)
                    for i = 1:ns
                        X(i,j) = ((2*i-1)/2)*delta(j);
                    end
                end
                
                P = zeros(ns,sum(data.n));
                P(:,1) = (1:ns)';
                if (mod(ns,2) == 0)
                    k = ns/2;
                else
                    k = (ns-1)/2;
                    P(k+1,:) = (k+1)*ones(1,sum(data.n));
                end
                
                for j = 2:sum(data.n)
                    P(1:k,j) = randperm(k)';
                    for i = 1:k
                        if (rand(1) <= 0.5)
                            P(ns+1-i,j) = ns+1-P(i,j);
                        else
                            P(ns+1-i,j) = P(i,j);
                            P(i,j) = ns+1-P(i,j);
                        end
                    end
                end
                
                xn = zeros(ns,sum(data.n));
                for j = 1:sum(data.n)
                    for i = 1:ns
                        xn(i,j) = X(P(i,j),j);
                    end
                end
                p_vec = (bsxfun(@plus,lb,bsxfun(@times,xn,(ub-lb))))';
            elseif strcmp(data.sampling_type,'RS')
                p_vec=[];
                for i=1:size(data.params_keys,2)
                    if data.params_keys(i)~=0
                        for j=1:data.n(i)
                            r=(ub(1)-lb(1)).*rand(ns,1)+lb(1);
                            p_vec=[p_vec;r'];
                            ub(1)=[];lb(1)=[];
                        end
                    end
                end
            else
                if ns > 2^sum(data.n)+1
                    error('There are not enough corners to create the requested number of starting points')
                end
                
                S=zeros(2^sum(data.n),sum(data.n)); %initialize sample site matrix
                for ii = 1:sum(data.n)
                    S(:,ii)=repmat([repmat(lb(ii),2^(sum(data.n)-ii),1);...
                        repmat(ub(ii),2^(sum(data.n)-ii),1)],2^(ii-1),1);
                end
                
                
                %if there are more corner points than desired number of starting  points,
                %randomly select corner points
                if size(S,1) >ns-1
                    r=randperm(size(S,1));
                    use=S(r(1:ns-1),:); %select random corners
                    samplesites=[use;lb+0.5*(ub-lb)]; %add center point
                else
                    samplesites=[S;lb+0.5*(ub-lb)];%add center point
                end
                p_vec=samplesites';
                
            end
        end
        function tvector=create_t_vector(data)
            % This function creates the equally spaced vector of "t_points"
            % between 't0' and 'Tend'
            tvector=linspace(data.t0,data.Tend,data.t_points);
        end
        function dy = ode_without_control(t,y,data)
            %the system of ODE for current example; without control
            %function
            dy=zeros(2,1);
            chunk=y(1)*(1-exp(-data.p(6)*y(1)))*y(2);
            dy(1)=data.p(1)*y(1)+data.p(2)*y(1)^2+data.p(5)*chunk;
            dy(2)=data.p(3)*y(2)+ data.p(4)*(y(2))^2+data.p(7)*data.p(5)*chunk;
        end
        function [Y1,Y2]=split_matrix(Y,data)
            % Since the matrix of snapshot Y created has both y1 and y2,
            % this function takes odd rows of snapshot matrix to make
            % matrix Y1 and even rows of snapshot matrix to make Y2
            for i=1:data.t_points
                Y1(i,:)=Y(2*i-1,:);
                Y2(i,:)=Y(2*i,:);
            end
        end
        function Y = make_matrix(system,data,i)
            %This function evaluates the system at t_points time instances
            %and stores it into a matrix Y
            for j=1:data.t_points
                t=data.tvector(1,j);
                response= deval(system, t);
                Y(2*j-1,i)=response(1);
                Y(2*j,i)= response(2);
            end
        end
        function [Y,Y1,Y2]=create_snapshots(data,p_train)
            % This function takes the vector of sampling points generated
            % using LHSdesign and calculates system repsponses for each
            % sampling point and stores it into a snapshot matrix Y, and
            % finally splits the Y matrix using previous function to get Y1
            % and Y2.
            
            % This if condition is applied to make this function applicable
            % not just for a whole training set, but also one training
            % example, such as when we use this function in the other
            % function create_approx_real.
            if size(p_train,1)== sum(data.n)
                mysample=transpose(p_train);
            else
                mysample=p_train;
            end
            for i=1:size(mysample,1)
                for k=1:data.t_points
                    t=data.tvector(1,k);
                    m=data.interpolate(t, mysample(i,:));
                    for j=1:size(data.params_keys,2)
                        if data.params_keys(j)~=0
                            data.p(j)=m(1);
                            m(1)=[];
                        end
                    end

                system = ode15s(@ode_without_control,data.tvector,...
                data.y0,data.ode_options,data);
                    y= deval(system, t);
                    Y(2*k-1,i)=y(1);
                    Y(2*k,i)= y(2);
                end
            end
            [Y1,Y2]=split_matrix(Y,data);
        end 
        function R2= Rsquared(data,approx,real)
            % The function to calculate coefficient of determination R2. It
            % takes two vectors, one with real/original responses and one
            % with approximated responses and calculates correlation
            % coefficient. It then squares it to get coefficient of
            % determination and takes the first value of the resulting
            % vector which is our R2. 
            R=corrcoef(approx,real);
            R2=R.*R;
            R2=R2(1,2);
        end
        function y = maxabs(data,approx,real)
            % The function calculates Maximum Absolute error. It
            % takes two vectors, one with real/original responses and one
            % with approximated responses and calculates error. It then
            % uses in-built max function to calculate the max of errors.  
            e=approx-real;
            y = max(abs(e(:)));
        end
        function y = relmaxabs(data,approx,real)
            % The function calculates Maximum Absolute error. It
            % takes two vectors, one with real/original responses and one
            % with approximated responses and calculates error. It then
            % uses in-built max function to calculate the max of errors.  
            e=(approx-real)/real;
            y = max(abs(e(:)));
        end
        function Ei=commulative_energy(data,s)
            % This function calculates the vector of commulative energies
            % for the singular values given in s.  
            s_transpose=s';
            s_squared= s_transpose.^2;
            Etotal = sum(s_squared(1:data.points));
            Ei=zeros(data.points,1);
            for i=1:data.points
                Ei(i,1)= sum(s_squared(1:i))/Etotal;
            end
        end
        function [B]=rbf(data,A,p_vec)
            % Interpolation by the use of RBFs (Linear, Cubic and Thin Plate
            % splines type. The input of this function is p_vect, which is vector
            % of sampling points and the output of this function is B, 
            % the basis of interpolation. The type of RBF is predefined
            % in the data structure, where L stands for linear spline ,C
            % for cubic spline and T for Thin Plate Spline
            % Normalization of P [0 1]
            for j=1:size(p_vec,1)
                minP(j,1)=min(p_vec(j,:));
                maxP(j,1)=max(p_vec(j,:));
                for i=1:size(p_vec,2)
                    x(j,i)=(p_vec(j,i)-minP(j))/(maxP(j)-minP(j));
                end
                
            end
            N=size(x,2);
            for i=1:N
                for j=1:N
                    if data.rbf_type == 'L' %Linear spline
                        G(i,j)=sum((x(:,i)-x(:,j)).^2).^0.5;
                    elseif data.rbf_type == 'C' %Cubic spline
                        G(i,j)=sum((x(:,i)-x(:,j)).^2).^1.5;
                    else   %Thin Plate Spline
                        G(i,j)=sum((x(:,i)-x(:,j)).^2)*log(sum((x(:,i)-x(:,j)).^2));
                    end
                end
            end
            Bt=inv(G')*A';
            B=Bt';
        end
        function [pX]=create_pX(data,p_trial)
            %This function creates pX that is then used to calculate the
            %matrix of coefficients G for RBF. 
            key_set={};value_set={};
            for i=1:size(data.params_keys,2)
                if data.params_keys(i) ==1
                    for j=1:data.n(i)
                        key_set{end+1}=int2str(i);
                        value_set{end+1}=[data.params_values(i,:)];
                    end
                end
            end
            M=containers.Map(key_set,value_set);
            for i= 1:sum(data.n)
                temporary_var= key_set(i);
                temporary_var=char(temporary_var);
                K=M(temporary_var);
                newentry=(p_trial(i)-K(1))/((K(2)-K(1)));
                pX(i,:)=newentry;
            end
        end
        function G=pod_G_vec(data,p_vec,pX)
            % Function that constructs G vector as function of given
            % parameters, pX comes from the function create pX and p_vec is
            % the vector of ns sample points
            N=size(p_vec,2); % The number of generated snapshots
            %Normalization of p
            for j=1:size(p_vec,1)
                minP(j,1)=min(p_vec(j,:));
                maxP(j,1)=max(p_vec(j,:));
                for i=1:size(p_vec,2)
                    x(j,i)=(p_vec(j,i)-minP(j))/(maxP(j)-minP(j));
                end
            end
            if data.rbf_type == 'L' %Linear
                gi=inline('(sum((x-y).^2).^0.5)');
            elseif data.rbf_type == 'C' %Cubic
                gi=inline('(sum((x-y).^2).^1.5)');
            else %Thin Plate Spline
                gi=inline('(sum((x-y).^2)).*log((sum((x-y).^2)))');
            end
            value=pX;
            for k=1:N
                G(k,1)=gi(value,x(:,k));
            end
        end
        function [approx_y1,real_y1,approx_y2,real_y2]= ...
            create_real_approx(data,phi,B,G)
            % This function takes the inputs phi (basis from SVD), B (Basis
            % from RBF) and G (coefficient vector of RBF) and creates an
            % approximate surface of surrogate model for an arbitrary value
            % of varying parameters, that are stored in p_trial. It also 
            % calculates the original system responses and then provides
            % four matrices in the end, two for approximated y1 and y2, and
            % two for real y1 and y2.
            approx_y=phi*B*G;
            [approx_y1,approx_y2]=split_matrix(approx_y,data);
            [~, real_y1, real_y2]=create_snapshots(data,data.p_trial);
        end
%% PLOTS 
 % This subsection is devoted to making plots for visualization of the
 % results that were obtained in the previous subsection. 
         function plot_actual_surface(data)
             %This function solves the ode system and present the plot of
             %solution of original ODE with y1 and y2 plotted on the same
             %graph.
             system = ode15s(@ode_without_control,[data.t0,data.Tend],...
                 data.y0,data.ode_options,data);
             y_t=deval(system,data.tvector);
             figure
             hold on
             plot(data.tvector,y_t(1,:))
             plot(data.tvector,y_t(2,:), 'LineStyle','--')
             legend({'y_1','y_2'})
             hold off
             xlabel('t'), ylabel('y'), title('Actual Surface')
         end
         function plot_rank_approximation(data,rank_required,snapshot_matrix)
             % This function plots the approximated surface for for all the
             % sample points for the rank provided as input. It calcualtes
             % SVD itself and approximates the surface using the given rank
             [u,s,v]=svd(snapshot_matrix);
             approx_y=u(:,1:rank_required)*s(1:rank_required,1: ...
                 rank_required)*v(:,1:rank_required)';
             [approx_y1,approx_y2]=split_matrix(approx_y,data);
             plot(data.tvector,approx_y1)
             plot(data.tvector,approx_y2)
             xlabel('t'), ylabel('y')
             title(['Rank ',num2str(rank_required),' approximation'])
             legend({'approx\_y_1','approx\_y_2'},'FontSize',6)
             R2=data.Rsquared(approx_y,snapshot_matrix);
             e=approx_y-snapshot_matrix;
             maxerror=data.maxabs(approx_y,snapshot_matrix);
             meanerror = mae(e);
             fprintf('for %d singular value/s : R-Squared= %4.2f ,Mean Absolute Error = %d and Maximum Absolute Error= %d\n',rank_required,R2,meanerror,maxerror)
         end
         function plot_singular_values(data,s)
             % This function plots the ratio of all singular values (stored
             % in vector s) with first singular value which is the largest
             s_diag=diag(s);
             ratio=s_diag./s_diag(1,1);
             x = linspace(1,data.points,data.points);
             plot(x,ratio,'-*', 'Color','b')
             xlabel('number'), ylabel('ratio \sigma_i/\sigma_1'),...
                 title('Ratio of Singular values of y')
         end
         function plot_commulative_energy(data,s)
             % This function plots the commulative energy of singular values
             % s, that was computed in the function "commulative energy"
             s=diag(s);
             x = linspace(1,data.points,data.points);
             Ei=commulative_energy(data,s);
             plot(x,Ei,'-x', 'Color','k')
             xlabel('number'), ylabel('E(i) of POD basis'),title('E(i)')
         end
         function plot_approx_real(data, phi, B,G)
             % Just like the function create_real_appprox, this function
             % calculates responses from both original and surrogate model
             % for a single sample point and plots them on the same graph
             % for comparision.
             [approx_y1,real_y1,approx_y2,real_y2]= ...
             create_real_approx(data,phi,B,G);
             figure
             grid on
             subplot(2,1,1)
             hold on
             tvector_final=linspace(data.t0,data.Tend,data.t_points);
             plot(tvector_final, approx_y1)
             plot(tvector_final, real_y1','LineStyle','--')
             legend({'approximated', 'Original'})
             title(['comparision of approximated and actual surface y_1 for Rank ',num2str(data.rank)])
             xlabel('t'), ylabel('y')
             hold off
             grid minor

             subplot(2,1,2)
             hold on
             tvector_final=linspace(data.t0,data.Tend,data.t_points);
             plot(tvector_final, approx_y2)
             plot(tvector_final, real_y2','LineStyle','--')
             legend({'approximated', 'Original'})
             title(['comparision of approximated and actual surface y_2 for Rank ',num2str(data.rank)])
             xlabel('t'), ylabel('y')
             hold off
             grid minor
         end
 %% METHODS FOR OPTIMIZATION
 % This subsection is devoted to organize functions related to optimization
 % problem
         function data = Problem(b)
             %This function sets up the optimization parameter b for
             %initializing th problem
          data.p_trial=[];
            %Adding a value of p_trial
            for i=1:size(data.params_keys,2)
                if data.params_keys(i)~=0
                    data.p_trial=[data.p_trial,repmat(data.p(i),1,data.n(i))];
                end
            end

%              if nargin == 0 % use default init. guess if not input is given
%                  a=str2double(data.key_set);
%                  for i=1:data.change_var
%                      k=a(i);
%                      b(i) = data.p(k);
%                  end
%              end
             data.tvector=create_t_vector(data);
         end
         function f = interpolate(data, t, v)
             %function for interpolation;
             %if there are more than one optimization parameters, this
             %function interpolates each parameter seperately and finally
             %produces a vector whose length is equal to number of changing
             % variables "change_var".
             f=[];
             for i=1:size(data.n,2) %loops to check for all values in data.n
                 if data.n(i)~=0  %If the value is non-zero and non-one, it 
                     %takes that many(=value) columns from v vector,
                     %interpolates and deletes those entries from v vector.
                     %Here it is assumed that all the values in vvector
                     %occur sequentially. 
                     if data.n(i)~=1
                         interval=linspace(data.t0,data.Tend,data.n(i));
                         z=v(1:data.n(i));
                         f(end+1)=interp1(interval,z,t,data.method);
                     else
                         f(end+1)=v(1);
                     end
                     v(1:data.n(i)) = [];
                 end
             end
             
         end
         function u = u(data, t)
             % This function interpolates the optimization parameter b at
             % every iteration of optimization algorithm without asking 
             % for b as input.
             u = data.interpolate(t, data.b);
         end
         function u0 = u0(data, t)
             % This function interpolates the optimization parameter b at
             % initail stage of optimization algorithm without asking 
             % for b0 as input.
             u0 = data.interpolate(t, data.b0);
         end
         function dy = ode_with_control(data,t,y)
             % This function represents ODE with control, where each
             % control parameter is replaced by its interpolated
             % counterpart before generating system responses. 
             p=data.p;
%              a=str2double(data.key_set);
%              for i=1:data.change_var
%                  k=a(i);
%                  m=data.u(t);
%                  p(k)=m(i);
%              end
                 m=data.u(t);
                 for i=1:size(data.params_keys,2)
                     if data.params_keys(i)~=0
                         p(i)=m(1);
                         m(1)=[];
                     end
                 end
             dy=zeros(2,1);
             chunk=y(1)*(1-exp(-p(6)*y(1)))*y(2);
             dy(1)=p(1)*y(1)+p(2)*y(1)^2+p(5)*chunk;
             dy(2)=p(3)*y(2)+ p(4)*(y(2))^2+p(7)*p(5)*chunk;
         end
         function setControl(data, b)
             %This function does not only help in setting the value of b,
             %but also clears cache from the previous calculations, so that
             %new value of criteria is calculated at each iteration; if it
             %is not done, the value of criteria will be retrieved from
             %cache at each step, without being updated and it will produce
             %wrong results in optimization. 
             %This function first checks if the new value of b is different
             %than before, and if so, it invalidates the cache, so that the
             %direct problem will be solved again. 
             if any(data.b ~= b)
                 data.b = b;
                 data.cache = {};
             end;
         end
         function [t, y] = direct(data)
             % This function obtains system responses through direct method
             % using ODE; it first checks if the system response is there
             % in the cache, if not, it calculates it and saves it in
             % cache. 
             if ~isempty(data.cache)
                 t = data.cache{1};
                 y = data.cache{2};
             else
                 [t, y] = ode15s(@data.ode_with_control, ...
                 data.tvector, data.y0);
                 data.cache = {t, y};
             end
         end
         function [t,y]=surrogate_response(data)
             % This function performs the necessary steps to get the system
             % responses using surrogate model. It uses some functions from
             % the first subsection.
             pX=create_pX(data,data.b);
             G=data.pod_G_vec(data.p_train,pX);
             surrogatemodel_response=data.phi*data.B*G;
             % Prepare surrogate model
             [s1,s2]=split_matrix(surrogatemodel_response,data);
             y=[s1,s2];
             t=transpose(data.tvector);
         end
         function psi0 = criteria(data)
             %This function calculates the value of optimization criteria
             %for both direct method and surrogate model.
             %response_type= 1 for surrogate model, 2 for direct method
             if data.response_type == 1
                 [t,y]= surrogate_response(data);
             elseif data.response_type == 2
                 [t, y] = data.direct();
             else
                 [t, y] = data.direct();
             end
             i = data.optY;
             psi0 = trapz(t,(y(:,i) - data.yd(i)).^2);
         end
         function psi0 = optTarget(data, b)
             %This function helps to set b, such that the optimization
             %criteria will change, as the value of optimization parameter
             %b changes in each iteration. 
             data.setControl(b);
             psi0 = data.criteria();
         end
         function psi1 = constraint(data)
             %This function calculates the value of optimization constraint
             %for both direct method and surrogate model.
             %response_type= 1 for surrogate model, 2 for direct method
             if data.response_type==1
                 [t,y]= surrogate_response(data);
             elseif data.response_type==2
                 [t, y] = data.direct();
             else
                 [t, y] = data.direct();
             end
             i = abs(data.optY-3); % 1->2 || 2->1
             yDiff = y(:,i) - data.yMax(i);
             psi1 = trapz(t,(abs(yDiff) + yDiff).^2);
         end
         function [c, ceq] = optConstraint(data, b)
            %This function helps to set b, such that the optimization
            %constraint will change, as the value of optimization parameter
            %b changes in each iteration.
             data.setControl(b);
             c = [];
             ceq = data.constraint();
         end
         function optimize(data, yConstraint, bConstraint)
             % This is the main function for optimization. It takes too
             % inputs yContraint (standing for optimization constraint
             % and bConstraint (standing for bilateral constraints), both
             % of them should be logical values 'True' or 'False', if they
             % should be considered, or not, respectively. \
             
             nonLinCon = [];
             
             %if the value of yConstraint is given to be true, the if 
             %condition sets up non-linea-contraint using the function of
             %optimization constraint
             if yConstraint
                 nonLinCon = @data.optConstraint;
             end

             bLower = [];
             bUpper = [];
                          
             %if the value of bConstraint is given to be true, the if 
             %condition sets up billateral-contraints using the equations
             %given in the problem statement

             %See description in p_vec function.
            for i=1:size(data.n,2)
                lb1(i)=data.params_values(i,1); %lowebounds
                ub1(i)=data.params_values(i,2); %upperbounds
            end
            
            if bConstraint
            for i=1:size(data.n,2)
                if data.n(i)~=0
                    bLower=[bLower,repmat(lb1(i),1,data.n(i))];
                    bUpper=[bUpper,repmat(ub1(i),1,data.n(i))];
                end
            end
            %since optimization tends to perform worse on boundary points,
            %we reduce the domain of our optimzation
            %let l be the current width of our domain and l_new be the
            %reduction parameter
            l=bUpper-bLower; l_new=0.1*l;
            bLower=bLower+l_new; bUpper=bUpper-l_new;
            end
             %options = optimset('LargeScale', 'on','MaxFunEvals', 3000);
             options = optimoptions( 'fmincon' , 'Algorithm' , 'sqp');
             
             %Note: if there is any kind of constraint, matlab package
             %fmincon is used for optimization, whereas, in absence of
             %constraints, fminsearch is used. 
             if yConstraint || bConstraint
                 b = fmincon(@data.optTarget,data.b0,[],[],[],[],bLower,...
                 bUpper,nonLinCon,options);
             else
                 b = fminsearch(@data.optTarget, data.b0, []);
             end
             data.setControl(b);
         end

    end 
end