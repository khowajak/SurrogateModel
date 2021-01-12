%% Plots

%% Plotting Obtained results 
figure
plot_actual_surface(data)
%% Plots of Singular Values
figure
subplot(2,1,1)
plot_singular_values(data,svd_s)
% Computation of E(i)
subplot(2,1,2)
plot_commulative_energy(data,svd_s)

%% Plot approximated surface
figure
for i=1:3
    rank=i;
    subplot(3,1,i)
    plot_rank_approximation(data,rank,snapshot_matrix)
end
%% Plot actual vs real for a trial point
plot_approx_real(data, phi, B,G)

%% Plot of Optimal Control Function
%     [t,y] = data.direct();
%     plot(t, data.u0(t)*0.9, '--');
%     plot(t, data.u0(t), ':c');
%     plot(t, data.u0(t)*1.1, '--');
%     plot(t, data.u(t));
%     title('Control function');
%     xlabel('t');
%     ylabel('u');
%     
% %%

figure
subplot(1,2,2)
interval= linspace(data.t0, data.Tend,2);
t=data.tvector;
Ut = interp1(interval, data.b, t);

[t,y]=data.direct();
hold on
%plot(t, U0*0.9,'--')
plot(t,Ut)
%plot(t, U0*1.1,'--')
hold off
grid on
title('Optimal control function')
ylabel('Control Function, u_1')
xlabel('time, t')

%% Plot error vector
figure
plot(linspace(1,count,count),epsilon_vec, '-*')
%plot(linspace(1,3,3),[e(1),e(3),e(4)], '-*')
title('RMAE for each iteration')
xlabel('iterations')
ylabel('\epsilon')
box on

%% Plot of Optimized surface
%[t,y_t]=data.surrogate_response();
[t,y_t]=data.direct();
figure
hold on
plot(data.tvector,y_t(:,1))
plot(data.tvector,y_t(:,2), 'LineStyle','--')
legend({'y_1','y_2'})
xlabel('t'), ylabel('y')
title('Optimized surface for constrained problem ')
hold off

%%
figure
plot(y_t(:,2),y_t(:,1));
%% Plot of Sampling Points
data.sampling_type='LHS';
lhs=create_p_samplingmethod(data,10);
data.sampling_type='SLHS';
slhs=create_p_samplingmethod(data,10);
data.sampling_type='RS';

figure
subplot(1,3,1)
scatter(lhs(1,:),lhs(2,:))
title('LHS')
xlabel('p_4'); ylabel('p_7')
subplot(1,3,2)
scatter(slhs(1,:),slhs(2,:))
title('SLHS')
xlabel('p_4'); ylabel('p_7')
subplot(1,3,3)
scatter(rs(1,:),rs(2,:))
xlabel('p_4'); ylabel('p_7')
title('RS')
rs=create_p_samplingmethod(data,10);