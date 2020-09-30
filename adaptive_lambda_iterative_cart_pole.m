clear
clc
close all

tStart = tic;


syms xi1 xi2 xi3 xi4 xi5 dt
m1 = 1;
m2 = 0.3;
l = 0.5;
g = 9.81;
F_real= [xi3;
         xi4;
         1/(m1+m2*(1-cos(xi2)^2))*(l*m2*sin(xi2)*xi4^2+xi5+m2*g*cos(xi2)*sin(xi2));
         -1/(l*m1+l*m2*(1-cos(xi2)^2))*(l*m2*cos(xi2)*sin(xi2)*(xi4)^2+xi5*cos(xi2)+(m1+m2)*g*sin(xi2));
         0];
F_real = simplify(F_real);

p = 5;

[Phi,Psi_p,JPhi] = compute_Phi_and_JPhi(p,F_real,[xi1 xi2 xi3 xi4 xi5],dt);
disp( 'Done!' );

tEnd = toc(tStart)

save('cart_pole_model.mat','Phi','Psi_p','JPhi')

%% Model:
%% Part 1: 
clear
clc
close all
load('cart_pole_model.mat','Phi','Psi_p','JPhi')
p=5;

tic
dt = 0.01;
T = 2;
iter_max = ceil(T/dt);

x0 = [0;0;0;0];
x_target = [1;pi;0;0];
u = 0.1*ones(iter_max,1);
% u = 0.1*rand(iter_max,1);



rejected = false;
lambda_coeff = 0.01;
continue_iterating = true;
while continue_iterating  % Iterating scheme: U --> trajactory [x(0)...x(N)] --> H --> U' 
    if rejected==false
    % Apply U to the real system   
    x = x0;
    x_traj = [];
    
    for iter = 1:iter_max
        % Prepare A,B to later calculate H
        R_big = JPhi(dt,x(1),x(2),x(3),x(4),u(iter));          
        R_store{iter} = R_big(1:4,1:4);     % A in Ax+Bu    
        B_store{iter} = R_big(1:4,5:5);     % B in Ax+Bu 
        
        % The actual trajectory using adaptive step size mechanism
        [~, x_trajJ_fine] = adaptive_taylor(p,Phi,Psi_p,[0 dt],[x;u(iter)]); 
        x = x_trajJ_fine(end,:)'; 
        x = x(1:4);
        x_traj = [x_traj, x];      
    end
    xn = x;
    
    cost_old = norm(xn-x_target);  
    if cost_old <= 0.01 
        cost_old
        continue_iterating = false;
    end
%     vav
    
    % Calculate H:
    H = B_store{1};
    for iter = 2:iter_max
        H = [R_store{iter}*H, B_store{iter}];  
    end
    end

    lambda = lambda_coeff*cost_old^1.5;
    du_proposal = -(H'*H+lambda*eye(iter_max))\(H'*(xn-x_target));
    u_proposal = u + du_proposal;
%     u = u_proposal;
    
%     display('here1')
    
    % simulate real system using proposed u
    x_sim = x0;
    for ite = 1:iter_max
        [~, x_trajJ_fine_sim] = adaptive_taylor(p,Phi,Psi_p,[0 dt],[x_sim;u_proposal(ite)]); 
        x_sim = x_trajJ_fine_sim(end,:)'; 
        x_sim = x_sim(1:4);
    end
%     x_sim
    cost_proposal_actual = norm(x_sim-x_target);
    
    if cost_proposal_actual < cost_old      
        lambda_coeff = 0.95*lambda_coeff;
        u = u_proposal;
        rejected = false;     % accept u_proposal due to cost benefits
    else                                    
        lambda_coeff = 1.2*lambda_coeff;
        rejected = true;      % reject u_proposal
    end
    
%     display('here2')
    [cost_old, lambda_coeff, rejected]
% avv

end
dt*norm(u)^2
step1_u = u;
disp( 'Done1!' );
toc
%% part 2:
tic
u = step1_u;
% dt = 0.01;
% T = 2;
% iter_max = ceil(T/dt);

for iters = 1:800
    % Apply U to the real system   
    x = x0;
    x_traj = [];
    
    for iter = 1:iter_max
        % Prepare A,B to later calculate H
        R_big = JPhi(dt,x(1),x(2),x(3),x(4),u(iter));          
        R_store{iter} = R_big(1:4,1:4);     % A in Ax+Bu    
        B_store{iter} = R_big(1:4,5:5);     % B in Ax+Bu 
        
        % The actual trajectory using adaptive step size mechanism
        [~, x_trajJ_fine] = adaptive_taylor(p,Phi,Psi_p,[0 dt],[x;u(iter)]); 
        x = x_trajJ_fine(end,:)'; 
        x = x(1:4);
        x_traj = [x_traj, x];      
    end
    xn = x;
    
    target = norm(xn-x_target);
    [iters target dt*norm(u)^2]
    
%     cost_old = norm(xn_real-x_target)^2;

    % Calculate H:
    H = B_store{1};
    for iter = 2:iter_max
        H = [R_store{iter}*H, B_store{iter}];  
    end

    % Lambda Adaptive Scheme
    mu = 5;
    options = optimset('display','off');
    du_proposal = quadprog((1+mu)*eye(iter_max,iter_max),u,[],[],H,x_target - xn,[],[],[],options);
%     u_proposal = u + du_proposal;
    u = u + du_proposal;
    
    if norm(du_proposal)<0.001
        break
    end
end
step2_u = u;
toc

%%
figure
hold on
% plot(x_traj_real(1,:))
% plot(x_traj_real(2,:))
plot(step1_u)
plot(step2_u)

%% Animation
x_traj2 = x_traj;

% end

close all

% L = 0.3;
L = 1;
H = 0.2; 
W = 0.4;
% 
% LB = min(x_traj2(1,:))-1;
% UB = max(x_traj2(1,:))+1;

for k = 1:1:length(x_traj2)
    
    clf    
    hold on
    grid on
    
    % plot cart
    rectangle('Position', [x_traj2(1,k)-0.5*W -0.5*H  W H],'FaceColor','b','EdgeColor','b')    
    
    % plot pole
    plot([x_traj2(1,k), x_traj2(1,k)+L*sin(x_traj2(2,k))],[0 -L*cos(x_traj2(2,k))], 'r','LineWidth',2.5);
    
    % plot mass
    scatter(x_traj2(1,k)+L*sin(x_traj2(2,k)), -L*cos(x_traj2(2,k)),100,'k','filled')
    
    axis(2*[-1 1 -1 1])
%     pause(0.1)
    drawnow;

end


