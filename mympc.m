% Constrained Model Predictive Control Function
% A: A matrix of the linear considered dynamic system
% B: B matrix of the linear considered dynamic system
% Q: status weight into the MPC cost function
% R: inputs weight into the MPC cost function
% S: final state weight (prediction horizon instant time) into the MPC cost function
% N: prediction horizon
% umin: inputs lower limit (scalar)
% umax: inputs upper limit (scalar)
% X: measured status at the current instant time
function u = mympc(A,B,Q,R,S,N,umin,umax,u_bar,x_sat_max, x_sat_min,xref,x)
    m=size(B,2);
    n=size(A,1);
    
    %% Q and R matrices for Open-Loop MPC (with Q1)
    Qsig = blkdiag(kron(eye(N-1),Q),S);
    Rsig = kron(eye(N),R);

    %% A and B matrices for the Open-Loop MPC
    % A matrix
    Asig = A;
    for i = 2:N
        Asig = [Asig; A^i];
    end

    % B matrix
    Bsig = [];
    temp = [];
    for i = 1:N
        temp = zeros(size(B,1)*(i-1),size(B,2));
        for j = 0:N-i
            temp = [temp; A^(j)*B];
        end
        Bsig = [Bsig temp];
    end

    %% H, F 
    H = Bsig'*Qsig*Bsig + Rsig;
    F = Asig'*Qsig*Bsig;
    ft = x'*F;
    f=ft';
    %input and status constraints definition
    lb = [repmat(umin-u_bar, N*m,1)];
    ub = [repmat(umax-u_bar, N*m,1)];
    
    xlim = 15;
    
    %% Define input and status constraint

    xl = repmat(x_sat_min-xref, N,1)
    xu = repmat(x_sat_max-xref, N,1)

    A_constraint = [Bsig; -Bsig];
    b_constraint = [xu - Asig * x; -(xl - Asig * x)];








    % xl = [repmat(x_sat_min-xref, N,1)]
    % xu = [repmat(x_sat_max-xref, N,1)];
    % 
    % u_bar_vec = repmat(u_bar, N, 1);
    % b_constraint = [ xu - Asig*x - Bsig*u_bar_vec;
    %                   xl + Asig*x + Bsig*u_bar_vec ];
    % 
    % A_constraint = [ Bsig;-Bsig ];

    % 
    % u_bar_vec = repmat(u_bar, N, 1);
    % b_constraint = [ xlim*ones(N*n,1) - Asig*x - Bsig*u_bar_vec;
    %                  xlim*ones(N*n,1) + Asig*x + Bsig*u_bar_vec ];
    % 
    % A_constraint = [ Bsig;-Bsig ];

    % xlim = 15;
    % 
    % A_constraint = [ Bsig; -Bsig ];
    % b_constraint = [ xlim*ones(N*n,1) - Asig*x;
    %              xlim*ones(N*n,1) + Asig*x ];

    options = optimset('Algorithm', 'interior-point-convex','Diagnostics','off', ...
        'Display','off');
    %solve the quadratic programming problem
    U = quadprog(H,f,A_constraint,b_constraint,[],[],lb,ub,[],options);
    %U = quadprog(H,f,[],[],[],[],lb,ub,[],options);
    %U = quadprog(H,f,[],[],[],[],[],[],[],options);
    %get the optimal input value (the receding horizon principle is applied)
    u = U(1:m);
end