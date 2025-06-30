% Constrained Model Predictive Control Function
% A: A matrix of the linear considered dynamic system
% B: B matrix of the linear considered dynamic system
% Q: status weight into the MPC cost function
% R: inputs weight into the MPC cost function
% S: final state weight (prediction horizon instant time) into the MPC cost function
% N: prediction horizon
% umin: inputs lower limit (scalar)
% umax: inputs upper limit (scalar)
% X: measured status at the current instant timet
function u = myvelocitympc(A,B,C,Q,R,S,N,umin,umax,u_bar,x_sat_max, x_sat_min,xref,input)
    m=size(B,2)
    n=size(A,1)
    p=size(C,1)

    %% Compute augmented state vector
    % e(k) = y⁰ - y(k)
    e = input(1:2);
    % δx(k) = x(k) - x(k-1)
    delta_x = input(3:4);
    %% errore qui sotto
    u_prev = input(5:6);

    % Augmented state: [δx(k); e(k)]
    x_aug = [delta_x; e];

      %% Point 7 matrix definition
    Av = [A, zeros(2,2);
      -C*A, eye(2)];
    Bv = [B;-C*B];

    %% Q and R matrices for Open-Loop MPC (with Q1)
    Qsig = blkdiag(kron(eye(N-1),Q),S);
    Rsig = kron(eye(N),R);

    %% Matrix for velocity form
    Q_aug = blkdiag(Q, eye(2));
    S_aug = blkdiag(S, eye(2));
    
    %% Compute prediction matrices (calligraphic notation)
    % Qσ matrix
    Qsig = blkdiag(kron(eye(N-1), Q_aug), S_aug);
    % Rσ matrix  
    Rsig = kron(eye(N), R);

    %% A and B matrices for the Open-Loop MPC
    % A matrix
    Asig = Av;
    for i = 2:N
        Asig = [Asig; Av^i];
    end

    % B matrix
    Bsig = [];
    temp = [];
    for i = 1:N
        temp = zeros(size(Bv,1)*(i-1),size(Bv,2))
        for j = 0:N-i
            temp = [temp; Av^(j)*Bv];
        end
        Bsig = [Bsig temp];
    end

    %% H, F 
    H = Bsig'*Qsig*Bsig + Rsig;
    F = Asig'*Qsig*Bsig;
    ft = x_aug'*F;
    f=ft';

    %% input and status constraints definition
    lb = [repmat(umin, N*m,1)];
    ub = [repmat(umax, N*m,1)];
    
    xlim = 15;
    
    x_sat_max_aug = [xlim;xlim;+inf;+inf];
    x_sat_min_aug = [-xlim;-xlim;-inf;-inf];

    xl = repmat(x_sat_min_aug, N,1)
    xu = repmat(x_sat_max_aug, N,1)

    A_constraint = [Bsig; -Bsig];
    b_constraint = [xu - Asig * x_aug; -(xl - Asig * x_aug)];


    options = optimset('Algorithm', 'interior-point-convex','Diagnostics','off', ...
        'Display','off');
    %solve the quadratic programming problem
    U = quadprog(H,f,A_constraint,b_constraint,[],[],lb,ub,[],options);
    %U = quadprog(H,f,[],[],[],[],lb,ub,[],options);
    %U = quadprog(H,f,[],[],[],[],[],[],[],options);
    %get the optimal input value (the receding horizon principle is applied)
    du = U(1:m);
    disp(size(Asig))
    u = du + u_prev;
end