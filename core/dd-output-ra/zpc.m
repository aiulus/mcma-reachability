function R = zpc()
    % Adapted from: https://github.com/aalanwar/Data-Driven-Predictive-Control/blob/main/ZPC.m
    %------------- BEGIN CODE --------------
    % dimension of x
    dim_x = 5;
    
    % System in cont time
    A = [-1 -4 0 0 0; 4 -1 0 0 0; 0 0 -3 1 0; 0 0 -1 -3 0; 0 0 0 0 -2];
    B_ss = ones(5,1);
    C = [1,0,0,0,0];
    D = 0;
    % define continuous time system
    sys_c = ss(A,B_ss,C,D);
    % convert to discrete system
    samplingtime = 0.05;
    sys_d = c2d(sys_c,samplingtime);
    %number of trajectories
    initpoints =100;
    %number of steps for each trajectory
    steps =5;
    %Total number of samples
    totalsamples = initpoints*steps;
    %% initial set and input
    %reference input
    uref = 8;
    %reference output
    ref = inv(eye(5)-sys_d.A)*sys_d.B*uref; 
    
    %output constraint
    y_lb = [-10;2;-10;-10;-10]; 
    y_ub = [10;10;10;10;10]; 
    intc = interval(y_lb,y_ub);
    
    
    %initial point
    y0 = [-2;4;3;-2.5;5.5];
    
    X0 = zonotope(y0, 25*eye(dim_x));
    U = zonotope(uref - 1, 19);  % 20 - 1 = 19
    
    %noise zontope W (modeling noise)
    %less noise
    wfac=0.01;
    %more noise
    %wfac=0.1;
    W = zonotope([zeros(dim_x,1),wfac*eye(dim_x)]);
    
    for i=1:size(W.generators,2)
        vec=W.G(:,i);
        GW{i}= [ vec,zeros(dim_x,totalsamples-1)];
        for j=1:totalsamples-1
            GW{j+i}= [GW{i+j-1}(:,2:end) GW{i+j-1}(:,1)];
        end
    end
    % matrix zonotpe of noise w (M_w)
    Wmatzono= matZonotope(zeros(dim_x,totalsamples),GW);
    
    %measurement noise
    %less measurement noise
    vfac = 0.002;
    %more measurement noise
    %vfac = 0.02;
    V = zonotope([zeros(dim_x,1),vfac*eye(dim_x)]);
    CV = zeros(dim_x,totalsamples);
    for i=1:size(V.generators,2)
        vec=V.G(:,i);
        GV{i}= [ vec,zeros(dim_x,totalsamples-1)];
        for j=1:totalsamples-1
            GV{j+i}= [GV{i+j-1}(:,2:end) GV{i+j-1}(:,1)];
        end
    end
    % matrix zonotpe of noise v (M_v)
    % Vmatzono= matZonotope(CV,GV);
    G_w = cat(3, GW{:});
    Wmatzono = matZonotope(zeros(dim_x, totalsamples), G_w);
    
    G_v = cat(3, GV{:});
    Vmatzono = matZonotope(zeros(dim_x, totalsamples), G_v);

    
    AV = sys_d.A*V;
    % matrix zonotpe of  Av (M_Av)
    VAmatzono = sys_d.A*Vmatzono;
    
    
    % randomly choose constant inputs for each step / sampling time
    for i=1:totalsamples
        uopt(i) = randPoint(U);
    end
    
    %generate data from different trajectories with noise
    x0 = X0.center;
    x(:,1) = x0;
    index=1;
    for j=1:dim_x:initpoints*dim_x
        x(j:j+dim_x-1,1) = randPoint(X0);
        x_v(j:j+dim_x-1,1) =  x(j:j+dim_x-1,1) + randPoint(V);
        
        for i=1:steps
            utraj(j,i) = uopt(index);
            x(j:j+dim_x-1,i+1) = sys_d.A*x(j:j+dim_x-1,i) + sys_d.B*uopt(index) + randPoint(W);
            x_v(j:j+dim_x-1,i+1) =  x(j:j+dim_x-1,i+1) + randPoint(V);
            index=index+1;
        end
    end
    
    
    %prepeare Y_+ Y_-
    index_0 =1;
    index_1 =1;
    for j=1:dim_x:initpoints*dim_x
        for i=2:steps+1
            x_meas_vec_1_v(:,index_1) = x_v(j:j+dim_x-1,i);
            x_meas_vec_1(:,index_1) = x(j:j+dim_x-1,i);
            index_1 = index_1 +1;
        end
        for i=1:steps
            u_mean_vec_0(:,index_0) = utraj(j,i);
            x_meas_vec_0(:,index_0) = x(j:j+dim_x-1,i);
            x_meas_vec_0_v(:,index_0) = x_v(j:j+dim_x-1,i);
            index_0 = index_0 +1;
        end
    end
    % U_data is U_-, Y_0T is Y_- , Y_1T is Y_+
    U_data = u_mean_vec_0(:,1:totalsamples); %same as u
    Y_0T = x_meas_vec_0_v(:,1:totalsamples);
    Y_1T = x_meas_vec_1_v(:,1:totalsamples);     
    
    
    % plot simulated trajectory
    figure;
    subplot(1,2,1); hold on; box on; plot(x(1,:),x(2,:),'b'); xlabel('x_1'); ylabel('x_2');
    subplot(1,2,2); hold on; box on; plot(x(3,:),x(4,:),'b'); xlabel('x_3'); ylabel('x_4');
    close;
    
    %prepare M_Sigma which is a set of [A B]
    AB = (Y_1T + -1* Vmatzono + -1*Wmatzono+VAmatzono)*pinv([Y_0T;U_data]);
    
    
    %double check if the true A B is part of M_Sigma
    intAB11 = intervalMatrix(AB);
    intAB1 = intAB11.int;
    intAB1.sup >= [sys_d.A,sys_d.B]
    intAB1.inf <= [sys_d.A,sys_d.B]
    % check the rank of the data
    data_rank = rank([Y_0T;U_data])
    
    N = 2; % Time horizon
    maxsteps = 80;

    y_t = zeros(dim_x, maxsteps);      % or NaNs if you prefer
    y_t_model = zeros(dim_x, maxsteps);
    YPred = zeros(dim_x, maxsteps);
    
    for timesteps = 1:maxsteps
        if timesteps == 1
            % set the initial output to y0
            y_t(:,timesteps) = y0;
            y_t_model(:,timesteps) = y0;
            YPred(:,1) = y0;
        end
        
    
        % sdpvar variables
        uopt = sdpvar(1*ones(1,N),ones(1,N));
        y = sdpvar(5*ones(1,N+1),ones(1,N+1));
        alpha_u = sdpvar(1,N);
        sinf = sdpvar(5*ones(1,N+1),ones(1,N+1));
        ssup = sdpvar(5*ones(1,N+1),ones(1,N+1));
        R={};
        R{1} = zonotope([y_t(:,timesteps)]);
        %set the first constraint as y_t = current y
        Constraints = [y_t(:,timesteps) == y{1}];%,...
    
        
        
        for i = 1:N
            %compute the reachable set for ZPC
            card_cen = [R{i}.center;uopt{i}];
            genleni = size(R{i}.generators,2);
            %card_zono = zonotope([card_cen,[R{i}.generators;zeros(1,genleni)]]);
            Z_center = card_cen;
            Z_generators = [R{i}.generators; zeros(1, genleni)];
            card_zono = zonotope(Z_center, Z_generators);

            ABcard = intervalMatrix(AB)* card_zono;
            %R{i+1} = zonotope([ABcard.center,[ABcard.generators,W.generators,V.generators,AV.generators]]);%AB * card_zono + W_sdp;
            R{i+1} = zonotope(ABcard.center, [ABcard.generators, W.generators, V.generators, AV.generators]);
            %convert R to interval
            %extract center
            c = R{i+1}.G(:,1);
            
            %determine left and right limit of the reahable set (convert to
            %interval)
            delta = sum(abs(R{i+1}.G),2) - abs(c);
            leftLimit{i} = c - delta;
            rightLimit{i} = c + delta;
            
            %specify the constraints
            Constraints = [Constraints,...
                uopt{i} == U.center + alpha_u(i) * U.generators,...
                y{i+1} - sinf{i} == leftLimit{i},...
                y{i+1} + ssup{i} == rightLimit{i},...
                y{i+1}   - sinf{i} >= intc.inf,...
                y{i+1}   + ssup{i} <= intc.sup,...
                sinf{i} >= zeros(5,1),...
                ssup{i} >= zeros(5,1),...
                alpha_u(i) <= 1 , ...
                alpha_u(i) >= -1, ...
                ];
        end
    end
end

