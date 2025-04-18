function [uPred, uPred_model, y_t, y_t_model, execTimeZPC, execTimeRMPC] = ...
    runZPC(sys_d, y0, intc, r_u, r_y, U, N, maxsteps, timestep_plot, Qy, Qu, W, V, AV, M_Sigma)

    execTimeZPC=[];
    execTimeRMPC=[];

    for k = 1:maxsteps
        if k == 1
            % set the initial output to y0
            y_t(:,k) = y0;
            y_t_model(:,k) = y0;
            YPred(:,1) = y0;
        end        
        
        % sdpvar variables
        u = sdpvar(1*ones(1,N),ones(1,N));
        y = sdpvar(5*ones(1,N+1),ones(1,N+1));
        alpha_u = sdpvar(1,N);
        sinf = sdpvar(5*ones(1,N+1),ones(1,N+1));
        ssup = sdpvar(5*ones(1,N+1),ones(1,N+1));
        R={};
        R{1} = zonotope([y_t(:,k)]);

        %set the first constraint as y_t = current y
        Constraints = [y_t(:,k) == y{1}]; %,...       
        
        for i = 1:N        
            genleni = size(R{i}.generators,2);

            placeholder_u = r_u;  % safe numeric value
            dummy_cen = [R{i}.center; placeholder_u];
            dummy_gen = [R{i}.generators; zeros(1, genleni)];
            card_zono = zonotope(dummy_cen, dummy_gen);
    
            
            ABcard = [sys_d.A , sys_d.B]* card_zono;
            R{i+1} = zonotope([ABcard.center,[ABcard.generators,W.generators,V.generators,AV.generators]]);%AB * card_zono + W_sdp;
           
    
    
            % convert R to interval        
            c = R{i+1}.c; % extract center
            delta = sum(abs(R{i+1}.G),2);
            
            delta = sum(abs(R{i+1}.Z),2) - abs(c);
            leftLimit{i} = c - delta; % determine left limit of the reahable set
            rightLimit{i} = c + delta; % determine right limit of the reahable set
            
            %specify the constraints
            Constraints = [Constraints,...
                u{i} == U.center + alpha_u(i) * U.generators,...
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
        
        % chose the cost of ZPC
        Cost=0;
        for i=1:N
            Cost = Cost + (y{i+1}-r_y)'*Qy*(y{i+1}-r_y)+ (u{i}-r_u)'*Qu*(u{i}-r_u);
        end
        %solve ZPC
        options = sdpsettings('verbose',0,'solver','mosek');
        tic
        Problem = optimize(Constraints,Cost,options)
        execTimeZPC=[execTimeZPC,toc];
        Objective = double(Cost);
        uPred(k) = double(u{1})
        YPred(:,k+1) = double(y{2});
        %%
    
        %% save for plotting
        Rplotall{k} = interval(zonotope(R{2}.c, R{2}.G));
        %%  ploting
        if timestep_plot == k
            for i =1:N+1
                R_N{i} = zonotope(R{i}.c, R{i}.G);
                R_N_int{i} = interval(R_N{i});
                y_N{i} =double(y{i});
                if i<N+1
                    u_N{i} =double(u{i});
                end
            end
        end
	    
	    
	    
        %% ZPC given the model (RMPC-zono)
        % Control					
        alpha_u = sdpvar(1,N);
        sinf = sdpvar(5*ones(1,N+1),ones(1,N+1));
        ssup = sdpvar(5*ones(1,N+1),ones(1,N+1));
        R={};
        R{1} = zonotope([y_t_model(:,k)]);
        u_model = sdpvar(1*ones(1,N),ones(1,N));
        y_model = sdpvar(5*ones(1,N+1),ones(1,N+1));
        Constraints = [y_t_model(:,k) == y_model{1}];
        for i = 1:N
            %compute the reachable set for ZPC
            genleni = size(R{i}.generators,2);
            %card_cen = [R{i}.center;u{i}];        
            %card_zono = zonotope([card_cen,[R{i}.generators;zeros(1,genleni)]]);
            % Use numeric placeholder for center in zonotope creation
            placeholder_u = r_u;  % safe numeric value
            dummy_cen = [R{i}.center; placeholder_u];
            dummy_gen = [R{i}.generators; zeros(1, genleni)];
            card_zono = zonotope(dummy_cen, dummy_gen);
    
            ABcard = intervalMatrix(M_Sigma)* card_zono;
            R{i+1} = zonotope([ABcard.center,[ABcard.generators,W.generators,V.generators,AV.generators]]);%AB * card_zono + W_sdp;
            %convert R to interval
            %extract center
            c = R{i+1}.Z(:,1);      
        
            delta = sum(abs(R{i+1}.Z),2) - abs(c);
            leftLimit{i} = c - delta;
            rightLimit{i} = c + delta;
            
            
            Constraints = [Constraints,...
                u_model{i} == U.center + alpha_u(i) * U.generators,...
                y_model{i+1} - sinf{i} == leftLimit{i},...
                y_model{i+1} + ssup{i} == rightLimit{i},...
                y_model{i+1}   - sinf{i} >= intc.inf,...
                y_model{i+1}   + ssup{i} <= intc.sup,...
                sinf{i} >= zeros(5,1),...
                ssup{i} >= zeros(5,1),...
                alpha_u(i) <= 1 , ...
                alpha_u(i) >= -1, ...
                ];
        end
        
        
        
        Cost_model=0;
        for i=1:N
            Cost_model = Cost_model + (y_model{i+1}-r_y)'*Qy*(y_model{i+1}-r_y)+ (u_model{i}-r_u)'*Qu*(u_model{i}-r_u);
        end
        options = sdpsettings('verbose',0,'solver','mosek');
        tic
        Problem = optimize(Constraints,Cost_model,options)
        execTimeRMPC=[execTimeRMPC,toc];
        Objective = double(Cost_model);
        uPred_model(k) = double(u_model{1});
        YPred_model(:,k+1) = double(y_model{2});
    
        
        % apply the optimal control input to the plant 
        w_point = randPoint(W);
        v_point = randPoint(V);
        y_t(:,k+1) = sys_d.A * y_t(:,k) + sys_d.B * uPred(k) + w_point +v_point - sys_d.A *v_point;
        y_t_model(:,k+1) = sys_d.A * y_t_model(:,k) + sys_d.B * uPred_model(k) + w_point +v_point - sys_d.A *v_point;
        
        yt2ref(k)= norm(y_t(:,k)-r_y,2)
        yt2ref_model(k)= norm(y_t_model(:,k)-r_y,2)
        halt = 1;
    end
end

