%% pure-pursuit QP + controller with reachablity set (no tube constraint) - 2D 

clc; clear; close all;

%% Parameters
ts_dt= 0.01; %control step, according to our system frequency (100 Hz)
fr_dt =0.04;%changing frame rate to be visual (25 fps)
T_time= 12;% total time
look_ahead_s = 0.001; % distance ahead of the closest point (m) ~1mm

v_max= 1; %0.577; % speed limit (m/s) - physical constraints - max for ur5e: 1m/s (magnitude in 3D)
a_max= 2.5; % acceleration limit (m/s^2) - max for ur5e: 2.5m/s^2 
N = ceil(T_time/ts_dt); %number of steps
e_tol_min= 0.001; % how much error is tolerable 
eps_p_bound = 0.0; % (m) position-noise bound per sample (for robustness-on, it should be >0)
eps_v_bound = 0.0; % (m/s) velocity-noise bound per sample


%% Reference path - random waypoint-based spline

n_wp = 6; % number of waypoints
wp = zeros(2,n_wp+1); % column 1 is the origin

theta = 2*pi*rand(1,n_wp); % random angles/directions in a circle (0, 2pi) - around the origin

%Scale positions to fit to real workspace - random targets inside an L_x m x L_y m rectangle
L_x = 0.5; % Lenght of the workspace in x axis - meter
L_y = 0.5; % Lenght of the workspace in y axis - meter

%converting directions into 2D coordinates
wp(1,2:end) = L_x*cos(theta); % x coordinates
wp(2,2:end) = L_y*sin(theta); % y coordinates

dLine = [0, vecnorm(diff(wp,1,2))]; % length of segments -  distance of difference between columns of wp. 0 for the first segment lenght 
sLine = cumsum(dLine); % cumulative sum of segment lenghts 
t_wp = sLine/sLine(end)*4*pi;  % Instead of real time, using cumulative distance along the path to define how far along each waypoint is, so: ~ (0,3pi)

% cubic spline through the points
ppx = spline(t_wp, wp(1,:)); % spline of t_wp and all columns in row 1 - x(t)
ppy = spline(t_wp, wp(2,:)); % y(t)



dppX = fnder(ppx,1);  % 1st derivative splines
dppY = fnder(ppy,1);
ddppX= fnder(ppx,2); % 2nd derivative splines
ddppY= fnder(ppy,2);


tau_grid = linspace(0, ppx.breaks(end), 150000);
p_ref  = [ ppval(ppx , tau_grid) ;
           ppval(ppy , tau_grid) ];%position at each time step
dp_ref = [ ppval(dppX, tau_grid) ;
           ppval(dppY, tau_grid) ];
ddp_ref= [ ppval(ddppX,tau_grid) ;
           ppval(ddppY,tau_grid) ];
v_ref  = vecnorm(dp_ref);


%% offline reachability checks
%using second derivative = 0 and finding optimums

sigma_off = (2*eps_p_bound)/ts_dt + eps_v_bound;  % same sigma as online
delta_list = []; bad_idx = [];

for i = 1:(numel(tau_grid)-1)
    % treat the reference as the current state
    p_k = p_ref(:,i);
    v_k = dp_ref(:,i);

    % same look-ahead rule as online
    vLA  = v_ref(i);
    sLA  = max(look_ahead_s, vLA*ts_dt);
    idxLA = i; s_acc = 0;
    while s_acc < sLA && idxLA < numel(tau_grid)
        s_acc = s_acc + norm(p_ref(:,idxLA+1) - p_ref(:,idxLA));
        idxLA = idxLA + 1;
    end
    p_LA_i = p_ref(:,idxLA);

    % one-step axial demand and robust margin
    r  = p_LA_i - p_k - ts_dt*v_k;
    ur = 2*norm(r)/ts_dt^2;
    di = ur - (a_max - sigma_off);

    delta_list(end+1) = di; 
    if di > 0, bad_idx(end+1) = i; end 
end

if isempty(bad_idx)
    fprintf('[OFFLINE] All samples satisfy |u_req| <= a_max - sigma (robust one-step feasible).\n');
else
    fprintf('[OFFLINE] %d samples exceed the robust bound; shrink look-ahead or cap local speed there.\n', numel(bad_idx));
end

%if overshoot <= 0
 %   fprintf('[OK]  Path is SAFE - max |a| = %.2f m/s^2 (limit = %.2f)\n', maxAbsA, a_max);
%else
 %   fprintf('[NOT OK] Path exceeds acceleration limit!\n');
  %  fprintf('Worst case |a| = %.2f m/s^2 (limit = %.2f)\n', maxAbsA, a_max);
   % fprintf('Overshoot = %.3f m/s^2 at t = %.3f s\n', overshoot, t_samp(violation_idx(1)));
%return; % Stop to skip online simulation if path is not safe
%end



%% Initial values
x = 0; % initial value, starting position 
y = 0;
vx = 0; % Initial value for velocity (zero velocity)
vy = 0;
t = (0:N-1)*ts_dt; %time vector

%% Plots
fig1 = figure;
hold on; grid on;
plot(p_ref(1,:),p_ref(2,:), 'k--','DisplayName','reference');
h_traj= animatedline('Color',[0.85 0 0],'LineWidth',1.4,'DisplayName','actual');
h_act = plot(x,y,'ks','MarkerSize',8,'LineWidth',1.5,'MarkerFaceColor','w','DisplayName','actual');
h_LA= plot(NaN,NaN,'ro','MarkerSize',8,'LineWidth',1.5,'DisplayName','look-ahead');
h_status = uicontrol(fig1,'Style','text','Units','normalized', 'Position',[0.45 0.85 0.15 0.08], 'FontSize',11,'FontWeight','bold','String','initialising','BackgroundColor',[0.9 0.9 0.9]);
legend('Location','best'); xlabel('x (m)'); ylabel('y (m)');
title('Pure pursuit with QP on random wp based spline - reach set')
h_Ctext = text(0.02, 0.98, '', 'Units','normalized', 'HorizontalAlignment','left','VerticalAlignment','top', 'Interpreter','none');


%% Simulation – random freeze object
interruptstartstep = randi([ceil(0.4*N) , ceil(0.6*N)]); % freeze begins 40–60 % through the run
interrupt_duration = 1;  % sec – how long the look-ahead is frozen
interrupt_steps_max= ceil(interrupt_duration/ts_dt);
freeze_end= interruptstartstep + interrupt_steps_max - 1; % last frozen step
frozen_pLA = nan(2,1);     % will hold the first frozen target


%% Initial arrays for logging
x_hist=zeros(1,N);  
y_hist=zeros(1,N);
vx_hist=zeros(1,N);
vy_hist=zeros(1,N);
ax_hist=zeros(1,N); 
ay_hist=zeros(1,N);
pLA_hist=zeros(2,N);
dpLA_hist=zeros(2,N);
delta_hist = zeros(1,N); 
inside_hist = false(1,N);
tol_req_hist= zeros(1,N);
sigma_allow_hist = zeros(1,N);   
eps_p_allow_hist = zeros(1,N);
eps_v_allow_hist = zeros(1,N);
% initialize for KKT-aligned C update 
ax = 0; ay = 0;  % needed for first-step a_star usage elsewhere
C_prev = 1e-3; %as initialize C  
C_hist = zeros(1,N);
C_KKT_hist = zeros(1,N);


%% the main loop
nextFrameT = 0;
reached_stop = true; stop_target = [NaN;NaN];
C_used_disp = 0; % default; inside reach-set the optimizer ignores C

for k = 1:N

    % Look ahead target - closest reference point
    d2 = (p_ref(1,:)-x).^2 + (p_ref(2,:)-y).^2; % calculating squared distance from the current positions
    [minDist, indexClosest] = min(d2);
    idxClosest = indexClosest;

    v_LA = v_ref(idxClosest);
    s_LA = max(look_ahead_s, v_LA*ts_dt);       % variable look-ahead distance = v_ref.dt - at least 1 mm
    idxLA = idxClosest;  s_acc = 0;
    while s_acc < s_LA && idxLA < numel(tau_grid)
            s_acc = s_acc + norm(p_ref(:,idxLA+1)-p_ref(:,idxLA));
            idxLA = idxLA + 1;
    end
    p_LA  = p_ref(:,idxLA);
    dp_LA = dp_ref(:,idxLA);
    
    freezeActive = (k >= interruptstartstep) && (k <= freeze_end);

    if freezeActive
        if k==interruptstartstep
            stop_target = [x; y];   % freeze at current pose
        end
        vx = 0; vy = 0; ax = 0; ay = 0;
        % no position integration; keep x,y fixed
        p_LA = stop_target; dp_LA = [0;0];
        reached_stop = true;
        sigma_allow = NaN; 
        eps_p_allow = NaN; 
        eps_v_allow = NaN; 
        delta = 0; 
    else


        % r-axis projection: r = p_LA - p - t_s v,  \hat e = r/||r||
        r_vec = p_LA - [x; y] - ts_dt*[vx; vy];
        r_norm = norm(r_vec);
        if r_norm < 1e-12
            e_axis = [1;0];   % arbitrary when r~~0
        else
            e_axis = r_vec / r_norm;
        end
        
        % scalar projections along \hat e
        p_rel  = dot([x; y],    e_axis);
        p_star = dot(p_LA,      e_axis);
        v_rel  = dot([vx; vy],  e_axis);
        
        % required axial accel (no ad-hoc noise injection)
        u_req = 2*(p_star - p_rel - v_rel*ts_dt) / ts_dt^2;

        pos_noise=0;%rand(1,1); % +- m arbitarily position noise bound
        vel_noise=0;%rand(1,1); % +- m/s vel noise bound
        eps_p= pos_noise*(2*rand - 1);  % noise in[-n,+n]
        eps_v= vel_noise*(2*rand - 1);
        p_rel_noisy= p_rel+ sign(p_star-p_rel)*eps_p;
        v_rel_noisy= v_rel+ sign(p_star-p_rel)*eps_v;        
        p_reach_low= p_rel_noisy+v_rel_noisy*ts_dt - 0.5*a_max*ts_dt^2;
        p_reach_high= p_rel_noisy+v_rel_noisy*ts_dt + 0.5*a_max*ts_dt^2;
        tol_req= max([p_star-p_reach_high , p_reach_low-p_star , 0]);
        e_tol_min= max(e_tol_min, tol_req);
    
        % robust one-step buffer:  sigma = 2*eps_p/t_s + eps_v
        sigma = (2*eps_p_bound)/ts_dt + eps_v_bound;
        
        % reachability margin (sufficient condition): delta <= 0 means reachable
        delta  = abs(u_req) - (a_max - sigma);
        inside = (delta <= 0);
        
        %how much more noise the system could tolerate in position or velocity and still be feasible?
        sigma_allow = max(a_max - abs(u_req), 0) ;    % reachability margin
        eps_p_allow = 0.5 * sigma_allow * ts_dt;  
        eps_v_allow = sigma_allow;
    
    
        e_p_vec = p_LA - ([x; y] + [vx; vy]*ts_dt); % this equals r_vec
        e_v_vec = dp_LA - [vx; vy]; % this equals d_v
        ep = norm(e_p_vec);
        ev = norm(e_v_vec);
        % QP (KKT - no tube constraint)
        C_used = C_prev; % use previous step's C to form the QP this step
        if inside
            % QP cost matrix J= ||p_LA-p||^2
            H = 2*((0.5*ts_dt^2)^2)*eye(2);
            d_f = -2*[(p_LA(1)-(x+vx*ts_dt))*0.5*ts_dt^2;
                      (p_LA(2)-(y+vy*ts_dt))*0.5*ts_dt^2];
        else
    

            % QP cost matrix J= ||p_LA-p||^2 + C.||v_LA - v||^2
            H   = 2 * ( ((0.5*ts_dt^2)^2) + C_used*(ts_dt^2) ) * eye(2);
            d_f = -2 * [ e_p_vec(1)*0.5*ts_dt^2 + C_used*e_v_vec(1)*ts_dt ;
                         e_p_vec(2)*0.5*ts_dt^2 + C_used*e_v_vec(2)*ts_dt ];
            C_used_disp = C_used;
        end 


        A = [ 1 0; -1 0; 0 1; 0 -1;  % accel box
            ts_dt 0; -ts_dt 0; 0 ts_dt; 0 -ts_dt]; % speed box
        b = [ a_max; a_max; a_max; a_max;
            v_max-vx; v_max+vx; v_max-vy; v_max+vy];
        C_used_disp = C_used;
    
        axay = quadprog(H,d_f,A,b,[],[],[],[],[],optimoptions("quadprog","Display","off"));
        
        if isempty(axay), axay=[0;0]; end
        ax = axay(1); 
        ay = axay(2);
    
        u_k = [ax; ay];                       % actual QP solution this step
        e_p_res = e_p_vec - 0.5*ts_dt^2*u_k;  % = r - 0.5 t_s^2 u_k
        e_v_res = e_v_vec - ts_dt*u_k;        % = d_v - t_s u_k
        
        eps_safe = 1e-9;                      % small denom safeguard
        C_raw = (ts_dt * norm(e_p_res)) / (2*norm(e_v_res) + eps_safe);
    
        
        % Bias rule: when outside the reach set (delta>0), prioritize position -> smaller C
        rho_out = 0.4; % 0<rho_out<=1 (smaller => more position priority when outside: stronger snap back)
        rho_in = 1.0; % inside: no bias
        bias  = (delta>0)*rho_out + (delta<=0)*rho_in;
        C_target = bias * C_raw;
        
        % Smooth + clamp for robustness as C for nex step
        beta = 0.65;  % follow KKT/bias how fast each step              
        Cmin = 1e-6; Cmax = 1e3;
        C_next = min(max((1-beta)*C_prev + beta*C_target, Cmin), Cmax);
        
        C_hist(k) = C_used_disp; % C actually used in THIS step's QP
        C_KKT_hist(k) = C_raw;  % raw KKT value (pre-bias), for diagnostics
        C_prev = C_next; % will be used to build H,d_f on the next step
    
        vx_old = vx;  vy_old = vy;
        vx = vx_old + ax*ts_dt;
        vy = vy_old + ay*ts_dt;
        x  = x + vx_old*ts_dt + 0.5*ax*ts_dt^2;
        y  = y + vy_old*ts_dt + 0.5*ay*ts_dt^2;


  %% live plot
        if t(k) >= nextFrameT
            addpoints(h_traj,x,y);
            set(h_act,'XData',x,'YData',y);
            set(h_LA,'XData',p_LA(1),'YData',p_LA(2));
            r= 0.5*a_max*ts_dt^2; % half-length on axis e
            reach_end1= p_LA-r*e_axis; % world-coords end-points
            reach_end2= p_LA+r*e_axis;
            set(h_Ctext, 'String', sprintf('C = %.3g', C_used_disp));
    
        
            if delta <= 0
                set(h_status,'String',' INSIDE ','BackgroundColor',[0 1 0]); % green
            else
                set(h_status,'String',' OUTSIDE ','BackgroundColor',[1 0 0]); % red
            end
                drawnow
                nextFrameT = nextFrameT + fr_dt;
        end

        eps_p_max = min(eps_p_allow_hist, [], 'omitnan');
        eps_v_max = min(eps_v_allow_hist, [], 'omitnan');
        max_eps_p_allow = max(eps_p_allow_hist, [], 'omitnan');
        max_eps_v_allow = max(eps_v_allow_hist, [], 'omitnan');
    end



    %% Log/history,
    x_hist(k)= x;
    y_hist(k)= y;
    vx_hist(k)= vx;
    vy_hist(k) = vy;
    ax_hist(k)= ax;
    ay_hist(k) = ay;
    speed_hist(k)= hypot(vx, vy);
    pLA_hist(:,k)=p_LA;  
    dpLA_hist(:,k)=dp_LA;
    tol_req_hist(k) = tol_req;;
    delta_hist(k)= delta;
    sigma_allow_hist(k) = sigma_allow;
    eps_p_allow_hist(k) = eps_p_allow;
    eps_v_allow_hist(k) = eps_v_allow;


end

eps_p_max = min(eps_p_allow_hist); %the maximum noise we can guarantee to tolerate at all steps
eps_v_max = min(eps_v_allow_hist);

max_eps_p_allow = max(eps_p_allow_hist);  
max_eps_v_allow = max(eps_v_allow_hist);

freeze = false(1,N); %masking freezing time
freeze(interruptstartstep:freeze_end) = true;
moving = ~freeze; 

t_moving = ts_dt * (0:sum(moving)-1);  % no freeze gap in time
x_moving = x_hist(moving);
y_moving = y_hist(moving);
vx_moving = vx_hist(moving);
vy_moving = vy_hist(moving);
v_moving = speed_hist(moving);

T_norm = (0:sum(moving)-1) ./ max(sum(moving)-1, 1);

%% Summary plots
figure;

% position
subplot(5,2,1); 
plot(t,x_hist,'b',t,y_hist,'r'); grid on
ylabel('pos (m)'); legend('x','y'); title('Position')

% speed & mag.
subplot(5,2,2);
%speed_hist = hypot(vx_hist,vy_hist);
plot(t,vx_hist,'b',t,vy_hist,'r', t, speed_hist,'k--'); grid on;
ylabel('speed (m/s)'); legend('vx','vy','v'); title('Speed');

% position error (distance to LA)
subplot(5,2,3)
pos_err = hypot(x_hist - pLA_hist(1,:), y_hist - pLA_hist(2,:));
plot(t,pos_err), grid on
ylabel('error (m)'); title('Position error (to LA)')

% Velocity error (difference to look-ahead velocity)
subplot(5,2,4)
vel_err = hypot(vx_hist - dpLA_hist(1,:), vy_hist - dpLA_hist(2,:));
plot(t,vel_err), grid on
ylabel('error (m/s)'); title('Velocity error (to LA)')

% magnitude of acc. 
subplot(5,2,5)
acc_hist = hypot(ax_hist, ay_hist);
plot(t,acc_hist), grid on
ylabel('m/s^2'); title('Acceleration')

subplot(5,2,6)
plot(t,delta_hist); hold on; grid on
plot(t,zeros(size(t)),'k:');      % zero line
ylabel('\delta = |u_{req}|-a_{max}');
title('Reachability margin  (below 0 = tube-safe)')

subplot(5,2,8)
plot(t, sigma_allow_hist); hold on; grid on
title('\sigma_{allow} per step'); ylabel('m/s^2');

subplot(5,2,7)
plot(t, tol_req_hist); hold on; grid on
title('e_{tol} history');

subplot(5,2,9)
plot(t, eps_p_allow_hist, 'b'); hold on;
yline(max_eps_p_allow, 'r--','DisplayName','max allowable');
xlabel('Time (s)'); ylabel('e_p^{max} (m)');
title('Allowable Position Noise Margin');
legend('actual','max'); grid on;

subplot(5,2,10)
plot(t, eps_v_allow_hist, 'r'); hold on;
yline(max_eps_v_allow, 'k--', 'DisplayName','max allowable');
xlabel('Time (s)'); ylabel('e_v^{max} (m/s)');
title('Allowable Velocity Noise Margin');
legend('actual','max'); grid on;


figure;
plot(t, C_hist,'b', t, C_KKT_hist, 'r'); grid on;
xlabel('Time (s)'); ylabel('C (gain)');
title('Adaptive gain C'); legend;

figure;
% position
subplot(2,1,1); grid on;
plot(t_moving, x_moving, 'b',t_moving, y_moving, 'r');
ylabel('pos (m)');
title('position vs time (freezing time removed)');
legend('x', 'y');

% velocity
subplot(2,1,2); grid on;
plot(t_moving, vx_moving, 'b', t_moving, vy_moving, 'r', t_moving, v_moving, 'k--');
xlabel('Time (s)'); ylabel('speed (m/s)'); legend('vx','vy','v'); title('Speed');


