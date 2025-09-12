function batch_qp_vs_pp()
% comparison: QP+reachability vs Pure Pursuit (PP)
% (r-axis, robust delta, KKT C-update)
clc; clear; close all;
%% 
Nruns = 1;     % number of random splines
Show_live  = true;  % live XY plot per run
T_grid = linspace(0,1,200);  % normalized moving-time grid for aggregation

%% aggregate holders (rows=runs, cols=normalized time)
ep_QP = nan(Nruns, numel(T_grid));  % ||p - p_LA|| (moving only, then resampled)
ev_QP = nan(Nruns, numel(T_grid));  % ||v - v_LA||
del_QP = nan(Nruns, numel(T_grid));  % delta(t)

ep_PP = nan(Nruns, numel(T_grid));
ev_PP = nan(Nruns, numel(T_grid));
del_PP = nan(Nruns, numel(T_grid));

mean_delta_QP = nan(Nruns,1);
mean_delta_PP = nan(Nruns,1);

fprintf('Running %d random-spline simulations (QP vs PP)...\n', Nruns);
for r = 1:Nruns
    o = simulate_one_run(Show_live);  % one path, both controllers

    moving = o.moving;
    if nnz(moving) < 5
        warning('Run %d has too few moving samples; skipping.', r);
        continue;
    end

    % magnitudes during moving time
    ep_m_QP = hypot(o.x_QP(moving) - o.pLA_QP(1,moving),  o.y_QP(moving) - o.pLA_QP(2,moving));
    ev_m_QP = hypot(o.vx_QP(moving) - o.dpLA_QP(1,moving), o.vy_QP(moving) - o.dpLA_QP(2,moving));
    del_m_QP = o.delta_QP(moving);

    ep_m_PP = hypot(o.x_PP(moving) - o.pLA_PP(1,moving),  o.y_PP(moving) - o.pLA_PP(2,moving));
    ev_m_PP = hypot(o.vx_PP(moving) - o.dpLA_PP(1,moving), o.vy_PP(moving) - o.dpLA_PP(2,moving));
    del_m_PP = o.delta_PP(moving);

    mean_delta_QP(r) = mean(del_m_QP, 'omitnan');
    mean_delta_PP(r) = mean(del_m_PP, 'omitnan');

    % normalize moving time to [0,1] and resample to T_grid
    t_m = o.t_moving;
    t_norm = (t_m - t_m(1)) / max(t_m(end) - t_m(1), eps);

    ep_QP(r,:) = interp1(t_norm, ep_m_QP,  T_grid, 'linear','extrap');
    ev_QP(r,:) = interp1(t_norm, ev_m_QP,  T_grid, 'linear','extrap');
    del_QP(r,:) = interp1(t_norm, del_m_QP, T_grid, 'linear','extrap');

    ep_PP(r,:) = interp1(t_norm, ep_m_PP,  T_grid, 'linear','extrap');
    ev_PP(r,:) = interp1(t_norm, ev_m_PP,  T_grid, 'linear','extrap');
    del_PP(r,:) = interp1(t_norm, del_m_PP, T_grid, 'linear','extrap');
end

%% Fig 1: RMS ep/ev vs normalized moving time
cb = [0.0000 0.4470 0.7410]; % blue (QP)
co = [0.8500 0.3250 0.0980]; % orange (PP)

ep_RMS_QP = sqrt(mean(ep_QP.^2, 1, 'omitnan'));
ev_RMS_QP = sqrt(mean(ev_QP.^2, 1, 'omitnan'));
ep_RMS_PP = sqrt(mean(ep_PP.^2, 1, 'omitnan'));
ev_RMS_PP = sqrt(mean(ev_PP.^2, 1, 'omitnan'));

fig1 = figure('Name','RMS Position and Velocity Errors (QP vs PP)', ...
              'Units','inches','Position',[0 0 5 3],'Color','w');
tiledlayout(2,1,'Padding','compact','TileSpacing','compact');

% top: position RMS
nexttile; hold on;
plot(T_grid, ep_RMS_QP, 'LineWidth',1.8, 'Color', cb);
plot(T_grid, ep_RMS_PP, 'LineWidth',1.8, 'Color', co);
grid on; ylabel('RMS Position Error (m)');
set(gca,'FontName','Times','FontSize',8,'LineWidth',0.8);
legend({'QP+reachability','Pure pursuit'},'Location','best');

% bottom: velocity RMS
nexttile; hold on;
plot(T_grid, ev_RMS_QP, 'LineWidth',1.8, 'Color', cb);
plot(T_grid, ev_RMS_PP, 'LineWidth',1.8, 'Color', co);
grid on; xlabel('Normalized moving time'); ylabel('RMS Velocity Error (m/s)');
set(gca,'FontName','Times','FontSize',8,'LineWidth',0.8);
legend({'QP+reachability','Pure pursuit'},'Location','best');

%% Fig 2: Histogram of per-run mean \delta (QP vs PP) — dynamic y-limits

XMIN = -200;
XMAX = 2000;
NBINS = 80;
edges = linspace(XMIN, XMAX, NBINS+1);
clip = @(x) min(max(x, XMIN), XMAX);
mdQ = clip(mean_delta_QP);
mdP = clip(mean_delta_PP);

fig2 = figure('Units','inches','Position',[0 0 5 3], 'Color','w','Renderer','painters');
tiledlayout(2,1,'Padding','compact','TileSpacing','compact');

%Top: QP
nexttile; hold on; box on;
hQ = histogram(mdQ,'BinEdges',edges,'FaceColor',cb,'EdgeColor','k','LineWidth',0.75);
xline(0,'k--','\delta=0', 'LabelHorizontalAlignment','left', 'LabelVerticalAlignment','bottom');
grid on; xlim([XMIN XMAX]); ylabel('Count'); title('QP + reachability');
% dynamic y-limit (with fallback if Values not available)
if isprop(hQ,'Values'), yMaxQ = max(hQ.Values); else, yMaxQ = max(histcounts(mdQ,edges)); end
yPadQ = max(1, ceil(0.06*yMaxQ));           % ~6% headroom, at least 1
ylim([0, yMaxQ + yPadQ]);

set(gca,'FontName','Times','FontSize',8,'LineWidth',0.8);

%Bottom: PP
nexttile; hold on; box on;
hP = histogram(mdP,'BinEdges',edges,'FaceColor',co,'EdgeColor','k','LineWidth',0.75);
xline(0,'k--','\delta=0', 'LabelHorizontalAlignment','left', 'LabelVerticalAlignment','bottom');
grid on; xlim([XMIN XMAX]); xlabel('Mean \delta (m/s^2)'); ylabel('Count'); title('Pure pursuit');
% dynamic y-limit (independent of the top)
if isprop(hP,'Values'), yMaxP = max(hP.Values); else, yMaxP = max(histcounts(mdP,edges)); end
yPadP = max(1, ceil(0.06*yMaxP));
ylim([0, yMaxP + yPadP]);

set(gca,'FontName','Times','FontSize',8,'LineWidth',0.8);


%% Fig 3: mean \delta(t) over runs (QP vs PP)
fig3 = figure('Units','inches','Position',[0 0 5 3],'Color','w');
tiledlayout(2,1,'Padding','compact','TileSpacing','compact');

nexttile; hold on;
plot(T_grid, mean(del_QP,1,'omitnan'), 'LineWidth',1.8, 'Color', cb);
yline(0,'k--'); grid on; ylabel('\delta (m/s^2)');
title('QP + reachability'); set(gca,'FontName','Times','FontSize',8,'LineWidth',0.8);

nexttile; hold on;
plot(T_grid, mean(del_PP,1,'omitnan'), 'LineWidth',1.8, 'Color', co);
yline(0,'k--'); grid on; xlabel('Normalized moving time'); ylabel('\delta (m/s^2)');
title('Pure pursuit'); set(gca,'FontName','Times','FontSize',8,'LineWidth',0.8);

%% Fig 4: Example run (QP) — position & speed (freeze removed) 
o_ex = simulate_one_run(false);  % no live plot
moving = o_ex.moving;
t_moving = o_ex.t_moving;
x_moving = o_ex.x_QP(moving);
y_moving = o_ex.y_QP(moving);
vx_moving= o_ex.vx_QP(moving);
vy_moving = o_ex.vy_QP(moving);
v_moving = hypot(vx_moving,vy_moving);
x_plan = o_ex.pLA_QP(1,moving);
y_plan = o_ex.pLA_QP(2,moving);

% freeze seam marker (time at which the 1s freeze was cut out)
freeze_start_idx = find(~o_ex.moving, 1, 'first');
if isempty(freeze_start_idx)
    t_cut = NaN;
else
    pre_m = sum(moving(1:freeze_start_idx-1));
    if pre_m > 0
        t_cut = t_moving(pre_m);
    else
        t_cut = 0;
    end

end

fig4 = figure('Units','inches','Position',[0 0 5 3],'Color','w');
tiledlayout(2,1,'Padding','compact','TileSpacing','compact');

% position
nexttile; hold on; grid on; box on;
plot(t_moving, x_moving, 'b-', 'LineWidth',1, 'DisplayName','x');
plot(t_moving, y_moving, 'r-', 'LineWidth',1, 'DisplayName','y');
plot(t_moving, x_plan,   'k--','LineWidth',0.8,'DisplayName','planned');
plot(t_moving, y_plan,   'k--','LineWidth',0.8,'HandleVisibility','off');
if ~isnan(t_cut), xline(t_cut,':','Color',[0.35 0.35 0.35], 'LineWidth',1.0, 'DisplayName','freeze'); end
ylabel('pos (m)'); legend('Location','best','NumColumns',2,'FontSize',8);
set(gca,'FontName','Times','FontSize',8,'LineWidth',0.8);

% velocity
nexttile; hold on; grid on; box on;
plot(t_moving, vx_moving, 'b-', 'LineWidth',1, 'DisplayName','v_x');
plot(t_moving, vy_moving, 'r-', 'LineWidth',1, 'DisplayName','v_y');
plot(t_moving, v_moving,  'k-', 'LineWidth',1, 'DisplayName','|v|');
if ~isnan(t_cut), xline(t_cut,':','Color',[0.35 0.35 0.35], 'LineWidth',1.0, 'DisplayName','freeze'); end
xlabel('Time (s)'); ylabel('speed (m/s)');
legend('Location','best','NumColumns',2,'FontSize',8);
set(gca,'FontName','Times','FontSize',8,'LineWidth',0.8);

end 



function o = simulate_one_run(Show_live)
% One random path; QP & PP track it with a single 1s freeze - Returns logs for aggregation.

%% Parameters
ts_dt= 0.01; fr_dt = 0.04; T_time= 12;
look_ahead_s = 0.001;
v_max= 1.0; a_max= 2.5;
N = ceil(T_time/ts_dt);
eps_p_bound = 0.0; eps_v_bound = 0.0;     % set >0 to use robust buffer
sigma = (2*eps_p_bound)/ts_dt + eps_v_bound;

% reference path: random waypoint spline
n_wp=6; L_x=0.5; L_y=0.5;
wp = zeros(2,n_wp+1);
theta = 2*pi*rand(1,n_wp);
wp(1,2:end) = L_x*cos(theta);
wp(2,2:end) = L_y*sin(theta);
dLine = [0, vecnorm(diff(wp,1,2))];
sLine = cumsum(dLine);
t_wp = sLine/sLine(end)*4*pi;

ppx  = spline(t_wp, wp(1,:));
ppy  = spline(t_wp, wp(2,:));
dppX = fnder(ppx,1);
dppY = fnder(ppy,1);

tau_grid = linspace(0, ppx.breaks(end), 150000);
p_ref  = [ ppval(ppx , tau_grid) ; ppval(ppy , tau_grid) ];
dp_ref = [ ppval(dppX, tau_grid) ; ppval(dppY, tau_grid) ];
v_ref  = vecnorm(dp_ref);

% initial states
x_QP=0; y_QP=0; vx_QP=0; vy_QP=0;
x_PP=0; y_PP=0; vx_PP=0; vy_PP=0;
t = (0:N-1)*ts_dt;

% live plot option
if Show_live
    cb = [0.0000 0.4470 0.7410]; % blue (QP)
    co = [0.8500 0.3250 0.0980]; % orange (PP)
    fig = figure('Units','inches','Position',[0 0 5 3],'Color','w'); hold on; grid on; box on; axis equal;
    plot(p_ref(1,:),p_ref(2,:), 'k--','DisplayName','reference');
    h_traj_QP = animatedline('Color',cb, 'LineWidth',1.8,'DisplayName','QP');
    h_traj_PP = animatedline('Color',co, 'LineWidth',1.8,'DisplayName','PP');
    h_act_QP  = plot(x_QP,y_QP,'s','MarkerSize',8,'LineWidth',1.2,'MarkerEdgeColor',cb,'MarkerFaceColor','w');
    h_act_PP  = plot(x_PP,y_PP,'o','MarkerSize',8,'LineWidth',1.2,'MarkerEdgeColor',co,'MarkerFaceColor','w');
    h_LA_QP   = plot(NaN,NaN,'o','MarkerSize',6,'LineWidth',1.2,'Color',cb);
    h_LA_PP   = plot(NaN,NaN,'^','MarkerSize',6,'LineWidth',1.2,'Color',co);
    legend('Location','best'); xlabel('x (m)'); ylabel('y (m)');
end

% freeze schedule
interruptstartstep = randi([ceil(0.4*N), ceil(0.6*N)]);
interrupt_duration = 1.0;  
freeze_end = interruptstartstep + ceil(interrupt_duration/ts_dt) - 1;
stop_target_QP = [NaN;NaN]; reached_stop_QP = true;
stop_target_PP = [NaN;NaN]; reached_stop_PP = true;

% logs
x_QP_hist=zeros(1,N); y_QP_hist=zeros(1,N);
vx_QP_hist=zeros(1,N); vy_QP_hist=zeros(1,N);
pLA_QP_hist=zeros(2,N); dpLA_QP_hist=zeros(2,N); delta_QP_hist=zeros(1,N);

x_PP_hist=zeros(1,N); y_PP_hist=zeros(1,N);
vx_PP_hist=zeros(1,N); vy_PP_hist=zeros(1,N);
pLA_PP_hist=zeros(2,N); dpLA_PP_hist=zeros(2,N); delta_PP_hist=zeros(1,N);

% C-gain
C_prev = 1e-3;  Cmin=1e-6; Cmax=1e3; rho_out=0.4; rho_in=1.0; beta=0.65;

nextFrameT = 0;

for k = 1:N
    %look-ahead for QP
    d2 = (p_ref(1,:)-x_QP).^2 + (p_ref(2,:)-y_QP).^2;
    [~, idxC] = min(d2);
    vLA  = v_ref(idxC);
    sLA  = max(look_ahead_s, vLA*ts_dt);
    idxLA = idxC; s_acc=0;
    while s_acc < sLA && idxLA < numel(tau_grid)
        s_acc = s_acc + norm(p_ref(:,idxLA+1)-p_ref(:,idxLA));
        idxLA = idxLA + 1;
    end
    p_LA_QP  = p_ref(:,idxLA);
    dp_LA_QP = dp_ref(:,idxLA);

    %look-ahead for PP (same policy)
    d2p = (p_ref(1,:)-x_PP).^2 + (p_ref(2,:)-y_PP).^2;
    [~, idxCp] = min(d2p);
    vLAp  = v_ref(idxCp);
    sLAp  = max(look_ahead_s, vLAp*ts_dt);
    idxLAp= idxCp; s_acc=0;
    while s_acc < sLAp && idxLAp < numel(tau_grid)
        s_acc = s_acc + norm(p_ref(:,idxLAp+1)-p_ref(:,idxLAp));
        idxLAp = idxLAp + 1;
    end
    p_LA_PP  = p_ref(:,idxLAp);
    dp_LA_PP = dp_ref(:,idxLAp);

    freezeActive = (k >= interruptstartstep) && (k <= freeze_end);

    %% QP dynamics
    if freezeActive
        if k==interruptstartstep
            v0 = [vx_QP;vy_QP]; v0n=norm(v0);
            if v0n<1e-5, stop_target_QP=[x_QP;y_QP];
            else, d_stop=v0n^2/(2*a_max); stop_target_QP=[x_QP;y_QP]+v0/v0n*d_stop; end
            reached_stop_QP=false;
        end
        p_LA_QP = stop_target_QP; dp_LA_QP = [0;0];
        if ~reached_stop_QP
            vnorm=norm([vx_QP;vy_QP]);
            if vnorm<1e-4, reached_stop_QP=true; ax_QP=0; ay_QP=0; vx_QP=0; vy_QP=0;
            else
                ax_QP=-a_max*vx_QP/max(vnorm,1e-9); ay_QP=-a_max*vy_QP/max(vnorm,1e-9);
                vx_QP=vx_QP+ax_QP*ts_dt; vy_QP=vy_QP+ay_QP*ts_dt;
                x_QP=x_QP+vx_QP*ts_dt+0.5*ax_QP*ts_dt^2;
                y_QP=y_QP+vy_QP*ts_dt+0.5*ay_QP*ts_dt^2;
            end
        else
            ax_QP=0; ay_QP=0;
        end
        delta_QP = 0;
    else
        % r-axis projection & u_req
        r_vec = p_LA_QP - [x_QP;y_QP] - ts_dt*[vx_QP;vy_QP];
        r_n = norm(r_vec); e_axis = (r_n<1e-12) * [1;0] + (r_n>=1e-12) * (r_vec/r_n);
        p_rel = dot([x_QP;y_QP], e_axis);
        p_star= dot(p_LA_QP,     e_axis);
        v_rel = dot([vx_QP;vy_QP],e_axis);
        u_req = 2*(p_star - p_rel - v_rel*ts_dt)/ts_dt^2;

        % reachability margin (sigma=0 here unless you set bounds above)
        delta_QP = abs(u_req) - (a_max - sigma);
        inside   = (delta_QP <= 0);

        % cost terms
        e_p_vec = p_LA_QP - ([x_QP;y_QP] + [vx_QP;vy_QP]*ts_dt);
        e_v_vec = dp_LA_QP - [vx_QP;vy_QP];

        C_used = C_prev;  % use last C to build this step's QP
        if inside
            H  = 2*((0.5*ts_dt^2)^2)*eye(2);
            df = -2*[(e_p_vec(1))*0.5*ts_dt^2; (e_p_vec(2))*0.5*ts_dt^2];
        else
            H  = 2*(((0.5*ts_dt^2)^2) + C_used*(ts_dt^2))*eye(2);
            df = -2*[ e_p_vec(1)*0.5*ts_dt^2 + C_used*e_v_vec(1)*ts_dt;
                      e_p_vec(2)*0.5*ts_dt^2 + C_used*e_v_vec(2)*ts_dt ];
        end

        % constraints: accel box + next-step speed box
        A = [ 1 0; -1 0; 0 1; 0 -1;   ts_dt 0; -ts_dt 0; 0 ts_dt; 0 -ts_dt];
        b = [ a_max; a_max; a_max; a_max;  v_max-vx_QP; v_max+vx_QP; v_max-vy_QP; v_max+vy_QP];

        axay = quadprog(H,df,A,b,[],[],[],[],[],optimoptions("quadprog","Display","off"));
        if isempty(axay), axay=[0;0]; end
        ax_QP=axay(1); ay_QP=axay(2);

        % integrate
        vx_old=vx_QP; vy_old=vy_QP;
        vx_QP=vx_old+ax_QP*ts_dt; vy_QP=vy_old+ay_QP*ts_dt;
        x_QP = x_QP + vx_old*ts_dt + 0.5*ax_QP*ts_dt^2;
        y_QP = y_QP + vy_old*ts_dt + 0.5*ay_QP*ts_dt^2;

        % KKT-inspired C update (for next step)
        u_k   = [ax_QP;ay_QP];
        epres = e_p_vec - 0.5*ts_dt^2*u_k;
        evres = e_v_vec - ts_dt*u_k;
        eps_safe = 1e-9;
        C_raw = (ts_dt * norm(epres)) / (2*norm(evres) + eps_safe);

        bias = (delta_QP>0)*rho_out + (delta_QP<=0)*rho_in;
        C_prev = min(max((1-beta)*C_prev + beta*(bias*C_raw), Cmin), Cmax);
    end

    %% PP dynamics
    if freezeActive
        if k==interruptstartstep
            v0=[vx_PP;vy_PP]; v0n=norm(v0);
            if v0n<1e-5, stop_target_PP=[x_PP;y_PP];
            else, d_stop=v0n^2/(2*a_max); stop_target_PP=[x_PP;y_PP]+v0/v0n*d_stop; end
            reached_stop_PP=false;
        end
        p_LA_PP = stop_target_PP; dp_LA_PP = [0;0];
        if ~reached_stop_PP
            vnorm=norm([vx_PP;vy_PP]);
            if vnorm<1e-4, reached_stop_PP=true; ax_PP=0; ay_PP=0; vx_PP=0; vy_PP=0;
            else
                ax_PP=-a_max*vx_PP/max(vnorm,1e-9); ay_PP=-a_max*vy_PP/max(vnorm,1e-9);
                vx_PP=vx_PP+ax_PP*ts_dt; vy_PP=vy_PP+ay_PP*ts_dt;
                x_PP=x_PP+vx_PP*ts_dt + 0.5*ax_PP*ts_dt^2;
                y_PP=y_PP+vy_PP*ts_dt + 0.5*ay_PP*ts_dt^2;
            end
        else
            ax_PP=0; ay_PP=0;
        end
        delta_PP = 0;
    else
        dir_PP = p_LA_PP - [x_PP;y_PP]; nd=norm(dir_PP);
        if nd==0, e_axis_PP=[1;0]; else, e_axis_PP=dir_PP/nd; end
        p_rel_PP = dot([x_PP;y_PP], e_axis_PP);
        p_star_PP= dot(p_LA_PP,     e_axis_PP);
        v_rel_PP = dot([vx_PP;vy_PP],e_axis_PP);
        u_req_PP = 2*(p_star_PP - p_rel_PP - v_rel_PP*ts_dt)/ts_dt^2;

        a_des_PP = u_req_PP * e_axis_PP;
        a_mag = norm(a_des_PP);
        if a_mag > a_max, a_des_PP = a_des_PP * (a_max / a_mag); end
        ax_PP = a_des_PP(1); ay_PP = a_des_PP(2);

        % next-step speed clipping to box
        vx_next = max(min(vx_PP + ax_PP*ts_dt, v_max), -v_max);
        vy_next = max(min(vy_PP + ay_PP*ts_dt, v_max), -v_max);

        x_PP = x_PP + vx_PP*ts_dt + 0.5*ax_PP*ts_dt^2;
        y_PP = y_PP + vy_PP*ts_dt + 0.5*ay_PP*ts_dt^2;
        vx_PP = vx_next; vy_PP = vy_next;

        delta_PP = abs(u_req_PP) - (a_max - sigma);
    end

    % live plot
    if Show_live && t(k) >= nextFrameT
        addpoints(h_traj_QP, x_QP, y_QP); set(h_act_QP,'XData',x_QP,'YData',y_QP);
        addpoints(h_traj_PP, x_PP, y_PP); set(h_act_PP,'XData',x_PP,'YData',y_PP);
        set(h_LA_QP,'XData',p_LA_QP(1),'YData',p_LA_QP(2));
        set(h_LA_PP,'XData',p_LA_PP(1),'YData',p_LA_PP(2));
        drawnow;
        nextFrameT = nextFrameT + fr_dt;
    end

    % logs
    x_QP_hist(k)=x_QP; y_QP_hist(k)=y_QP; vx_QP_hist(k)=vx_QP; vy_QP_hist(k)=vy_QP;
    pLA_QP_hist(:,k)=p_LA_QP; dpLA_QP_hist(:,k)=dp_LA_QP; delta_QP_hist(k)=delta_QP;

    x_PP_hist(k)=x_PP; y_PP_hist(k)=y_PP; vx_PP_hist(k)=vx_PP; vy_PP_hist(k)=vy_PP;
    pLA_PP_hist(:,k)=p_LA_PP; dpLA_PP_hist(:,k)=dp_LA_PP; delta_PP_hist(k)=delta_PP;
end

% moving mask & gapless time (freeze removed)
freeze = false(1,N); freeze(interruptstartstep:freeze_end) = true;
moving = ~freeze; 
t_moving = ts_dt * (0:sum(moving)-1);

% outputs
o.moving   = moving;
o.t_moving = t_moving;

o.x_QP = x_QP_hist; o.y_QP = y_QP_hist; o.vx_QP=vx_QP_hist; o.vy_QP=vy_QP_hist;
o.pLA_QP = pLA_QP_hist; o.dpLA_QP = dpLA_QP_hist; o.delta_QP = delta_QP_hist;

o.x_PP = x_PP_hist; o.y_PP = y_PP_hist; o.vx_PP=vx_PP_hist; o.vy_PP=vy_PP_hist;
o.pLA_PP = pLA_PP_hist; o.dpLA_PP = dpLA_PP_hist; o.delta_PP = delta_PP_hist;
end

function [bw, L] = nice_bins(all_vals, shrink)
% Freedman-Diaconis bin width - Returns bin width bw and a symmetric half-range L so you can do: edges = -L:bw:L

    if nargin < 2 || isempty(shrink), shrink = 1; end   % 3 makes 3× more bins
    all_vals = all_vals(isfinite(all_vals));
    if isempty(all_vals)
        bw = 0.1; L = 1.0; 
        return;
    end

    % FD suggested width
    bw_fd = 2*iqr(all_vals) / max(numel(all_vals)^(1/3), eps);
    if ~isfinite(bw_fd) || bw_fd <= 0
        bw_fd = range(all_vals)/12;   % fallback
    end

    tier  = 10.^floor(log10(max(bw_fd, realmin)));
    cands = [1 2 5].*tier;
    [~,ix] = min(abs(cands - bw_fd));
    bw = cands(ix);

    bw = bw / shrink;

    % Symmetric extent around zero (so a bar always sits at 0)
    L = ceil(max(abs(all_vals))/bw)*bw;
end

