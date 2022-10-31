

clear all

constants.g = 10;

constants.m_1 = 1;
constants.l_1 = 2;
constants.l_MCD_1 = 1;
constants.I_1 = 1/3 .* constants.m_1 .* constants.l_MCD_1.^2;

constants.m_2 = 1;
constants.l_2 = 2;
constants.l_MCD_2 = 1;
constants.I_2 = 1/3 .* constants.m_2 .* constants.l_MCD_2.^2;

g = constants.g;

m_1 = constants.m_1;
l_1 = constants.l_1;
l_MCD_1 = constants.l_MCD_1;
I_1 = constants.I_1;

m_2 = constants.m_2;
l_2 = constants.l_2;
l_MCD_2 = constants.l_MCD_2;
I_2 = constants.I_2;

func_tau_1 = @(t,q,varargin) 0 + zeros( size( t ) );
func_tau_2 = @(t,q,varargin) 1 + zeros( size( t ) );

func_tau = @(t,q,varargin) [ func_tau_1(t,q,varargin{:}); func_tau_2(t,q,varargin{:}); ];

ode_fun = @(t,q) ddt_ode45(t,q, constants, func_tau );
event_fun = @(t,q) event_one_swing(t, q, constants, func_tau );

ode_option = odeset('Events',event_fun);


t = 0:1e-3:10;
t = t(:);

th_1_0 = 1/2*pi + deg2rad(90);
omega_1_0 = 0;

th_2_0 = deg2rad( 0 );
omega_2_0 = 0;

q0 = [th_1_0, omega_1_0, th_2_0, omega_2_0]';

% [time, q, te, qe, ie] = ode45(ode_fun, t, q0, ode_option);
[time, q] = ode45(ode_fun, t, q0);

trimming_index = ~isnan(time);
time = time( trimming_index );
q = q( trimming_index, : );

th_1 = q(:,1);
omega_1 = q(:,2);

th_2 = q(:,3);
omega_2 = q(:,4);

%%
p_CoM_whole = find_p_CoM_whole(l_1,l_MCD_1,l_MCD_2,m_1,m_2,th_1,th_2);
zeros_array = zeros( size( p_CoM_whole(:,1) ) );
NaN_array = NaN( size( p_CoM_whole(:,1) ) );
p_top = [0,0] + zeros_array;
p_bottom_1 = p_top + l_1 * [cos( th_1 ), sin( th_1 )];
p_bottom_2 = p_bottom_1 + l_2 * [cos( th_1 + th_2 ), sin( th_1 + th_2 )];


%% anime

xArray = [p_top(:,1), p_bottom_1(:,1), p_bottom_2(:,1)];
yArray = [p_top(:,2), p_bottom_1(:,2), p_bottom_2(:,2)];
zArray = zeros( size( xArray ) );

anime_fig = SimplestAnime_exported(time, xArray, yArray, zArray);
ax_range = [min( [xArray, yArray, zArray], [], 'all' ), max( [xArray, yArray, zArray], [], 'all' )];
anime_fig.axAnime.XLim = ax_range;
anime_fig.axAnime.YLim = ax_range;
anime_fig.axAnime.ZLim = ax_range;
view(anime_fig.axAnime, [0,0,1])
grid(anime_fig.axAnime, 'on')













