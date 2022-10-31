clear all
tic

%% prepare variables

syms g real

syms m_1 real
syms I_1 real
syms l_1 l_MCD_1 real

syms th_1_pre(t)
syms th_1 omega_1 alpha_1 real
syms tau_1 real

syms m_2 real
syms I_2 real
syms l_2 l_MCD_2 real

syms th_2_pre(t)
syms th_2 omega_2 alpha_2 real
syms tau_2 real

syms_pre_tobereplaced = [
    th_1_pre, diff(th_1_pre, t), diff(th_1_pre, t, t), ...
    th_2_pre, diff(th_2_pre, t), diff(th_2_pre, t, t), ...
    ];
syms_pre_replacing = [
    th_1, omega_1, alpha_1, ...
    th_2, omega_2, alpha_2, ...
    ];
subs_pre = @(input) subs( input, syms_pre_tobereplaced, syms_pre_replacing );

syms x_1(t) y_1(t)
syms_restraint_tobereplaced = [
    x_1, y_1, ...
    ];
syms_restraint_replacing = [
    0, 0, ...
    ];
syms_restraint_replacing = sym( syms_restraint_replacing );
subs_restraing = @(input) subs( input, syms_restraint_tobereplaced, syms_restraint_replacing );

%% model definition
% define bottom point for each segment. 
p_top = formula( ...
    [x_1, y_1] ...
    );
p_bottom_1 = formula( ...
    p_top + l_1 * [ cos( th_1_pre ), sin( th_1_pre ) ] ...
    );
p_bottom_2 = formula( ...
    p_bottom_1 + l_2 * [ cos( th_1_pre + th_2_pre ), sin( th_1_pre + th_2_pre ) ] ...
    );

% calculate each segemnt`s CoM and the whole CoM
p_CoM_1 = formula( ...
    p_top + l_MCD_1 * [ cos( th_1_pre ), sin( th_1_pre ) ] ...
    );
p_CoM_2 = formula( ...
    p_bottom_1 + l_MCD_2 * [ cos( th_1_pre + th_2_pre ), sin( th_1_pre + th_2_pre ) ] ...
    );
p_CoM_whole = formula( ...
    ( m_1 * p_CoM_1 + m_2 * p_CoM_2 ) / ( m_1 + m_2 ) ...
    );

% calculate each segemnt`s CoM velocity and the whole CoM velocity
v_top = formula( diff( p_top, t) );
v_CoM_1 = formula( diff( p_CoM_1, t) );
v_CoM_2 = formula( diff( p_CoM_2, t) );
v_CoM_whole = formula( diff( p_CoM_whole, t) );

% calculate each segemnt`s CoM velocity and the whole CoM acceleration
a_top = formula( diff( v_top, t) );
a_CoM_1 = formula( diff( v_CoM_1, t) );
a_CoM_2 = formula( diff( v_CoM_2, t) );
a_CoM_whole = formula( diff( v_CoM_whole, t) );

%% prepare lagrangian & other physical quantities

% calculate differentiation of momentum.
% later this should be equal to [0, -g] + ext_F
momentum = formula( m_1 * v_CoM_1 + m_2 * v_CoM_2 );
d_momentum = formula( diff( momentum, t ) );

% calculate differentiation of angular momentum.
% later this should be equalt to tau_1.
% first, calculate angular momentum related to Inertial
AM_omega = I_1 * diff( th_1_pre, t ) + I_2 * diff( th_1_pre + th_2_pre, t );

% second, calculate angular momentum related to CoM position and velocity
AM_CoM_each_CoM = 0 ...
    + m_1 * cross( [ p_CoM_1 - p_CoM_whole, 0 ], [ v_CoM_1 - v_CoM_whole, 0 ] ) ...
    + m_2 * cross( [ p_CoM_2 - p_CoM_whole, 0 ], [ v_CoM_2 - v_CoM_whole, 0 ] ) ...
    ;
AM_CoM_each_CoM = formula( AM_CoM_each_CoM );
AM_CoM_each_CoM = AM_CoM_each_CoM(3);
AM_O_each_CoM = 0 ...
    + m_1 * cross( [ p_CoM_1 - p_top, 0 ], [ v_CoM_1 - v_top, 0 ] ) ...
    + m_2 * cross( [ p_CoM_2 - p_top, 0 ], [ v_CoM_2 - v_top, 0 ] ) ...
    ;
AM_O_each_CoM = formula( AM_O_each_CoM );
AM_O_each_CoM = AM_O_each_CoM(3);

% third, calculate angular momentum related to CoM position and velocity
AM_CoM = AM_CoM_each_CoM + AM_omega;
d_AM_CoM = formula( diff( AM_CoM , t ) );
AM_O = AM_O_each_CoM + AM_omega;
d_AM_O = formula( diff( AM_O , t ) );

% define kinetic energy
T = 0 ...
    + 1/2 * m_1 * (v_CoM_1 * v_CoM_1') ...
    + 1/2 * I_1 * diff( th_1_pre, t ).^2 ...
    + 1/2 * m_2 * (v_CoM_2 * v_CoM_2') ...
    + 1/2 * I_2 * diff( th_1_pre + th_2_pre, t ).^2 ...
    ;

% define potential energy
U = ( m_1 + m_2 ) * p_CoM_whole(2) * g;

% define KineticEnergy and KinmaticEnergy
KE = formula( T + U );
KME = formula( T );

% calculate Inertia around CoM
Inertia_CoM = 0 ...
    + m_1 * ( (p_CoM_1 - p_CoM_whole) * (p_CoM_1 - p_CoM_whole)' )...
    + m_2 * ( (p_CoM_2 - p_CoM_whole) * (p_CoM_2 - p_CoM_whole)' )...
    + I_1 ...
    + I_2 ...
    ;
Inertia_CoM = formula( Inertia_CoM );

% calculate length between O and CoM
L_CoM = sqrt( (p_top - p_CoM_whole) * (p_top - p_CoM_whole)' );
L_CoM = formula( L_CoM );
dL_CoM = formula( diff( L_CoM, t ) );
ddL_CoM = formula( diff( dL_CoM, t ) );

% define lagrangian
L = T - U;

% calc F by differetiating with connected_xb and connected_yb
Fx = -functionalDerivative(L, x_1);
Fy = -functionalDerivative(L, y_1);

% replace tmporary variable into actual variable
L = subs_restraing( L );
Fx = subs_restraing( Fx );
Fy = subs_restraing( Fy );
momentum = subs_restraing( momentum );
d_momentum = subs_restraing( d_momentum );
AM_CoM = subs_restraing( AM_CoM );
d_AM_CoM = subs_restraing( d_AM_CoM );
AM_O = subs_restraing( AM_O );
d_AM_O = subs_restraing( d_AM_O );
p_top = subs_restraing( p_top );
p_CoM_1 = subs_restraing( p_CoM_1 );
p_CoM_2 = subs_restraing( p_CoM_2 );
p_CoM_whole = subs_restraing( p_CoM_whole );
v_top = subs_restraing( v_top );
v_CoM_1 = subs_restraing( v_CoM_1 );
v_CoM_2 = subs_restraing( v_CoM_2 );
v_CoM_whole = subs_restraing( v_CoM_whole );
a_top = subs_restraing( a_top );
a_CoM_1 = subs_restraing( a_CoM_1 );
a_CoM_2 = subs_restraing( a_CoM_2 );
a_CoM_whole = subs_restraing( a_CoM_whole );
KE = subs_restraing( KE );
KME = subs_restraing( KME );
Inertia_CoM = subs_restraing( Inertia_CoM );
L_CoM = subs_restraing( L_CoM );
dL_CoM = subs_restraing( dL_CoM );
ddL_CoM = subs_restraing( ddL_CoM );

%% prepare equations to solve

% calc lagrangian equation for each variable
% all the equations are written with variables like diff(th_1_pre(t), t, t)
th_1_eq = -functionalDerivative( L, th_1_pre ) == tau_1;
th_2_eq = -functionalDerivative( L, th_2_pre ) == tau_2;

equations = [
    th_1_eq;
    th_2_eq;
    ];

% because diff(th_1_pre(t), t, t) is not capable of solving in MATLAB,
% replace them t-separeted variables by subs_pre
% after this procedure, equations are written in form of 
% alpha_1

equations = subs_pre( equations );
Fx = subs_pre( Fx );
Fy = subs_pre( Fy );
momentum = subs_pre( momentum );
d_momentum = subs_pre( d_momentum );
AM_CoM = subs_pre( AM_CoM );
d_AM_CoM = subs_pre( d_AM_CoM );
AM_O = subs_pre( AM_O );
d_AM_O = subs_pre( d_AM_O );
p_top = subs_pre( p_top );
p_CoM_whole = subs_pre( p_CoM_whole );
v_top = subs_pre( v_top );
v_CoM_whole = subs_pre( v_CoM_whole );
a_top = subs_pre( a_top );
a_CoM_whole = subs_pre( a_CoM_whole );
KE = subs_pre( KE );
KME = subs_pre( KME );
Inertia_CoM = subs_pre( Inertia_CoM );
L_CoM = subs_pre( L_CoM );
dL_CoM = subs_pre( dL_CoM );
ddL_CoM = subs_pre( ddL_CoM );

%% solve the equations

% solve the equation by linear argebra
% because all the two differentiated variables are linear in the equations
% solved output X is equalt to each variables
variables = [alpha_1, alpha_2]';
[A, B] = equationsToMatrix(equations, variables);
if det(A) == 0
    error('det(A) == 0')
end
X = inv(A)*B;

subs_ans = @(input) subs( input, variables, X );

% verify equations by check d_momentum and d_angular_momentum is equal to designated variable
% by replacing variables into solved output X,
% the formula are now written in form of tau_1
d_momentum = subs_ans( d_momentum );
d_AM_CoM = subs_ans( d_AM_CoM );
d_AM_O = subs_ans( d_AM_O );

% replace two differentiated variables into solved output X
% now Fx and Fy are written in form of tau_1
Fx = subs_ans( Fx );
Fy = subs_ans( Fy );

% define ext_F and ext_N, which are external force and external momentum
% respectively.
ext_F = formula( [Fx, Fy] );
ext_N = cross( [p_top, 0] - [p_CoM_whole, 0], [ext_F, 0] );
ext_N = formula( ext_N );
ext_N = ext_N(3);

% verify d_momentum = ext_F, d_angular_momentum = ext_N
simplify( d_momentum - ext_F )
simplify( d_AM_CoM - ext_N )
simplify( d_AM_O - ( -(m_1 + m_2)*g*p_CoM_whole(1) ) )

%% solve the equations for [tau_1, tau_2]'
variables = [tau_1, tau_2]';
[A, B] = equationsToMatrix(equations, variables);
if det(A) == 0
    error('det(A) == 0')
end
X_tau = inv(A)*B;

%% solve the equations for [alpha_1, tau_2]'
variables = [alpha_1, tau_2]';
[A, B] = equationsToMatrix(equations, variables);
if det(A) == 0
    error('det(A) == 0')
end
X_only_tau_2 = inv(A)*B;

%% output functions

parallel.defaultClusterProfile('local');
c = parcluster();

% output X
job = createJob(c);
createTask(job, @matlabFunction, 1, ...
    { ...
        X(1), X(2), ...
        'file', 'find_acc.m', ...
        'outputs', {'alpha_1', 'alpha_2'} ...
    });
submit(job)
job.Tasks

% output p_CoM_whole
job = createJob(c);
createTask(job, @matlabFunction, 1, ...
    { ...
        p_CoM_whole, ...
        'file', 'find_p_CoM_whole.m', ...
        'outputs', {'p_CoM_whole'} ...
    });

submit(job)
job.Tasks
job = createJob(c);
createTask(job, @matlabFunction, 1, ...
    { ...
        p_CoM_whole(1), ...
        'file', 'find_px_CoM_whole.m', ...
        'outputs', {'px_CoM_whole'} ...
    });
submit(job)

job.Tasks
job = createJob(c);
createTask(job, @matlabFunction, 1, ...
    { ...
        p_CoM_whole(2), ...
        'file', 'find_py_CoM_whole.m', ...
        'outputs', {'py_CoM_whole'} ...
    });
submit(job)
job.Tasks

% output v_CoM_whole
job = createJob(c);
createTask(job, @matlabFunction, 1, ...
    { ...
        v_CoM_whole, ...
        'file', 'find_v_CoM_whole.m', ...
        'outputs', {'v_CoM_whole'} ...
    });
submit(job)
job.Tasks

job = createJob(c);
createTask(job, @matlabFunction, 1, ...
    { ...
        v_CoM_whole(1), ...
        'file', 'find_vx_CoM_whole.m', ...
        'outputs', {'vx_CoM_whole'} ...
    });
submit(job)
job.Tasks

job = createJob(c);
createTask(job, @matlabFunction, 1, ...
    { ...
        v_CoM_whole(2), ...
        'file', 'find_vy_CoM_whole.m', ...
        'outputs', {'vy_CoM_whole'} ...
    });
submit(job)
job.Tasks

% output Fx
Fx = formula( Fx );
job = createJob(c);
createTask(job, @matlabFunction, 1, ...
    { ...
        Fx, ...
        'file', 'find_Fx.m', ...
        'outputs', {'Fx'} ...
    });
submit(job)
job.Tasks

% output Fy
Fy = formula( Fy );
job = createJob(c);
createTask(job, @matlabFunction, 1, ...
    { ...
        Fy, ...
        'file', 'find_Fy.m', ...
        'outputs', {'Fy'} ...
    });
submit(job)
job.Tasks

% output v_top
job = createJob(c);
createTask(job, @matlabFunction, 1, ...
    { ...
        v_top, ...
        'file', 'find_v_top.m', ...
        'outputs', {'v_top'} ...
    });
submit(job)
job.Tasks

% output a_top
a_top = formula( a_top );
job = createJob(c);
createTask(job, @matlabFunction, 1, ...
    { ...
        a_top, ...
        'file', 'find_a_top.m', ...
        'outputs', {'a_top'} ...
    });
submit(job)
job.Tasks

% output a_CoM_whole
a_CoM_whole = formula( a_CoM_whole );
job = createJob(c);
createTask(job, @matlabFunction, 1, ...
    { ...
        a_CoM_whole, ...
        'file', 'find_a_CoM_whole.m', ...
        'outputs', {'a_CoM_whole'} ...
    });
submit(job)
job.Tasks

% output X_tau
job = createJob(c);
createTask(job, @matlabFunction, 1, ...
    { ...
        X_tau(1), X_tau(2), ...
        'file', 'find_tau.m', ...
        'outputs', {'tau_1', 'tau_2'} ...
    });
submit(job)
job.Tasks

% output X_tau
job = createJob(c);
createTask(job, @matlabFunction, 1, ...
    { ...
        X_only_tau_2(1), X_only_tau_2(2), ...
        'file', 'find_alpha_1_tau_2.m', ...
        'outputs', {'alpha_1', 'tau_2'} ...
    });
submit(job)
job.Tasks

% output angular_momentum
job = createJob(c);
createTask(job, @matlabFunction, 1, ...
    { ...
        AM_CoM, ...
        'file', 'find_AM_CoM.m', ...
        'outputs', {'AM_CoM'} ...
    });
submit(job)
job.Tasks

% output angular_momentum
job = createJob(c);
createTask(job, @matlabFunction, 1, ...
    { ...
        AM_O, ...
        'file', 'find_AM_O.m', ...
        'outputs', {'AM_O'} ...
    });
submit(job)
job.Tasks

% output KME and KE
job = createJob(c);
createTask(job, @matlabFunction, 1, ...
    { ...
        KE, KME, ...
        'file', 'find_KE_KME.m', ...
        'outputs', {'KE', 'KME'} ...
    });
submit(job)
job.Tasks

% output Inertia_CoM
job = createJob(c);
createTask(job, @matlabFunction, 1, ...
    { ...
        Inertia_CoM, ...
        'file', 'find_Inertia_CoM.m', ...
        'outputs', {'Inertia_CoM'} ...
    });
submit(job)
job.Tasks

% output L_CoM, dL_CoM, ddL_CoM
ddL_CoM = subs_ans( ddL_CoM );
job = createJob(c);
createTask(job, @matlabFunction, 1, ...
    { ...
        L_CoM, ...
        'file', 'find_L_CoM.m', ...
        'outputs', {'L_CoM'} ...
    });
submit(job)
job.Tasks
job = createJob(c);
createTask(job, @matlabFunction, 1, ...
    { ...
        dL_CoM, ...
        'file', 'find_dL_CoM.m', ...
        'outputs', {'dL_CoM'} ...
    });
submit(job)
job.Tasks
job = createJob(c);
createTask(job, @matlabFunction, 1, ...
    { ...
        ddL_CoM, ...
        'file', 'find_ddL_CoM.m', ...
        'outputs', {'ddL_CoM'} ...
    });
submit(job)
job.Tasks