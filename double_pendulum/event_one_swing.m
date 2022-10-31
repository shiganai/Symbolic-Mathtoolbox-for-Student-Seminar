function [value, isterminal, direction] = event_one_swing(t,q,constants,func_tau)

g = constants.g;

m_1 = constants.m_1;
l_1 = constants.l_1;
l_MCD_1 = constants.l_MCD_1;
I_1 = constants.I_1;

m_2 = constants.m_2;
l_2 = constants.l_2;
l_MCD_2 = constants.l_MCD_2;
I_2 = constants.I_2;

tau_vec = func_tau(t,q);
tau_1 = tau_vec(1,:,:);
tau_2 = tau_vec(2,:,:);

th_1 = q(1,:,:);
omega_1 = q(2,:,:);

th_2 = q(3,:,:);
omega_2 = q(4,:,:);

pop_num = size( th_1, 2 );
event_value_num = 2;
value = zeros( event_value_num, pop_num );
isterminal = zeros( size( value ) );
direction = zeros( size( value ) );
event_index = 0;

py_CoM_whole = find_py_CoM_whole(l_1,l_MCD_1,l_MCD_2,m_1,m_2,th_1,th_2);

vy_CoM_whole = find_vy_CoM_whole(l_1,l_MCD_1,l_MCD_2,m_1,m_2,omega_1,omega_2,th_1,th_2);

event_index = event_index + 1;
value( event_index,: ) = py_CoM_whole;
direction( event_index,: ) = 1;
isterminal( event_index, th_1 > 3/2*pi ) = 1;

event_index = event_index + 1;
value( event_index,: ) = vy_CoM_whole;
direction( event_index,: ) = -1;
isterminal( event_index, th_1 > 3/2*pi ) = 1;

event_index = event_index + 1;
value( event_index,: ) = abs( th_2 ) - 1/2*pi;
direction( event_index,: ) = 0;
isterminal( event_index, : ) = 1;

end

