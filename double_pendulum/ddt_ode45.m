function dotq = ddt_ode45(t, q, constants, func_tau)

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

[alpha_1,alpha_2] = find_acc(I_1,I_2,g,l_1,l_MCD_1,l_MCD_2,m_1,m_2,omega_1,omega_2,tau_1,tau_2,th_1,th_2);

%% output

dotq = [omega_1; alpha_1; omega_2; alpha_2];

end

