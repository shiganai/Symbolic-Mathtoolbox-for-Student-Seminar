function [alpha_1,alpha_2] = find_acc(I_1,I_2,g,l_1,l_MCD_1,l_MCD_2,m_1,m_2,omega_1,omega_2,tau_1,tau_2,th_1,th_2)
%FIND_ACC
%    [ALPHA_1,ALPHA_2] = FIND_ACC(I_1,I_2,G,L_1,l_MCD_1,l_MCD_2,M_1,M_2,OMEGA_1,OMEGA_2,TAU_1,TAU_2,TH_1,TH_2)

%    This function was generated by the Symbolic Math Toolbox version 9.0.
%    01-Nov-2022 02:12:40

t2 = cos(th_1);
t3 = sin(th_1);
t4 = th_1+th_2;
t5 = l_1.^2;
t6 = l_MCD_1.^2;
t7 = l_MCD_2.^2;
t8 = m_2.^2;
t9 = omega_1.^2;
t10 = omega_2.^2;
t11 = I_1.*I_2;
t16 = -tau_1;
t12 = t2.^2;
t13 = t3.^2;
t14 = cos(t4);
t15 = sin(t4);
t19 = g.*l_1.*m_2.*t2;
t20 = g.*l_MCD_1.*m_1.*t2;
t17 = t14.^2;
t18 = t15.^2;
t21 = g.*l_MCD_2.*m_2.*t14;
t22 = I_2.*m_2.*t5.*t12;
t23 = I_2.*m_1.*t6.*t12;
t24 = l_1.*l_MCD_2.*m_2.*t2.*t14;
t25 = I_2.*m_2.*t5.*t13;
t26 = I_2.*m_1.*t6.*t13;
t27 = l_1.*l_MCD_2.*m_2.*t3.*t15;
t33 = l_1.*l_MCD_2.*m_2.*t3.*t9.*t14;
t34 = l_1.*l_MCD_2.*m_2.*t2.*t9.*t15;
t35 = l_1.*l_MCD_2.*m_2.*t3.*t10.*t14;
t36 = l_1.*l_MCD_2.*m_2.*t2.*t10.*t15;
t37 = l_1.*l_MCD_2.*m_2.*omega_1.*omega_2.*t3.*t14.*2.0;
t38 = l_1.*l_MCD_2.*m_2.*omega_1.*omega_2.*t2.*t15.*2.0;
t48 = t2.*t3.*t5.*t7.*t8.*t14.*t15.*2.0;
t28 = -t21;
t29 = m_2.*t7.*t17;
t30 = m_2.*t7.*t18;
t39 = -t38;
t40 = -t34;
t41 = -t36;
t46 = t5.*t7.*t8.*t13.*t17;
t47 = t5.*t7.*t8.*t12.*t18;
t49 = -t48;
t31 = I_1.*t29;
t32 = I_1.*t30;
t42 = m_1.*t6.*t12.*t29;
t43 = m_1.*t6.*t13.*t29;
t44 = m_1.*t6.*t12.*t30;
t45 = m_1.*t6.*t13.*t30;
t50 = t28+t33+t40+tau_2;
t51 = I_2+t24+t27+t29+t30;
t52 = t16+t19+t20+t21+t35+t37+t39+t41;
t53 = t11+t22+t23+t25+t26+t31+t32+t42+t43+t44+t45+t46+t47+t49;
t54 = 1.0./t53;
alpha_1 = t51.*t54.*(t21-t33+t34-tau_2)-t52.*t54.*(I_2+t29+t30);
if nargout > 1
    alpha_2 = -t54.*(t21-t33+t34-tau_2).*(I_1+t24+t27+t51+m_1.*t6.*t12+m_2.*t5.*t12+m_1.*t6.*t13+m_2.*t5.*t13)+t51.*t52.*t54;
end
