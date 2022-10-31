function [KE,KME] = find_KE_KME(I_1,I_2,g,l_1,l_MCD_1,l_MCD_2,m_1,m_2,omega_1,omega_2,th_1,th_2)
%find_KE_KME
%    [KE,KME] = find_KE_KME(I_1,I_2,G,L_1,l_MCD_1,l_MCD_2,M_1,M_2,OMEGA_1,OMEGA_2,TH_1,TH_2)

%    This function was generated by the Symbolic Math Toolbox version 9.0.
%    01-Nov-2022 02:13:01

t2 = cos(th_1);
t3 = sin(th_1);
t4 = omega_1+omega_2;
t5 = th_1+th_2;
t6 = l_MCD_1.^2;
t7 = omega_1.^2;
t8 = t2.^2;
t9 = t3.^2;
t10 = cos(t5);
t11 = sin(t5);
t12 = l_1.*omega_1.*t2;
t13 = l_1.*omega_1.*t3;
t14 = t4.^2;
t15 = (I_1.*t7)./2.0;
t16 = l_MCD_2.*t4.*t10;
t17 = l_MCD_2.*t4.*t11;
t18 = (I_2.*t14)./2.0;
t19 = t6.*t7.*t8;
t20 = t6.*t7.*t9;
t21 = t12+t16;
t22 = t13+t17;
t25 = t19+t20;
t23 = t21.^2;
t24 = t22.^2;
t26 = (m_1.*t25)./2.0;
t27 = t23+t24;
t28 = (m_2.*t27)./2.0;
KE = t15+t18+t26+t28+g.*(m_2.*(l_1.*t3+l_MCD_2.*t11)+l_MCD_1.*m_1.*t3);
if nargout > 1
    KME = t15+t18+t26+t28;
end
