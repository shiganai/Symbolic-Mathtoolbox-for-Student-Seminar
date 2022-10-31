function p_CoM_whole = find_p_CoM_whole(l_1,l_MCD_1,l_MCD_2,m_1,m_2,th_1,th_2)
%find_p_CoM_whole
%    p_CoM_whole = find_p_CoM_whole(L_1,l_MCD_1,l_MCD_2,M_1,M_2,TH_1,TH_2)

%    This function was generated by the Symbolic Math Toolbox version 9.0.
%    01-Nov-2022 02:12:41

t2 = cos(th_1);
t3 = sin(th_1);
t4 = m_1+m_2;
t5 = th_1+th_2;
t6 = 1.0./t4;
p_CoM_whole = [t6.*(m_2.*(l_1.*t2+l_MCD_2.*cos(t5))+l_MCD_1.*m_1.*t2),t6.*(m_2.*(l_1.*t3+l_MCD_2.*sin(t5))+l_MCD_1.*m_1.*t3)];
