function [ ydot ] = odefun5( t,y )

    load('dhc2_vars.mat')
    
    Ix = 5368.39; Iy = 6.92893e3; Iz = 11158.75;
    Ixy = 0; Ixz = 1.1764e2; Iyz = 0;
    m = 2288.231; g = 9.81;
    
    Delv = -9.3083e-002; 
    Dail = 9.6242e-003; 
    Drud = -4.9242e-002; 
    Dflap = 0;
    n = 1800; 
    Pz = 20; 
    
    p = y(1); q = y(2); r = y(3); phi = y(4); theta = y(5); psi = y(6);
    u = y(7); v = y(8); w = y(9); xe = y(10); ye = y(11); z = y(12);
    
    V = sqrt(u^2 + v^2 + w^2);
    alpha = atan2(w,u);
    beta = asin(v/V);
    
    rho = 1.225*exp(-g*(z)/287.05/(288-0.0065*(z)));
    P = 0.7355*(-326.5+(0.00412*(Pz+7.4)*(n+2010)+(408-0.0965*n)*(1-rho/1.225)));
    dpt = 0.08696+191.18*(P*2/rho/V^3);
    qdyn = 0.5*rho*V^2;
    
    Cx = Cx0 + Cx_alpha*alpha + Cx_alpha2*alpha^2 + Cx_alpha3*alpha^3 + Cx_q*q*beaver_c/V...
         + Cx_dr*Drud + Cx_df*Dflap + Cx_df_alpha*alpha*Dflap...
         + Cx_dpt*dpt + Cx_dpt2_alpha*dpt^2*alpha;
     
    Cy = Cy0 + Cy_beta*beta + Cy_p*p*beaver_b/2/V + Cy_r*r*beaver_b/2/V...
         + Cy_da*Dail + Cy_dr*Drud + Cy_dr_alpha*alpha*Drud;
     
    Cz = Cz0 + Cz_alpha*alpha + Cz_alpha3*alpha^3 + Cz_q*q*beaver_c/V + Cz_de*Delv...
         + Cz_de_beta*Delv*beta^2 + Cz_df*Dflap...
         + Cz_df_alpha*alpha*Dflap + Cz_dpt*dpt;
     
    Cl = Cl0 + Cl_beta*beta + Cl_p*p*beaver_b/2/V + Cl_r*r*beaver_b/2/V...
         + Cl_da*Dail+Cl_dr*Drud...
         + Cl_da_alpha*alpha*Dail + Cl_alpha2_dpt*alpha^2*dpt;

    Cm = Cm0 + Cm_alpha*alpha + Cm_alpha2*alpha^2 + Cm_q*q*beaver_c/V...
         + Cm_de*Delv + Cm_beta2*beta^2 + Cm_r*r*beaver_b/2/V + Cm_df*Dflap...
         + Cm_dpt*dpt;

    Cn = Cn0 + Cn_beta*beta + Cn_p*p*beaver_b/2/V + Cn_r*r*beaver_b/2/V + Cn_da...
        *Dail + Cn_dr*Drud + Cn_q*q*beaver_c/V + Cn_beta3*beta^3 + Cn_dpt3*dpt^3;

    Fx = Cx*qdyn*beaver_S;
    Fy = Cy*qdyn*beaver_S;
    Fz = Cz*qdyn*beaver_S;
    L = Cl*qdyn*beaver_S*beaver_b;
    M = Cm*qdyn*beaver_S*beaver_c;
    N = Cn*qdyn*beaver_S*beaver_b;
    
    ydot(1) = (Iz*L + Ixz*N + Ixz^2*q*r - Iz^2*q*r - Ixz*Iy*p*q + 2*Ixz*Iz*p*q + Iy*Iz*q*r)/(Ixz^2 + Ix*Iz);
    ydot(2) = (M - Ixz*p^2 + Ixz*r^2 - Ix*q*r + Iz*q*r)/Iy;
    ydot(3) = -(Ixz*L - Ix*N + Ixz^2*p*q + Ix*Iy*p*q - Ix*Iz*p*q - Ix*Ixz*q*r + Ixz*Iy*q*r - Ixz*Iz*q*r)/(Ixz^2 + Ix*Iz);
 
    euler = [[1, sin(y(4))*tan(y(5)), cos(y(4))*tan(y(5))];...
             [0, cos(y(4)), -sin(y(4))];...
             [0, sin(y(4))*(1/cos(y(5))), cos(y(4))*(1/cos(y(5)))]] * ...
             [y(1); y(2); y(3)];
  
    ydot(4) = euler(1);
    ydot(5) = euler(2);
    ydot(6) = euler(3);
  
    ydot(7) = (Fx - g*m*sin(theta) - m*q*w + m*r*v)/m;
    ydot(8) = (Fy + m*p*w - m*r*u + g*m*cos(theta)*sin(phi))/m;
    ydot(9) = (Fz - m*p*v + m*q*u + g*m*cos(phi)*cos(theta))/m;
    
    C = [[cos(y(5))*cos(y(6)), sin(y(4))*sin(y(5))*cos(y(6))-(cos(y(4))*sin(y(6))), cos(y(4))*sin(y(5))*cos(y(6))+(sin(y(4))*sin(y(6)))];
         [cos(y(5))*sin(y(6)), sin(y(4))*sin(y(5))*sin(y(6))+(cos(y(4))*cos(y(6))), cos(y(4))*sin(y(5))*sin(y(6))-(sin(y(4))*cos(y(6)))];
         [-sin(y(5)), sin(y(4))*cos(y(5)), cos(y(4))*cos(y(5))]];
    
    Cmat = C * [y(7); y(8); y(9)];
    
    ydot(10) = Cmat(1);
    ydot(11) = Cmat(2);
    ydot(12) = Cmat(3);
    ydot = [ydot(1); ydot(2); ydot(3); ydot(4); ydot(5); ydot(6);...
            ydot(7); ydot(8); ydot(9); ydot(10); ydot(11); ydot(12)];
end








