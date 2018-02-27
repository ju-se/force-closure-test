function fc = is_forceclosure(H, THETA, PHI)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constants
% PHI MUST NOT BE ZERO!!!!!!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%H = 0.5;
%THETA = pi/4;
%PHI = 0.1;

R = 0.2; % radius of the object
D = 0.3; % the distance between A and B along the axis of the object
CW_CONE_A = 0.2; % the size of the wrench cone at A: the angle between the axis and the surface 
CW_CONE_B = 0.2; % the size of the wrench cone at B: the angle between the axis and the surface
CW_CONE_G = 0.7; % the size of the wrench cone at G: the angle between the axis and the surface

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Coordinates of A, B, G and their contact normals before rotating about y-axis
% relative to frame {g} attached to G
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
aprime_g = [R*cos(PHI) R*sin(PHI) 0.0]' - [R 0.0 0.0]';
bprime_g = [R*cos(PHI+pi) R*sin(PHI+pi) 0.0]' - [R 0.0 0.0]';
a_g = aprime_g + [0.0 0.0 H+D]'; % coordinate of A
b_g = bprime_g + [0.0 0.0 H]'; % coordinate of B
g_g = [0.0 0.0 0.0]'; % coordinate of G
na_g = bprime_g - aprime_g; 
na_g = na_g / norm(na_g); % contact normal at A
nb_g = aprime_g - bprime_g;
nb_g = nb_g / norm(nb_g); % contact normal at B
ng_g = [0.0 0.0 1.0]'; % contact normal at G

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Coordinates after rotating about y-axis of G frame by pi/2-theta
% relative to frame {g} attached to G
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ry_g = [cos(pi/2-THETA) 0.0 sin(pi/2-THETA); 0.0 1.0 0.0; -sin(pi/2-THETA) 0.0 cos(pi/2-THETA)];
a_g = ry_g * a_g;
b_g = ry_g * b_g;
na_g = ry_g * na_g;
nb_g = ry_g * nb_g;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Normal vector of the plane of ABG
% relative to frame {g} attached to G
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nabg_g = cross(a_g-g_g, b_g-g_g);
nabg_g = nabg_g / norm(nabg_g);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Angles between the contact normals at A, B, G and the plane of ABG
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ang_nA_ABG = abs(pi/2 - acos(dot(nabg_g, na_g)));
ang_nB_ABG = abs(pi/2 - acos(dot(nabg_g, nb_g)));
ang_nG_ABG = abs(pi/2 - acos(dot(nabg_g, ng_g)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check if the angles are less than the cone angles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ang_nA_ABG < CW_CONE_A
    is_ABGinA = true;
else
    is_ABGinA = false;
end
if ang_nB_ABG < CW_CONE_B
    is_ABGinB = true;
else
    is_ABGinB = false;
end
if ang_nG_ABG < CW_CONE_G
    is_ABGinG = true;
else
    is_ABGinG = false;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check if the cones intersect ABG. If so, return the wrenches on ABG
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
z = is_ABGinA * is_ABGinB * is_ABGinG;

if z == false
    fc = false;
else
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Find the wrenches that are the intersection of the cone at G and ABG
    % relative to frame {g} attached to G
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    eqnG = [1.0+(nabg_g(2)/nabg_g(1))^2, 2.0*nabg_g(2)*nabg_g(3)/(nabg_g(1)^2), (nabg_g(3)/nabg_g(1))^2-(tan(CW_CONE_G))^2];
    wGy_g = roots(eqnG);
    wG1_g = [-(nabg_g(2)/nabg_g(1))*wGy_g(1)-(nabg_g(3)/nabg_g(1)); wGy_g(1); 1.0];
    wG1_g = wG1_g/norm(wG1_g);
    wG2_g = [-(nabg_g(2)/nabg_g(1))*wGy_g(2)-(nabg_g(3)/nabg_g(1)); wGy_g(2); 1.0];
    wG2_g = wG2_g/norm(wG2_g);
    %check_1 = dot(w_G_1, n_ABG);
    %check_2 = dot(w_G_2, n_ABG);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Find the wrenches that are the intersection of the cone at A and ABG
    % relative to frame {a} attached to A
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    r_ga = [[cos(THETA); 0.0; sin(THETA)], cross(na_g, [cos(THETA); 0.0; sin(THETA)]), na_g];
    h_ga = [r_ga, a_g];
    h_ga = [h_ga; [0.0 0.0 0.0 1.0]];
    h_ag = [r_ga', -r_ga'*a_g];
    h_ag = [h_ag; [0.0 0.0 0.0 1.0]];
    %check = h_ag*[[cos(THETA); 0.0; sin(THETA)]; 0.0];
    a_a = h_ag * [a_g; 1.0];
    a_a = a_a(1:3);
    b_a = h_ag * [b_g; 1.0];
    b_a = b_a(1:3);
    g_a = h_ag * [g_g; 1.0];
    g_a = g_a(1:3);
    nabg_a = cross(a_a-g_a, b_a-g_a);
    nabg_a = nabg_a / norm(nabg_a);
    %check = h_ag * [nabg_g; 0.0];
    eqnA = [1.0+(nabg_a(2)/nabg_a(1))^2, 2.0*nabg_a(2)*nabg_a(3)/(nabg_a(1)^2), (nabg_a(3)/nabg_a(1))^2-(tan(CW_CONE_A))^2];
    wAy_a = roots(eqnA);
    wA1_a = [-(nabg_a(2)/nabg_a(1))*wAy_a(1)-(nabg_a(3)/nabg_a(1)); wAy_a(1); 1.0];
    wA1_a = wA1_a/norm(wA1_a);
    wA2_a = [-(nabg_a(2)/nabg_a(1))*wAy_a(2)-(nabg_a(3)/nabg_a(1)); wAy_a(2); 1.0];
    wA2_a = wA2_a/norm(wA2_a);
    %check_1 = dot(wA1_a, nabg_a);
    %check_2 = dot(wA2_a, nabg_a);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Find the wrenches that are the intersection of the cone at B and ABG
    % relative to frame {b} attached to B
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    r_gb = [[cos(THETA); 0.0; sin(THETA)], cross(nb_g, [cos(THETA); 0.0; sin(THETA)]), nb_g];
    h_gb = [r_gb, b_g];
    h_gb = [h_gb; [0.0 0.0 0.0 1.0]];
    h_bg = [r_gb', -r_gb'*b_g];
    h_bg = [h_bg; [0.0 0.0 0.0 1.0]];
    %check = h_ag*[[cos(THETA); 0.0; sin(THETA)]; 0.0];
    a_b = h_bg * [a_g; 1.0];
    a_b = a_b(1:3);
    b_b = h_bg * [b_g; 1.0];
    b_b = b_b(1:3);
    g_b = h_bg * [g_g; 1.0];
    g_b = g_b(1:3);
    nabg_b = cross(a_b-g_b, b_b-g_b);
    nabg_b = nabg_b / norm(nabg_b);
    %check = h_bg * [nabg_g; 0.0];
    eqnB = [1.0+(nabg_b(2)/nabg_b(1))^2, 2.0*nabg_b(2)*nabg_b(3)/(nabg_b(1)^2), (nabg_b(3)/nabg_b(1))^2-(tan(CW_CONE_B))^2];
    wBy_b = roots(eqnB);
    wB1_b = [-(nabg_b(2)/nabg_b(1))*wBy_b(1)-(nabg_b(3)/nabg_b(1)); wBy_b(1); 1.0];
    wB1_b = wB1_b/norm(wB1_b);
    wB2_b = [-(nabg_b(2)/nabg_b(1))*wBy_b(2)-(nabg_b(3)/nabg_b(1)); wBy_b(2); 1.0];
    wB2_b = wB2_b/norm(wB2_b);
    %check_1 = dot(wB1_b, nabg_b);
    %check_2 = dot(wB2_b, nabg_b);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % All the wrenches in G frame
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    wA1_g = h_ga * [wA1_a; 0.0];
    wA1_g = wA1_g(1:3);
    wA2_g = h_ga * [wA2_a; 0.0];
    wA2_g = wA2_g(1:3);
    wB1_g = h_gb * [wB1_b; 0.0];
    wB1_g = wB1_g(1:3);
    wB2_g = h_gb * [wB2_b; 0.0];
    wB2_g = wB2_g(1:3);
    %check1 = dot(wG1_g, nabg_g);
    %check2 = dot(wG2_g, nabg_g);
    %check3 = dot(wA1_g, nabg_g);
    %check4 = dot(wA2_g, nabg_g);
    %check5 = dot(wB1_g, nabg_g);
    %check6 = dot(wB2_g, nabg_g);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % All the wrenches projected on the plane of ABG
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    xabg_g = (a_g - g_g)/norm(a_g - g_g); % the x-axis on ABG: GA
    yabg_g = cross(nabg_g, xabg_g); % the y-axis on ABG
    %ax_abg = dot(a_g, xabg_g);
    %ay_abg = dot(a_g, yabg_g);
    %bx_abg = dot(b_g, xabg_g);
    %by_abg = dot(b_g, yabg_g);
    wrenches(:,1) = [dot(wG1_g, xabg_g); dot(wG1_g, yabg_g); 0.0]; % from wG1_g
    wrenches(:,2) = [dot(wG2_g, xabg_g); dot(wG2_g, yabg_g); 0.0]; % from wG2_g
    wrenches(:,3) = [dot(wA1_g, xabg_g); dot(wA1_g, yabg_g); dot(cross(a_g, wA1_g), nabg_g)]; % from wA1_g
    wrenches(:,4) = [dot(wA2_g, xabg_g); dot(wA2_g, yabg_g); dot(cross(a_g, wA2_g), nabg_g)]; % from wA2_g
    wrenches(:,5) = [dot(wB1_g, xabg_g); dot(wB1_g, yabg_g); dot(cross(b_g, wB1_g), nabg_g)]; % from wB1_g
    wrenches(:,6) = [dot(wB2_g, xabg_g); dot(wB2_g, yabg_g); dot(cross(b_g, wB1_g), nabg_g)]; % from wB2_g
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Linear program: check force-closure on the plane of ABG
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if rank(wrenches) < 3
        fc = false;
    else
        [xx,fval,exitflag,output] = linprog([1 1 1 1 1 1],[],[],wrenches,[0.0; 0.0; 0.0],[1.0, 1.0, 1.0, 1.0, 1.0, 1.0],[]);
        if exitflag == 1
            fc = true;
        else
            fc = false;
        end
    end
end