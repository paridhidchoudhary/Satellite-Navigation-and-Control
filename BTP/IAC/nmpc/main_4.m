% NMPC for Spacecraft - Component Isolation Debug
clear all; close all; clc;
import casadi.*

% Create CasADi test function to check cross product
disp('========== TEST 1: Cross Product ==========');
a = MX.sym('a', 3);
b = MX.sym('b', 3);

% Define cross product as inline calculation
c_inline = [a(2)*b(3) - a(3)*b(2);
            a(3)*b(1) - a(1)*b(3);
            a(1)*b(2) - a(2)*b(1)];

% Create a function
cross_test = Function('cross_test', {a, b}, {c_inline});

% Test with numerical values
a_val = [1; 2; 3];
b_val = [4; 5; 6];
c_val = cross_test(a_val, b_val);

disp('Cross product test:');
disp(['a = [', num2str(a_val'), ']']);
disp(['b = [', num2str(b_val'), ']']);
disp(['a × b = [', num2str(full(c_val)'), ']']);
disp(['MATLAB cross check: [', num2str(cross(a_val, b_val)'), ']']);


% Test 2: Check rotation matrix skew-symmetric operation
disp('========== TEST 2: Rotation Matrix Dynamics ==========');
R = MX.sym('R', 3, 3);
omega = MX.sym('omega', 3);

omega_skew = [0, -omega(3), omega(2);
              omega(3), 0, -omega(1);
              -omega(2), omega(1), 0];
          
R_dot = R * omega_skew;

R_dot_func = Function('R_dot_func', {R, omega}, {R_dot});

% Test with numerical values
R_val = eye(3);
omega_val = [0.1; 0.2; 0.3];
R_dot_val = R_dot_func(R_val, omega_val);

disp('Rotation matrix dynamics test:');
disp('R = eye(3)'); 
disp(['omega = [', num2str(omega_val'), ']']);
disp('R_dot = ');
disp(full(R_dot_val));

% Test 3: Angular velocity dynamics (simplified)
disp('========== TEST 3: Angular Velocity Dynamics ==========');

% Constants
Jc = diag([6.0833, 1.5, 6.0833]) * 1e3;
Jt = diag([300, 250, 350]);

% Symbolic variables
R_ct = MX.sym('R_ct', 3, 3);
omega_ct = MX.sym('omega_ct', 3);
omega_cce = MX.sym('omega_cce', 3);
omega_tte = MX.sym('omega_tte', 3);
Tau = MX.sym('Tau', 3);

% Angular velocity dynamics
Jc_inv = inv(Jc); % Using explicit inverse instead of backslash
Jt_inv = inv(Jt);

term1 = Jc_inv * cross_casadi(Jc*omega_cce, omega_cce);
term2 = Jc_inv * Tau;
term3 = -Jt_inv * cross_casadi(Jt*omega_tte, omega_tte);
term4 = cross_casadi(omega_ct, R_ct * omega_cce);

omega_dot = term1 + term2 + term3 + term4;

omega_dot_func = Function('omega_dot_func', {R_ct, omega_ct, omega_cce, omega_tte, Tau}, ...
                          {omega_dot, term1, term2, term3, term4});

% Test with numerical values
R_ct_val = eye(3);
omega_ct_val = [-0.0046; -0.0046; -0.0046];
omega_cce_val = [0; 0; 0];
omega_tte_val = [0.0046; 0.0046; 0.0046];
Tau_val = [0; 0; 0];

[omega_dot_val, term1_val, term2_val, term3_val, term4_val] = omega_dot_func(R_ct_val, omega_ct_val, ...
                                                                omega_cce_val, omega_tte_val, Tau_val);

disp('Angular velocity dynamics test with zero inputs:');
disp(['omega_dot = [', num2str(full(omega_dot_val)'), ']']);
disp('Terms:');
disp(['term1 (Jc\(cross(Jc*omega_cce, omega_cce))): [', num2str(full(term1_val)'), ']']);
disp(['term2 (Jc\Tau): [', num2str(full(term2_val)'), ']']);
disp(['term3 (-Jt\(cross(Jt*omega_tte, omega_tte))): [', num2str(full(term3_val)'), ']']);
disp(['term4 (cross(omega_ct, R_ct*omega_cce)): [', num2str(full(term4_val)'), ']']);

% Test 4: Position dynamics (simplified)
disp('========== TEST 4: Position Dynamics ==========');

% Constants
m_c = 4000;
r_pc_c = [0; 0; -1.75];

% Symbolic variables
F = MX.sym('F', 3);
omega_ct = MX.sym('omega_ct', 3);
omega_dot = MX.sym('omega_dot', 3);
R_ct = MX.sym('R_ct', 3, 3);

% Position dynamics
rot_pc = R_ct * r_pc_c;
term1 = F / m_c;
term2 = cross_casadi(omega_dot, rot_pc);
term3 = cross_casadi(omega_ct, cross_casadi(omega_ct, rot_pc));
rho_ddot = term1 + term2 + term3;

rho_ddot_func = Function('rho_ddot_func', {F, omega_ct, omega_dot, R_ct}, ...
                         {rho_ddot, term1, term2, term3});

% Test with numerical values
F_val = [0; 0; 0];
omega_dot_val = [0; 0; 0];
R_ct_val = eye(3);
omega_ct_val = [-0.0046; -0.0046; -0.0046];

[rho_ddot_val, term1_val, term2_val, term3_val] = rho_ddot_func(F_val, omega_ct_val, ...
                                                                 omega_dot_val, R_ct_val);

disp('Position dynamics test with zero force:');
disp(['rho_ddot = [', num2str(full(rho_ddot_val)'), ']']);
disp('Terms:');
disp(['term1 (F/m_c): [', num2str(full(term1_val)'), ']']);
disp(['term2 (cross(omega_dot, rot_pc)): [', num2str(full(term2_val)'), ']']);
disp(['term3 (cross(omega_ct, cross(omega_ct, rot_pc))): [', num2str(full(term3_val)'), ']']);

% Test 5: Simple NLP problem to check if IPOPT works
disp('========== TEST 5: Simple NLP Test ==========');

opti = Opti();

% Define simple problem
x = opti.variable(2);
opti.minimize(x(1)^2 + x(2)^2);
opti.subject_to(x(1)^2 + x(2) >= 1);

% Solver options
opts = struct;
opts.ipopt.print_level = 5;
opti.solver('ipopt', opts);

try
    sol = opti.solve();
    disp('Simple NLP solved successfully!');
    disp(['Solution: [', num2str(sol.value(x)'), ']']);
catch ME
    disp(['Simple NLP failed: ', ME.message]);
end

% Test 6: Try a very simplified version of our actual problem
disp('========== TEST 6: Simplified NMPC First Step ==========');

opti = Opti();

% Constants
m_c = 4000;
fmax = 20;

% Just one-step prediction
X = opti.variable(6); % Only position and velocity [rho; rhodot]
U = opti.variable(3); % Just a single 3D force vector

% Initial state
rho_0 = [-1.1394; 113.0847; -7.6204];
rhodot_0 = [0.0081; 0.0039; 1.3490];
x0 = [rho_0; rhodot_0];

% Dynamics (simple double integrator)
Ts = 60;
A = [eye(3), Ts*eye(3); zeros(3), eye(3)];
B = [0.5*Ts^2*eye(3)/m_c; Ts*eye(3)/m_c];

% Constraint: Initial state
opti.subject_to(X == A*x0 + B*U);

% Constraint: Force magnitude
opti.subject_to(U'*U <= fmax^2);

% Cost function: Minimize position error
opti.minimize(X(1:3)'*X(1:3));

% Solver options
opts = struct;
opts.ipopt.print_level = 5;
opti.solver('ipopt', opts);

try
    sol = opti.solve();
    disp('Simple spacecraft NMPC step solved successfully!');
    disp(['Position: [', num2str(sol.value(X(1:3))'), ']']);
    disp(['Velocity: [', num2str(sol.value(X(4:6))'), ']']);
    disp(['Force: [', num2str(sol.value(U)'), ']']);
catch ME
    disp(['Simple spacecraft NMPC step failed: ', ME.message]);
end

disp('========== ALL TESTS COMPLETE ==========');


% Define helper function for subsequent tests
function c = cross_casadi(a, b)
    c = [a(2)*b(3) - a(3)*b(2);
         a(3)*b(1) - a(1)*b(3);
         a(1)*b(2) - a(2)*b(1)];
end
