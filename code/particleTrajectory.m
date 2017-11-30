function [t,y] = particleTrajectory(Ctensor,Dtensor,eulerAngle, maxTime)
%particleTrajectory   particle trajectory via ICEP
%
%   [t,y] = particleOrbit(Resistance, eulerAngle, maxTime) computes the 
%   dynamical trajectory of a rigid particle in a viscous shear flow. The 
%   particle is characterized by its hydrodynamic resistance tensors input
%   in the structure Resistance. The initial orientation of the
%   particle is characterized by the Euler angle vector [phi,theta,psi]
%   using the (3,1,3) convention.  Flow is directed in the x-direction
%   with unit gradients in the y-direction, ux = y.  The dynamics are
%   integrated in time from t = 0 to t = maxTime.
%
%   The rigid body dynamics are performed using quaternions as detailed in
%   J. Diebel, Representing attitude: Euler angles, unit quaternions, and 
%   rotation vectors. Matrix 58, 1-35 (2006).  
%   
%   Kyle Bishop, June 20, 2015
% 


%% Initial Position and Orientation
% position is the 3D vector specifying the center of the particle in the
% world coordinate system.
position = [0,0,0]; 

% the initial orientation is specifed by eulerAngle [phi,theta,psi] using
% the 3,1,3 convention. A vector z in the world coordinate can be
% mapped to one in the body-fixed coordinate z' as 
% 
%    z' = R3(phi) * R1(theta) * R3(psi) * z
%
% where the R are rotation matrices (see Diebel section 2.4).
quaternion = eulerToQuaternion(eulerAngle);


%% Applied field, E
E = [0 0 1];
EE = reshape(E' * E, 9, 1);

 
%% Integrate particle dynamics
% initial condition for ode113
y0 = [position, quaternion];

% integration options
options = odeset('RelTol',1e-10,'AbsTol',1e-10);

% integrate
[t,y] = ode113(@dynamics, [0,maxTime], y0, options, Ctensor, Dtensor, EE);

% check quaternion normalization
error = abs(1 - sqrt((y(end,4)^2 + y(end,5)^2 + y(end,6)^2 + y(end,7)^2)));
if error > 10*1e-8
    fprintf('Error is %e\n',error);
end



function dydt = dynamics(t, y, Ctensor,Dtensor, EE)
%% equations of motion

% Preallocate
position = y(1:3);
quaternion = y(4:7);
dydt = zeros(7,1);

% Compute rotation matrix
Rq = quaternionToRotationMatrix(quaternion);

% Rotate resistance tensors into the world frame
C = tensorTransform3(Rq, Ctensor);
D = tensorTransform3(Rq, Dtensor);

% Compute particle velocity
U = reshape(C,3,9) * EE;
W = reshape(D,3,9) * EE;

% Compute position rates (translation velocity)
dydt(1:3) = U;   

% Diebel, equations 156 & 109
Qmatrix = [-quaternion(2), -quaternion(3), -quaternion(4); ...
            quaternion(1),  quaternion(4), -quaternion(3); ...
           -quaternion(4),  quaternion(1),  quaternion(2); ...
            quaternion(3), -quaternion(2),  quaternion(1)];
         
% Compute quaternion rates (angular velocity)
dydt(4:7) = 0.5 * Qmatrix * W;  


function Tw = tensorTransform2(R,Tb)
%% Transform Diadics: Tw_ij = Tb_pq * R_pi * R_qj
%
% Tb is a second order tensor in the body-fixed coordinates 
% R is a transformation matrix that maps from the world coordinates to the
% body-fxed coordinates as zb = R * zw
%
% Consider the following equation among vectors & tensors in the body-fixed
% coordinates: yb = Tb * xb.  To convert to world coordinates, we write
% that R * yw = Tb * R * xw or, equivalently, yw = (Rinv * Tb * R) * xw.
% Here, the tensor in parenthesis is the tensor Tw in the world
% coordinates.  Note also that Rinv = R' for transformation matrices.

Tw = R' * Tb * R;


function Tw = tensorTransform3(R,Tb)
%% Transform Triadics: Tw_ijk = Tb_pqr * R_pi * R_qj * R_rk
% 
% Tb is a second order tensor in the body-fixed coordinates 
% R is a rotation matrix that maps from the world coordinates to the
% body-fixed coordinates as zb = R * zw
%
% Consider the following equation among vectors & tensors in the body-fixed
% coordinates: yb = Tb * xb (where Tb is a 3rd order tensor, xb is a 2nd 
% order tensor, and yb is a vector).  To convert to world coordinates, we 
% write that R * yw = Tb * R * xw * R' or, equivalently,
%   
%   yw_i = R'_ip * Tb_pqr * (R_qj * xw_jk * R'_kr)
%   yw_i = (Tb_pqr * R_pi * R_qj * R_rk) * xw_jk
%

Tijr = zeros(3,3,3);
Tijr(:,:,1) = R' * Tb(:,:,1) * R; % T_ijr, r = 1
Tijr(:,:,2) = R' * Tb(:,:,2) * R; % T_ijr, r = 2
Tijr(:,:,3) = R' * Tb(:,:,3) * R; % T_ijr, r = 3
Tw = reshape( reshape(Tijr,9,3) * R, 3, 3, 3);

