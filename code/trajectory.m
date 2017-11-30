function sol = trajectory(Ctensor,Dtensor,eulerAngle, maxTime)
%trajectory   particle trajectory via ICEP
%
%   sol = particleOrbit(Resistance, eulerAngle, maxTime) computes the 
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
sol = ode113(@dynamics, [0,maxTime], y0, options, Ctensor, Dtensor, EE);

% check quaternion normalization
error = abs(1 - sqrt((sol.y(4,end)^2 + sol.y(5,end)^2 + sol.y(6,end)^2 + sol.y(7,end)^2)));
if error > 10*1e-8
    fprintf('Error is %e\n',error);
end

[dydt, U, W] = dynamics(maxTime, sol.y(:,end), Ctensor,Dtensor, EE);
sol.U = U;
sol.W = W;

function [dydt, U, W] = dynamics(t, y, Ctensor,Dtensor, EE)
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


