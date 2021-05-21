% INTRO_DRIVER
%   Driver file for the example ligand binding model used as sample code
%   for modeling systems of ODEs in MATLAB
%
%   FUNCTIONS:
%       intro_main
%       intro_eqns

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Species

m.A         = 1;
m.B         = 2;
m.RA        = 3;
m.RB        = 4;
m.CoR       = 5;
m.A_RA      = 6;
m.B_RB      = 7;
m.A_RA_CoR  = 8;
m.B_RB_CoR  = 9;
m.Acl       = 10;
m.Bcl       = 11;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameters

% Parameters required for the simulation
p.V             = 0.1;  % mL
p.qA            = 2e-7; % nM/s
p.qB            = 1e-7; % nM/s
p.kclA          = 3e-6; % s-1
p.kclB          = 1e-6; % s-1
p.kon_A_RA      = 3e-4; % nM-1 s-1
p.koff_A_RA     = 1e-2; % s-1
p.kon_B_RB      = 8e-4; % nM-1 s-1
p.koff_B_RB     = 5e-3; % s-1
p.kon_A_RA_CoR  = 2e-4; % nM-1 s-1
p.koff_A_RA_CoR = 2e-2; % s-1
p.kon_B_RB_CoR  = 1e-4; % nM-1 s-1
p.koff_B_RB_CoR = 3e-2; % s-1

% Conditions for the simulation
cells = 1e5;                                        % cells
receptors = 1e4;                                    % # receptors/cell
avogadro = 6.022e23;                                % molecules/mole

% Alpha parameter to use in ligand equations for correct units
p.alpha = cells / p.V / avogadro * 1e12;
% Units = nM/(# receptors/cell)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run ODE solver for different cases

% Set time span for simulation
tspan = [0 3600];

% Specify initial conditions
y0 = zeros(1,11);
y0(m.A)   = 1e-2; % nM
y0(m.B)   = 5e-3; % nM
y0(m.RA)  = receptors * p.alpha; % nM
y0(m.RB)  = receptors * p.alpha; % nM
y0(m.CoR) = receptors * p.alpha; % nM

% Run ODE model
[T,Y] = intro_main(tspan,y0,m,p);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


















