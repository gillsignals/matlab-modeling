function [T,Y] = intro_main(tspan,y0,m,p)

% INTRO_MAIN [T,Y] = intro_main(tspan,y0,m,p)
%   Solves the system of ODEs given in the intro_eqns file based on
%   the parameters given. Inputs are the parameters and the length of
%   time; outputs are the time and concentration matrices from the ODE
%   solver
%
%   INPUT:
%       times = vector of times to use for the beginning and end of the
%       system, given in seconds; if multiple periods, give the beginning,
%       end, and every change point
%       y0 = initial values for the model, input as a vector
%       m = structure containing the molecules included in the model and
%       their assigned numbers, input as a structure, correct naming for
%       species is ligand name or receptor complex name (with
%       components separated by underscores)
%       p = structure containing parameters for the model, input as a
%       structure, correct naming for parameters is "kon_" or "koff_",
%       followed by the ligand name or receptor complex name
%       (with components separated by underscores), ending with the
%       species being bound to; e.g., kon_A_RA_CoR
%   OUTPUT:
%       T = vector of time values corresponding to the concentrations in Y,
%       given in seconds
%       Y = matrix of concentration values output from the ODEs; each
%       column corresponds to an equation in order, each row is a time
%       point

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ODE Solver

% Set the ODE solver options
options = odeset('AbsTol', 1e-12,'RelTol', 1e-8, 'NonNegative', 1:11, ...
    'InitialStep', 1e-2);

% Run the ODE solver on the equations in the eqns file
[T, Y] = ode45(@intro_eqns, tspan, y0, options, m, p);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Mole Balance

% Check mole balance (in nmol)
balanceA = y0(m.A)*p.V + y0(m.A_RA)*p.V + y0(m.A_RA_CoR)*p.V + ...
    p.qA*T*p.V - Y(:,m.A)*p.V - Y(:,m.A_RA)*p.V - Y(:,m.A_RA_CoR)*p.V - ...
    Y(:,m.Acl)*p.V;
balanceB = y0(m.B)*p.V + y0(m.B_RB)*p.V + y0(m.B_RB_CoR)*p.V + ...
    p.qB*T*p.V - Y(:,m.B)*p.V - Y(:,m.B_RB)*p.V - Y(:,m.B_RB_CoR)*p.V - ...
    Y(:,m.Bcl)*p.V;

balanceRA = y0(m.RA)*p.V + y0(m.A_RA)*p.V + y0(m.A_RA_CoR)*p.V - ...
    Y(:,m.RA)*p.V - Y(:,m.A_RA)*p.V - Y(:,m.A_RA_CoR)*p.V;
balanceRB = y0(m.RB)*p.V + y0(m.B_RB)*p.V + y0(m.B_RB_CoR)*p.V - ...
    Y(:,m.RB)*p.V - Y(:,m.B_RB)*p.V - Y(:,m.B_RB_CoR)*p.V;

balanceCoR = y0(m.CoR)*p.V + y0(m.A_RA_CoR)*p.V + y0(m.B_RB_CoR)*p.V - ...
    Y(:,m.CoR)*p.V - Y(:,m.A_RA_CoR)*p.V - Y(:,m.B_RB_CoR)*p.V;

Species = fieldnames(m);
Species = Species(1:9);
MaxBalance = [max(abs(balanceA)); max(abs(balanceB)); ...
    max(abs(balanceRA)); max(abs(balanceRB)); max(abs(balanceCoR))];

% Specify tolerance for mole balance check
tol = 0.005;

% Print error message and table if the moles are not balanced
if any(MaxBalance > tol)
    disp('The moles of the components are not balanced.')
    balanceTable = table(Species,MaxBalance);
    disp(balanceTable)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%










