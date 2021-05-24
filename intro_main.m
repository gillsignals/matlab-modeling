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

% AbsTol = absolute error tolerance; if the solution is smaller than this
% value, the solver will not try to get correct digits for this value

% RelTol = relative error tolerance; the error tolerance relative to the
% magnitude of each solution component (maximum acceptable error = maximum
% of RelTol*abs(y(i)) and AbsTol)

% NonNegative = specifies which outputs are required to be positive values;
% since this model is a biological system, all 11 equations are required to
% be positive values

% InitialStep = Suggested initial step size (upper bound, MATLAB might take
% a smaller step than this)

% Run the ODE solver on the equations in the eqns file
[T, Y] = ode45(@intro_eqns, tspan, y0, options, m, p);

% ode45 is the default ODE solver in MATLAB and should work for many
% applications; if ode45 is particularly slow, the system might be "stiff",
% and you can try ode15s (first choice when ode45 is inefficient)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Mole Balance

% Always need to check that the moles of each component are balanced -
% Input + Production - Output - Consumption - Still in System = 0

% Convert all concentrations to moles (or mass) and calculate all of the
% components of the mass balance - what was in the system at time 0, what
% has been added to the system since time 0, what has exited the system
% since time 0, if anything has been converted since time 0, and what
% remains in the system at the final time

% m is present in all of the code files to refer to the species numbers,
% can see "intro_species.m" for the full list (e.g., m.B = 2, so y0(m.B)
% will grab the 2nd value of y0, which is the initial concentration of B in
% the system)

% y0 = initial concentrations
% Y = output from the ODE solver (rows = time, columns = species in same
% order as m)

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

% Create a table to display if the moles are not balanced
Species = fieldnames(m);        % Grabs all of the species names from m
Species = Species(1:5);         % Select the ones participating in balance
MaxBalance = [max(abs(balanceA)); max(abs(balanceB)); ...
    max(abs(balanceRA)); max(abs(balanceRB)); max(abs(balanceCoR))];

% Specify tolerance for mole balance check
tol = 1e-13;

% Tolerance depends on the values of your outputs - my smallest outputs
% were in the range of 1e-10, so I selected a tolerance 1000 times smaller
% than my smallest outputs

% Print error message and table if the moles are not balanced
if any(MaxBalance > tol)
    disp('The moles of the components are not balanced.')
    balanceTable = table(Species,MaxBalance);
    disp(balanceTable)
end

% It is helpful to plot your mass balance versus time when you are first
% running the simulations, but you will ultimately want an automated way to
% check the mass balance, especially when you are running several
% simulations in a loop

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%










