function dydt = intro_eqns(~,y,m,p)

% INTRO_EQNS [dydt] = intro_eqns(t,y,m,p)
%   Sets up the system of ordinary differential equations for the binding
%   of ligands A and B to the receptors A and B (RA and RB) and co-receptor
%   (CoR). Used as an input for an ODE solver
%
%   INPUT:
%       t = time, only present because it is necessary for the ODE solver,
%       represented by a tilde (~)
%       y = concentration of each component (in the order of the species 
%       structure), input as a vector, input needs to be the initial 
%       concentrations at the starting time
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
%       dydt = matrix of concentration values at each time point; columns
%       are in the order that the equations are defined in (same as species
%       structure order)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Equation Snippets

% Idea is that the full differential equations can be built from
% combinations of equation snippets

% Each snippet represents one binding complex with both on-binding and
% off-binding reactions included

% y() = the concentration of the molecules, m.species has the number that
% refers to each molecule (e.g., m.RA = 3, so y(m.RA) will select the 3rd 
% value of y(), corresponding to RA); "intro_species.m" has the full list)

% The p structure contains all of the parameters, "intro_parameters.m" has
% the full list and their default values

% Binding of each ligand to its respective receptor
v_A_RA = p.kon_A_RA * y(m.A) * y(m.RA) - p.koff_A_RA * y(m.A_RA);
v_B_RB = p.kon_B_RB * y(m.B) * y(m.RB) - p.koff_B_RB * y(m.B_RB);

% Binding of each ligand-receptor complex to the co-receptor
v_A_RA_CoR = p.kon_A_RA_CoR * y(m.A_RA) * y(m.CoR) - ...
    p.koff_A_RA_CoR * y(m.A_RA_CoR);
v_B_RB_CoR = p.kon_B_RB_CoR * y(m.B_RB) * y(m.CoR) - ...
    p.koff_B_RB_CoR * y(m.B_RB_CoR);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Equations

% dydt() = the derivative of the molecules' concentrations with respect to
% time, same idea as y() above

% Need one equation for each species in the model (both bound and free),
% also need one equation for each species that is cleared from the model
% (necessary for the mole balance)

% Each equation is built from the processes that the molecule can undergo -
% production, binding, clearance - and all of the binding processes are
% expressed with the equation snippets built above

% Signs need to reflect what is happening to the molecule - all of the
% equation snippets are positive for the product of the binding and
% negative for the components of the binding (e.g., v_A_RA needs to be
% negative for A and RA and positive for A_RA)

% Need to initialize the output vector to the same length as the number of
% equations
dydt = zeros(11,1);

% Change in ligand concentration - includes production of ligand (q),
% clearance of ligand (kcl * concentration), and binding of ligand to
% receptor
dydt(m.A) = p.qA - v_A_RA - p.kclA * y(m.A);
dydt(m.B) = p.qB - v_B_RB - p.kclB * y(m.B);

% Change in receptor concentrations - include all binding processes that
% the receptor or co-receptor can participate in
dydt(m.RA) = -v_A_RA;
dydt(m.RB) = -v_B_RB;

dydt(m.CoR) = -v_A_RA_CoR - v_B_RB_CoR;

% Change in complex concentations - formation of complex is positive,
% consumption of complex is negative
dydt(m.A_RA) = v_A_RA - v_A_RA_CoR;
dydt(m.B_RB) = v_B_RB - v_B_RB_CoR;

dydt(m.A_RA_CoR) = v_A_RA_CoR;
dydt(m.B_RB_CoR) = v_B_RB_CoR;

% Clearance of ligand from system
dydt(m.Acl) = p.kclA * y(m.A);
dydt(m.Bcl) = p.kclB * y(m.B);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

















