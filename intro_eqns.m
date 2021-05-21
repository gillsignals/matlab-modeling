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

v_A_RA = p.kon_A_RA * y(m.A) * y(m.RA) - p.koff_A_RA * y(m.A_RA);
v_B_RB = p.kon_B_RB * y(m.B) * y(m.RB) - p.koff_B_RB * y(m.B_RB);

v_A_RA_CoR = p.kon_A_RA_CoR * y(m.A_RA) * y(m.CoR) - ...
    p.koff_A_RA_CoR * y(m.A_RA_CoR);
v_B_RB_CoR = p.kon_B_RB_CoR * y(m.B_RB) * y(m.CoR) - ...
    p.koff_B_RB_CoR * y(m.B_RB_CoR);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Equations

dydt = zeros(11,1);

dydt(m.A) = p.qA - v_A_RA - p.kclA * y(m.A);
dydt(m.B) = p.qB - v_B_RB - p.kclB * y(m.B);

dydt(m.RA) = -v_A_RA;
dydt(m.RB) = -v_B_RB;

dydt(m.CoR) = -v_A_RA_CoR - v_B_RB_CoR;

dydt(m.A_RA) = v_A_RA - v_A_RA_CoR;
dydt(m.B_RB) = v_B_RB - v_B_RB_CoR;

dydt(m.A_RA_CoR) = v_A_RA_CoR;
dydt(m.B_RB_CoR) = v_B_RB_CoR;

dydt(m.Acl) = p.kclA * y(m.A);
dydt(m.Bcl) = p.kclB * y(m.B);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

















