% INTRO_DRIVER
%   Driver file for the example ligand binding model used as sample code
%   for modeling systems of ODEs in MATLAB
%
%   FUNCTIONS:
%       intro_main
%       intro_eqns

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% User-Specified Options

% I generally have saving figures, displaying figures, etc. be options that
% I can toggle depending on what I am trying to do with the simulation -
% all of these values are checked with if statements later in the code

% "visibleQ" is helpful for when you want to save several plots each
% iteration and don't need to look at them

% "soundQ" is helpful for slow simulations (like global sensitivity) so you
% can work on other things and know when the code is done

simulation  = 'equal';      % Which case to run with switch-case
visibleQ    = 0;            % Option to display figures; 1 = display
saveQ       = 0;            % Option to save figures; 1 = save, 0 = nothing
outputQ     = 0;            % Option to save data output; 1 = save
closeQ      = 0;            % Option to close all figures; 1 = close
soundQ      = 0;            % Option to play sound when done running

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Species

% Setting up the numbers for each species to refer to throughout the code
% (same as the reference file, "intro_species.m")

% For more complicated models, you can write a reusable function to create
% the molecule structure

% Species numbered list for use in equations and initial conditions
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

% Setting up the parameter structure to store all of the values needed for
% the model; all parameters are "dummy" values - this is a toy model
% loosely based on cytokines binding to receptors, and the parameter values
% are not based on a real system

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
p.alpha = cells / p.V / avogadro * 1e12;        % 1e12 converts from mol to nmol
% Units = nM/(# receptors/cell)

% I used receptors in terms of # receptors/cell, so I needed to convert the
% units on the receptor concentrations to nM to match the units in the rest
% of the model

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run ODE solver for different cases

% Switch-case is very useful for running different simulations - functions
% similarly to "if-else", but it is more streamlined for variables
% containing a single string; can set which simulation to run in the
% "User-Defined Options" section of the code

switch simulation
    case 'equal'    % Similar concentrations of ligand and receptor
        
        % Set time span for simulation (in seconds)
        tspan = [0 3600];
        
        % ode45 requires a start time and an end time, formatted as a
        % vector; if you want outputs at specific time points, you can
        % provide those values (e.g., tspan = 0:30:3600 will give outputs
        % every 30 seconds from 0 seconds to 3600 seconds)
        
        % Specify initial conditions
        y0 = zeros(1,11); % Initialize vector
        y0(m.A)   = 1e-2; % nM
        y0(m.B)   = 5e-3; % nM
        y0(m.RA)  = receptors * p.alpha; % nM
        y0(m.RB)  = receptors * p.alpha; % nM
        y0(m.CoR) = receptors * p.alpha; % nM
        
        % Need to provide the concentrations at the start time for all
        % molecules of the model (the zeros function provides zeros for all
        % of the molecules that aren't present at time 0)
        
        % Run ODE solver in the main function file
        [T,Y] = intro_main(tspan,y0,m,p);
        
        % Set plot title for specific case
        plot_title = '[Ligand] = [Receptor]';
        
    case 'receptor' % More receptor than ligand initially
        
        % Set time span for simulation
        tspan = [0 3600];
        
        % Specify initial conditions
        y0 = zeros(1,11);
        y0(m.A)   = 1e-2; % nM
        y0(m.B)   = 5e-3; % nM
        y0(m.RA)  = receptors * p.alpha * 100; % nM
        y0(m.RB)  = receptors * p.alpha * 100; % nM
        y0(m.CoR) = receptors * p.alpha * 100; % nM
        
        % Run ODE model
        [T,Y] = intro_main(tspan,y0,m,p);
        
        % Plot title
        plot_title = '[Ligand] < [Receptor]';
        
    case 'ligand'   % More ligand than receptor initially
        
        % Set time span for simulation
        tspan = [0 3600];
        
        % Specify initial conditions
        y0 = zeros(1,11);
        y0(m.A)   = 1e-2 * 100; % nM
        y0(m.B)   = 5e-3 * 100; % nM
        y0(m.RA)  = receptors * p.alpha; % nM
        y0(m.RB)  = receptors * p.alpha; % nM
        y0(m.CoR) = receptors * p.alpha; % nM
        
        % Run ODE model
        [T,Y] = intro_main(tspan,y0,m,p);
        
        % Plot title
        plot_title = '[Ligand] > [Receptor]';
        
    case 'localsens'   % Local sensitivity analysis
        
        % In local sensitivity analysis, we vary each parameter by the same
        % percentage, and then we examine the percent change in output
        % divided by the percent change in parameter
        
        % The output can be whatever we are interested in - for this model,
        % I examined how the AUCs of the different complexes were
        % sensitive to the model parameters
        
        % AUC can be calculated with the trapz function - calculates the
        % area under the curve with the trapezoid rule;
        % syntax: "AUC = trapz(x,y)"
        
        % General idea is to first get the output(s) for the base parameter
        % values ("base case") and then vary each parameter one at a time
        % and compare the output(s) to the base case
        
        % ----- Base case -----
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
        [T0,Y0] = intro_main(tspan,y0,m,p);
        
        % Calculate base AUC - output of the "base case"
        aucA_RA0 = trapz(T0,Y0(:,m.A_RA));
        aucB_RB0 = trapz(T0,Y0(:,m.B_RB));
        aucA_RA_CoR0 = trapz(T0,Y0(:,m.A_RA_CoR));
        aucB_RB_CoR0 = trapz(T0,Y0(:,m.B_RB_CoR));
        
        % Plot title
        plot_title = 'Local Sensitivity';
        
        % ----- Local Sensitivity -----
        % Parameters to vary
        params = fieldnames(p);  % Outputs parameter names as a cell array
        params = params(~strcmp(params, 'alpha'));
        % Not varying alpha because it is a unit conversion factor, not a
        % model parameter
        
        % Store initial parameter values for resetting after each change
        p0 = p;
        
        % Set percentage to vary parameters by
        delta = 0.05;
        
        % Initialize output matrix - important for speed of for loop
        sens = zeros(4,length(params));
        % Each row corresponds to one complex of interest (A-RA, B-RB,
        % A-RA-CoR, and B-RB-CoR); each column corresponds to a parameter
        % being varied
        
        % Vary each parameter, one at a time
        for i = 1:length(params)
            
            % Reset to initial values after each iteration
            p = p0;     
            
            % Calculate new value from original value and specified delta
            p.(params{i}) = p0.(params{i}) * (1+delta);
            
            % Alpha parameter to use in ligand equations for correct units
            % Have to recalculate when p.V changes
            p.alpha = cells / p.V / avogadro * 1e12;
            % Units = nM/(# receptors/cell)
            
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
            
            % Calculate new AUC after changing parameter value
            aucA_RA = trapz(T,Y(:,m.A_RA));
            aucB_RB = trapz(T,Y(:,m.B_RB));
            aucA_RA_CoR = trapz(T,Y(:,m.A_RA_CoR));
            aucB_RB_CoR = trapz(T,Y(:,m.B_RB_CoR));
            
            % Percent change in AUC over percent change in parameter
            sens(1,i) = ((aucA_RA - aucA_RA0)/aucA_RA0)/delta;
            sens(2,i) = ((aucB_RB - aucB_RB0)/aucB_RB0)/delta;
            sens(3,i) = ((aucA_RA_CoR - aucA_RA_CoR0)/...
                aucA_RA_CoR0)/delta;
            sens(4,i) = ((aucB_RB_CoR - aucB_RB_CoR0)/...
                aucB_RB_CoR0)/delta;
            
        end
        
    case 'globalsens'   % Global sensitivity analysis
        
        % In global sensitvity, we vary a few parameters or model inputs
        % over a large range to see how the model output changes across a
        % range of inputs
        
        % The output can be whatever we are interested in - for this model,
        % I examined how the AUCs of the different complexes were
        % sensitive to the initial ligand and receptor concentrations
        
        % AUC can be calculated with the trapz function - calculates the
        % area under the curve with the trapezoid rule;
        % syntax: "AUC = trapz(x,y)"
        
        % We do not need a "base case" for this analysis since we are
        % looking at the outputs over a large range of possible input
        % values
        
        % Set up values for simulation
        ligands = logspace(0, -6, 121);     % nM
        receptors = logspace(1, 7, 121);    % # receptors/cell
        % Will iterate over all possible combinations of these two input
        % vectors
        
        % Initialize output array - important for speed of loop
        sens = zeros(length(ligands), length(receptors), 4);
        % Rows correspond to changing ligand concentrations, columns
        % correspond to changing receptor concentrations, third dimension
        % (like a stack of matrices) corresponds to complex of interest
        % (A-RA, B-RB, A-RA-CoR, B-RB-CoR)
        
        % Set plot title for specific case
        plot_title = 'Global Sensitivity';
        
        % Vary each parameter - loop over each input vector
        for l = 1:length(ligands)
            for r = 1:length(receptors)
                
                % Set time span for simulation
                tspan = [0 3600];
                
                % Specify initial conditions from input vectors
                y0 = zeros(1,11);
                y0(m.A)   = ligands(l); % nM
                y0(m.B)   = ligands(l)/2; % nM
                y0(m.RA)  = receptors(r) * p.alpha; % nM
                y0(m.RB)  = receptors(r) * p.alpha; % nM
                y0(m.CoR) = receptors(r) * p.alpha; % nM
                
                % Run ODE model
                [T,Y] = intro_main(tspan,y0,m,p);
                
                % Calculate new AUC
                sens(l,r,1) = trapz(T,Y(:,m.A_RA));
                sens(l,r,2) = trapz(T,Y(:,m.B_RB));
                sens(l,r,3) = trapz(T,Y(:,m.A_RA_CoR));
                sens(l,r,4) = trapz(T,Y(:,m.B_RB_CoR));
                
            end
        end 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plots

% Toggle plot visibility based on user input
if visibleQ == 1; fig1 = figure; else; fig1 = figure('visible', 'off'); end

% Subplot generates one figure with multiple panels - arguments are rows,
% columns, and which figure you are working on
subplot(2,2,1)  % 2 rows, 2 columns, panel 1

% T/60 to convert time to minutes, 'LineWidth' of 2 makes the lines thicker
plot(T/60,Y(:,m.A),T/60,Y(:,m.B), 'LineWidth', 2)

% Legend will be labeled in order, 'location' allows the legend to be moved
legend('Ligand A', 'Ligand B', 'Location', 'east')
title('Free Ligand')
xlabel('Time (min)')
ylabel('Concentration (nM)')

% gca is for "get current axes" - necessary for some plot modifications,
% like changing font size
ax1 = gca;
ax1.FontSize = 12;

subplot(2,2,2)  % 2 rows, 2 columns, panel 1
plot(T/60,Y(:,m.RA),T/60,Y(:,m.RB), 'LineWidth', 2)
legend('Receptor A', 'Receptor B', 'Location', 'east')
title('Free Receptor')
xlabel('Time (min)')
ylabel('Concentration (nM)')
ax2 = gca;
ax2.FontSize = 12;

subplot(2,2,3)  % 2 rows, 2 columns, panel 3
plot(T/60,Y(:,m.A_RA),T/60,Y(:,m.B_RB), 'LineWidth', 2)
legend('A-RA Complex', 'B-RB Complex', 'Location', 'east')
title({'Ligand-Receptor', 'Complex'})
xlabel('Time (min)')
ylabel('Concentration (nM)')
ax3 = gca;
ax3.FontSize = 12;

subplot(2,2,4)  % 2 rows, 2 columns, panel 4
plot(T/60,Y(:,m.A_RA_CoR),T/60,Y(:,m.B_RB_CoR), 'LineWidth', 2)
legend('A-RA-CoR Complex', 'B-RB-CoR Complex', 'Location', 'east')
title({'Ligand-Receptor-Co-Receptor', 'Complex'})
xlabel('Time (min)')
ylabel('Concentration (nM)')
ax4 = gca;
ax4.FontSize = 12;

% sgtitle adds a main title to a figure with subplots
sgtitle(plot_title, 'FontSize', 16, 'FontWeight', 'bold')

% If you set figure as a variable (like "fig1 = figure;", done in the
% visibleQ if statement in this code), you can refer to it later to make
% modifications to the full figure, like the size of the figure
fig1.PaperUnits = 'inches';
fig1.PaperPosition = [0 0 10 7]; % [left bottom width height]

% sprintf builds a string from text and given variables, very handy for
% dynamic naming of figure and output files; "%s" is the flag for string
% variables, "%d" is flag for decimal variables, etc. 

% strcat concatenates multiple strings into one string, no spacing

% Save plot based on user input
if saveQ == 1
    
    % Create file name from simulation case
    filename = sprintf('Plots/SampleModel_%s', simulation);
    
    % Save figure as a .png file with 300 dpi resolution
    print(filename, '-dpng', '-r300')
end

% Save output variables based on user input
if outputQ == 1
    
    % Create part of file name from simulation
    folder = sprintf('Output/SampleModel_%s_', simulation);
    
    % Select which variables to save based on current simulation
    if strcmp(simulation, 'localsens') == 1 || ...
            strcmp(simulation, 'globalsens') == 1
        
        % Create full file name and save output as .mat file
        save(strcat(folder,'sens'), 'sens')
    else
        % Create full file name and save output as .mat file
        save(strcat(folder,'T'), 'T')
        save(strcat(folder,'Y'), 'Y')
    end
end

% Close any open figures based on user input
if closeQ == 1
    close all
end

% Play sound when code is finished running
% Useful for slow code, like the global sensitivity analysis
if soundQ == 1
    load gong
    sound(y/16,Fs*2)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%














