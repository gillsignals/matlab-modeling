% INTRO_DRIVER
%   Driver file for the example ligand binding model used as sample code
%   for modeling systems of ODEs in MATLAB
%
%   FUNCTIONS:
%       intro_main
%       intro_eqns

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% User-Specified Options

simulation  = 'globalsens';   % Which case to run
visibleQ    = 0;            % Option to display figures; 1 = display
saveQ       = 0;            % Option to save figures; 1 = save, 0 = nothing
outputQ     = 1;            % Option to save data output; 1 = save
closeQ      = 1;            % Option to close all figures; 1 = close
soundQ      = 1;            % Option to play sound when done running

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

switch simulation
    case 'equal'    % Similar concentrations of ligand and receptor
        
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
        
        % Plot title
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
        
        % Calculate base AUC
        aucA_RA0 = trapz(T0,Y0(:,m.A_RA));
        aucB_RB0 = trapz(T0,Y0(:,m.B_RB));
        aucA_RA_CoR0 = trapz(T0,Y0(:,m.A_RA_CoR));
        aucB_RB_CoR0 = trapz(T0,Y0(:,m.B_RB_CoR));
        
        % Plot title
        plot_title = 'Local Sensitivity';
        
        % ----- Local Sensitivity -----
        % Parameters to vary
        params = fieldnames(p);
        params = params(~strcmp(params, 'alpha'));
        
        % Store initial parameter values
        p0 = p;
        
        % Set percentage to vary parameters by
        delta = 0.05;
        
        % Initialize output matrix
        sens = zeros(4,length(params));
        
        % Vary each parameter, one at a time
        for i = 1:length(params)
            p = p0;
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
            
            % Calculate new AUC
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
        
        % Set up values for simulation
        ligands = logspace(0, -6, 121);
        receptors = logspace(1, 7, 121);
        
        % Initialize output array
        sens = zeros(length(ligands), length(receptors), 4);
        
        % Vary each parameter
        for l = 1:length(ligands)
            for r = 1:length(receptors)
                
                % Set time span for simulation
                tspan = [0 3600];
                
                % Specify initial conditions
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

% Toggle plot visibility
if visibleQ == 1; fig1 = figure; else; fig1 = figure('visible', 'off'); end

subplot(2,2,1)
plot(T/60,Y(:,m.A),T/60,Y(:,m.B), 'LineWidth', 2)
legend('Ligand A', 'Ligand B', 'Location', 'east')
title('Free Ligand')
xlabel('Time (min)')
ylabel('Concentration (nM)')
ax1 = gca;
ax1.FontSize = 12;

subplot(2,2,2)
plot(T/60,Y(:,m.RA),T/60,Y(:,m.RB), 'LineWidth', 2)
legend('Receptor A', 'Receptor B', 'Location', 'east')
title('Free Receptor')
xlabel('Time (min)')
ylabel('Concentration (nM)')
ax2 = gca;
ax2.FontSize = 12;

subplot(2,2,3)
plot(T/60,Y(:,m.A_RA),T/60,Y(:,m.B_RB), 'LineWidth', 2)
legend('A-RA Complex', 'B-RB Complex', 'Location', 'east')
title({'Ligand-Receptor', 'Complex'})
xlabel('Time (min)')
ylabel('Concentration (nM)')
ax3 = gca;
ax3.FontSize = 12;

subplot(2,2,4)
plot(T/60,Y(:,m.A_RA_CoR),T/60,Y(:,m.B_RB_CoR), 'LineWidth', 2)
legend('A-RA-CoR Complex', 'B-RB-CoR Complex', 'Location', 'east')
title({'Ligand-Receptor-Co-Receptor', 'Complex'})
xlabel('Time (min)')
ylabel('Concentration (nM)')
ax4 = gca;
ax4.FontSize = 12;

sgtitle(plot_title, 'FontSize', 16, 'FontWeight', 'bold')
fig1.PaperUnits = 'inches';
fig1.PaperPosition = [0 0 10 6]; % [left bottom width height]

% Save plot
if saveQ == 1
    filename = sprintf('Plots/SampleModel_%s', simulation);
    print(filename, '-dpng', '-r300')
end

% Save output variables
if outputQ == 1
    folder = sprintf('Output/SampleModel_%s_', simulation);
    if strcmp(simulation, 'localsens') == 1 || ...
            strcmp(simulation, 'globalsens') == 1
        save(strcat(folder,'sens'), 'sens')
    else
        save(strcat(folder,'T'), 'T')
        save(strcat(folder,'Y'), 'Y')
    end
end

% Close open figures
if closeQ == 1
    close all
end

% Play sound when code is finished running
if soundQ == 1
    load gong
    sound(y/16,Fs*2)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%














