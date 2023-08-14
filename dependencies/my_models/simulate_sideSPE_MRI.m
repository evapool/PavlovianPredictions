function [VV, RPE, MB_V, SPE, UNC, Vcs, Tp]= simulate_sideSPE_MRI (alpha, eta, task_data, displayPlot, Vcs, Tp, v0)

%--------------------------------------------------------------------------
% input : free parameter alpha (learing rate for reward PE) and eta (learning
% rate for state PE), and information about the task (trial
% sequence: CS, outcome presentation order for the participant)
% output: timeseries of the values and the prediction order

%--------------------------------------------------------------------------
% last modified nov 2021

% get the specific outcome
sp_outcome = (task_data.US_identity*10) + task_data.US_side;

% define how the reward is tracked based on the model
task_data.US_identity(task_data.US_identity >0) = 1; % if there is a reward or not independenlty of the side
outcome = task_data.US_identity;


% initialise value and prediction error vectors to store the data.
VV   = nan (length(outcome),1);
RPE  = nan (length(outcome),1);
MB_V = nan (length(outcome),1);
SPE  = nan (length(outcome),1);
UNC  = nan (length(outcome),1);

for t = 1:length(outcome)
    
    CS       = task_data.CSname (t); %which CS is displayed
    r        = outcome (t); % reward
    sr       = sp_outcome(t); % specific code for the outcome
    
    switch CS
        
        case 31
            
            % MF
            dv = r-Vcs.AL; % compute prediction error
            VV (t) = Vcs.AL; % store value in timeserie
            RPE (t) = dv;    % store prediction error in timeserie
            
            % MB
            % most likly transition probability
            [T, idx] = max(structfun(@(x)max(x(:)),Tp.AL));
            % what is the value of the most likely subsequent status?
            tmp  = fieldnames(Tp.AL);
            name = char(tmp(idx));
            
            %V    = Vcs.(name);
            R     = 1;
            
            MB_V (t) = T * R;   % MB value is the transition probability multiplied by the value of the outcome found the s'
            
            UNC (t)  = 1 - T ; % 1 - max T
            
            % update MF
            Vcs.AL = Vcs.AL+ alpha * dv;% update CSpL values
            
            
            % update MB based on specific US
            switch sr
                
                case 31
                    sdv      = 1 - Tp.AL.AL;
                    SPE  (t) = sdv ;                 % store state prediction error
                    Tp.AL.AL = Tp.AL.AL + eta * sdv; % update transition probability
                    
                    % update transition probability for states not arrived in
                    Tp.AL.AR = Tp.AL.AR * (1-eta);
                    Tp.AL.BL = Tp.AL.BL * (1-eta);
                    Tp.AL.BR = Tp.AL.BR * (1-eta);
                    Tp.AL.m  = Tp.AL.m  * (1-eta);
                    
                case 32
                    sdv      = 1 - Tp.AL.AR;
                    SPE  (t) = sdv ;                 % store state prediction error
                    Tp.AL.AR = Tp.AL.AR + eta * sdv; % update transition probability
                    
                    % update transition probability for states not arrived in
                    Tp.AL.AL = Tp.AL.AL * (1-eta);
                    Tp.AL.BL = Tp.AL.BL * (1-eta);
                    Tp.AL.BR = Tp.AL.BR * (1-eta);
                    Tp.AL.m  = Tp.AL.m  * (1-eta);
                    
                case 41
                    sdv      = 1 - Tp.AL.BL;
                    SPE  (t) = sdv ;                 % store state prediction error
                    Tp.AL.BL = Tp.AL.BL + eta * sdv; % update transition probability
                    
                    % update transition probability for states not arrived in
                    Tp.AL.AL = Tp.AL.AL * (1-eta);
                    Tp.AL.AR = Tp.AL.AR * (1-eta);
                    Tp.AL.BR = Tp.AL.BR * (1-eta);
                    Tp.AL.m  = Tp.AL.m  * (1-eta);
                    
                case 42
                    sdv      = 1 - Tp.AL.BR;
                    SPE  (t) = sdv ;                 % store state prediction error
                    Tp.AL.BR = Tp.AL.BR + eta * sdv; % update transition probability
                    
                    % update transition probability for states not arrived in
                    Tp.AL.AL = Tp.AL.AL * (1-eta);
                    Tp.AL.AR = Tp.AL.AR * (1-eta);
                    Tp.AL.BL = Tp.AL.BL * (1-eta);
                    Tp.AL.m  = Tp.AL.m  * (1-eta);
                    
                case 0
                    sdv      = 1 - Tp.AL.m;
                    SPE  (t) = sdv ;                 % store state prediction error
                    Tp.AL.m = Tp.AL.m + eta * sdv;   % update transition probability
                    
                    % update transition probability for states not arrived in
                    Tp.AL.AL = Tp.AL.AL * (1-eta);
                    Tp.AL.AR = Tp.AL.AR * (1-eta);
                    Tp.AL.BL = Tp.AL.BL * (1-eta);
                    Tp.AL.BR = Tp.AL.BR * (1-eta);
                    
                    
            end
            
        case 32
            
            % MF
            dv = r-Vcs.AR;% compute prediction error
            VV (t) = Vcs.AR;% store value in timeserie
            RPE (t) = dv;% store prediction error in timeserie
            
            % MB
            % most likly transition probability
            [T, idx] = max(structfun(@(x)max(x(:)),Tp.AR));
            % what is the value of the most likely subsequent status?
            tmp  = fieldnames(Tp.AR);
            name = char(tmp(idx));
            
            %V    = Vcs.(name);
            R     = 1;
            
            MB_V (t) = T * R;   % MB value is the transition probability multiplied by the value of the outcome found the s'
            
            UNC (t)  = 1 - T ; % 1 - max T
            
            % update MF
            Vcs.AR = Vcs.AR + alpha*dv;% update CSpR values
            
            % update MB based on specific US
            switch sr
                
                case 31
                    sdv      = 1 - Tp.AR.AL;
                    SPE  (t) = sdv ;                 % store state prediction error
                    Tp.AR.AL = Tp.AR.AL + eta * sdv; % update transition probability
                    
                    % update transition probability for states not arrived in
                    Tp.AR.AR = Tp.AR.AR * (1-eta);
                    Tp.AR.BL = Tp.AR.BL * (1-eta);
                    Tp.AR.BR = Tp.AR.BR * (1-eta);
                    Tp.AR.m  = Tp.AR.m  * (1-eta);
                    
                case 32
                    sdv      = 1 - Tp.AR.AR;
                    SPE  (t) = sdv ;                 % store state prediction error
                    Tp.AR.AR = Tp.AR.AR + eta * sdv; % update transition probability
                    
                    % update transition probability for states not arrived in
                    Tp.AR.AL = Tp.AR.AL * (1-eta);
                    Tp.AR.BL = Tp.AR.BL * (1-eta);
                    Tp.AR.BR = Tp.AR.BR * (1-eta);
                    Tp.AR.m  = Tp.AR.m  * (1-eta);
                    
                case 41
                    sdv      = 1 - Tp.AR.BL;
                    SPE  (t) = sdv ;                 % store state prediction error
                    Tp.AR.BL = Tp.AR.BL + eta * sdv; % update transition probability
                    
                    % update transition probability for states not arrived in
                    Tp.AR.AL = Tp.AR.AL * (1-eta);
                    Tp.AR.AR = Tp.AR.AR * (1-eta);
                    Tp.AR.BR = Tp.AR.BR * (1-eta);
                    Tp.AR.m  = Tp.AR.m  * (1-eta);
                    
                case 42
                    sdv      = 1 - Tp.AR.BR;
                    SPE  (t) = sdv ;                 % store state prediction error
                    Tp.AR.BR = Tp.AR.BR + eta * sdv; % update transition probability
                    
                    % update transition probability for states not arrived in
                    Tp.AR.AL = Tp.AR.AL * (1-eta);
                    Tp.AR.AR = Tp.AR.AR * (1-eta);
                    Tp.AR.BL = Tp.AR.BL * (1-eta);
                    Tp.AR.m  = Tp.AR.m  * (1-eta);
                    
                case 0
                    sdv      = 1 - Tp.AR.m;
                    SPE  (t) = sdv ;                 % store state prediction error
                    Tp.AR.m = Tp.AR.m + eta * sdv; % update transition probability
                    
                    % update transition probability for states not arrived in
                    Tp.AR.AL = Tp.AR.AL * (1-eta);
                    Tp.AR.AR = Tp.AR.AR * (1-eta);
                    Tp.AR.BL = Tp.AR.BL * (1-eta);
                    Tp.AR.BR = Tp.AR.BR * (1-eta);
                    
            end
            
            
        case 41
            
            % MF
            
            dv = r-Vcs.BL;% compute prediction error
            VV (t) = Vcs.BL;% store value in timeserie
            RPE (t) = dv;% store prediction error in timeserie
            
            % MB
            % most likly transition probability
            [T, idx] = max(structfun(@(x)max(x(:)),Tp.BL));
            % what is the value of the most likely subsequent status?
            tmp  = fieldnames(Tp.BL);
            name = char(tmp(idx));
            
            %V    = Vcs.(name);
            R     = 1;
            
            MB_V (t) = T * R;   % MB value is the transition probability multiplied by the value of the outcome found the s'
            
            UNC (t)  = 1 - T ; % 1 - max T
            
            % update MF
            Vcs.BL = Vcs.BL + alpha*dv;% update CSpR values
            
            % update MB based on specific US
            switch sr
                
                case 31
                    sdv      = 1 - Tp.BL.AL;
                    SPE  (t) = sdv ;                 % store state prediction error
                    Tp.BL.AL = Tp.BL.AL + eta * sdv; % update transition probability
                    
                    % update transition probability for states not arrived in
                    Tp.BL.AR = Tp.BL.AR * (1-eta);
                    Tp.BL.BL = Tp.BL.BL * (1-eta);
                    Tp.BL.BR = Tp.BL.BR * (1-eta);
                    Tp.BL.m  = Tp.BL.m  * (1-eta);
                    
                case 32
                    sdv      = 1 - Tp.BL.AR;
                    SPE  (t) = sdv ;                 % store state prediction error
                    Tp.BL.AR = Tp.BL.AR + eta * sdv; % update transition probability
                    
                    % update transition probability for states not arrived in
                    Tp.BL.AL = Tp.BL.AL * (1-eta);
                    Tp.BL.BL = Tp.BL.BL * (1-eta);
                    Tp.BL.BR = Tp.BL.BR * (1-eta);
                    Tp.BL.m  = Tp.BL.m  * (1-eta);
                    
                case 41
                    sdv      = 1 - Tp.BL.BL;
                    SPE  (t) = sdv ;                 % store state prediction error
                    Tp.BL.BL = Tp.BL.BL + eta * sdv; % update transition probability
                    
                    % update transition probability for states not arrived in
                    Tp.BL.AL = Tp.BL.AL * (1-eta);
                    Tp.BL.AR = Tp.BL.AR * (1-eta);
                    Tp.BL.BR = Tp.BL.BR * (1-eta);
                    Tp.BL.m  = Tp.BL.m  * (1-eta);
                    
                case 42
                    sdv      = 1 - Tp.BL.BR;
                    SPE  (t) = sdv ;                 % store state prediction error
                    Tp.BL.BR = Tp.BL.BR + eta * sdv; % update transition probability
                    
                    % update transition probability for states not arrived in
                    Tp.BL.AL = Tp.BL.AL * (1-eta);
                    Tp.BL.AR = Tp.BL.AR * (1-eta);
                    Tp.BL.BL = Tp.BL.BL * (1-eta);
                    Tp.BL.m  = Tp.BL.m  * (1-eta);
                    
                case 0
                    sdv      = 1 - Tp.BL.m;
                    SPE  (t) = sdv ;                 % store state prediction error
                    Tp.BL.m = Tp.BL.m + eta * sdv; % update transition probability
                    
                    % update transition probability for states not arrived in
                    Tp.BL.AL = Tp.BL.AL * (1-eta);
                    Tp.BL.AR = Tp.BL.AR * (1-eta);
                    Tp.BL.BL = Tp.BL.BL * (1-eta);
                    Tp.BL.BR  = Tp.BL.BR  * (1-eta);
            end
            
        case 42
            
            % MF
            dv = r-Vcs.BR;% compute prediction error
            VV (t) = Vcs.BR;% store value in timeserie
            RPE (t) = dv;% store prediction error in timeserie
            
            % MB
            % most likly transition probability
            [T, idx] = max(structfun(@(x)max(x(:)),Tp.BR));
            % what is the value of the most likely subsequent status?
            tmp  = fieldnames(Tp.BR);
            name = char(tmp(idx));
            
            %V    = Vcs.(name);

            R     = 1;
            
            MB_V (t) = T * R;   % MB value is the transition probability multiplied by the value of the outcome found the s'
            UNC (t)  = 1 - T ; % 1 - max T
            
            % update MF
            Vcs.BR = Vcs.BR + alpha*dv;% update CSpR values
            
            % update MB based on specific US
            switch sr
                
                case 31
                    sdv      = 1 - Tp.BR.AL;
                    SPE  (t) = sdv ;                 % store state prediction error
                    Tp.BR.AL = Tp.BR.AL + eta * sdv; % update transition probability
                    
                    % update transition probability for states not arrived in
                    Tp.BR.AR = Tp.BR.AR * (1-eta);
                    Tp.BR.BL = Tp.BR.BL * (1-eta);
                    Tp.BR.BR = Tp.BR.BR * (1-eta);
                    Tp.BR.m  = Tp.BR.m  * (1-eta);
                    
                case 32
                    sdv      = 1 - Tp.BR.AR;
                    SPE  (t) = sdv ;                 % store state prediction error
                    Tp.BR.AR = Tp.BR.AR + eta * sdv; % update transition probability
                    
                    % update transition probability for states not arrived in
                    Tp.BR.AL = Tp.BR.AL * (1-eta);
                    Tp.BR.BL = Tp.BR.BL * (1-eta);
                    Tp.BR.BR = Tp.BR.BR * (1-eta);
                    Tp.BR.m  = Tp.BR.m  * (1-eta);
                    
                case 41
                    sdv      = 1 - Tp.BR.BL;
                    SPE  (t) = sdv ;                 % store state prediction error
                    Tp.BR.BL = Tp.BR.BL + eta * sdv; % update transition probability
                    
                    % update transition probability for states not arrived in
                    Tp.BR.AL = Tp.BR.AL * (1-eta);
                    Tp.BR.AR = Tp.BR.AR * (1-eta);
                    Tp.BR.BR = Tp.BR.BR * (1-eta);
                    Tp.BR.m  = Tp.BR.m  * (1-eta);
                    
                case 42
                    sdv      = 1 - Tp.BR.BR;
                    SPE  (t) = sdv ;                 % store state prediction error
                    Tp.BR.BR = Tp.BR.BR + eta * sdv; % update transition probability
                    
                    % update transition probability for states not arrived in
                    Tp.BR.AL = Tp.BR.AL * (1-eta);
                    Tp.BR.AR = Tp.BR.AR * (1-eta);
                    Tp.BR.BL = Tp.BR.BL * (1-eta);
                    Tp.BR.m  = Tp.BR.m  * (1-eta);
                    
                    
                case 0
                    sdv      = 1 - Tp.BR.m;
                    SPE  (t) = sdv ;                 % store state prediction error
                    Tp.BR.m = Tp.BR.m + eta * sdv; % update transition probability
                    
                    % update transition probability for states not arrived in
                    Tp.BR.AL = Tp.BR.AL * (1-eta);
                    Tp.BR.AR = Tp.BR.AR * (1-eta);
                    Tp.BR.BL = Tp.BR.BL * (1-eta);
                    Tp.BR.BR = Tp.BR.BR * (1-eta);
                    
            end
            
        case 50
            
            % MF
            dv = r-Vcs.m; % compute prediction error
            VV (t) = Vcs.m; % store value in timeserie
            RPE (t) = dv;% store prediction error in timeserie
            
            % MB
            % most likly transition probability
            [T, idx] = max(structfun(@(x)max(x(:)),Tp.m));
            % what is the value of the most likely subsequent status?
            tmp  = fieldnames(Tp.m);
            name = char(tmp(idx));
            
            %V    = Vcs.(name);

            R     = 0;
            
            MB_V (t) = T * R;   % MB value is the transition probability multiplied by the value of the outcome found the s'
            UNC (t)  = 1 - T ; % 1 - max T
            
            % update MF
            Vcs.m = Vcs.m + alpha*dv;% update CSpR values
            
            % update MB based on specific US
            switch sr
                
                case 31
                    sdv      = 1 - Tp.m.AL;
                    SPE  (t) = sdv ;               % store state prediction error
                    Tp.m.AL = Tp.m.AL + eta * sdv; % update transition probability
                    
                    % update transition probability for states not arrived in
                    Tp.m.AR = Tp.m.AR * (1-eta);
                    Tp.m.BL = Tp.m.BL * (1-eta);
                    Tp.m.BR = Tp.m.BR * (1-eta);
                    Tp.m.m  = Tp.m.m  * (1-eta);
                    
                case 32
                    sdv      = 1 - Tp.m.AR;
                    SPE  (t) = sdv ;               % store state prediction error
                    Tp.m.AR = Tp.m.AR + eta * sdv; % update transition probability
                    
                    % update transition probability for states not arrived in
                    Tp.m.AL = Tp.m.AL * (1-eta);
                    Tp.m.BL = Tp.m.BL * (1-eta);
                    Tp.m.BR = Tp.m.BR * (1-eta);
                    Tp.m.m  = Tp.m.m  * (1-eta);
                    
                case 41
                    sdv      = 1 - Tp.m.BL;
                    SPE  (t) = sdv ;               % store state prediction error
                    Tp.m.BL = Tp.m.BL + eta * sdv; % update transition probability
                    
                    % update transition probability for states not arrived in
                    Tp.m.AR = Tp.m.AR * (1-eta);
                    Tp.m.AL = Tp.m.AL * (1-eta);
                    Tp.m.BR = Tp.m.BR * (1-eta);
                    Tp.m.m  = Tp.m.m  * (1-eta);
                    
                    
                case 42
                    sdv      = 1 - Tp.m.BR;
                    SPE  (t) = sdv ;               % store state prediction error
                    Tp.m.BR = Tp.m.BR + eta * sdv; % update transition probability
                    
                    % update transition probability for states not arrived in
                    Tp.m.AR = Tp.m.AR * (1-eta);
                    Tp.m.AL = Tp.m.AL * (1-eta);
                    Tp.m.BL = Tp.m.BL * (1-eta);
                    Tp.m.m  = Tp.m.m  * (1-eta);
                    
                case 0
                    sdv      = 1 - Tp.m.m;
                    SPE  (t) = sdv ;             % store state prediction error
                    Tp.m.m = Tp.m.m + eta * sdv; % update transition probability
                    
                    % update transition probability for states not arrived in
                    Tp.m.AR = Tp.m.AR * (1-eta);
                    Tp.m.AL = Tp.m.AL * (1-eta);
                    Tp.m.BL = Tp.m.BL * (1-eta);
                    Tp.m.BR = Tp.m.BR * (1-eta);
            end
            
            
    end
    
end

%---------------------------------------------------------------------
if displayPlot
    
    task_data.CSname(task_data.CSname >45) = 0;
    task_data.CSname(task_data.CSname >0) = 1;
    outcome(outcome ==0) = nan;
    outcome(outcome ==1) = 1.1;
    
    figure
    %suptitle (mdl_name)
    
    subplot(3,1,1)
    value = plot (VV, '-b');
    set(value, 'LineWidth', 2);
    hold
    mb_val = plot(MB_V,'-r');
    set(mb_val, 'LineWidth', 2);
    
    cs = plot (task_data.CSname,'-ok');
    Outcome = plot (outcome, 'r*');
    set(cs, 'LineWidth', 1);
    ylim ([-0.5 1.5]);
    xlabel ('trials');
    ylabel ('Expected Values');
    legend([value mb_val cs Outcome],'MF-val','MB-val','CS', 'Outcome');
    
    subplot(3,1,2)
    VD = plot (RPE, '-b');
    set(VD, 'LineWidth', 2);
    hold
    VD2 = plot(SPE,'-r');
    set(VD2, 'LineWidth', 2);
    cs = plot (task_data.CSname,'-ok');
    Outcome = plot (outcome, 'r*');
    set(cs, 'LineWidth', 1);
    ylim ([-1.5 1.5]);
    xlabel ('trials');
    ylabel ('Prediction Errors');
    legend([VD VD2 cs Outcome],'Reward Prediction Errors','State Prediction Error','CS', 'Outcome');
    
    subplot(3, 1, 3)
    VD = plot(UNC, '-b');
    set(VD, 'LineWidth', 2');
    hold
    cs = plot (task_data.CSname,'-ok');
    Outcome = plot (outcome, 'r*');
    set(cs, 'LineWidth', 1);
    ylim ([-0.5 1.5]);
    xlabel ('trials');
    ylabel ('Uncertainty');
    legend([VD cs Outcome],'Uncertainty','CS', 'Outcome');
    
end
%--------------------------------------------------------------------------
% input : free parameter alpha (learing rate for reward PE) and eta (learning
% rate for state PE), and information about the task (trial
% sequence: CS, outcome presentation order for the participant)
% output: timeseries of the values and the prediction order

%--------------------------------------------------------------------------
% last modified nov 2021

% get the specific outcome
sp_outcome = (task_data.US_identity*10) + task_data.US_side;

% define how the reward is tracked based on the model
task_data.US_identity(task_data.US_identity >0) = 1; % if there is a reward or not independenlty of the side
outcome = task_data.US_identity;


% initialise value and prediction error vectors to store the data.
VV   = nan (length(outcome),1);
RPE  = nan (length(outcome),1);
MB_V = nan (length(outcome),1);
SPE  = nan (length(outcome),1);
UNC  = nan (length(outcome),1);

for t = 1:length(outcome)
    
    CS       = task_data.CSname (t); %which CS is displayed
    r        = outcome (t); % reward
    sr       = sp_outcome(t); % specific code for the outcome
    
    switch CS
        
        case 31
            
            % MF
            dv = r-Vcs.AL; % compute prediction error
            VV (t) = Vcs.AL; % store value in timeserie
            RPE (t) = dv;    % store prediction error in timeserie
            
            % MB
            % most likly transition probability
            [T, idx] = max(structfun(@(x)max(x(:)),Tp.AL));
            % what is the value of the most likely subsequent status?
            tmp  = fieldnames(Tp.AL);
            name = char(tmp(idx));
            
            %V    = Vcs.(name);
            R     = 1;
            
            MB_V (t) = T * R;   % MB value is the transition probability multiplied by the value of the outcome found the s'
            
            UNC (t)  = 1 - T ; % 1 - max T
            
            % update MF
            Vcs.AL = Vcs.AL+ alpha * dv;% update CSpL values
            
            
            % update MB based on specific US
            switch sr
                
                case 31
                    sdv      = 1 - Tp.AL.AL;
                    SPE  (t) = sdv ;                 % store state prediction error
                    Tp.AL.AL = Tp.AL.AL + eta * sdv; % update transition probability
                    
                    % update transition probability for states not arrived in
                    Tp.AL.AR = Tp.AL.AR * (1-eta);
                    Tp.AL.BL = Tp.AL.BL * (1-eta);
                    Tp.AL.BR = Tp.AL.BR * (1-eta);
                    Tp.AL.m  = Tp.AL.m  * (1-eta);
                    
                case 32
                    sdv      = 1 - Tp.AL.AR;
                    SPE  (t) = sdv ;                 % store state prediction error
                    Tp.AL.AR = Tp.AL.AR + eta * sdv; % update transition probability
                    
                    % update transition probability for states not arrived in
                    Tp.AL.AL = Tp.AL.AL * (1-eta);
                    Tp.AL.BL = Tp.AL.BL * (1-eta);
                    Tp.AL.BR = Tp.AL.BR * (1-eta);
                    Tp.AL.m  = Tp.AL.m  * (1-eta);
                    
                case 41
                    sdv      = 1 - Tp.AL.BL;
                    SPE  (t) = sdv ;                 % store state prediction error
                    Tp.AL.BL = Tp.AL.BL + eta * sdv; % update transition probability
                    
                    % update transition probability for states not arrived in
                    Tp.AL.AL = Tp.AL.AL * (1-eta);
                    Tp.AL.AR = Tp.AL.AR * (1-eta);
                    Tp.AL.BR = Tp.AL.BR * (1-eta);
                    Tp.AL.m  = Tp.AL.m  * (1-eta);
                    
                case 42
                    sdv      = 1 - Tp.AL.BR;
                    SPE  (t) = sdv ;                 % store state prediction error
                    Tp.AL.BR = Tp.AL.BR + eta * sdv; % update transition probability
                    
                    % update transition probability for states not arrived in
                    Tp.AL.AL = Tp.AL.AL * (1-eta);
                    Tp.AL.AR = Tp.AL.AR * (1-eta);
                    Tp.AL.BL = Tp.AL.BL * (1-eta);
                    Tp.AL.m  = Tp.AL.m  * (1-eta);
                    
                case 0
                    sdv      = 1 - Tp.AL.m;
                    SPE  (t) = sdv ;                 % store state prediction error
                    Tp.AL.m = Tp.AL.m + eta * sdv;   % update transition probability
                    
                    % update transition probability for states not arrived in
                    Tp.AL.AL = Tp.AL.AL * (1-eta);
                    Tp.AL.AR = Tp.AL.AR * (1-eta);
                    Tp.AL.BL = Tp.AL.BL * (1-eta);
                    Tp.AL.BR = Tp.AL.BR * (1-eta);
                    
                    
            end
            
        case 32
            
            % MF
            dv = r-Vcs.AR;% compute prediction error
            VV (t) = Vcs.AR;% store value in timeserie
            RPE (t) = dv;% store prediction error in timeserie
            
            % MB
            % most likly transition probability
            [T, idx] = max(structfun(@(x)max(x(:)),Tp.AR));
            % what is the value of the most likely subsequent status?
            tmp  = fieldnames(Tp.AR);
            name = char(tmp(idx));
            
            %V    = Vcs.(name);
            R     = 1;
            
            MB_V (t) = T * R;   % MB value is the transition probability multiplied by the value of the outcome found the s'
            
            UNC (t)  = 1 - T ; % 1 - max T
            
            % update MF
            Vcs.AR = Vcs.AR + alpha*dv;% update CSpR values
            
            % update MB based on specific US
            switch sr
                
                case 31
                    sdv      = 1 - Tp.AR.AL;
                    SPE  (t) = sdv ;                 % store state prediction error
                    Tp.AR.AL = Tp.AR.AL + eta * sdv; % update transition probability
                    
                    % update transition probability for states not arrived in
                    Tp.AR.AR = Tp.AR.AR * (1-eta);
                    Tp.AR.BL = Tp.AR.BL * (1-eta);
                    Tp.AR.BR = Tp.AR.BR * (1-eta);
                    Tp.AR.m  = Tp.AR.m  * (1-eta);
                    
                case 32
                    sdv      = 1 - Tp.AR.AR;
                    SPE  (t) = sdv ;                 % store state prediction error
                    Tp.AR.AR = Tp.AR.AR + eta * sdv; % update transition probability
                    
                    % update transition probability for states not arrived in
                    Tp.AR.AL = Tp.AR.AL * (1-eta);
                    Tp.AR.BL = Tp.AR.BL * (1-eta);
                    Tp.AR.BR = Tp.AR.BR * (1-eta);
                    Tp.AR.m  = Tp.AR.m  * (1-eta);
                    
                case 41
                    sdv      = 1 - Tp.AR.BL;
                    SPE  (t) = sdv ;                 % store state prediction error
                    Tp.AR.BL = Tp.AR.BL + eta * sdv; % update transition probability
                    
                    % update transition probability for states not arrived in
                    Tp.AR.AL = Tp.AR.AL * (1-eta);
                    Tp.AR.AR = Tp.AR.AR * (1-eta);
                    Tp.AR.BR = Tp.AR.BR * (1-eta);
                    Tp.AR.m  = Tp.AR.m  * (1-eta);
                    
                case 42
                    sdv      = 1 - Tp.AR.BR;
                    SPE  (t) = sdv ;                 % store state prediction error
                    Tp.AR.BR = Tp.AR.BR + eta * sdv; % update transition probability
                    
                    % update transition probability for states not arrived in
                    Tp.AR.AL = Tp.AR.AL * (1-eta);
                    Tp.AR.AR = Tp.AR.AR * (1-eta);
                    Tp.AR.BL = Tp.AR.BL * (1-eta);
                    Tp.AR.m  = Tp.AR.m  * (1-eta);
                    
                case 0
                    sdv      = 1 - Tp.AR.m;
                    SPE  (t) = sdv ;                 % store state prediction error
                    Tp.AR.m = Tp.AR.m + eta * sdv; % update transition probability
                    
                    % update transition probability for states not arrived in
                    Tp.AR.AL = Tp.AR.AL * (1-eta);
                    Tp.AR.AR = Tp.AR.AR * (1-eta);
                    Tp.AR.BL = Tp.AR.BL * (1-eta);
                    Tp.AR.BR = Tp.AR.BR * (1-eta);
                    
            end
            
            
        case 41
            
            % MF
            
            dv = r-Vcs.BL;% compute prediction error
            VV (t) = Vcs.BL;% store value in timeserie
            RPE (t) = dv;% store prediction error in timeserie
            
            % MB
            % most likly transition probability
            [T, idx] = max(structfun(@(x)max(x(:)),Tp.BL));
            % what is the value of the most likely subsequent status?
            tmp  = fieldnames(Tp.BL);
            name = char(tmp(idx));
            
            %V    = Vcs.(name);
            R     = 1;
            
            MB_V (t) = T * R;   % MB value is the transition probability multiplied by the value of the outcome found the s'
            
            UNC (t)  = 1 - T ; % 1 - max T
            
            % update MF
            Vcs.BL = Vcs.BL + alpha*dv;% update CSpR values
            
            % update MB based on specific US
            switch sr
                
                case 31
                    sdv      = 1 - Tp.BL.AL;
                    SPE  (t) = sdv ;                 % store state prediction error
                    Tp.BL.AL = Tp.BL.AL + eta * sdv; % update transition probability
                    
                    % update transition probability for states not arrived in
                    Tp.BL.AR = Tp.BL.AR * (1-eta);
                    Tp.BL.BL = Tp.BL.BL * (1-eta);
                    Tp.BL.BR = Tp.BL.BR * (1-eta);
                    Tp.BL.m  = Tp.BL.m  * (1-eta);
                    
                case 32
                    sdv      = 1 - Tp.BL.AR;
                    SPE  (t) = sdv ;                 % store state prediction error
                    Tp.BL.AR = Tp.BL.AR + eta * sdv; % update transition probability
                    
                    % update transition probability for states not arrived in
                    Tp.BL.AL = Tp.BL.AL * (1-eta);
                    Tp.BL.BL = Tp.BL.BL * (1-eta);
                    Tp.BL.BR = Tp.BL.BR * (1-eta);
                    Tp.BL.m  = Tp.BL.m  * (1-eta);
                    
                case 41
                    sdv      = 1 - Tp.BL.BL;
                    SPE  (t) = sdv ;                 % store state prediction error
                    Tp.BL.BL = Tp.BL.BL + eta * sdv; % update transition probability
                    
                    % update transition probability for states not arrived in
                    Tp.BL.AL = Tp.BL.AL * (1-eta);
                    Tp.BL.AR = Tp.BL.AR * (1-eta);
                    Tp.BL.BR = Tp.BL.BR * (1-eta);
                    Tp.BL.m  = Tp.BL.m  * (1-eta);
                    
                case 42
                    sdv      = 1 - Tp.BL.BR;
                    SPE  (t) = sdv ;                 % store state prediction error
                    Tp.BL.BR = Tp.BL.BR + eta * sdv; % update transition probability
                    
                    % update transition probability for states not arrived in
                    Tp.BL.AL = Tp.BL.AL * (1-eta);
                    Tp.BL.AR = Tp.BL.AR * (1-eta);
                    Tp.BL.BL = Tp.BL.BL * (1-eta);
                    Tp.BL.m  = Tp.BL.m  * (1-eta);
                    
                case 0
                    sdv      = 1 - Tp.BL.m;
                    SPE  (t) = sdv ;                 % store state prediction error
                    Tp.BL.m = Tp.BL.m + eta * sdv; % update transition probability
                    
                    % update transition probability for states not arrived in
                    Tp.BL.AL = Tp.BL.AL * (1-eta);
                    Tp.BL.AR = Tp.BL.AR * (1-eta);
                    Tp.BL.BL = Tp.BL.BL * (1-eta);
                    Tp.BL.BR  = Tp.BL.BR  * (1-eta);
            end
            
        case 42
            
            % MF
            dv = r-Vcs.BR;% compute prediction error
            VV (t) = Vcs.BR;% store value in timeserie
            RPE (t) = dv;% store prediction error in timeserie
            
            % MB
            % most likly transition probability
            [T, idx] = max(structfun(@(x)max(x(:)),Tp.BR));
            % what is the value of the most likely subsequent status?
            tmp  = fieldnames(Tp.BR);
            name = char(tmp(idx));
            
            %V    = Vcs.(name);

            R     = 1;
            
            MB_V (t) = T * R;   % MB value is the transition probability multiplied by the value of the outcome found the s'
            UNC (t)  = 1 - T ; % 1 - max T
            
            % update MF
            Vcs.BR = Vcs.BR + alpha*dv;% update CSpR values
            
            % update MB based on specific US
            switch sr
                
                case 31
                    sdv      = 1 - Tp.BR.AL;
                    SPE  (t) = sdv ;                 % store state prediction error
                    Tp.BR.AL = Tp.BR.AL + eta * sdv; % update transition probability
                    
                    % update transition probability for states not arrived in
                    Tp.BR.AR = Tp.BR.AR * (1-eta);
                    Tp.BR.BL = Tp.BR.BL * (1-eta);
                    Tp.BR.BR = Tp.BR.BR * (1-eta);
                    Tp.BR.m  = Tp.BR.m  * (1-eta);
                    
                case 32
                    sdv      = 1 - Tp.BR.AR;
                    SPE  (t) = sdv ;                 % store state prediction error
                    Tp.BR.AR = Tp.BR.AR + eta * sdv; % update transition probability
                    
                    % update transition probability for states not arrived in
                    Tp.BR.AL = Tp.BR.AL * (1-eta);
                    Tp.BR.BL = Tp.BR.BL * (1-eta);
                    Tp.BR.BR = Tp.BR.BR * (1-eta);
                    Tp.BR.m  = Tp.BR.m  * (1-eta);
                    
                case 41
                    sdv      = 1 - Tp.BR.BL;
                    SPE  (t) = sdv ;                 % store state prediction error
                    Tp.BR.BL = Tp.BR.BL + eta * sdv; % update transition probability
                    
                    % update transition probability for states not arrived in
                    Tp.BR.AL = Tp.BR.AL * (1-eta);
                    Tp.BR.AR = Tp.BR.AR * (1-eta);
                    Tp.BR.BR = Tp.BR.BR * (1-eta);
                    Tp.BR.m  = Tp.BR.m  * (1-eta);
                    
                case 42
                    sdv      = 1 - Tp.BR.BR;
                    SPE  (t) = sdv ;                 % store state prediction error
                    Tp.BR.BR = Tp.BR.BR + eta * sdv; % update transition probability
                    
                    % update transition probability for states not arrived in
                    Tp.BR.AL = Tp.BR.AL * (1-eta);
                    Tp.BR.AR = Tp.BR.AR * (1-eta);
                    Tp.BR.BL = Tp.BR.BL * (1-eta);
                    Tp.BR.m  = Tp.BR.m  * (1-eta);
                    
                    
                case 0
                    sdv      = 1 - Tp.BR.m;
                    SPE  (t) = sdv ;                 % store state prediction error
                    Tp.BR.m = Tp.BR.m + eta * sdv; % update transition probability
                    
                    % update transition probability for states not arrived in
                    Tp.BR.AL = Tp.BR.AL * (1-eta);
                    Tp.BR.AR = Tp.BR.AR * (1-eta);
                    Tp.BR.BL = Tp.BR.BL * (1-eta);
                    Tp.BR.BR = Tp.BR.BR * (1-eta);
                    
            end
            
        case 50
            
            % MF
            dv = r-Vcs.m; % compute prediction error
            VV (t) = Vcs.m; % store value in timeserie
            RPE (t) = dv;% store prediction error in timeserie
            
            % MB
            % most likly transition probability
            [T, idx] = max(structfun(@(x)max(x(:)),Tp.m));
            % what is the value of the most likely subsequent status?
            tmp  = fieldnames(Tp.m);
            name = char(tmp(idx));
            
            %V    = Vcs.(name);

            R     = 0;
            
            MB_V (t) = T * R;   % MB value is the transition probability multiplied by the value of the outcome found the s'
            UNC (t)  = 1 - T ; % 1 - max T
            
            % update MF
            Vcs.m = Vcs.m + alpha*dv;% update CSpR values
            
            % update MB based on specific US
            switch sr
                
                case 31
                    sdv      = 1 - Tp.m.AL;
                    SPE  (t) = sdv ;               % store state prediction error
                    Tp.m.AL = Tp.m.AL + eta * sdv; % update transition probability
                    
                    % update transition probability for states not arrived in
                    Tp.m.AR = Tp.m.AR * (1-eta);
                    Tp.m.BL = Tp.m.BL * (1-eta);
                    Tp.m.BR = Tp.m.BR * (1-eta);
                    Tp.m.m  = Tp.m.m  * (1-eta);
                    
                case 32
                    sdv      = 1 - Tp.m.AR;
                    SPE  (t) = sdv ;               % store state prediction error
                    Tp.m.AR = Tp.m.AR + eta * sdv; % update transition probability
                    
                    % update transition probability for states not arrived in
                    Tp.m.AL = Tp.m.AL * (1-eta);
                    Tp.m.BL = Tp.m.BL * (1-eta);
                    Tp.m.BR = Tp.m.BR * (1-eta);
                    Tp.m.m  = Tp.m.m  * (1-eta);
                    
                case 41
                    sdv      = 1 - Tp.m.BL;
                    SPE  (t) = sdv ;               % store state prediction error
                    Tp.m.BL = Tp.m.BL + eta * sdv; % update transition probability
                    
                    % update transition probability for states not arrived in
                    Tp.m.AR = Tp.m.AR * (1-eta);
                    Tp.m.AL = Tp.m.AL * (1-eta);
                    Tp.m.BR = Tp.m.BR * (1-eta);
                    Tp.m.m  = Tp.m.m  * (1-eta);
                    
                    
                case 42
                    sdv      = 1 - Tp.m.BR;
                    SPE  (t) = sdv ;               % store state prediction error
                    Tp.m.BR = Tp.m.BR + eta * sdv; % update transition probability
                    
                    % update transition probability for states not arrived in
                    Tp.m.AR = Tp.m.AR * (1-eta);
                    Tp.m.AL = Tp.m.AL * (1-eta);
                    Tp.m.BL = Tp.m.BL * (1-eta);
                    Tp.m.m  = Tp.m.m  * (1-eta);
                    
                case 0
                    sdv      = 1 - Tp.m.m;
                    SPE  (t) = sdv ;             % store state prediction error
                    Tp.m.m = Tp.m.m + eta * sdv; % update transition probability
                    
                    % update transition probability for states not arrived in
                    Tp.m.AR = Tp.m.AR * (1-eta);
                    Tp.m.AL = Tp.m.AL * (1-eta);
                    Tp.m.BL = Tp.m.BL * (1-eta);
                    Tp.m.BR = Tp.m.BR * (1-eta);
            end
            
            
    end
    
end

%---------------------------------------------------------------------
if displayPlot
    
    task_data.CSname(task_data.CSname >45) = 0;
    task_data.CSname(task_data.CSname >0) = 1;
    outcome(outcome ==0) = nan;
    outcome(outcome ==1) = 1.1;
    
    figure
    %suptitle (mdl_name)
    
    subplot(3,1,1)
    value = plot (VV, '-b');
    set(value, 'LineWidth', 2);
    hold
    mb_val = plot(MB_V,'-r');
    set(mb_val, 'LineWidth', 2);
    
    cs = plot (task_data.CSname,'-ok');
    Outcome = plot (outcome, 'r*');
    set(cs, 'LineWidth', 1);
    ylim ([-0.5 1.5]);
    xlabel ('trials');
    ylabel ('Expected Values');
    legend([value mb_val cs Outcome],'MF-val','MB-val','CS', 'Outcome');
    
    subplot(3,1,2)
    VD = plot (RPE, '-b');
    set(VD, 'LineWidth', 2);
    hold
    VD2 = plot(SPE,'-r');
    set(VD2, 'LineWidth', 2);
    cs = plot (task_data.CSname,'-ok');
    Outcome = plot (outcome, 'r*');
    set(cs, 'LineWidth', 1);
    ylim ([-1.5 1.5]);
    xlabel ('trials');
    ylabel ('Prediction Errors');
    legend([VD VD2 cs Outcome],'Reward Prediction Errors','State Prediction Error','CS', 'Outcome');
    
    subplot(3, 1, 3)
    VD = plot(UNC, '-b');
    set(VD, 'LineWidth', 2');
    hold
    cs = plot (task_data.CSname,'-ok');
    Outcome = plot (outcome, 'r*');
    set(cs, 'LineWidth', 1);
    ylim ([-0.5 1.5]);
    xlabel ('trials');
    ylabel ('Uncertainty');
    legend([VD cs Outcome],'Uncertainty','CS', 'Outcome');
    
end