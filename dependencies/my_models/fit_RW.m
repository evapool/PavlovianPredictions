function [loglik]= fit_RW(x, task_data)

alpha = x(1);

% we fix eta since is irrelvant to the RW model
eta   = 0.5;

% init value
v0 = 0.5;
Vcs.AL = v0; Vcs.AR = v0; Vcs.BL = v0; Vcs.BR = v0; Vcs.m  = v0;


% init trans probabilities
Tp0  = 0.2;
Tp.AL.AL = Tp0; Tp.AL.BL = Tp0; Tp.AL.AR = Tp0; Tp.AL.BR = Tp0; Tp.AL.m = Tp0;
Tp.AR.AL = Tp0; Tp.AR.BL = Tp0; Tp.AR.AR = Tp0; Tp.AR.BR = Tp0; Tp.AR.m = Tp0;
Tp.BR.AL = Tp0; Tp.BR.BL = Tp0; Tp.BR.AR = Tp0; Tp.BR.BR = Tp0; Tp.BR.m = Tp0;
Tp.BL.AL = Tp0; Tp.BL.BL = Tp0; Tp.BL.AR = Tp0; Tp.BL.BR = Tp0; Tp.BL.m = Tp0;
Tp.m.AL  = Tp0; Tp.m.BL  = Tp0; Tp.m.AR  = Tp0; Tp.m.BR  = Tp0; Tp.m.m = Tp0;


% fit 
[VV, ~, ~, ~, ~, ~, ~] = simulate_RPE_SPE_MRI(alpha, eta, task_data, 0, Vcs, Tp, v0);
mdl = fitlm(VV,task_data.fitbehavior); 
loglik = mdl.LogLikelihood;