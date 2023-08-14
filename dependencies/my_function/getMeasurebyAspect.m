function [measureByCondition] = getMeasurebyAspect (measure, data, aspect, devalued_outcome)

% last modified oct 2016

% aspect = 10 = identity (the identity is coded in decimal); 1 = side (side
% is coded in the unit)

% in case the experiment was stopped before the end
data.CSname = data.CSname (1:length(measure));


if aspect == 10 % in case we want the measure to be grouped by outcome identity
    
    if nargin == 4 % if a devaluated outcome is specified: 3 = sweet 4 = salty
        
        if devalued_outcome == 3
            
            measureByCondition.deval = measure(data.CSname==31 | data.CSname==32);
            measureByCondition.val = measure (data.CSname == 41 | data.CSname == 42 );
            measureByCondition.CSm = measure (data.CSname == 50);
            
        elseif devalued_outcome == 4
            
            measureByCondition.val = measure(data.CSname==31 | data.CSname==32);
            measureByCondition.deval = measure (data.CSname == 41 | data.CSname == 42);
            measureByCondition.CSm = measure (data.CSname == 50);
            
        end
        
    else
        
        measureByCondition.CSA = measure(data.CSname ==31 | data.CSname == 32);
        measureByCondition.CSB = measure (data.CSname == 41 | data.CSname == 42);
        measureByCondition.CSm = measure (data.CSname == 50);
        
    end
    
elseif aspect == 1 % in case we want the measure to be grouped by lateratlity
    
    measureByCondition.CSL = measure(data.CSname ==31 | data.CSname == 41);
    measureByCondition.CSR = measure (data.CSname == 32 | data.CSname == 42);
    measureByCondition.CSm = measure (data.CSname == 50);
    
end


end