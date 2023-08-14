function [] = betas_mat2txt(ana_name, ana_dir, sub)

% create text file with colons: ID, & Betas
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


cd (ana_dir.roi)

%% ROI NAMES
roi_list = dir(['*' ana_name '.mat']);

for i = 1:length(roi_list)
    load(roi_list(i).name);
    
    Rdatabase(:,i,:) = result;
    
    % now we need to rename the name of the colons in a meaningful way
    roi_name = roi_list(i).name(1:end-(4 + length(ana_name) + 1 + length('betas') + 1));
    
    
    for ii = 1:size(result,2)
        
        con_name = result(1,ii);
        
        label = [char(con_name) '_in_' char(roi_name)];
        
        Rdatabase{1,i,ii} = label;
        
    end
 
end

% concatenate along the third dimention to get a 2d matrix
tmp = reshape(Rdatabase,size(Rdatabase,1),[])';
Rdatabase = tmp';

% add participant id
ID    = cell(length(sub.list)+1,1);
ID{1} = 'ID';
ID(2:end) = cellstr(sub.list);

fdatabase =[ID, Rdatabase];

%cd (out_dir)
% Convert cell to a table and use first row as variable names

% make title compatible
titles = strrep(fdatabase(1,:),'.','_');
titles = strrep(titles ,'-','_');

T = cell2table(fdatabase(2:end,:),'VariableNames',titles);
 
% Write the table to a CSV file
%writetable(T, ['extracted_betas_' ROI_name '.txt'],'Delimiter','\t');
writetable(T, ['extracted_betas_' ana_name '.txt'],'Delimiter','\t');
%writetable(T, ['extracted_betas_' con_name_orig '_via_' con_name '.txt'],'Delimiter','\t');
