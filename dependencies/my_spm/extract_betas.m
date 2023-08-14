function [] = extract_betas (ana_name, ana_dir, con, sub)

% last modified on May 2019

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% | intialize spm 
spm('defaults','fmri');
spm_jobman('initcfg');


% | get the list of the ROI from which we want to extract the betas
roi_list = char(spm_select('FPList', ana_dir.roi, ['^'  '.*' '.nii ']));


% | go to the data directory
cd(ana_dir.data)


% | loop across ROIs first
for r=1:size(roi_list,1)
    result = con.names;
    maskName = roi_list(r,:);
    v_mask = spm_vol(maskName); %extract voxels from ROIs as mask
    outputFile = [v_mask.fname(1:end-4) '_betas_' ana_name '.mat'];
    
    if exist(outputFile,'file')==0 %extract betas only if file doesn't exist    
        Y = spm_read_vols(v_mask);
        roi_volume_mm = length(find(Y > 0))*abs(det(v_mask.mat));
        clear Y;
        
        %loop across contrasts
        for c = 1:length(con.list)
            conName = con.names(c);
            % List of files to extract data from
            for s0 = 1:length(sub.list)  % select con image from every subject
                conDir = fullfile( ana_dir.data, ['sub-' sub.list(s0,:) '_' char(con.list(c))]); %add Model directory if necessary
                images(s0,:) = conDir;
            end
            numImages = size(images,1);

            fprintf(1,[char(conName) ': Memory mapping images...\n']);
            for i=1:numImages
                v{i} = spm_vol(images(i,:));
                % Verify images are in same space
                if ~isequal(v{i}.dim(1:3),v{1}.dim(1:3))
                error('Images must have same dimensions.')
                end
                % Verify orientation/position are the same
                if ~isequal(v{i}.mat,v{1}.mat)
                error('Images must have same orientation/position.')
                end
            end
            fprintf(1,'done.\n');

            [Y, XYZmm] = spm_read_vols(v{1});
            clear Y

            XYZmask = inv(v_mask.mat)*([XYZmm; ones(1,size(XYZmm,2))]);
            ind = find(spm_sample_vol(v_mask, XYZmask(1,:), XYZmask(2,:), XYZmask(3,:),0) > 0);

            ind = ind(:)';
            vals = [];
            for j=1:length(v)
                Y = spm_read_vols(v{j});
                % rows:  images.  cols:  voxels
                vals = [vals; Y(ind)];
            end

            % Consider only values which are finite for all images
            vals = vals(:, find(all(isfinite(vals))));
            intersection_volume_mm = size(vals,2)*abs(det(v{1}.mat));
            if intersection_volume_mm==0
                error('No voxels in ROI are in-brain for all images to be sampled.');
            end
            m = mean(vals, 2);

            result(2:length(v)+1,c) = num2cell(m); 

            save(outputFile,'result');
           

        end
    else
        disp('********************************************************* ATTENTION: It is likely that files with betas already exist, so they are not extracted.')
 
    end
         
end

end