function [t_FC, FC_rest, FC_movie] = mlm_FC(EFC_rest, EFC_movie, num_sub, roi_lh, roi_rh)
    %  Inputs:
    %    FC_rest   [360 x 360 x n_sub]  resting-state FC matrices
    %    FC_movie  [360 x 360 x n_sub]  movie-watching FC matrices
    %    n_sub     scalar               number of subjects
    %    roi_lh    scalar               left-hemisphere seed parcel index
    %    roi_rh    scalar               right-hemisphere seed parcel index
    %
    %  Outputs:
    %    MRD           [360 x 1]  t-statistic (movie > rest) for each parcel
    %    FC_rest_mean  [360 x 1]  group-average rest FC from seed
    %    FC_movie_mean [360 x 1]  group-average movie FC from seed    
    %%% SETTINGS
    addpath(genpath('C:\MICA\tools\BrainStat-0.4.2/brainstat_matlab/'))
    area_lh = setdiff(1:180, roi_lh); % Excludes roi_lh from left hemisphere range
    area_rh = setdiff(181:360, roi_rh); % Excludes roi_rh from right hemisphere range

    % Initialize storage matrices
    FC_rest_v1_all = zeros(360, num_sub);
    FC_movie_v1_all = zeros(360, num_sub);

    % Get FC in movie and rest for each subject
    for ii = 1:num_sub
        stream = EFC_rest(:,:,ii);
        FC_rest_v1_all(area_lh, ii) = stream(roi_lh, area_lh)';
        FC_rest_v1_all(area_rh, ii) = stream(roi_rh, area_rh)';

        stream = EFC_movie(:,:,ii);
        FC_movie_v1_all(area_lh, ii) = stream(roi_lh, area_lh)';
        FC_movie_v1_all(area_rh, ii) = stream(roi_rh, area_rh)';
    end

    % Concatenate the results for further analysis
    fc = [FC_rest_v1_all, FC_movie_v1_all]; % 360 * 2*num_sub

    %%% Group Averages for Outputs
    FC_rest = mean(FC_rest_v1_all, 2);
    FC_movie = mean(FC_movie_v1_all, 2);

    %%% STATISTICAL ANALYSIS
    % Generate terms for mixed linear model
    term_sub1 = repmat({'rest'}, num_sub, 1);
    term_sub2 = repmat({'movie'}, num_sub, 1);
    state = [term_sub1; term_sub2];

    sub1 = arrayfun(@(x) {['sub', num2str(x)]}, 1:num_sub)';
    subj = [sub1; sub1];

    % Define and fit mixed model using BrainStat
    term_state = FixedEffect(state);
    term_subject = MixedEffect(subj);
    model_mixed = term_state + term_subject;
    contrast_state = (state == "movie") - (state == "rest");

    slm_mixed = SLM(model_mixed, contrast_state, 'correction', 'fdr', 'two_tailed', true);

    % Process FC for t
    fc_358 = fc';
    fc_358(:, [roi_lh roi_rh]) = [];
    slm_mixed.fit(fc_358);
    t_FC = zeros(360, 1);
    t_FC([area_lh, area_rh]) = slm_mixed.t;

end