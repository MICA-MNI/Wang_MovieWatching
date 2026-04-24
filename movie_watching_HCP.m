%% ========================================================================
%  MAIN ANALYSIS SCRIPT
%  Movie vs. Rest Functional Connectivity (MRD) Analysis
%  
%  REQUIREMENTS (add to MATLAB path before running):
%    - BrainSpace toolbox:  https://github.com/MICA-MNI/BrainSpace
%    - gifti toolbox:       https://www.artefact.tk/software/matlab/gifti/
%    - spin_permutations:   included in BrainSpace or OSF data package
%    - plotSurfaceROIBoundary: https://github.com/StuartJO/plotSurfaceROIBoundary
%    - MATLAB Statistics and Machine Learning Toolbox
%
%  DATA: Download from [OSF/Zenodo link] and place in the `data/` folder.
%
%  USAGE:
%    1. Set `base_path` below to the folder containing this script and `data/`
%    2. Run each section (%%)'s sequentially, or run all at once
%
%  AUTHORS: [Your names]
%  DATE:    [Year]
%  VERSION: 1.0
% ========================================================================
 
%% ---- 0. CONFIGURATION: SET YOUR PATH HERE (only edit this section) ----
clear; clc;
% >>> EDIT THIS LINE ONLY <<<
base_path = '/path/to/your/project/folder';   % e.g. '/home/user/MRD_project'
 
% Derived paths (do not edit below)
data_path   = fullfile(base_path, 'data');
cmap_path   = fullfile(base_path, 'colormap');
output_path = fullfile(base_path, 'outputs');
if ~exist(output_path, 'dir'); mkdir(output_path); end
 
% Add required toolboxes (edit paths as needed, or use `addpath` in startup.m)
% addpath(genpath('/path/to/BrainSpace'));
% addpath(genpath('/path/to/gifti'));
% addpath(genpath('/path/to/plotSurfaceROIBoundary'));
 
fprintf('Base path set to: %s\n', base_path);
fprintf('Output path: %s\n', output_path);
 
 
%% ---- 1. SETTINGS & DATA LOADING ----------------------------------------
num_parc = 360;   % number of brain regions (Glasser-360 parcellation)
num_sub  = 93;    % number of subjects
 
% --- Subject list ---
sub_list_file = fullfile(data_path, 'sublist_hcp.csv');
assert(exist(sub_list_file, 'file') == 2, ...
    'Subject list not found: %s', sub_list_file);
sub_list = readtable(sub_list_file);
sub_list = char(string(sub_list.Var1));
 
% --- Colormaps ---
cmap_names = {'buda','cork','tokyo','nuuk','vik','batlow','imola',...
              'bamako','hawaii','lajolla','roma'};
for c = 1:numel(cmap_names)
    cmap_file = fullfile(cmap_path, [cmap_names{c} '.mat']);
    assert(exist(cmap_file,'file')==2, 'Colormap file missing: %s', cmap_file);
    load(cmap_file, cmap_names{c});
end
 
% --- Network labels ---
load(fullfile(data_path, 'label_12n.mat'),                  'label_12n');
load(fullfile(data_path, 'label_mesulam_glasser360.mat'),   'label_mesulam_glasser');
load(fullfile(data_path, 'color12.mat'),                    'color12');
 
% --- Surface & parcellation ---
surf_file = fullfile(data_path, 'fsaverage.midthickness_mni_32k_fs_LR.mat');
assert(exist(surf_file,'file')==2, 'Surface file missing: %s', surf_file);
load(surf_file, 'G');
c69 = G; clear G;
[surf_l32k, surf_r32k] = split_surfaces(c69);
 
parc_glasser = readmatrix(fullfile(data_path, 'glasser-360_conte69.csv'));
% Remove duplicate index of parcel 180 on right sphere
parc_180 = find(parc_glasser == 180);
parc_glasser(parc_180(59:121)) = 0;
 
% --- Flat surfaces ---
glasser_flat_lh = gifti(fullfile(data_path, 'S1200.L.flat.32k_fs_LR.surf.gii'));
glasser_flat_rh = gifti(fullfile(data_path, 'S1200.R.flat.32k_fs_LR.surf.gii'));
 
glasser_flat_boundry.vertices  = [glasser_flat_lh.vertices; glasser_flat_rh.vertices];
glasser_flat_boundry.faces     = [glasser_flat_lh.faces;   glasser_flat_rh.faces + 32492];
glasser_flat_boundry_lh.vertices = glasser_flat_lh.vertices;
glasser_flat_boundry_lh.faces    = glasser_flat_lh.faces;
 
% --- Functional connectivity matrices (93 subjects) ---
% Shape: [360 x 360 x 93]  (full session average)
% Shape: [360 x 360 x 93 x 4]  (per-session)
fc_files = {'FC_rest_93.mat','FC_movie_93.mat','FC_rest_93_4.mat','FC_movie_93_4.mat'};
fc_vars  = {'FC_rest_93',    'FC_movie_93',    'FC_rest_93_4',    'FC_movie_93_4'};
for f = 1:numel(fc_files)
    fc_file = fullfile(data_path, fc_files{f});
    assert(exist(fc_file,'file')==2, 'FC file missing: %s', fc_file);
    load(fc_file, fc_vars{f});
end
 
fprintf('All data loaded successfully.\n');
 
 
%% ---- 2. FIGURE 1A: Visualise per-session and mean FC matrices ----------
%
%  Input:  FC_rest_93_4, FC_movie_93_4  [360 x 360 x 93 x 4]
%          label_12n  [360 x 1]  — network label for each parcel
%  Output: figures (optionally saved to output_path)
%
%  Modify `id_sub` to plot a different subject (1 … 93).
 
id_sub = 92;   % <<< change subject index here (1-indexed)
 
% Sort parcels by network label for display
[~, sorted_idx] = sort(label_12n);
 
for ii = 1:4
    restFC  = FC_rest_93_4(:,:,id_sub,ii);
    movieFC = FC_movie_93_4(:,:,id_sub,ii);
 
    fig = figure('Color','w');
    imagesc(restFC(sorted_idx,sorted_idx)); axis image;
    colorbar('ticklength',0); colormap(imola);
    title(sprintf('Rest FC — Subject %d, Session %d', id_sub, ii));
    % saveas(fig, fullfile(output_path, sprintf('restFC_sub%d_ses%d.png', id_sub, ii)));
 
    fig = figure('Color','w');
    imagesc(movieFC(sorted_idx,sorted_idx)); axis image;
    colorbar('ticklength',0); colormap(imola);
    title(sprintf('Movie FC — Subject %d, Session %d', id_sub, ii));
    % saveas(fig, fullfile(output_path, sprintf('movieFC_sub%d_ses%d.png', id_sub, ii)));
end
 
% Mean FC across 4 sessions
rest_FC_subj  = mean(FC_rest_93_4(:,:,id_sub,:), 4);
movie_FC_subj = mean(FC_movie_93_4(:,:,id_sub,:), 4);
 
fig = figure('Color','w');
imagesc(rest_FC_subj(sorted_idx,sorted_idx)); axis image;
colorbar('ticklength',0); colormap(imola);
title(sprintf('Mean Rest FC — Subject %d', id_sub));
% saveas(fig, fullfile(output_path, sprintf('meanRestFC_sub%d.png', id_sub)));
 
fig = figure('Color','w');
imagesc(movie_FC_subj(sorted_idx,sorted_idx)); axis image;
colorbar('ticklength',0); colormap(imola);
title(sprintf('Mean Movie FC — Subject %d', id_sub));
% saveas(fig, fullfile(output_path, sprintf('meanMovieFC_sub%d.png', id_sub)));
 
% Network colour bar
label_sorted = label_12n(sorted_idx);
label_colors = color12(label_sorted, :);   % 360 x 3
figure('Color','w');
image(permute(label_colors, [1 3 2]));
set(gca, 'YTick',[], 'XTick',[]);
title('Network colour bar (sorted)');
 
 
%% ---- 3. IDENTIFY REPRESENTATIVE ROI PER NETWORK -----------------------
%
%  For each of the 12 networks, the parcel whose resting-state FC profile
%  has the smallest mean cosine distance to all other parcels in that
%  network is chosen as the representative ROI.
%
%  Output: roi_12network  [12 x 1]  — parcel index of representative ROI
%          Brain surface maps saved to figures
 
num_roi      = 12;
restFC_group = mean(FC_rest_93, 3);    % [360 x 360] group-average rest FC
movieFC_group= mean(FC_movie_93, 3);   % [360 x 360] group-average movie FC
 
roi_12network = zeros(12, 1);
for i = 1:12
    idx_temp    = find(label_12n == i);
    FC_movie_sub= movieFC_group(idx_temp, idx_temp);
    FC_rest_sub = restFC_group(idx_temp, idx_temp);
 
    zmovie = squareform(pdist(FC_movie_sub', 'cosine'));
    zrest  = squareform(pdist(FC_rest_sub',  'cosine'));
 
    csdist = mean(zrest, 1);
    roi_12network(i) = idx_temp(csdist == min(csdist));
end
 
% --- Visualise 12 seed ROIs (manually defined for figure labels) ---
% NOTE: These indices are fixed for the Glasser-360 / HCP-MMP parcellation.
%       Left-hemisphere parcel i corresponds to right-hemisphere parcel i+180.
roi_12 = zeros(360, num_roi);
seed_lh = [1; 20; 51; 59; 50; 25; 83; 24; 35; 141; 122; 110];
seed_rh = seed_lh + 180;
for i = 1:num_roi
    roi_12(seed_lh(i), i) = 1;
    roi_12(seed_rh(i), i) = 1;
end
 
roi_labels_1to5  = {'V1','LO1 V2','S1 SMN','a24pr CIN','MIP dAN'};
roi_labels_6to10 = {'PSL LANG','p9-46v FPN','A1','31pv DMN','TPOJ3 PMM'};
roi_labels_11to12= {'PeEc VMM','Pir (LMN)'};
 
B = plot_hemispheres(roi_12(:,1:5), {surf_l32k,surf_r32k}, ...
    'parcellation', parc_glasser, 'labeltext', roi_labels_1to5);
colormap([0.8 0.8 0.8; vik]);
 
B = plot_hemispheres(roi_12(:,6:10), {surf_l32k,surf_r32k}, ...
    'parcellation', parc_glasser, 'labeltext', roi_labels_6to10);
colormap([0.8 0.8 0.8; vik]);
 
B = plot_hemispheres(roi_12(:,11:12), {surf_l32k,surf_r32k}, ...
    'parcellation', parc_glasser, 'labeltext', roi_labels_11to12);
colormap([0.8 0.8 0.8; vik]);
 
for ii = 1:num_roi
    figure;
    plotSurfaceROIBoundary(glasser_flat_boundry_lh, parc_glasser(1:32492), ...
        roi_12(1:180,ii), 'faces', [0.8 0.8 0.8; vik], [0, 1]);
    title(sprintf('ROI %d flat map', ii));
end
 
 
%% ---- 4. FIGURES 1B & 2: MRD and FC for 12 seed ROIs ------------------
%
%  For each seed ROI (bilateral), compute:
%    - Mean Rest FC  profile across subjects
%    - Mean Movie FC profile across subjects
%    - MRD (Movie–Rest Difference t-statistic) via mlm_FC()
%    - Spearman correlation between rest and movie FC across parcels
%    - Spin-test p-value
%
%  Input:  FC_rest_93, FC_movie_93  [360 x 360 x 93]
%  Output: MRD_12roi, FC_rest, FC_movie  [360 x 12]
%          rho/p/pspin  [12 x 1]
%          scatter plots saved to output_path
 
num_roi   = 12;
idx_12roi = [1,181; 20,200; 51,231; 59,239; 50,230; 25,205; ...
             83,263; 24,204; 35,215; 141,321; 122,302; 110,290];
 
MRD_12roi        = zeros(num_parc, num_roi);
FC_rest          = zeros(num_parc, num_roi);
FC_movie         = zeros(num_parc, num_roi);
rho_rest_movie   = zeros(num_roi, 1);
p_rest_movie     = zeros(num_roi, 1);
pspin_rest_movie = zeros(num_roi, 1);
 
[sphere_lh, sphere_rh] = load_conte69('spheres');
 
for ii = 1:num_roi
    % --- Compute MRD and group-average FC ---
    [MRD_12roi(:,ii), FC_rest(:,ii), FC_movie(:,ii)] = ...
        mlm_FC(FC_rest_93, FC_movie_93, num_sub, idx_12roi(ii,1), idx_12roi(ii,2));
 
    % --- Scatter: rest vs movie FC ---
    fig = figure('Color','w');
    scatter(FC_rest(:,ii), FC_movie(:,ii), 30, ...
        'MarkerEdgeColor', color12(ii,:), ...
        'MarkerFaceColor', color12(ii,:), ...
        'MarkerFaceAlpha', 0.6, 'LineWidth', 1);
    minVal = min([FC_rest(:,ii); FC_movie(:,ii)]) - 0.05;
    maxVal = max([FC_rest(:,ii); FC_movie(:,ii)]) + 0.05;
    xlim([minVal, maxVal]); ylim([minVal, maxVal]);
    xlabel('Rest FC'); ylabel('Movie FC');
    title(sprintf('ROI %d: rest vs movie FC', ii));
    set(gca, 'FontSize', 20); grid off;
    % saveas(fig, fullfile(output_path, sprintf('scatter_restVSmovieFC_roi%d.png', ii)));
 
    % --- Spearman correlation ---
    [rho_rest_movie(ii), p_rest_movie(ii)] = ...
        corr(FC_rest(:,ii), FC_movie(:,ii), 'type', 'spearman');
 
    % --- Spin test ---
    data_c69    = mica_parcelData2surfData([0; FC_rest(:,ii)],  c69, parc_glasser)';
    FCmovie_c69 = mica_parcelData2surfData([0; FC_movie(:,ii)], c69, parc_glasser)';
 
    n_permutations = 1000;
    y_rand = spin_permutations( ...
        {data_c69(1:32492), data_c69(32493:64984)}, ...
        {sphere_lh, sphere_rh}, n_permutations, 'random_state', 0);
 
    FCrest_rotated   = squeeze([y_rand{1}(:,1,:); y_rand{2}(:,1,:)]);
    r_rand           = corr(FCmovie_c69, FCrest_rotated, ...
                            'rows','pairwise','type','spearman');
    pspin_rest_movie(ii) = mean(rho_rest_movie(ii) < r_rand);
end
 
% Display summary table
fprintf('\n--- Rest vs Movie FC correlation (12 ROIs) ---\n');
fprintf('%-6s %-8s %-8s %-10s\n','ROI','rho','p','p_spin');
for ii = 1:num_roi
    fprintf('%-6d %-8.3f %-8.4f %-10.4f\n', ...
        ii, rho_rest_movie(ii), p_rest_movie(ii), pspin_rest_movie(ii));
end
 
 
%% ---- 5. FIGURE 2A: Subject-level FC maps and Mesulam boxplot ----------
%
%  (a) Show rest/movie FC maps for 3 example subjects
%  (b) Boxplot of MRD across Mesulam hierarchy levels (3 seed ROIs)
%
%  NOTE: Subject indices below are 1-indexed. Subjects 61–63 are used as
%        examples; change as desired. Both rest and movie use the SAME
%        subject set for a fair comparison.
 
% --- (a) Example subject FC maps ---
roi_lh    = 24;           % seed parcel (A1, left hemisphere)
subj_ids  = [61, 62, 63]; % <<< change subject indices here (1-indexed, max 93)
 
temp_rest  = FC_rest_93(:, roi_lh, subj_ids);   % [360 x 1 x 3] → [360 x 3]
temp_rest  = squeeze(temp_rest);
temp_rest(181:360, :) = 0;                       % show left hemisphere only
temp_rest(roi_lh, :)  = -Inf;                    % mask seed parcel
 
A = plot_hemispheres(temp_rest, {surf_l32k,surf_r32k}, ...
    'parcellation', parc_glasser, ...
    'labeltext', arrayfun(@(s) sprintf('restFC sub%d',s), subj_ids, 'UniformOutput',false));
A.colormaps([0.8 0.8 0.8; imola]);
A.colorlimits([-0.062, 1.156]);
 
temp_movie = FC_movie_93(:, roi_lh, subj_ids);
temp_movie = squeeze(temp_movie);
temp_movie(181:360, :) = 0;
temp_movie(roi_lh, :)  = -Inf;
 
A = plot_hemispheres(temp_movie, {surf_l32k,surf_r32k}, ...
    'parcellation', parc_glasser, ...
    'labeltext', arrayfun(@(s) sprintf('movieFC sub%d',s), subj_ids, 'UniformOutput',false));
A.colormaps([0.8 0.8 0.8; imola]);
A.colorlimits([-0.062, 1.156]);
 
% --- (b) MRD boxplot across Mesulam levels ---
% Seeds: V1 (roi 1), MIP (roi 5), 31pv (roi 9)
roi_names_3   = {'V1', 'MIP', '31pv'};
roi_idx_3     = [1, 5, 9];
roi_colors_3  = color12(roi_idx_3, :);
MRD_3roi      = MRD_12roi(:, roi_idx_3);   % [360 x 3]
 
% Mesulam level strings and plotting order
mesulam_str = strings(size(label_mesulam_glasser));
mesulam_str(label_mesulam_glasser == 4) = "Idio";
mesulam_str(label_mesulam_glasser == 3) = "Uni";
mesulam_str(label_mesulam_glasser == 2) = "Hetero";
mesulam_str(label_mesulam_glasser == 1) = "Paralim";
hierarchy_order = ["Idio","Uni","Hetero","Paralim"];
 
num_parcel = size(MRD_3roi, 1);
data_long        = [];
group_hierarchy  = [];
group_roi        = [];
 
for i = 1:3
    data_long       = [data_long;       MRD_3roi(:, i)];
    group_hierarchy = [group_hierarchy; mesulam_str];
    group_roi       = [group_roi;       repmat(roi_names_3(i), num_parcel, 1)];
end
 
group_hierarchy = categorical(group_hierarchy, hierarchy_order, 'Ordinal', true);
group_roi       = categorical(group_roi, roi_names_3, 'Ordinal', true);
 
figure('Color','w','Position',[100 100 670 500]);
hold on;
boxplot(data_long, {group_hierarchy, group_roi}, ...
    'FactorSeparator',1, 'Colors','k', 'Symbol','o', 'Widths',0.6, 'Labels',[]);
 
h = findobj(gca, 'Tag','Box');
for j = 1:length(h)
    roi_idx = mod(length(h) - j, 3) + 1;
    patch(get(h(j),'XData'), get(h(j),'YData'), roi_colors_3(roi_idx,:), ...
        'FaceAlpha',0.5, 'EdgeColor','k');
end
 
xticks(1.5:3:12); xticklabels(hierarchy_order);
ylabel('MRD','FontSize',18); set(gca,'FontSize',16);
 
legend_colors_leg = 0.7 * roi_colors_3 + 0.3;
legend_handles    = gobjects(3, 1);
for i = 1:3
    legend_handles(i) = plot(NaN, NaN, 's', ...
        'MarkerFaceColor', legend_colors_leg(i,:), ...
        'MarkerEdgeColor', legend_colors_leg(i,:), 'MarkerSize',10);
end
legend(legend_handles, roi_names_3, 'Location','northeastoutside','FontSize',14);
 
% Mesulam level brain maps
mesulam_view = double([label_mesulam_glasser==4, label_mesulam_glasser==3, ...
                       label_mesulam_glasser==2, label_mesulam_glasser==1]);
A = plot_hemispheres(mesulam_view, {surf_l32k,surf_r32k}, ...
    'parcellation', parc_glasser, 'labeltext',{'Idio','Uni','Hetero','Paralim'});
A.colormaps([1 1 1; 0.3 0.3 0.3]);
 
 
%% ---- 6. FIGURES 3A & 3B: GD and MRD–GD association (12 seed ROIs) ----
%
%  Computes partial correlations between geodesic distance (GD) and MRD
%  controlling for activation difference (movie minus rest mean FC).
%
%  Input:  gd_glasser_7subjs.mat  — parcel-level GD [180 x 180 x 7]
%          HCP_MPC_15subj.mat     — MPC matrices [360 x 360 x 15]
%  Output: scatter plots (partial effect of GD on MRD)
%          slopes  [12 x 1]
 
num_roi   = 12;
num_parc  = 360;
idx_12roi = [1,181; 20,200; 51,231; 59,239; 50,230; 25,205; ...
             83,263; 24,204; 35,215; 141,321; 122,302; 110,290];
 
load(fullfile(data_path,'gd_glasser_7subjs.mat'),  'parcel_gd');
load(fullfile(data_path,'HCP_MPC_15subj.mat'),     'MPC_360_all');
 
restFC_group  = mean(FC_rest_93,  3);
movieFC_group = mean(FC_movie_93, 3);
gd_group      = mean(parcel_gd,   3);   % [180 x 180] group-average GD
mpc_group     = mean(MPC_360_all, 3);
 
actvt_rest  = mean(restFC_group,  2);
actvt_movie = mean(movieFC_group, 2);
 
num_parcel = 180;
 
% Pre-allocate
gd_lh_12roi  = zeros(num_parcel-1, num_roi);
mpc_lh_12roi = zeros(num_parcel-1, num_roi);
cov_lh_12roi = zeros(num_parcel-1, num_roi);
mrd_lh_12roi = zeros(num_parcel-1, num_roi);
mrd_12roi    = zeros(num_parc,     num_roi);
 
for ii = 1:num_roi
    % MRD (full cortex)
    [mrd_12roi(:,ii), ~, ~] = mlm_FC(FC_rest_93, FC_movie_93, num_sub, ...
                                      idx_12roi(ii,1), idx_12roi(ii,2));
    % Left-hemisphere MRD without seed parcel
    temp_lh = mrd_12roi(1:num_parcel, ii);
    temp_lh(idx_12roi(ii,1)) = [];
    mrd_lh_12roi(:,ii) = temp_lh;
 
    % GD (left hemisphere, remove seed)
    temp_gd = gd_group(1:num_parcel, idx_12roi(ii,1));
    temp_gd(idx_12roi(ii,1)) = [];
    gd_lh_12roi(:,ii) = temp_gd;
 
    % MPC (left hemisphere, remove seed)
    temp_mpc = mpc_group(1:num_parcel, idx_12roi(ii,1));
    temp_mpc(idx_12roi(ii,1)) = [];
    mpc_lh_12roi(:,ii) = temp_mpc;
 
    % Covariate: activation difference (left hemisphere, remove seed)
    cov_FC = actvt_movie(1:num_parcel) - actvt_rest(1:num_parcel);
    cov_FC(idx_12roi(ii,1)) = [];
    cov_lh_12roi(:,ii) = cov_FC;
end
 
% --- Partial correlation scatter plots: GD vs MRD ---
roi_names_12 = {'V1','LO1','S1','a24pr','MIP','PSL', ...
                'p9-46v','A1','31pv','TPOJ3','PeEc','Pir'};
slopes_12roi = zeros(num_roi, 1);
 
for i = 1:num_roi
    x   = gd_lh_12roi(:,i);
    y   = mrd_lh_12roi(:,i);
    cov = cov_lh_12roi(:,i);
 
    % Partial out covariate
    [~,~,x_resid] = regress(x, [ones(size(cov)), cov]);
    [~,~,y_resid] = regress(y, [ones(size(cov)), cov]);
 
    p_fit         = polyfit(x_resid, y_resid, 1);
    y_fit         = polyval(p_fit, x_resid);
    slopes_12roi(i) = p_fit(1);
 
    fig = figure('Color','w');
    scatter(x_resid, y_resid, 30, ...
        'MarkerEdgeColor', color12(i,:), ...
        'MarkerFaceColor', color12(i,:), 'MarkerFaceAlpha',0.6);
    hold on;
    plot(x_resid, y_fit, '-', 'Color', color12(i,:), 'LineWidth',2);
    buffer = 0.05;
    xlim([min(x_resid)-buffer*range(x_resid), max(x_resid)+buffer*range(x_resid)]);
    ylim([min(y_resid)-buffer*range(y_resid), max(y_resid)+buffer*range(y_resid)]);
    xlabel('GD (residual)','FontSize',16);
    ylabel('MRD (residual)','FontSize',16);
    title(sprintf('Partial: MRD vs GD — %s', roi_names_12{i}),'FontSize',14);
    set(gca,'FontSize',20); grid off;
    % saveas(fig, fullfile(output_path, sprintf('partial_MRD_GD_%s.tiff', roi_names_12{i})));
end
 
 
%% ---- 7. FIGURE 3C: Slopes for all 180 LH parcels ----------------------
%
%  For every left-hemisphere parcel, regress MRD on GD (+ covariate) and
%  store the GD slope.  Plot results as boxplots per network and per
%  Mesulam level, and on the brain surface.
%
%  Output: slopes_180  [180 x 1]
 
num_parcel   = 180;
slopes_180   = zeros(num_parcel, 1);
gd_lh_180    = zeros(num_parcel-1, num_parcel);
mrd_lh_180   = zeros(num_parcel-1, num_parcel);
cov_lh_180   = zeros(num_parcel-1, num_parcel);
 
actvt_rest  = mean(mean(FC_rest_93,  2), 3);   % recompute if not in workspace
actvt_movie = mean(mean(FC_movie_93, 2), 3);
 
for ii = 1:num_parcel
    % GD
    temp_gd = gd_group(1:num_parcel, ii); temp_gd(ii) = [];
    gd_lh_180(:,ii) = temp_gd;
 
    % MRD
    [temp_region, ~, ~] = mlm_FC(FC_rest_93, FC_movie_93, num_sub, ii, ii+180);
    temp_region = temp_region(1:num_parcel); temp_region(ii) = [];
    mrd_lh_180(:,ii) = temp_region;
 
    % Covariate
    cov_FC = actvt_movie(1:num_parcel) - actvt_rest(1:num_parcel);
    cov_FC(ii) = [];
    cov_lh_180(:,ii) = cov_FC;
 
    % OLS: MRD ~ intercept + GD + covariate
    X = [ones(num_parcel-1,1), gd_lh_180(:,ii), cov_lh_180(:,ii)];
    b = X \ mrd_lh_180(:,ii);
    slopes_180(ii) = b(2);   % slope for GD
end
 
% --- Boxplot by network ---
num_networks   = 12;
label_12n_lh   = label_12n(1:180);

tmp_data = arrayfun(@(n) slopes_180(label_12n_lh==n), 1:num_networks, 'UniformOutput', false);
grouped_data_n = vertcat(tmp_data{:});

tmp_idx = arrayfun(@(n) n*ones(sum(label_12n_lh==n),1), 1:num_networks, 'UniformOutput', false);
group_idx_n = vertcat(tmp_idx{:});
 
figure('Color','w'); hold on;
boxplot(grouped_data_n, group_idx_n, 'Colors','k', 'Symbol','o', 'Widths',0.5);
h = findobj(gca,'Tag','Box');
for j = 1:num_networks
    patch(get(h(num_networks-j+1),'XData'), get(h(num_networks-j+1),'YData'), ...
          color12(j,:), 'FaceAlpha',0.6, 'EdgeColor','k');
end
set(gca,'XTick',1:num_networks,'XTickLabel',1:num_networks,'FontSize',14);
ylabel('Slope (GD on MRD)','FontSize',16); xlabel('Network','FontSize',16);
ylim padded; box off; grid off;
 
% --- Boxplot by Mesulam level ---
mesulam_labels = label_mesulam_glasser(1:180);
mesulam_order  = [4, 3, 2, 1];   % Idio, Uni, Hetero, Paralim
data_mes       = arrayfun(@(o) slopes_180(mesulam_labels==o), ...
                           mesulam_order, 'UniformOutput',false);
grouped_mes = vertcat(data_mes{:});   % data_mes 本身已是 cell，这行本身没问题

tmp_mes = arrayfun(@(n) n*ones(length(data_mes{n}),1), 1:4, 'UniformOutput', false);
group_idx_mes = vertcat(tmp_mes{:});

figure('Color','w'); hold on;
boxplot(grouped_mes, group_idx_mes, 'Colors','k', 'Symbol','o', 'Widths',0.5);
h = findobj(gca,'Tag','Box');
for j = 1:4
    patch(get(h(4-j+1),'XData'), get(h(4-j+1),'YData'), ...
          [0.4 0.4 0.4], 'FaceAlpha',0.6, 'EdgeColor','k');
end
set(gca,'XTick',1:4,'XTickLabel',{'Idio','Uni','Hetero','Paralim'},'FontSize',14);
ylabel('Slope (GD on MRD)','FontSize',16);
xlabel('Mesulam Hierarchy Level','FontSize',16);
ylim padded; box off; grid off;
 
% --- Brain surface map of slopes ---
B = plot_hemispheres([slopes_180; slopes_180], {surf_l32k,surf_r32k}, ...
    'parcellation',parc_glasser,'labeltext',{'slope'});
colormap([0.8 0.8 0.8; autumn]); B.colorlimits([-0.099, 0.046]);
 
figure; plotSurfaceROIBoundary(glasser_flat_boundry_lh, parc_glasser(1:32492), ...
    slopes_180, 'faces', [0.8 0.8 0.8; autumn], [], [-0.099, 0.046]);
title('Slopes — flat map');
 
 
%% ---- 8. FIGURE 3D: Slope vs FC / MPC gradients -----------------------
%
%  Correlate the per-parcel GD slope with the first FC gradient (FCG1)
%  and first MPC gradient (MPCG1) from an independent dataset (MICs).
%
%  Input:  fsLR_5k_YWNE.mat, FCGS_100subjs_*.mat, MPCGS_100subjs_*.mat,
%          index_glasser_fsLR5k.mat
 
load(fullfile(data_path,'fsLR_5k_YWNE.mat'),            'coord','tri','mask');
C69_10k    = struct('coord', coord, 'tri', int32(tri));
mask_crtx  = mask;
[surf_l, surf_r] = split_surfaces(C69_10k);
num_parcel  = 360;
 
load(fullfile(data_path,'FCGS_100subjs_ses1_FWHM3_group.mat'),     'FCGS_all_group');
load(fullfile(data_path,'MPCGS_100subjs_ses1_FWHM3_YW_group.mat'), 'MPCGS_all_group');
load(fullfile(data_path,'index_glasser_fsLR5k.mat'),                'index_glasser_fsLR5k');
 
label_parcel360 = index_glasser_fsLR5k';
 
% Helper: map vertex-wise gradient to Glasser-360 parcels
map_gradient_to_glasser = @(G1) arrayfun(@(ii) ...
    mean(G1(label_parcel360 == ii)), 1:num_parcel)';
 
% --- FCG1 ---
FCG_group_align = mean(FCGS_all_group, 3);
FCG1 = FCG_group_align(:,1);
FCG1 = FCG1 / max(abs(FCG1));
 
gs_all = zeros(size(mask_crtx,1), 1);
gs_all(mask_crtx) = FCG1;
FCG1_glasser = map_gradient_to_glasser(gs_all);
 
B = plot_hemispheres(FCG1_glasser, {surf_l32k,surf_r32k}, ...
    'parcellation',parc_glasser,'labeltext',{'FCG1'});
colormap([0.8 0.8 0.8; imola]); B.colorlimits([-0.70, 0.90]);
 
figure('Color','w');
scatter(FCG1_glasser(1:180), slopes_180, 'filled', ...
    'MarkerFaceColor',[0.1 0.1 0.1],'MarkerFaceAlpha',0.6);
xlabel('FCG1','FontSize',16); ylabel('Slope','FontSize',16);
set(gca,'FontSize',15);
[rho_slope_FCG1, p_slope_FCG1] = corr(FCG1_glasser(1:180), slopes_180, 'type','spearman');
fprintf('FCG1 vs slope: rho = %.3f, p = %.4f\n', rho_slope_FCG1, p_slope_FCG1);
 
% --- MPCG1 ---
MPCG_group_align = mean(MPCGS_all_group, 3);
MPCG1 = MPCG_group_align(:,1);
MPCG1 = MPCG1 / max(abs(MPCG1));
 
gs_all = zeros(size(mask_crtx,1), 1);
gs_all(mask_crtx) = MPCG1;
MPCG1_glasser = map_gradient_to_glasser(gs_all);
 
B = plot_hemispheres(MPCG1_glasser, {surf_l32k,surf_r32k}, ...
    'parcellation',parc_glasser,'labeltext',{'MPCG1'});
 
figure('Color','w');
scatter(MPCG1_glasser(1:180), slopes_180, 'filled', ...
    'MarkerFaceColor',[0.1 0.1 0.1],'MarkerFaceAlpha',0.6);
xlabel('MPCG1','FontSize',16); ylabel('Slope','FontSize',16);
set(gca,'FontSize',15);
[rho_slope_MPCG1, p_slope_MPCG1] = corr(MPCG1_glasser(1:180), slopes_180, 'type','spearman');
fprintf('MPCG1 vs slope: rho = %.3f, p = %.4f\n', rho_slope_MPCG1, p_slope_MPCG1);
 
 
%% ---- 9. SUPPLEMENTARY: Reorder networks by FCG1 and re-plot boxplot ---
%
%  Recompute FCG1 from the HCP rest FC group average (in-dataset gradient),
%  reorder the 12 networks by their mean FCG1 value, and redraw the
%  slope boxplot with colours matched to FCG1-ordered network position.
 
gm = GradientMaps();
gm = gm.fit(mean(FC_rest_93, 3));
FCG1_hcp = gm.gradients{1}(:,1);
 
A = plot_hemispheres(FCG1_hcp, {surf_l32k,surf_r32k}, ...
    'parcellation',parc_glasser,'labeltext',{'FCG1 (HCP)'});
 
% Mean FCG1 per network → sort ascending
mean_FCG1_per_net = arrayfun(@(i) mean(FCG1_hcp(label_12n==i)), 1:12)';
[~, idx_reorder]  = sort(mean_FCG1_per_net, 'ascend');
 
fprintf('Networks reordered by mean FCG1 (low → high): ');
fprintf('%d ', idx_reorder); fprintf('\n');
 
% Colours for the reordered networks
legend_colors_ro = color12(idx_reorder, :) .^ (1/1.2);  % contrast boost
legend_colors_ro = min(max(legend_colors_ro, 0), 1);
legend_colors_ro(9,:)  = [1, 0.4, 0.6];   % manual tweak: VMM
legend_colors_ro(10,:) = [1, 1,   0  ];   % manual tweak: LAN
 
label_12n_lh = label_12n(1:180);
data_reordered = arrayfun(@(i) slopes_180(label_12n_lh==idx_reorder(i)), ...
                           1:12, 'UniformOutput',false);
grouped_ro   = vertcat(data_reordered{:});
tmp_ro       = arrayfun(@(n) n*ones(length(data_reordered{n}),1), 1:12, 'UniformOutput', false);
group_idx_ro = vertcat(tmp_ro{:});
 
figure('Color','w'); hold on;
boxplot(grouped_ro, group_idx_ro, 'Colors','k', 'Symbol','o', 'Widths',0.5);
h = findobj(gca,'Tag','Box');
for j = 1:12
    patch(get(h(12-j+1),'XData'), get(h(12-j+1),'YData'), ...
          legend_colors_ro(j,:), 'FaceAlpha',0.8,'EdgeColor','k','LineWidth',1.2);
end
set(gca,'XTick',1:12,'XTickLabel',1:12,'FontSize',14);
xlabel('Network (ordered by mean FCG1)','FontSize',16);
ylabel('Slope','FontSize',16);
ylim padded; box off; grid off;
 
% Colour legend bar
figure('Color','w'); hold on;
for i = 1:12
    y0 = (i-1)*0.07;
    patch([0,0.9,0.9,0], [y0, y0, y0+0.05, y0+0.05], ...
          legend_colors_ro(i,:), 'EdgeColor','k','LineWidth',1.0);
    text(0.95, y0+0.025, num2str(i), 'FontSize',12, ...
         'HorizontalAlignment','left','VerticalAlignment','middle');
end
axis off; xlim([0 1]); ylim([0, 12*0.07]);
 
% Brain maps for manually adjusted networks (VMM = 11, LAN = 6)
for net_id = [11, 6]
    ID_ROI = double(label_12n_lh == net_id);
    A = plot_hemispheres([ID_ROI; zeros(180,1)], {surf_l32k,surf_r32k}, ...
        'parcellation',parc_glasser,'labeltext',{sprintf('Network %d',net_id)});
    if net_id == 11
        colormap([0.9 0.9 0.9; 1, 0.4, 0.6]);
    else
        colormap([0.9 0.9 0.9; 1, 1, 0]);
    end
end
 
fprintf('\nAnalysis complete. Figures saved to: %s\n', output_path);
