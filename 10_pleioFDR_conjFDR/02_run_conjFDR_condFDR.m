%% 13b_run_pleiofdr.m
%  Run conjFDR and condFDR for TRAIT_B × {TRAIT_A, TRAIT_C, TRAIT_D, TRAIT_E}
%  using pleiotropyFDR (Andreassen et al., 2013; Smeland et al., 2020)
%
%  Prerequisites:
%    1. MATLAB (>= 2015b)
%    2. Reference file:
%         wget https://precimed.s3-eu-west-1.amazonaws.com/pleiofdr/ref9545380_1kgPhase3eur_LDr2p1.mat
%       Place in: ~/pleiofdr_ref/
%    3. .mat sumstats created by 13a_format_sumstats_pleiofdr.py
%
%  Run from MATLAB command line:
%    cd ${PROJECT_ROOT} pair (A x B)/analysis
%    run('13b_run_pleiofdr.m')
%  Or from terminal:
%    matlab -nodisplay -nosplash -r "run('13b_run_pleiofdr.m'); exit"

clear; clc;

PLEIOFDR_DIR = '${PROJECT_ROOT}';
addpath(PLEIOFDR_DIR);

REF_FILE  = '${PROJECT_ROOT}';
MAT_DIR   = '${PROJECT_ROOT}';
OUT_BASE  = '${PROJECT_ROOT}';

if ~exist(REF_FILE, 'file')
    error(['Reference file not found: ' REF_FILE]);
end

ASD_MAT = fullfile(MAT_DIR, 'TRAIT_B.mat');

% Pairs: {gut_disease_label, .mat_file, display_name}
PAIRS = {
    'IBS',         fullfile(MAT_DIR, 'TRAIT_A.mat'),         'ASD_vs_IBS';
    'de Lange_IBD', fullfile(MAT_DIR, 'de Lange_IBD.mat'), 'ASD_vs_IBD_deLange';
    'de Lange_CD',  fullfile(MAT_DIR, 'de Lange_CD.mat'),  'ASD_vs_CD_deLange';
    'de Lange_UC',  fullfile(MAT_DIR, 'de Lange_UC.mat'),  'ASD_vs_UC_deLange';
    'Liu_IBD',     fullfile(MAT_DIR, 'Liu_IBD.mat'),     'ASD_vs_IBD_Liu';
    'Liu_CD',      fullfile(MAT_DIR, 'Liu_CD.mat'),      'ASD_vs_CD_Liu';
    'Liu_UC',      fullfile(MAT_DIR, 'Liu_UC.mat'),      'ASD_vs_UC_Liu';
};

% MHC region to exclude (hg19)
EXCLUDE_MHC = '6 26000000 34000000';

for i = 1:size(PAIRS, 1)
    gut_label  = PAIRS{i, 1};
    gut_mat    = PAIRS{i, 2};
    out_label  = PAIRS{i, 3};

    if ~exist(gut_mat, 'file')
        fprintf('WARNING: %s not found, skipping.\n', gut_mat);
        continue;
    end

    fprintf('\n====== %s ======\n', out_label);

    %% --- conjFDR: SNPs jointly associated with TRAIT_B AND gut disease ---
    conj_out = fullfile(OUT_BASE, [out_label '_conjFDR']);
    if ~exist(conj_out, 'dir'), mkdir(conj_out); end

    config_conj = TextConfig(fullfile(PLEIOFDR_DIR, 'config_template.txt'));
    config_conj.set('ref_file',    REF_FILE);
    config_conj.set('trait_folder','');
    config_conj.set('trait_file1', ASD_MAT);
    config_conj.set('trait_name1', 'ASD');
    config_conj.set('trait_file2', gut_mat);
    config_conj.set('trait_name2', gut_label);
    config_conj.set('out_dir',     conj_out);
    config_conj.set('stat_type',   'conjfdr');
    config_conj.set('fdr_thresh',  '0.05');
    config_conj.set('randprune_n', '20');
    config_conj.set('exclude_chr_pos', EXCLUDE_MHC);

    try
        pleiotropy_analysis(config_conj);
    catch ME
        fprintf('ERROR in conjFDR %s: %s\n', out_label, ME.message);
    end

    %% --- condFDR (TRAIT_B | gut): SNPs enriched for TRAIT_B given gut-disease signal ---
    cond_out = fullfile(OUT_BASE, [out_label '_condFDR_ASD_given_gut']);
    if ~exist(cond_out, 'dir'), mkdir(cond_out); end

    config_cond = TextConfig(fullfile(PLEIOFDR_DIR, 'config_template.txt'));
    config_cond.set('ref_file',    REF_FILE);
    config_cond.set('trait_folder','');
    config_cond.set('trait_file1', ASD_MAT);
    config_cond.set('trait_name1', 'ASD');
    config_cond.set('trait_file2', gut_mat);
    config_cond.set('trait_name2', gut_label);
    config_cond.set('out_dir',     cond_out);
    config_cond.set('stat_type',   'condfdr');
    config_cond.set('fdr_thresh',  '0.01');
    config_cond.set('randprune_n', '20');
    config_cond.set('exclude_chr_pos', EXCLUDE_MHC);

    try
        pleiotropy_analysis(config_cond);
    catch ME
        fprintf('ERROR in condFDR %s: %s\n', out_label, ME.message);
    end

    fprintf('Done: %s\n', out_label);
end

fprintf('\nAll pleiotropyFDR runs complete.\nResults in: %s\n', OUT_BASE);
fprintf('Next: run 13c_analyze_pleiofdr.R to parse and visualize results.\n');
