clear;clc;
addpath(genpath('C:/Users/nnu02/Documents/MATLAB/spm12'));
addpath('C:/Users/nnu02/Documents/MATLAB/fieldtrip-20220827');
ft_defaults

subjects = ["04" "05" "06" "07" "08" "09" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "24" "25" "26" "27" "28" "29"];

folder = 'E:/Gian/GG_SensAtt_Prediction/02Data';
results_first_level = '00Results1stLevel';
results_second_level = strcat(folder, '/2nd_level');

if ~isfolder(results_second_level)
    mkdir(results_second_level)
end

spm eeg

%% SECOND LEVEL STATISTICS: TOUCH vs noTOUCH vs High/Equal/Low PROB vs MOVE / noMOVE
matlabbatch{1}.spm.stats.factorial_design.dir = {results_second_level};
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(1).name = 'subject';
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(1).dept = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(1).variance = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(1).gmsca = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(1).ancova = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(2).name = 'Prob';
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(2).dept = 1;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(2).variance = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(2).gmsca = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(2).ancova = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(3).name = 'Move';
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(3).dept = 1;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(3).variance = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(3).gmsca = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(3).ancova = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(4).name = 'Touch';
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(4).dept = 1;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(4).variance = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(4).gmsca = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(4).ancova = 0;
for i = 1:length(subjects)
first_level_folder = strcat(folder, '/ID', ID, '/021stLevel/');
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fsuball.fsubject(i).scans = {
                                                                                  strcat(first_level_folder, results_first_level, '/beta_0001.nii,1')
                                                                                  strcat(first_level_folder, results_first_level, '/beta_0002.nii,1')
                                                                                  strcat(first_level_folder, results_first_level, '/beta_0003.nii,1')
                                                                                  strcat(first_level_folder, results_first_level, '/beta_0004.nii,1')
                                                                                  strcat(first_level_folder, results_first_level, '/beta_0005.nii,1')
                                                                                  strcat(first_level_folder, results_first_level, '/beta_0006.nii,1')
                                                                                  strcat(first_level_folder, results_first_level, '/beta_0007.nii,1')
                                                                                  strcat(first_level_folder, results_first_level, '/beta_0008.nii,1')
                                                                                  strcat(first_level_folder, results_first_level, '/beta_0009.nii,1')
                                                                                  strcat(first_level_folder, results_first_level, '/beta_0010.nii,1')
                                                                                  strcat(first_level_folder, results_first_level, '/beta_0011.nii,1')
                                                                                  strcat(first_level_folder, results_first_level, '/beta_0012.nii,1')
                                                                                  };
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fsuball.fsubject(i).conds = [1 1 1
                                                                                  2 1 1
                                                                                  3 1 1
                                                                                  1 2 1
                                                                                  2 2 1
                                                                                  3 2 1
                                                                                  1 1 2
                                                                                  2 1 2
                                                                                  3 1 2
                                                                                  1 2 2
                                                                                  2 2 2
                                                                                  3 2 2];
end
matlabbatch{1}.spm.stats.factorial_design.des.fblock.maininters{1}.inter.fnums = [2
                                                                                  3
                                                                                  4];
%matlabbatch{1}.spm.stats.factorial_design.des.fblock.maininters{2}.fmain.fnum = 4;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.maininters{2}.fmain.fnum = 1;
matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
Con=spm_make_contrasts([3 2 2]);
factors_numbers = ["1"; "2"; "3"];
factors_names = ["Prob"; "Movement"; "Touch"];
for i = 2:length(Con)
    matlabbatch{3}.spm.stats.con.consess{i-1}.fcon.name = replace(Con(i).name, factors_numbers, factors_names);
    matlabbatch{3}.spm.stats.con.consess{i-1}.fcon.weights = Con(i).c;
    matlabbatch{3}.spm.stats.con.consess{i-1}.fcon.sessrep = 'none';
end
%then add parametric effect
matlabbatch{3}.spm.stats.con.consess{8}.tcon.name = 'Parametric effect';
matlabbatch{3}.spm.stats.con.consess{8}.tcon.weights = [-3 3 -2 2 -1 1 1 -1 2 -2 3 -3];
matlabbatch{3}.spm.stats.con.consess{8}.tcon.sessrep = 'none';

matlabbatch{3}.spm.stats.con.delete = 0;

spm_jobman('run', matlabbatch);
clear matlabbatch