function [] = addContainingDirAndSubDir()

% Add Folder and Its Subfolders to Search Path

    here = mfilename('fullpath');
    [path, ~, ~] = fileparts(here);
    addpath(genpath(path));
end