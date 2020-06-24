function root_path = solverRootPath()
%SOLVERROOTPATH Returns root directory of polynomial matlab package

file_path  = fileparts(which('solverRootPath'));
file_dirs  = regexp(file_path, filesep, 'split');
root_path  = strjoin(file_dirs(1:end-1), filesep);
end