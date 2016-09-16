function [eigs,eigvs] = get_eigenvalues_eigenvectors(filenameop2)

%% Function that reads and outputs the eigenfrequencies and the
% eigenvectors from an op2 (PARAM,POST,-1)

% import pyNastran op2/bdf classes

import py.pyNastran.op2.op2.OP2 % import OP2 class

% Instantiate an OP2 class
op2_results = OP2();


%% Read the op2

op2_results.read_op2(filenameop2)

% Save eigenvector structure of a particular SUBCASE

subcase = 1;

eigenvector_struct = op2_results.eigenvectors{subcase}; % In MATLAB curly braces are needed to access to dictionaries

% Convert EIGENFREQUENCIES from list to MATLAB array   

eigrs = cell2mat(cell(eigenvector_struct.eigrs)); % NASTRAN eigenvalues real
eigis = cell2mat(cell(eigenvector_struct.eigis)); % NASTRAN eigenvalues imag

eigs = eigrs+eigis*1i; % tot = real + imag

% Convert EIGENVECTOR ndarray into MATLAB array (GRID,DISPLACEMENT,MODE_NUM)

eigvs = ndarray2mat(eigenvector_struct.data);

end






