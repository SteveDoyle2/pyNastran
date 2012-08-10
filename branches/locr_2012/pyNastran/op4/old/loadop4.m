% [a,b,c,...] = loadop4(filename [, nSkip ])
%
% Load matrices stored in a Nastran .op4 file into matlab variables.
% All .op4 file formats, storage schemes, and matrix types are supported.
% Matrices saved with either sparse format (BIGMAT=TRUE or BIGMAT=FALSE)
% in the .op4 file are loaded as matlab sparse matrices.  For binary files, 
% the code will automatically swap bytes if a little-endian machine reads
% files created on a big-endian machine and vice versa.
%
% Input:
%     filename    The name of the .op4 file.
%     nSkip       [optional] The number of matrices to skip over before
%                 loading the first matrix.
% Output:
%     a,b,c,...   The loaded matrices.
%
%  See also:      saveop4
