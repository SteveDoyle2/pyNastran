clear all
close all
clc

% Load eigenvalues/eigenvectors from op2 through pyNastran
% SOL145 performed with PK method, RESVEC = NO

[eigs,eigvs] = get_eigenvalues_eigenvectors('BAH_Plane.op2');

eigs = transpose(eigs);

% Load Flutter data from bulk data
[Density,Velocity,Mach,Nmodes] = get_flutter_bdf('BAH_Plane.bdf');
Nsweep_par = numel(Velocity);
Ncases_op2_check = Nsweep_par*Nmodes;
Ndofs = 6; % T1 T2 T3 R1 R2 R3 NASTRAN DOFs
Ngrids = size(eigvs,1); % number of GRIDs can be read from the eigenvectors

tmpEigvs = NaN*zeros(Nsweep_par,Ngrids*Ndofs,Nmodes);

for k = 1:Nsweep_par

    % Eigenvalues
    modeidx_nth_par = 1+Nmodes*(k-1):k*Nmodes;
    Eigs(:,k) = eigs(modeidx_nth_par); %(2D matrix MODES x NPAR array)

    % Eigenvectors
    tmp = squeeze(eigvs(:,:,modeidx_nth_par));
    for m = 1:Nmodes
        for g = 1:Ngrids
            tmpEigvs(k,(Ndofs*(g-1)+1:Ndofs*g),m) = tmp(g,:,m); % (3D matrix NPAR x [NGRIDSxNDOFS(6)] x NMODES)
        end
    end

end


% Creation of data structure for PoleTrkr

data.param = Velocity; % (1D array NPAR x 1)
data.Eigs = Eigs; %(2D matrix MODES x NPAR array)
data.Eigvs = tmpEigvs; % (3D matrix NPAR x [NGRIDSxNDOFS(6)] x NMODES)
