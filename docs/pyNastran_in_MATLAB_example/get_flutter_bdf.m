function [Density,Velocity,Mach,Nmodes] = get_flutter_bdf(filenamebdf)
%get_flutter_bdf Reads bdf and provides Flutter condition from FLFACT

import py.pyNastran.bdf.bdf.BDF % import BDF class
import py.numpy.asarray % import as array function (convert list into a ndarray)

% Instantiate an BDF class
bdf = BDF();

%% Read the BDF

bdf.read_bdf(filenamebdf);

% NMODES (is valid only if RESVEC = NO, that is no residual augmentation)

EIGRL_CARDS_ID = ndarray2mat(asarray(py.list(bdf.methods.keys())));
RESVEC_tuple = bdf.case_control_deck.get_subcase_parameter(0,'RESVEC');
RESVEC_cell = RESVEC_tuple.cell;
RESVEC_RH = char(RESVEC_tuple{1}); % RH side of the RESVEC case command

if numel(EIGRL_CARDS_ID) ~= 1
    
    errordlg('If bdf contains more than 1 EIGRL cards, the software does not work. (Things are much more complicated)');
    return
elseif strcmp(RESVEC_RH,'NO')~=1
        errordlg('If bdf does not contain RESVEC = NO commands, the software does not work. (Things are much more complicated)');
    return
end

Nmodes = double(bdf.Method(EIGRL_CARDS_ID).nd);


% FLUTTER SETUP

FLUTTER_CARDS_ID = ndarray2mat(asarray(py.list(bdf.flutters.keys())));

if numel(FLUTTER_CARDS_ID) ~= 1
    
    errordlg('If bdf contains more than 1 FLUTTER cards, the software does not work. (Things are much more complicated)');
    return
    
end


% DENSITY
Density_fact = ndarray2mat(bdf.FLFACT(1).factors);
% Density_fact = [1 2];

% MACH
Mach_fact = ndarray2mat(bdf.FLFACT(2).factors);
% Mach_fact = [1 2];



% VELOCITY
Velocity_fact = ndarray2mat(bdf.FLFACT(3).factors);
Velocity_fact(Velocity_fact>0) = [];
Velocity_fact = -1*Velocity_fact; % In the op2 NASTRAN gives only the eigenvalues/eigenvector associated to negative velocity (sic!)



if strcmp(char(bdf.flutters{FLUTTER_CARDS_ID}.method),'PK')
    
    tmpdensity = repmat(Density_fact(:),numel(Velocity_fact)*numel(Mach_fact),1);
    Density = reshape(tmpdensity,numel(Density_fact)*numel(Velocity_fact)*numel(Mach_fact),1);
    tmpvelocity = repmat(Velocity_fact(:),numel(Density_fact),numel(Mach_fact));
    Velocity = reshape(tmpvelocity',numel(Density_fact)*numel(Velocity_fact)*numel(Mach_fact),1);
    tmpmach = repmat(Mach_fact(:),1,numel(Density_fact)*numel(Velocity_fact));
    Mach = reshape(tmpmach',numel(Density_fact)*numel(Velocity_fact)*numel(Mach_fact),1);
    
else
    strcmp(char(bdf.flutters{FLUTTER_CARDS_ID}.method),'PKNL');
    Density = Density_fact;
    Velocity = Velocity_fact;
    Mach = Mach_fact;
end

end

