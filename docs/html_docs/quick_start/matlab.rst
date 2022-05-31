=============================
Calling pyNastran from Matlab
=============================

 * pyNastran also supports Matlab through the Matlab/Python interface.  [Information about setting up Matlab with Python can be found here.](http://www.mathworks.com/help/matlab/matlab-engine-for-python.html?s_tid=gn_loc_drop)

Note about Speed
================
There are two ways to pull large data from Python to Nastran.

   1.  Use the Matlab-Python Interface
   2.  Do an Matlab call to Python, dump your OP2 results matrices using to hdf5
       (using h5py) and load them into and load them Matlab.  It's recommended
       that you don't use scipy's MAT reader as it seems to be buggy, not to
       mention that hdf5 has replaced the MAT format in Matlab.

Intuitively, it seems to that Option #1 should be faster, but for large problems,
that doesn't seem to be the case.  Then again, Option #1, would be probably be
better for any geometry related operation.  In other words, test it.


Working around Matlab's oddities
================================
Replace the base redirectstdout.m file (that for my installation is located in the following folder):

    C:\Program Files\MATLAB\MATLAB Production Server\R2015a\toolbox\matlab\external\interfaces\python\+python\+internal\redirectstdout.m

with this file:

    https://github.com/SteveDoyle2/pyNastran/tree/main/docs/pyNastran_in_MATLAB_example/redirectstdout/redirectstdout.m

Also, instead of imports like:

    >>> import py.pyNastran.op2.op2.OP2

use:

    >>> import py.pyNastran.op2.op2.OP2.*
    >>> clear
    >>> import py.pyNastran.op2.op2.OP2.*
    >>> clear import
    >>> import py.pyNastran.op2.op2.OP2

If you don't need all this insanity, please post and say what you did.

Example 1 - BDF
===============
This example demonstrates how to call the BDF class and extract velocity, machs, and densities (FLFACT cards) from a SOL 145 deck

.. code-block:: matlab

    ..

        >>> rho, velocity, mach, nmodes = get_flutter_bdf("model_145.bdf");

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



Example 2 - OP2
===============
This example demonstrates how to call the OP2 class and extract the eigenvectors.


.. code-block:: matlab

    ..

        >>> eigs, eigvs = get_eigenvalues_eigenvectors("model_145.op2");

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






