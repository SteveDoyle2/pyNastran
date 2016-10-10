
Matlab
======

 * pyNastran also supports Matlab through the Matlab/Python interface.  [Information about setting up Matlab with Python can be found here.](http://www.mathworks.com/help/matlab/matlab-engine-for-python.html?s_tid=gn_loc_drop)


Example 1 - BDF
===============
Not done...

Example 2 - OP2
===============
This example demonstrates how to call the OP2 class and extract the eigenvectors.


.. code-block:: matlab

    ..

        >>> from pyNastran.bdf.bdf import BDF
        >>> model = BDF()
        >>> model.read_bdf(bdf_filename)

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






