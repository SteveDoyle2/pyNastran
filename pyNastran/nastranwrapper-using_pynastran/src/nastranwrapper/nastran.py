"""``nastran.py`` defines NastranComponent.
  
"""
import time
import os.path
import inspect
from tempfile import mkdtemp, gettempdir
from shutil import rmtree

import logging

from openmdao.lib.components.external_code import ExternalCode

from openmdao.lib.datatypes.api import Str, Bool, List

from pyNastran.bdf.bdf import BDF
from pyNastran.op2.op2 import OP2
from pyNastran.f06.f06 import F06
from pyNastran.f06.f06 import FatalError

class NastranComponent(ExternalCode):
    """All Nastran-capable components should be subclasses of NastranComponent.

    Your subclass
    must specify how to handle the input and output variables to NastranComponent
    by specifying nastran-specific attributes on the traits. All of these
    attributes are described further in the :ref:`NastranComponent` docs.

    Note: This component does nothing with ``external_files``. If you want to deal with
    that, then do so in your subclass.
    """

    nastran_filename = Str(iotype="in", desc="Input filename with \
                                              placeholder variables.")
    nastran_command = Str(iotype="in", desc="Location of nastran \
                                             executable.")
    nastran_command_args = List(Str, iotype="in",
                                desc="Arguments to the nastran command.")

    output_filename = Str(iotype="out", desc="Output filename.")

    delete_tmp_files = Bool(True, iotype="in", desc="Should I delete \
                            the temporary files?")

    output_tempdir_dir = Str(gettempdir(), iotype="in", desc="Directory in which \
                                            to put the output temp dir.")


    keep_first_iteration = Bool(True, iotype="in", desc="If I am \
    deleting temporary files, should I keep the first one?")

    keep_last_iteration = Bool(True, iotype="in", desc="If I am \
    deleting temporary files, should I keep the last one?")

    log_timings = Bool(False, iotype="in", desc="If True, print out timings for each section" )

    def __init__(self):
        super(NastranComponent, self).__init__()

        # references to Nastran input and output objects
        # All of these objects get their data from Nastran
        # I/O files
        self.bdf = None
        self.f06 = None
        self.op2 = None
        
        # These variables are just to keep track of what we've
        # deleted if you select keep_first_iteration or keep_last_iteration
        self._seen_first_iteration = False
        self._last_seen_iteration = ""

        self.t0 = None # for timing


    def execute(self):
        """Runs the NastranComponent.

        We are overiding ExternalCode's execute function.
        The steps are:
           1. Get a list of the input variables
           2. Read the BDF file into pyNastran's BDF object
           3. Using the info collected in Step #1, replace some
                of the values in the BDF object
           4. Call the update_hook method. Subclasses can override that
               to do processing of the BDF file before it is written out
               again
           5. Write the modified BDF file
           6. Run Nastran
           7. Read the OP2 file ( and the F06 if needed to get Nastran run error info )
           8. Read the results from the OP2 file and set output variables for this
                Component

        RuntimeError
            The component relies on ExternalCode which can throw all
            sorts of RuntimeError-like exceptions (RunStopped,
            RunInterrupted also included).
            
        Filesystem-type Errors
            NastranComponent makes a temporary directory with mkdtemp
            in the temp module. If that fails, the error just
            propagates up.


        While there are no explicit parameters or return values for this
        function, it gets all the input it needs from the design
        variables that are connected to the subclass of NastranComponent.
        This should be described pretty well in the :ref:`documentation<NastranComponent>`.

        """

        # all of these are {"traitname" : trait}
        smart_replacements = {}
        output_variables = {}
        grid_outputs = {}

        for name, trait in self.traits().iteritems():
            if trait.iotype == "in":
                if trait.nastran_card and trait.nastran_id and trait.nastran_field:
                    smart_replacements[name] = trait
                elif trait.nastran_card or trait.nastran_id or trait.nastran_fieldnum:
                    raise RuntimeError("You specified at least one of " + \
                                    "nastran_card, nastran_id, and " + \
                                    "nastran_fieldnum, but you did " + \
                                    "not specify all of them. You " + \
                                    "most probably mistyped.")

            elif trait.iotype == "out":

                # if we want to supply a function that will parse
                # out the wanted information from the output object
                if trait.nastran_func:
                    output_variables[name] = trait

                # this is the grid method of accessing. We have to
                # specify a header, id, and column and
                # the output variable will be set to that value
                elif trait.nastran_header and trait.nastran_constraints :
                    grid_outputs[name] = trait
                elif trait.nastran_header or trait.nastran_constraints:
                    raise RuntimeError("You specified at least one of " + \
                                    "nastran_header and nastran_constraints"+\
                                    ", but you " + \
                                    "did not specify all them. You " + \
                                    "most probably mistyped")

        # do our work in a tmp dir
        tmpdir = mkdtemp(dir = self.output_tempdir_dir)
        tmppath = os.path.join(tmpdir, "input.bdf")

        pyNastran_get_card_methods = {
            'PSHELL': 'Property',
            'PROD': 'Property',
            'FORCE': 'Load',
            'MAT1': 'Material',
            }

        ########## Read BDF ##########
        self.timing_section( "Read BDF" )

        self.bdf = BDF(debug=False,log=logging.getLogger() )
        self.bdf.readBDF(self.nastran_filename,xref=True)

        ########## Modify BDF ##########
        self.timing_section( "Modify BDF" )
        for name, trait in smart_replacements.iteritems():
            value = getattr(self, name)
            nastran_id = int( trait.nastran_id )
            get_method = getattr( self.bdf, pyNastran_get_card_methods[ trait.nastran_card ] )
            # some of these methods have an extra arg for error reporting
            args = inspect.getargspec(get_method).args
            if 'msg' in args:
                nastran_item = get_method( nastran_id, 'dummy msg' )
            else:
                nastran_item = get_method( nastran_id )

            if trait.nastran_card == 'FORCE' :
                nastran_item = nastran_item[0]

            setattr(nastran_item, trait.nastran_field, value)

        ########## update hook ##########
        self.update_hook()
 
        ########## write modified BDF ##########
        self.timing_section( "Write modified BDF" )
         #self.bdf.write_bdf(tmppath)
        self.bdf.write_bdf(tmppath,precision='double',size=16)

        ########## Run Nastran via subprocess ##########
        self.timing_section( "Run Nastran" )
        self.output_filename = os.path.join(tmpdir, "input.out")
        print self.output_filename # perhaps this should be logged, or something

        # Then we run the nastran file
        if self.nastran_command == 'python':  # True when using fake_nastran.py
            self.command = [self.nastran_command,
                            self.nastran_command_args[0], tmppath]
            self.command.extend(self.nastran_command_args[1:])
        else:
            self.command = [self.nastran_command, tmppath]
            self.command.extend(self.nastran_command_args)
        self.command.extend(["batch=no", "out=" + tmpdir, "dbs=" + tmpdir])

        # This calls ExternalCode's execute which will run
        super(NastranComponent, self).execute()

        ########## read OP2 ##########
        self.timing_section( "Read OP2" )
        op2_filename = self.output_filename[:-4] + '.op2'
        f06_filename = self.output_filename
        self.op2 = OP2(op2_filename, debug=False,log=None)
        #self.op2.make_op2_debug = True   # can create a HUGE file that slows things down a lot

        if os.path.exists(op2_filename):
            try:
                self.op2.read_op2()  # doesn't tell you what the error message is
            except FatalError:
                try:
                    self.f06 = F06(f06_filename,debug=False)  # debug True makes it slow
                    self.f06.read_f06()
                except FatalError as err:
                    raise RuntimeError('Nastran fatal error:' + str( err ) )
        elif os.path.exists(f06_filename):
            try:
                self.f06 = F06(f06_filename,debug=False)  # debug True makes it slow
                self.f06.read_f06()  # this will stop with a FatalError with the proper FATAL message
            except FatalError as err:
                raise RuntimeError('Nastran fatal error:' + str( err ) )
        else:
            raise RuntimeError('nastran fatal error' )

        ########## get the outputs using pyNastran ##########
        self.timing_section( "Set outputs using data from OP2" )

        from pyNastran.utils import object_attributes
        for name, trait in grid_outputs.iteritems():
            if trait.nastran_header == 'displacements' :
                var = getattr(self.op2, trait.nastran_header)
                case = var[trait.nastran_subcase]
                for key, eid in trait.nastran_constraints.iteritems():  # e.g., "translations", 18
                    if not hasattr(case, key):
                        #op2.displacements[isubcase=1] doesn't have an attribute "translation".
                        #valid attributes are:  ["translations", "rotations"]
                        msg = "op2.%s[isubcase=%i] does not have an attribute '%s'. " % (trait.nastran_header, trait.nastran_subcase, key )
                        msg += "Valid attributes are: %s" % str(['%s' % att for att in object_attributes(case)])
                        raise KeyError(msg)
                    results_data = getattr(case, key)  # get the translations
                    if trait.nastran_time_step_freq_mode:
                        disp = results_data[trait.nastran_time_step_freq_mode][eid]
                    else:  # "transient" result
                        disp = results_data[eid] # get the specific ID
                    setattr( self, name, disp[trait.nastran_dof] )
            else:
                raise RuntimeError("The Nastran header, %s, is not supported yet" % trait.nastran_header )

        # displacement_columns = ['T1','T2','T3']
        # for name, trait in grid_outputs.iteritems():
        #     table = trait.nastran_table
        #     subcase = trait.nastran_subcase
        #     nastran_id = trait.nastran_id
        #     column = trait.nastran_column
        #     if table == "displacement vector" :
        #         ixyz = displacement_columns.index( column )
        #         setattr(self, name, self.op2.displacements[subcase].translations[nastran_id][ixyz])
        
        for output_name, output_trait in output_variables.iteritems():
            # We run trait.nastran_func on op2 to get output values
            if output_trait.nastran_args:
                setattr(self, output_name,
                        output_trait.nastran_func(self.op2,**output_trait.nastran_args))
            else:
                setattr(self, output_name,
                        output_trait.nastran_func(self.op2))

        self.timing_section( "" )

        # get rid of our tmp dir
        tmpdir_to_delete = ""
        if self.delete_tmp_files:
            if self.keep_first_iteration:
                if not self._seen_first_iteration:
                    self._seen_first_iteration = True
                else:
                    if self.keep_last_iteration: # keep both
                        tmpdir_to_delete = self._last_seen_iteration
                        self._last_seen_iteration = tmpdir
                    else: # just keep first
                        tmpdir_to_delete = tmpdir
            else:
                if self.keep_last_iteration: # only keep last
                    tmpdir_to_delete = self._last_seen_iteration
                    self._last_seen_iteration = tmpdir
                else: # don't keep anything
                    tmpdir_to_delete = tmpdir

            if tmpdir_to_delete:
                rmtree(tmpdir_to_delete)

    def update_hook(self):
        """
        Override this method if you want to do some processing of the
        BDF file before it is written out again
        """
        pass

    def timing_section(self,section_name=None):
        '''Used to log the timings of the sections'''

        if not self.log_timings:
            return

        if self.t0 == None : # we know this is the first section
            self.t0 = time.time()
            self.current_section = section_name
        else:
            t1 = time.time()
            print self.current_section, "time", t1 - self.t0
            self.t0 = time.time()
            self.current_section = section_name
