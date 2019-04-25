import os

#from ??? import write_include

def write_rpy_file(fname_model, jobname, fname_user_subroutine=None, include=None):
    if fname_user_subroutine is None:
        fname_user_subroutine = ''

    model_name = os.path.basename(fname_model).rsplit('.', 1)[0]

    nprocessors = 1

    lines = ['from abaqus import *\n',
             'from abaqusConstants import *\n',
             'session.Viewport(name="Viewport: 1", origin=(0.0, 0.0), width=300.0, height=100.0)\n'
             'from caeModules import *\n',
             'from driverUtils import executeOnCaeStartup\n',
             'executeOnCaeStartup()\n',
             'Mdb()\n',
             'mdb.Model(name="%s",modelType=STANDARD_EXPLICIT)\n' % (model_name),
             'a = mdb.models["%s"].rootAssembly\n' % (model_name),
             'mdb.models["%s"].setValues(noPartsInputFile=ON)\n' % (model_name),
             'import job\n',
             'mdb.models["%s"].keywordBlock.synchVersions(storeNodesAndElements=False)\n' % (model_name)]
    if include is not None:
        include_lines = 'r"""\n'
        #for fname in include:
            #include_lines += write_include(fname)
        include_lines += '"""'
        lines.append('mdb.models["%s"].keywordBlock.insert(0, %s)\n' % (model_name, include_lines))

    indent = '        '
    lines.extend(['mdb.Job(name="%s", model="%s",\n' % (jobname, model_name),
                  '%sdescription="", type=ANALYSIS, atTime=None, waitMinutes=0, waitHours=0,\n' % (indent),
                  '%squeue=None, memory=90, memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True,\n' % (indent),
                  '%sexplicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=ON,\n' % (indent),
                  '%smodelPrint=ON, contactPrint=ON, historyPrint=ON,\n' % (indent),
                  '%suserSubroutine=r"%s",\n' % (indent, fname_user_subroutine),
                  '%sscratch="", resultsFormat=ODB, parallelizationMethodExplicit=DOMAIN,\n' % (indent),
                  '%snumDomains=1, activateLoadBalancing=False, multiprocessingMode=DEFAULT,\n' % (indent),
                  '%snumCpus=%i)\n' % (indent, nprocessors),
                  'mdb.jobs["%s"].submit(consistencyChecking=OFF)\n' % (jobname),
                  'import sys\n',
                  'sys.exit()\n'])
    return lines

def cmd_line():
    from docopt import docopt
    msg = 'Usage:\n'
    msg += 'spawn_abaqus MODEL_INP JOBNAME [--rpy RPY_FILENAME] [--user USER]\n'
    #msg += 'spawn_abaqus MODEL_INP JOBNAME [--rpy RPY_FILENAME] [--user USER] [--include INC]...\n'
    msg += '\n'
    msg += 'Positonal Arguments:\n'
    msg += '  MODEL_INP  the abaqus.inp file\n'
    msg += '  JOBNAME    the job name\n'
    msg += '\n'
    msg += 'Options:\n'
    msg += '  --user USER         defines a user subroutine\n'
    msg += '  --rpy RPY_FILENAME  defines the rpy file to write\n'
    from pyNastran import version
    data = docopt(msg, version=version)
    model = data['MODEL_INP']
    jobname = data['JOBNAME']
    user = data['--user']
    rpy_filename = data['--rpy']
    if rpy_filename is None:
        rpy_filename = jobname + '.rpy'
    lines = write_rpy_file(model, jobname, fname_user_subroutine=None, include=None)
    with open(rpy_filename, 'w') as rpy:
        rpy.writelines(lines)
    print('rpy_filename =', rpy_filename)

if __name__ == '__main__':  #pragma: no cover
    cmd_line()
