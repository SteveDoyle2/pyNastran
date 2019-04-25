import os
import psutil
import gc
import objgraph
import random

from pyNastran.op2.op2 import read_op2,OP2

# we need to check the memory usage
is_linux = None
is_memory = True
#is_memory = False
collect = False
if is_memory:
    try:  # pragma: no cover
        if os.name == 'nt':  # windows
            windows_flag = True
            is_linux = False
            import wmi
            comp = wmi.WMI()

            """Functions for getting memory usage of Windows processes."""
            __all__ = ['get_current_process', 'get_memory_info', 'get_memory_usage']

            import ctypes
            from ctypes import wintypes

            GetCurrentProcess = ctypes.windll.kernel32.GetCurrentProcess
            GetCurrentProcess.argtypes = []
            GetCurrentProcess.restype = wintypes.HANDLE

            SIZE_T = ctypes.c_size_t

            class PROCESS_MEMORY_COUNTERS_EX(ctypes.Structure):
                """
                Windows memory tool that has the same interface as
                `resource` on Linux/Mac.
                """
                _fields_ = [
                    ('cb', wintypes.DWORD),
                    ('PageFaultCount', wintypes.DWORD),
                    ('PeakWorkingSetSize', SIZE_T),
                    ('WorkingSetSize', SIZE_T),
                    ('QuotaPeakPagedPoolUsage', SIZE_T),
                    ('QuotaPagedPoolUsage', SIZE_T),
                    ('QuotaPeakNonPagedPoolUsage', SIZE_T),
                    ('QuotaNonPagedPoolUsage', SIZE_T),
                    ('PagefileUsage', SIZE_T),
                    ('PeakPagefileUsage', SIZE_T),
                    ('PrivateUsage', SIZE_T),
                ]

            GetProcessMemoryInfo = ctypes.windll.psapi.GetProcessMemoryInfo
            GetProcessMemoryInfo.argtypes = [
                wintypes.HANDLE,
                ctypes.POINTER(PROCESS_MEMORY_COUNTERS_EX),
                wintypes.DWORD,
            ]
            GetProcessMemoryInfo.restype = wintypes.BOOL

            def get_current_process():
                """Return handle to current process."""
                return GetCurrentProcess()

            def get_memory_info(process=None):
                """Return Win32 process memory counters structure as a dict."""
                if process is None:
                    process = get_current_process()
                counters = PROCESS_MEMORY_COUNTERS_EX()
                ret = GetProcessMemoryInfo(process, ctypes.byref(counters),
                                           ctypes.sizeof(counters))
                if not ret:
                    raise ctypes.WinError()
                info = dict((name, getattr(counters, name))
                            for name, _ in counters._fields_)
                return info

            def get_memory_usage(process=None):
                """Return this process's memory usage in bytes."""
                info = get_memory_info(process=process)
                return info['PrivateUsage']

        elif os.name in ['posix', 'mac']:  # linux/mac
            import resource
            windows_flag = False
            is_linux = True
        else:
            raise NotImplementedError('os.name=%r and must be nt, posix, mac' % os.name)
    except:
        raise
        is_memory = False

#--------------------------------------------------------------------------------------------

check_memory = True
op2_filename = 'bwb_saero.op2'
exclude_results = ['cquad4_composite_stress', 'ctria3_composite_stress',
                   'cquad4_composite_strain', 'ctria3_composite_strain'] # 'displacements',
exclude_results = ['stress', 'strain', 'grid_point_forces', 'spc_forces'] # 'displacements',
ntimes = 200
m_init = []
m_final = []
for i in range(ntimes):
    if check_memory and 1:
        if is_memory:  # pragma: no cover
            if is_linux: # linux
                kb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
            else: # windows
                kb = get_memory_usage() / 1024
            mb = kb / 1024.
            print("Memory usage start: %s (KB); %.2f (MB)" % (kb, mb))
        else:
            raise RuntimeError('wmi (for Windows) or resource (for Linux/Mac) cannot be found')
        m_init.append(mb)

    #m0 = psutil.Process(os.getpid()).memory_info().rss/1024**3

    #objgraph.show_growth(limit=3)
    if 1:
        model = read_op2(op2_filename,
                       build_dataframe=False, debug=False,
                       exclude_results=exclude_results,
                       skip_undefined_matrices=True)
    else:
        model = OP2()
        model._cleanup_words()

    if 0:
        model.remove_empty_results()
        del model.case_control_deck
        del model.matrices
        for name in model.object_methods():
            try:
                delattr(model, name)
            except AttributeError:
                pass
        objgraph.show_growth()
        objgraph.show_refs(model, too_many=100, filename='sample-graph.png')
        #objgraph.show_backrefs(model, too_many=1000, filename='sample-backref-graph.png')

        # bad...
        #objgraph.show_chain(objgraph.find_backref_chain(model), filename='chain.png')
        objgraph.show_chain(
            objgraph.find_backref_chain(
                random.choice(objgraph.by_type('dict')),
                objgraph.is_proper_module),
            filename='chain.png')

        roots = objgraph.get_leaking_objects()
        objgraph.show_most_common_types(objects=roots)
        objgraph.show_refs(roots[:3], refcounts=True, too_many=100, filename='roots.png')
        break
    del model

    if collect:
        print('Unreachable objects:', gc.collect())
        asdf
    #del gc.garbage[:]
    #print('Unreachable objects:', gc.collect())

    #--------------------------------------------------------------
    if is_memory and check_memory:  # pragma: no cover
        if is_linux: # linux
            kb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        else: # windows
            kb = get_memory_usage() / 1024
        mb = kb / 1024.
        m_final.append(mb)
        print("Memory usage   end: %s (KB); %.2f (MB)" % (kb, mb))

    #m1 = psutil.Process(os.getpid()).memory_info().rss/1024**3
    #print(m0, m1)
    #print(f'{m0:6.4f}|{m1:6.4f}|{id(op2)}')

if len(m_init) > 1:
    import matplotlib.pyplot as plt
    plt.plot(m_init, '-o', label='init')
    #plt.plot(m_final, '-o', label='final')
    plt.legend()
    plt.grid(True)
    plt.show()
