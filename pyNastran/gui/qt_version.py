found_gui = False

try:
    import PyQt5
    qt_version = 5
    found_gui = True
except ImportError:
    pass

if not found_gui:
    try:
        import PyQt4
        qt_version = 4
        found_gui = True
    except:
        pass


#if not found_gui:
    #try:
        #import PySide
        #qt_version = 'pyside'
        #found_gui = True
    #except:
        #pass

#if not found_gui:
    #raise ImportError('PyQt4, PySide, or PyQt5 is required')
if not found_gui:
    raise ImportError('PyQt4 or PyQt5 is required')

# required to make a pretty console
try:
    import pygments
    is_pygments = True
except ImportError:
    is_pygments = False
