try:
    import PyQt5
    qt_version = 5
except ImportError:
    try:
        import PyQt4
        qt_version = 4
    except:
        raise ImportError('PyQt4 or PyQt5 is required')

try:
    import pygments
    is_pygments = True
except ImportError:
    pass
