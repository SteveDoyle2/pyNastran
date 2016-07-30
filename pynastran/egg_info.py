import importlib

def egg_info(package_name):
    try:
        my_module = importlib.import_module(package_name)
    except ImportError:
        return 'NA', 'NA'

    path = my_module.__path__[0]

    version = 'NA'
    if package_name == 'vtk':
        version = my_module.VTK_VERSION
    else:
        version = my_module.__version__
    return path, version

def batch_egg_info(package_name):
    path, version = egg_info(package_name)
    print('%s (%s): %s' % (package_name, version, path))

def main():
    batch_egg_info('numpy')
    batch_egg_info('scipy')
    batch_egg_info('matplotlib')
    batch_egg_info('vtk')
    batch_egg_info('setuptools')
    batch_egg_info('pyNastran')
    batch_egg_info('pandas')
    #batch_egg_info('pyXML')
    batch_egg_info('pip')

if __name__ == '__main__':
    main()