import os
IS_DEV = (
    'TRAVIS' in os.environ or
    'APPVEYOR' in os.environ or
    'READTHEDOCS' in os.environ
)
#IS_TRAVIS = 'TRAVIS' in os.environ
#IS_RTD = 'READTHEDOCS' in os.environ
