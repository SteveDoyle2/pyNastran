import os
import yaml
#import pyNastran
dirname = os.path.dirname(__file__)
yml_filename = os.path.join(dirname, 'main.yml')

with open(yml_filename, 'r') as yml_file:
    try:
        print(yaml.safe_load(yml_file))
    except yaml.YAMLError as exc:
        raise
        print(exc)

asdf
