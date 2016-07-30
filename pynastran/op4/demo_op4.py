#from __future__ import print_function
#import os
#import sys
#import time
#from traceback import print_exc
#from six import iteritems

#import numpy as np

#from docopt import docopt

import pynastran
from pynastran.op4 import OP4


def run_op4(op4_filename, write_op4=True, debug=True,
			stop_on_failure=False):
	print('***debug=%s' % debug)
	assert '.op4' in op4_filename.lower(), 'op4_filename=%s is not an OP4' % op4_filename
	is_passed = False
	stop_on_failure = True
	delete_op4 = True

	#debug = True
	try:
		op4 = OP4(debug=debug)
		op4._new = True
		matrices = op4.read_op4(op4_filename)

		if 0:
			matrices2 = op4.read_op4(op4_filename)

			print(matrices)
			print('matrices =', matrices.keys())

			assert list(sorted(matrices.keys())) == list(sorted(matrices2.keys()))
			for key, (form, matrix) in sorted(iteritems(matrices)):
				form2, matrix2 = matrices2[key]
				assert form == form2
				delta = matrix - matrix2
				assert np.array_equal(matrix, matrix2), 'delta=\n%s' % delta

		if write_op4:
			model = os.path.splitext(op4_filename)[0]
			op4.write_op4(model+'.test_op4_ascii.op4', matrices, is_binary=False)
			op4.write_op4(model+'.test_op4_binary.op4', matrices, is_binary=True)
			if delete_op4:
				try:
					os.remove(model+'.test_op4_ascii.op4')
					os.remove(model+'.test_op4_binary.op4')
				except:
					pass

		del op4
		is_passed = True
	except KeyboardInterrupt:
		sys.stdout.flush()
		print_exc(file=sys.stdout)
		sys.stderr.write('**file=%s\n' % op4_filename)
		sys.exit('keyboard stop...')
	#except RuntimeError: # the op2 is bad, not my fault
	#	is_passed = True
	#	if stop_on_failure:
	#		raise
	#	else:
	#		is_passed = True

	except IOError: # missing file
		if stop_on_failure:
			raise
	#except AssertionError:
	#	is_passed = True
	#except RuntimeError:
	#	is_passed = True
	except SystemExit:
		#print_exc(file=sys.stdout)
		#sys.exit('stopping on sys.exit')
		raise
	#except NameError:  # variable isnt defined
	#	if stop_on_failure:
	#		raise
	#	else:
	#		is_passed = True
	#except IndexError:
	#	is_passed = True
	except SyntaxError: #Param Parse
		if stop_on_failure:
			raise
		is_passed = True
	except:
		#print e
		if stop_on_failure:
			raise
		else:
			print_exc(file=sys.stdout)
			is_passed = False
	return is_passed


def main():
	ver = str(pyNastran.__version__)

	msg = "Usage:\n"

	# all
	# release
	# current
	msg += "test_op4 [-o] [-d] OP4_FILENAME\n"
	msg += "  test_op4 -h | --help\n"
	msg += "  test_op4 -v | --version\n"
	msg += "\n"
	msg += "Tests to see if an OP4 will work with pyNastran %s.\n" % ver
	msg += "\n"
	msg += "Positional Arguments:\n"
	msg += "  OP4_FILENAME		 Path to OP4 file\n"
	msg += "\n"
	msg += "Options:\n"
	msg += "  -d, --debug		  Developer Debug (default=False)\n"
	msg += "  -o, --write_op4	  Writes the op2 to fem.test_op4.op4 (default=True)\n"
	msg += "  -h, --help		   Show this help message and exit\n"
	msg += "  -v, --version		Show program's version number and exit\n"

	if len(sys.argv) == 1:
		sys.exit(msg)

	data = docopt(msg, version=ver)
	#print("data", data)

	for key, value in sorted(iteritems(data)):
		print("%-12s = %r" % (key.strip('--'), value))

	t0 = time.time()
	run_op4(
		data['OP4_FILENAME'],
		write_op4=data['--write_op4'],
		debug=data['--debug'],
	)
	print("dt = %f" % (time.time() - t0))


if __name__ == '__main__':  # op4
	main()
