
def run_lots_of_files(files, write_op4=True,
                      debug=True, save_cases=True, skip_files=None,
                      stop_on_failure=False, nstart=0, nstop=1000000000):
    if skip_files is None:
        skip_files = []
    n = ''
    isubcases = []
    failed_cases = []
    nfailed = 0
    ntotal = 0
    npassed = 0
    t0 = time.time()
    for (i, op4file) in enumerate(files[nstart:nstop], nstart):  # 149
        base_name = os.path.basename(op4file)
        #if baseName not in skipFiles and not base_name.startswith('acms') and i not in nSkip:
        if base_name not in skip_files and '#' not in op4file:
            print("%"*80)
            print('file=%s\n' % op4file)
            n = '%s ' % i
            sys.stderr.write('%sfile=%s\n' %(n, op4file))
            ntotal += 1
            is_passed = run_op4(op4file,
                                debug=debug,
                                stop_on_failure=stop_on_failure) # True/False
            if not is_passed:
                sys.stderr.write('**file=%s\n' % op4file)
                failed_cases.append(op4file)
                nfailed += 1
            else:
                npassed += 1
            #sys.exit('end of test...test_op4.py')

    if save_cases:
        with open('failed_cases.in', 'wb') as failed_file:
            for op4file in failed_cases:
                failed_file.write('%s\n' % op4file)

    seconds = time.time()-t0
    minutes = seconds/60.
    print("dt = %s seconds = %s minutes" % (seconds, minutes))

    msg = '-----done with all models %s/%s=%.2f%%  nFailed=%s-----' %(
        npassed, ntotal, 100.*npassed/float(ntotal), ntotal-npassed)
    print(msg)
    sys.exit(msg)
