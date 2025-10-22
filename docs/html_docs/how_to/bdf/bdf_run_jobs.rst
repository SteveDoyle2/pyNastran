BDF Run Jobs
============

``bdf run_jobs`` lets you run jobs.


```bash
PS C:\work> bdf run_jobs --help

usage: run_jobs [-x EXE] [-c] [-r] [--skip SKIP [SKIP ...]] [-a ALL] [--args ARGS] [--infile INFILE | --outfile OUTFILE] 
                [--test] [--debug] [-v] bdf_dirname_filename [bdf_dirname_filename ...] [-h]

positional arguments:
  bdf_dirname_filename  path to Nastran filename/directory

options:
  -h, --help            show this help message and exit
  -x EXE, --exe EXE     path to Nastran executable
  -c, --cleanup         cleanup the junk output files (log, f04, plt)
  -r, --recursive       recursively search for directories
  --skip SKIP [SKIP ...]
                        dont process specific files
  -a ALL, --all ALL     dont skip files that have an op2
  --args ARGS           additional arguments
  --infile INFILE       run only files listed in the file; overwrites bdf_dirname_filename
  --outfile OUTFILE     write the files listed to be run
  --test                skip run the jobs
  --debug               more debugging
  -v, --version         show program's version number and exit
```
