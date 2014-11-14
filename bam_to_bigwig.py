#!/usr/bin/env python
'''Convert your alignment files in BAM format into coverage files in bigWig format

Usage:
    bam_to_bigwig.py BAM_FILE [BAM_FILE ...]
    [-o, --bigwig_filename=<output file name>
     -t, --tempfile
     -k, --keep-tempfile
     -s, --ignore-secondary
     -q, --ignore-qc-fail
     -d, --ignore-optical-pcr-duplicate
     -u, --ignore-supplementary]

--bigwig_filename: if not given will save to <BAM file prefix>.bigwig.
--tempfile: If this is set, will use tempfile library to generate file names instead of using <BAM file prefix>.wig and <BAM file prefix>.sizes.
--keep-tempfile: If this is set, will not delete the intermediate files <BAM file prefix>.wig and <BAM file prefix>.sizes.
--ignore-secondary: If this is set, rsem-bam2wig will ignore alignments with the "secondary alignment" flag bit 0x100 set.
--ignore-qc-fail: If this is set, rsem-bam2wig will ignore alignments with the "not passing quality controls" flag bit 0x200 set.
--ignore-optical-pcr-duplicate: If this is set, rsem-bam2wig will ignore alignments with the "PCR or optical duplicate" flag bit 0x400 set.
--ignore-supplementary: If this is set, rsem-bam2wig will ignore alignments with the "supplementary alignment" flag bit 0x800 set.
'''
'''
Han Lin, 2014
hotdogee[at]gmail[dot]com
Biomedical Electronics and Bioinformatics, National Taiwan University
National Agricultural Library, ARS, USDA

Dependencies:
    pysam (https://github.com/pysam-developers/pysam)
    wigToBigWig from UCSC (http://hgdownload.cse.ucsc.edu/admin/exe/)
    rsem-bam2wig from RSEM (http://deweylab.biostat.wisc.edu/rsem/)

Installation:
* pysam
  1. Install build dependencies
    * CentOS: yum install python-devel zlib-devel
    * Ubuntu: apt-get install python-dev zlib1g-dev
  2. Install using pip
    * pip install pysam
* wigToBigWig
  1. Download binary: http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/wigToBigWig
  2. Make sure wigToBigWig is in a directory in your $PATH
    * cp wigToBigWig /usr/local/bin
* rsem-bam2wig
  1. Download the source code: http://deweylab.biostat.wisc.edu/rsem/
  2. Compile
    * make
  3. Make sure rsem-bam2wig is in a directory in your $PATH
    * cp rsem-bam2wig /usr/local/bin

Changelog:
v1.3:
* Check write permissions before calling rsem-bam2wig and wigToBigWig
* More informative error messages
* Hide output when running dependency checking 
v1.2:
* Added dependency checking
* Better documentation
v1.1:
* Print progress
v1.0:
* Fixed usage documentation
'''

__version__ = '1.2'

import os
import sys
import tempfile
import errno
from subprocess import call, check_call, Popen, PIPE, STDOUT
from optparse import OptionParser
from contextlib import contextmanager
import logging
log = logging.getLogger(__name__)

unix_signals_tsv = """1	SIGHUP	Exit	Hangup
2	SIGINT	Exit	Interrupt
3	SIGQUIT	Core	Quit
4	SIGILL	Core	Illegal Instruction
5	SIGTRAP	Core	Trace/Breakpoint Trap
6	SIGABRT	Core	Abort
7	SIGEMT	Core	Emulation Trap
8	SIGFPE	Core	Arithmetic Exception
9	SIGKILL	Exit	Killed
10	SIGBUS	Core	Bus Error
11	SIGSEGV	Core	Segmentation Fault
12	SIGSYS	Core	Bad System Call
13	SIGPIPE	Exit	Broken Pipe
14	SIGALRM	Exit	Alarm Clock
15	SIGTERM	Exit	Terminated
16	SIGUSR1	Exit	User Signal 1
17	SIGUSR2	Exit	User Signal 2
18	SIGCHLD	Ignore	Child Status
19	SIGPWR	Ignore	Power Fail/Restart
20	SIGWINCH	Ignore	Window Size Change
21	SIGURG	Ignore	Urgent Socket Condition
22	SIGPOLL	Ignore	Socket I/O Possible
23	SIGSTOP	Stop	Stopped (signal)
24	SIGTSTP	Stop	Stopped (user)
25	SIGCONT	Ignore	Continued
26	SIGTTIN	Stop	Stopped (tty input)
27	SIGTTOU	Stop	Stopped (tty output)
28	SIGVTALRM	Exit	Virtual Timer Expired
29	SIGPROF	Exit	Profiling Timer Expired
30	SIGXCPU	Core	CPU time limit exceeded
31	SIGXFSZ	Core	File size limit exceeded
32	SIGWAITING	Ignore	All LWPs blocked
33	SIGLWP	Ignore	Virtual Interprocessor Interrupt for Threads Library
34	SIGAIO	Ignore	Asynchronous I/O"""
unix_signals = dict([(-int(t[0]), t) for t in [line.split('\t') for line in unix_signals_tsv.split('\n')]])

def query_yes_no(question, default="yes"):
    """Ask a yes/no question via raw_input() and return their answer.

    "question" is a string that is presented to the user.
    "default" is the presumed answer if the user just hits <Enter>.
        It must be "yes" (the default), "no" or None (meaning
        an answer is required of the user).

    The "answer" return value is one of "yes" or "no".
    """
    valid = {"yes":True,   "y":True,  "ye":True,
             "no":False,     "n":False}
    if default == None:
        prompt = " [y/n] "
    elif default == "yes":
        prompt = " [Y/n] "
    elif default == "no":
        prompt = " [y/N] "
    else:
        raise ValueError("invalid default answer: '%s'" % default)

    while True:
        sys.stdout.write(question + prompt)
        choice = raw_input().lower()
        if default is not None and choice == '':
            return valid[default]
        elif choice in valid:
            return valid[choice]
        else:
            sys.stdout.write("Please respond with 'yes' or 'no' "\
                             "(or 'y' or 'n').\n")

def bam_to_wig(bam_filename, wig_filename, ignore_secondary=False, ignore_qc_fail=False, ignore_optical_pcr_duplicate=False, ignore_supplementary=False):
    file_prefix = os.path.splitext(bam_filename)[0]
    if not wig_filename:
        wig_filename = '%s.wig' % file_prefix
    # check if we have permissions to write wig_filename if file exists or permissions to write its directory to create it
    if os.path.isfile(wig_filename):
        if not os.access(wig_filename, os.W_OK):
            log.critical('Write permission denied: %s' % (os.path.abspath(wig_filename)))
            return False
    elif not os.access(os.path.dirname(os.path.abspath(wig_filename)), os.W_OK):
        log.critical('Write permission denied: %s' % (os.path.dirname(os.path.abspath(wig_filename))))
        return False
    cl = ['rsem-bam2wig', bam_filename, wig_filename, file_prefix, '--no-fractional-weight']
    if ignore_secondary:
        cl.append('--ignore-secondary')
    if ignore_qc_fail:
        cl.append('--ignore-qc-fail')
    if ignore_optical_pcr_duplicate:
        cl.append('--ignore-optical-pcr-duplicate')
    if ignore_supplementary:
        cl.append('--ignore-supplementary')
    p = Popen(cl)
    p.communicate()
    rc = p.returncode
    if rc != 0:
        if rc < 0:
            try:
                log.critical('"%s" was terminated by signal %d: %s, %s, %s' % (' '.join(cl), -rc, unix_signals[rc][1], unix_signals[rc][2], unix_signals[rc][3]))
            except KeyError:
                log.critical('"%s" was terminated by signal %d' % (' '.join(cl), -rc))
        else:
            log.critical('"%s" terminated with non-zero return code: %d' % (' '.join(cl), rc))
    return rc == 0, rc

@contextmanager
def indexed_bam(bam_filename):
    import pysam
    if not os.path.exists(bam_filename + ".bai"):
        pysam.index(bam_filename)
    sam_reader = pysam.Samfile(bam_filename, "rb")
    yield sam_reader
    sam_reader.close()

def get_sizes_from_bam(bam_filename):
    with indexed_bam(bam_filename) as work_bam:
        chr_sizes = zip(work_bam.references, work_bam.lengths)
    return chr_sizes

def bam_to_sizes(bam_filename, chr_sizes_filename):
    chr_sizes = get_sizes_from_bam(bam_filename)
    try:
        fp = open(chr_sizes_filename, "wb")
    except IOError as e:
        if e.errno == errno.EACCES:
            log.critical('Write permission denied: %s' % (os.path.abspath(chr_sizes_filename)))
            return False
        # Not a permission error.
        raise
    else:
        with fp:
            for chrom, size in chr_sizes:
                fp.write("%s\t%s\n" % (chrom, size))
    return True

def wig_to_bigwig(wig_filename, bigwig_filename, chr_sizes_filename):
    file_prefix = os.path.splitext(wig_filename)[0]
    if not bigwig_filename:
        bigwig_filename = '%s.bigwig' % file_prefix
    # check if we have permissions to write bigwig_filename if file exists or permissions to write its directory to create it
    if os.path.isfile(bigwig_filename):
        if not os.access(bigwig_filename, os.W_OK):
            log.critical('Write permission denied: %s' % (os.path.abspath(bigwig_filename)))
            return False
    elif not os.access(os.path.dirname(os.path.abspath(bigwig_filename)), os.W_OK):
        log.critical('Write permission denied: %s' % (os.path.dirname(os.path.abspath(bigwig_filename))))
        return False
    cl = ['wigToBigWig', wig_filename, chr_sizes_filename, bigwig_filename]
    p = Popen(cl)
    p.communicate()
    rc = p.returncode
    if rc != 0:
        if rc < 0:
            try:
                log.critical('"%s" was terminated by signal %d: %s, %s, %s' % (' '.join(cl), -rc, unix_signals[rc][1], unix_signals[rc][2], unix_signals[rc][3]))
            except KeyError:
                log.critical('"%s" was terminated by signal %d' % (' '.join(cl), -rc))
        else:
            log.critical('"%s" terminated with non-zero return code: %d' % (' '.join(cl), rc))
    return rc == 0, rc

def main(bam_filename, bigwig_filename=None, use_tempfile=False, keep_tempfile=False, ignore_secondary=False, ignore_qc_fail=False, ignore_optical_pcr_duplicate=False, ignore_supplementary=False):
        #config = {'program': {'ucsc_bigwig' : 'wigToBigWig', 'rsem-bam2wig' : 'rsem-bam2wig'}}
    file_prefix = os.path.splitext(bam_filename)[0]
    if bigwig_filename is None:
        bigwig_filename = '%s.bigwig' % file_prefix
    if os.path.abspath(bam_filename) == os.path.abspath(bigwig_filename):
        sys.stderr.write('Bad arguments, input and output files are the same.\n')
        sys.exit(1)
    if os.path.exists(bigwig_filename) and os.path.getsize(bigwig_filename) > 0:
        if not query_yes_no('Output file exists and not empty, overwrite?', default='no'):
            sys.exit(1)
    
    # run rsem-bam2wig
    wig_filename = '%s.wig' % file_prefix
    if use_tempfile:
        #Use a temp file to avoid any possiblity of not having write permission
        wig_filename = tempfile.NamedTemporaryFile(delete=False).name
    sys.stdout.write('Building wig file: %s ...\n' % wig_filename)
    if not bam_to_wig(bam_filename, wig_filename, ignore_secondary, ignore_qc_fail, ignore_optical_pcr_duplicate, ignore_supplementary):
        sys.exit(1)
    sys.stdout.write('Done\n')
    
    # generate sizes file
    chr_sizes_filename = '%s.sizes' % file_prefix
    if use_tempfile:
        #Use a temp file to avoid any possiblity of not having write permission
        chr_sizes_filename = tempfile.NamedTemporaryFile(delete=False).name
    sys.stdout.write('Building sizes file: %s ...' % chr_sizes_filename)
    if not bam_to_sizes(bam_filename, chr_sizes_filename):
        sys.stdout.write(' Failed\n')
        sys.exit(1)
    sys.stdout.write(' Done\n')
    
    # run wigToBigWig
    sys.stdout.write('Building bigwig file: %s ...' % bigwig_filename)
    if not wig_to_bigwig(wig_filename, bigwig_filename, chr_sizes_filename):
        sys.stdout.write(' Failed\n')
        sys.exit(1)
    sys.stdout.write(' Done\n')
    
    # remove temp files
    if not keep_tempfile:
        os.remove(chr_sizes_filename)
        os.remove(wig_filename)
        
def dependences_exist():
    """The script requires:
    pysam (https://github.com/pysam-developers/pysam)
    wigToBigWig from UCSC (http://hgdownload.cse.ucsc.edu/admin/exe/)
    rsem-bam2wig from RSEM (http://deweylab.biostat.wisc.edu/rsem/)
    """
    try:
        from subprocess import DEVNULL # py3k
    except ImportError:
        import os
        DEVNULL = open(os.devnull, 'wb')
    all_good = True
    try:
        import pysam
    except ImportError:
        all_good = False
        log.critical('Missing dependency: pysam (https://github.com/pysam-developers/pysam)')
    try:
        call(['wigToBigWig'], stdout=DEVNULL, stderr=STDOUT)
    except OSError:
        all_good = False
        log.critical('Missing dependency: wigToBigWig from UCSC (http://hgdownload.cse.ucsc.edu/admin/exe/)')
    try:
        call(['rsem-bam2wig'], stdout=DEVNULL, stderr=STDOUT)
    except OSError:
        all_good = False
        log.critical('Missing dependency: rsem-bam2wig from RSEM (http://deweylab.biostat.wisc.edu/rsem/)')
    return all_good

if __name__ == '__main__':
    log.setLevel(logging.DEBUG)
    lh = logging.StreamHandler()
    lh.setFormatter(logging.Formatter('\033[01;31m%(levelname)-8s\033[0;0m %(message)s'))
    log.addHandler(lh)
    parser = OptionParser()
    parser.add_option('-o', '--bigwig_filename', dest='bigwig_filename')
    parser.add_option('-t', '--tempfile', dest='use_tempfile', action='store_true', default=False)
    parser.add_option('-k', '--keeptemp', dest='keep_tempfile', action='store_true', default=False)
    parser.add_option('-s', '--ignore-secondary', dest='ignore_secondary', action='store_true', default=False)
    parser.add_option('-q', '--ignore-qc-fail', dest='ignore_qc_fail', action='store_true', default=False)
    parser.add_option('-d', '--ignore-optical-pcr-duplicate', dest='ignore_optical_pcr_duplicate', action='store_true', default=False)
    parser.add_option('-u', '--ignore-supplementary', dest='ignore_supplementary', action='store_true', default=False)
    (options, args) = parser.parse_args()
    if len(args) == 0:
        print __doc__
        dependences_exist()
        sys.exit()
    if not dependences_exist():
        sys.exit()
    kwargs = dict(
        bigwig_filename=options.bigwig_filename,
        use_tempfile=options.use_tempfile,
        keep_tempfile=options.keep_tempfile,
        ignore_secondary=options.ignore_secondary,
        ignore_qc_fail=options.ignore_qc_fail,
        ignore_optical_pcr_duplicate=options.ignore_optical_pcr_duplicate,
        ignore_supplementary=options.ignore_supplementary)
    try:
        for bam_filename in args:
            main(bam_filename, **kwargs)
    except KeyboardInterrupt:
        sys.stdout.write('\n')
        log.info('Aborted by keyboard interrupt')