#!/usr/bin/python2
import argparse
import fcntl
import signal
import os
import select
import subprocess
import threading
import time

from linear_combinations_worker import *

WORKER_PATH='./linear_combinations_worker.py'

worker_process = None
"""The (singular, for now) worker process that actually finds results"""

last_heartbeat = 0
"""The last time the worker process gave any indication that it was alive"""

worker_args = None
"""Arguments to pass to workers (same as for this program)"""

def preexec_discard_signals():
    """Discard any signals so that the timing gets sent to this process,
    which can actually handle it.

    """
    os.setpgrp()

def get_next_output(out):
    """Given a file-ish object, read from it until there is no more to
    immediately read, then return the result.

    """
    result = ''
    incremental = ''
    readable = []
    try:
        while not readable:
            readable, w, e = select.select( [ out ], [], [], 0.2 )
    except select.error as e:
        _sigint_handler()
    try:
        while True:
            incremental = os.read(out.fileno(), 1024)
            if incremental == '':
                return result
            result = result + incremental
    except OSError:
        return result

def kickoff_worker():
    """Initializes a worker process, waits for it to print out some
    preliminary information, then returns that process.

    """

    new_process = subprocess.Popen(['python2',
                                    WORKER_PATH,
                                    '--spanfile', worker_args.spanfile,
                                    '--output', worker_args.output,
                                    '--precision', str(worker_args.precision),
                                    '-n', str(worker_args.n),
                                    '--lindep-precision', str(worker_args.lindep_precision)],
                                   bufsize = 1,
                                   shell = False,
                                   stdin = subprocess.PIPE,
                                   stdout = subprocess.PIPE,
                                   stderr = subprocess.STDOUT,
                                   preexec_fn = preexec_discard_signals)

    fcntl.fcntl(new_process.stdout, fcntl.F_SETFL,
                fcntl.fcntl(new_process.stdout, fcntl.F_GETFL) | os.O_NONBLOCK)

    output = ''
    while True:
        output = output + get_next_output(new_process.stdout)
        for line in output.split(os.linesep):
            if line == LCC_READY_TO_RECEIVE_DATA_STRINGS:
                last_heartbeat = time.time()
                return new_process
            elif (line == LCC_REFUSE_TO_CLOBBER
                  or line == LCC_CANT_DEAL_WITH_OUTPUT_FILE):
                print('Output path provided must be either a nonexistent '
                      + 'file, an empty file, or a file with the correct '
                      + 'header')
                sys.stdout.flush()
                sys.exit(1)

        time.sleep(1)

def _sigint_handler(signum, frame):
    if worker_process is not None:
        worker_process.terminate()
    try:
        print('Terminating.')
        sys.stdout.flush()
    except IOError:
        pass
    sys.exit(0)

class TimeoutExpired:
    pass

def _raise_timeout(signum, frame):
    raise TimeoutExpired

if __name__ == '__main__':
    signal.signal(signal.SIGINT, _sigint_handler)
    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument('--spanfile', '-s',
                            default = 'spans.csv',
                            help = 'Path of master list of spans for reading')
    arg_parser.add_argument('--output', '-o',
                            default = 'linear_combinations.csv',
                            help = 'Output path (appending handled well)')
    arg_parser.add_argument('-n',
                            default = '2',
                            metavar = 'n',
                            help = 'For considering SL(n,C)')
    arg_parser.add_argument('--precision', '-p',
                            type = int,
                            default = 500,
                            help = 'Precision to set PARI to')
    arg_parser.add_argument('--lindep-precision', '-l',
                            type = int,
                            default = 16,
                            help = 'Precision to use when calling lindep()')
    arg_parser.add_argument('--timeout', '-t',
                            type = int,
                            default = 120,
                            help = 'After this much time (in seconds) '
                            + 'without progress, kill worker')
    worker_args = arg_parser.parse_args()

    print('Starting up')
    sys.stdout.flush()

    # Can't figure out a nice way to pass in ranges of manifolds...
    manifolds = OrientableCuspedCensus[1152:]
    # manifolds = OrientableCuspedCensus[145:148]

    worker_process = kickoff_worker()
    signal.signal(signal.SIGALRM, _raise_timeout)
    for m in manifolds:
        print(str(m) + ' : Computing')
        sys.stdout.flush()
        num_tries = 0
        output = ''
        finished_this_manifold = False
        try:
            if worker_process.poll():
                worker_process = kickoff_worker()

            nice_str = (str(m) + os.linesep).encode(encoding = 'UTF-8')
            worker_process.stdin.write(nice_str)
            worker_process.stdin.flush()

            while not finished_this_manifold:
                if worker_process.poll():
                    num_tries += 1
                    if num_tries > 3:
                        break
                    print(str(m) + ' : Worker process crashed. Retrying ('
                          + str(num_tries) + ' / 3)')
                    sys.stdout.flush()
                    worker_process = kickoff_worker()
                    nice_str = (str(m) + os.linesep).encode(encoding = 'UTF-8')
                    worker_process.stdin.write(nice_str)
                    worker_process.stdin.flush()


                # Output is currently the tail end of whatever the
                # last response. We want to preserve this, so split it
                # into lines and only delete lines for which something
                # makes sense
                output = output + get_next_output(worker_process.stdout)
                outputs = output.split(os.linesep)
                idx = 0
                while True:
                    if idx >= len(outputs):
                        break
                    if outputs[idx] == LCC_MAKING_PROGRESS:
                        signal.alarm(0)
                        signal.alarm(worker_args.timeout)
                        outputs = outputs[(idx + 1):]
                        idx = 0
                    elif outputs[idx] == LCC_WANT_MORE_DATA:
                        signal.alarm(0)
                        signal.alarm(worker_args.timeout)
                        finished_this_manifold = True
                        break
                    elif outputs[idx] == LCC_ENCOUNTERED_ERROR:
                        print(str(m) + ' : Encountered error.')
                        sys.stdout.flush()
                        outputs = outputs[(idx + 1):]
                        idx = 0
                    elif outputs[idx] == LCC_CANT_DEAL_WITH_OUTPUT_FILE:
                        print(str(m) + ' : Some kind of problem with output')
                        sys.stdout.flush()
                        sys.exit(1)
                    else:
                        idx += 1

                output = os.linesep.join(outputs)

        except TimeoutExpired:
            print(str(m) + ' : Caused worker to timeout')
            sys.stdout.flush()
            try:
                worker_process.terminate()
            except:
                pass
            worker_process = kickoff_worker()
