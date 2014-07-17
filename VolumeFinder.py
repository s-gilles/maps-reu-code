#!/usr/bin/python

from __future__ import print_function
from snappy import *
from cypari import *
from ManifoldIterators import *
import Queue
import code
import errno
import fcntl
import itertools
import os
import pdb
import re
import signal
import subprocess
import sys
import threading
import time
import traceback

# `configuration' variables between runs/machines
CENSUS_CHUNK_SIZE = 50
SNAP_PATH='/usr/local/bin/snap'
TRIG_PATH='./triangulations/linkexteriors'

SIG_FINISH = 1
SIG_DIE = 2
SIG_MERGE = 3

ACT_DISTRUBUTE_WORK = 10
ACT_COLLECT = 11
ACT_COLLECT_THEN_DIE = 12
ACT_DIE = 14
main_action = ACT_DISTRUBUTE_WORK

RE_INV_TRACE_FIELD_NOT_FOUND = re.compile('.*Invariant trace field not found.*', re.DOTALL)
RE_FUNC_REQ_GROUP = re.compile('.*Function requires a group.*', re.DOTALL)
RE_ERROR = re.compile('.*Error.*', re.DOTALL)
RE_TRACE_FIELD = re.compile('.*Invariant trace field: ([-+*x0-9^ ]+) \[[0-9]+, *([0-9]+)\] [-0-9]+ R\([-0-9]+\) = ([-+*0-9.I]+).*', re.DOTALL)

SOL_TYPE_STRINGS = ['not_attempted', 'geometric', 'nongeometric', 'flat', 'degenerate', 'unrecognized', 'none_found']

pari.set_real_precision(100)

# Paired with snap process creation to prevent signal propagation
def preexec_discard_signals():
    os.setpgrp()

# Each snap process is managed by one python thread.
_snap_process = [ ]

# A thread-safe queue for distributing manifolds to python threads, that work
# may be done on them with snap
_ready_manifolds = Queue.Queue()

# A thread-safe queue for distributing messages from the python threads to the
# main thread. e.g. Main thread sends SIG_MERGE to all workers, wishes to know
# when the merge actually occurs.
_worker_to_main_messages = Queue.Queue()

# A list of all known hyperbolic volumes, organized by shape. The keys of this
# dictionary are polynomials defining shape fields (e.g. x^2 - x + 1), and the
# elements are lists of dictionaries. THESE dictionaries hav keys of hyperbolic
# volumes arising from that shape field, and elements of lists all known
# manifolds that exhibit that volume.
#
# e.g. _full_list['x^3 - x^2 - 3'] = {'-0.43 - 1.19*I' => [L12n2018, ...],
#                                     '-0.43 + 1.19*I' => L11n371, ...],
#                                     ...
#                                    }
_full_list = dict()
_full_list_lock = threading.Lock()

_discriminants = dict()
_discriminants_lock = threading.Lock()

def write_dict_to_output(output_filename = 'output.csv',  first_time = True, separator = ';'):
    global _full_list, _full_list_lock
    with _full_list_lock:
        if first_time:
            f = open(output_filename, 'w')
            f.write('Name' + separator +
                    'InvariantTraceField' + separator +
                    'Root' + separator +
                    'NumberOfComplexPlaces' + separator +
                    'Volume' + separator +
                    'InvariantTraceFieldDegree' + separator +
                    'SolutionType' + separator +
                    'Disc' + separator +
                    'DiscFactors' + separator +
                    'Tetrahedra\n')
        else:
            f = open(output_filename, 'a')

        with _discriminants_lock:
            for poly,data in sorted(_full_list.items()):
                dm = re.match('x\^([0-9]+).*', poly)
                deg = '0'
                if dm is not None:
                    deg = dm.group(1)

                if poly in _discriminants:
                    disc_str, disc_fact_str = _discriminants[poly]
                else:
                    disc = pari(poly).nfdisc()
                    disc_str = str(disc)
                    disc_fact_str = ''
                    try:
                        for p, e in disc.factor().mattranspose():
                            if e == 1:
                                disc_fact_str = disc_fact_str + str(p) + '*'
                            else:
                                disc_fact_str = disc_fact_str + str(p) + '^' + str(e) + '*'

                        if disc_fact_str == '':
                            disc_fact_str = '1'
                        else:
                            disc_fact_str = disc_fact_str[:-1]
                    except ValueError:
                        disc_fact_str = disc_str

                    _discriminants[poly] = (disc_str, disc_fact_str)

                for vol, l in sorted(data.items()):
                    for rec in l:
                        #unpack tuples
                        m = rec[0]
                        ncp = rec[1]
                        root = rec[2]
                        sol_type = rec[3]
                        f.write('"' + str(m) + '"' + separator)
                        f.write('"' + poly + '"' + separator)
                        f.write('"' + root + '"' + separator)
                        f.write('"' + str(ncp) + '"' + separator)
                        f.write('"' + vol + '"' + separator)
                        f.write('"' + deg + '"' + separator)
                        f.write('"' + sol_type + '"' + separator)
                        f.write('"' + disc_str + '"' + separator)
                        f.write('"' + disc_fact_str + '"' + separator)
                        f.write('"' + str(m.num_tetrahedra()) + '"\n')
        _full_list = dict()
        f.close()

def drain(out):
    result = ''
    incremental = ''
    while True:
        try:
            incremental = os.read(out.fileno(), 1024)
            if incremental == '':
                return result
            result = result + incremental
        except OSError:
            return result

# Initializes a snap process, waits for it to print out some preliminary
# information, then returns that snap process. If the snap process fails at
# all, this method will not return in a timely fashion. Therefore, all threads
# should perform
#
#     _snap_process[my_unique_thread_id] = kickoff_snap()
#
# as their first action. After a suitable waiting period (say 0.1s), those
# threads corresponding to id's for which _snap_process[id] is None are
# definitely hung, and may be restarted
def kickoff_snap(temp_file_dir):
    cprocess = subprocess.Popen([SNAP_PATH],
                                bufsize = 1,
                                shell = False,
                                stdin = subprocess.PIPE,
                                stdout = subprocess.PIPE,
                                stderr = subprocess.STDOUT,
                                preexec_fn = preexec_discard_signals)

    fcntl.fcntl(cprocess.stdout, fcntl.F_SETFL,
                fcntl.fcntl(cprocess.stdout, fcntl.F_GETFL) | os.O_NONBLOCK)

    drain(cprocess.stdout)

    # No output is expected after these
    send_cmd(cprocess, 'set precision 10\n')
    send_cmd(cprocess, 'set digits_printed 100 f\n')
    send_cmd(cprocess, 'set degree 8\n')
    send_cmd(cprocess, 'set simplify false\n')
    send_cmd(cprocess, 'set path ' + temp_file_dir + '\n')
    return cprocess

# Emulates typing string into a terminal running process. Note that since snap
# displays a prompt, entering '\n' causes snap's stdout to immediately echo the
# prompt plus string. This must be eaten, otherwise the process will hang.
def send_cmd(process, string):
    process.stdin.write((string + '\n').encode(encoding = 'UTF-8'))
    process.stdin.flush()
    return drain(process.stdout)

# Merge a dictionary to _full_list. This should be called in turn by each worker
# thread once the list of manifolds to analyze has been exhausted
def merge_up_dict(local_dict):
    with _full_list_lock:
        for polynomial_str, local_polydict in local_dict.items():
            fulls_polydict = _full_list.setdefault(polynomial_str, dict())
            for vol, manifold_records in local_polydict.items():
                fulls_manifolds = fulls_polydict.setdefault(vol, list())
                fulls_manifolds.extend(manifold_records)

# The action that a worker thread takes.  Simply read from the Queue and
# perform computations.
def compute_shape_fields(idx, temp_file_dir):
    global SIG_FINISH, SIG_DIE, SIG_MERGE
    global _snap_process
    global _worker_to_main_messages
    global SOL_TYPE_STRINGS
    local_dict = dict()
    fname = 'tmp_' + str(os.getpid()) + '_' + str(idx) + '.trig'
    full_fname = temp_file_dir + fname
    snap_output = ''
    while True:

        # This could either be the manifold itself, or the name of a file
        # directly. It should be pretty easy to switch over once the full
        # triangulation files are saved - I'm just using temp files until they
        # get uploaded.
        manifold, sig = _ready_manifolds.get()

        if sig is SIG_FINISH:
            merge_up_dict(local_dict)
            _ready_manifolds.task_done()
            break
        elif sig is SIG_DIE:
            _ready_manifolds.task_done()
            break
        elif sig is SIG_MERGE:
            merge_up_dict(local_dict)
            local_dict = dict()
            _ready_manifolds.task_done()
            _worker_to_main_messages.put((SIG_MERGE,))
            continue


        sol_type = SOL_TYPE_STRINGS[-1]
        sol_enum = int(manifold.solution_type(enum = True))
        try:
            sol_type = SOL_TYPE_STRINGS[sol_enum]
        except:
            pass

        # Bail if not 'geometric' or 'non-geometric' or 'flat'
        if sol_enum != 1 and sol_enum != 2 and sol_enum != 3:
            _ready_manifolds.task_done()
            continue

        # If non-geometric, spend a bit of time trying to make it geometric
        randomize_attempts = 0
        while manifold.solution_type(enum = True) == 2 and randomize_attempts < 16:
            manifold.randomize()
            randomize_attempts += 1


        manifold.save(full_fname)
        while True:
            try:
                send_cmd(_snap_process[idx], 'read file ' + fname)
                if _snap_process[idx].poll() is not None:
                    raise IOError

                snap_output = send_cmd(_snap_process[idx], 'compute invariant_trace_field')
                if _snap_process[idx].poll() is not None:
                    raise IOError

                # due to O_NONBLOCK and drain(), IOError doesn't show up until
                # the call AFTER the killer
                break
            except IOError:
                print(str(manifold) + ' crashed snap in ' + str(idx) + '!')
                try:
                    _snap_process[idx].terminate()
                except:
                    pass
                _snap_process[idx] = kickoff_snap(temp_file_dir)
                continue

        times_retried = 0
        while True:
            snap_output = snap_output + drain(_snap_process[idx].stdout)
            if RE_INV_TRACE_FIELD_NOT_FOUND.match(snap_output):
                break
            elif RE_FUNC_REQ_GROUP.match(snap_output):
                break
            elif RE_ERROR.match(snap_output):
                print(str(manifold) + ' crashed snap in ' + str(idx) + '!')
                try:
                    _snap_process[idx].terminate()
                except:
                    pass
                while True:
                    try:
                        _snap_process[idx] = kickoff_snap(temp_file_dir)
                        break
                    except Exception as e:
                        print(str(idx) + ' attempting to recover [' + str(e) + '] ...')
                        time.sleep(1.5)
                print(str(idx) + ' recovered')
                break

            trace_match = RE_TRACE_FIELD.match(snap_output)
            if trace_match is not None:
                try:
                    vol = str(manifold.high_precision().volume())
                except:
                    print(str(manifold) + ' crashed pari while determining volume in ' + str(idx) + '!')
                    break
                polynomial = trace_match.group(1).strip()
                dm = re.match('x\^([0-9]+).*', polynomial)
                degree = 0
                ncp = int(trace_match.group(2).strip())
                root = trace_match.group(3).strip()

                if dm is not None:
                    degree = int(dm.group(1))
                if degree <= 8:
                    by_poly = local_dict.setdefault(polynomial, dict())
                    by_volume = by_poly.setdefault(vol, list())
                    # All the information to be sent back from the threads is packed in a tuple:
                    by_volume.append((manifold,ncp,root, sol_type))
                break
            time.sleep(0.5)
            try:
                snap_output = snap_output + send_cmd(_snap_process[idx], '')
            except IOError:
                # Force the next RE_ERROR.match to trigger
                snap_output = 'Error'

            # This happens extremely rarely, but occasionally snap
            # will simply hang.  E.g. computing invariant trace field
            # for 10^3_3(-1,8)(1,3)(1,3)
            times_retried += 1
            if times_retried > 20:
                snap_output = 'Error'

        _ready_manifolds.task_done()


    if _snap_process[idx] is not None:
        send_cmd(_snap_process[idx], 'quit')
    _snap_process[idx] = None

    return

# Wrapper around compute_shape_fields, in case any one-time startup/shutdown
# code needs to be applied.
def worker_action(idx, temp_file_dir= '/tmp/'):
    _snap_process[idx] = kickoff_snap(temp_file_dir)
    compute_shape_fields(idx, temp_file_dir)
    return

def double_sigint_handler(signum, frame):
    print()
    sys.exit(0)

def sigint_handler(signum, frame):
    global main_action
    main_action = ACT_COLLECT_THEN_DIE
    print('\nDying. Please wait for worker threads to terminate. (^C again REALLY kills)')
    signal.signal(signal.SIGINT, double_sigint_handler)

def sigusr2_handler(sig, frame):
    global main_action
    main_action = ACT_COLLECT
    print('Writing out current progress. Please wait.')

def sigusr1_handler(sig, frame):
    id2name = dict([(th.ident, th.name) for th in threading.enumerate()])
    code = []
    for threadId, stack in sys._current_frames().items():
        code.append("\n# Thread: %s(%d)" % (id2name.get(threadId,""), threadId))
        for filename, lineno, name, line in traceback.extract_stack(stack):
            code.append('File: "%s", line %d, in %s' % (filename, lineno, name))
            if line:
                code.append("  %s" % (line.strip()))
    print("\n".join(code))

    d={'_frame':frame}         # Allow access to frame object.
    d.update(frame.f_globals)  # Unless shadowed by global
    d.update(frame.f_locals)

def begin_collection(iterator, output_filename = 'output.csv', thread_num = 12,
                     install_signal_handlers = True, is_appending = False, csv_separator = ';'):
    """Call this, given a batch iterator, to exhaust that batch iterator and
store the result to output_filename.  Example:

  beginCollection(BatchIterator(TorusBundleIterator), 50)

will set up the default thread state.

Optional parameters:

  output_filename (default = 'output.csv'). A filename to which to write output.

  install_signal_handlers (default = True). If set, the program will install signal handlers for the following:
    SIGINT will attempt to gracefully stop threads, gather their work, and close the program
    SIGINT twice will close the program
    SIGUSR1 will give a stack trace, without otherwise affecting operations
    SIGUSR2 will write out partial progress (since last write) to output_filename so that it can be inspected. Note thta if this happens, the program will switch to appending mode (see below) to avoid erasing its own results.  This is highly useful.

  thread_num (default = 12). Use this many worker threads.  It is recommended to be set to the number of logical cores on the system.

  is_appending (default = False).  If set to True, the program will assume that output_filename already contains the results of some previous run, and will not overwrite them.  If set to False (the default!), the program will completely overwrite output_filename with its own results.

  csv_separator (default = ';').  Specifies the separator used in writing out csv files."""
    global CENSUS_CHUNK_SIZE
    global SNAP_PATH
    global TRIG_PATH

    global SIG_FINISH
    global SIG_DIE
    global SIG_MERGE

    global ACT_DISTRUBUTE_WORK
    global ACT_COLLECT
    global ACT_COLLECT_THEN_DIE
    global ACT_DIE
    global main_action
    global _snap_process

    pari.set_real_precision(100)

    # Fiddle about with waiting for workers to startup
    print('Initializing...')
    worker_threads = list()
    _snap_process = list()
    for i in range(0, thread_num):
        _snap_process.append(None)

    for i in range(0, thread_num):
        new_thread = threading.Thread(group = None, target = worker_action, args = (i,))
        new_thread.daemon = True
        new_thread.start()
        while _snap_process[i] is None:
            time.sleep(1) #LMOD
        worker_threads.append(new_thread)

    have_written_out_already = is_appending

    # To trigger these, find the PID and issue
    #
    #     kill -s SIGUSR1 $PID
    #
    # or such.
    if install_signal_handlers:
        signal.signal(signal.SIGINT, sigint_handler)
        signal.signal(signal.SIGUSR1, sigusr1_handler)
        signal.signal(signal.SIGUSR2, sigusr2_handler)

    print('Working...')

    main_action = ACT_DISTRUBUTE_WORK
    next_batch = iterator.next_batch()
    last_pulled = next_batch[-1]
    while True:
        if main_action is ACT_DIE:
            with _ready_manifolds.mutex:
                _ready_manifolds.queue.clear()
            for i in range(0, thread_num):
                _ready_manifolds.put((None, SIG_DIE))
            break
        elif main_action is ACT_COLLECT:
            _worker_to_main_messages.queue.clear()
            for i in range(0, thread_num):
                _ready_manifolds.put((None, SIG_MERGE))
            for i in range(0, thread_num):
                _worker_to_main_messages.get()
            write_dict_to_output(output_filename, first_time = not have_written_out_already,
                                 separator = csv_separator)
            have_written_out_already = True
            print('Wrote out current progress.')
            try:
                print('That was for manifolds up to (possibly not including): ' + str(last_pulled))
            except:
                pass
            main_action = ACT_DISTRUBUTE_WORK
            continue
        elif main_action is ACT_COLLECT_THEN_DIE:
            for i in range(0, thread_num):
                _ready_manifolds.put((None, SIG_FINISH))
            break
        elif main_action is ACT_DISTRUBUTE_WORK:
            try:
                last_pulled = next_batch[-1]
                for m in next_batch:
                    _ready_manifolds.put((m, None))
                next_batch = iterator.next_batch()
                main_action = ACT_COLLECT
            except StopIteration:
                print('Done.')
                main_action = ACT_COLLECT_THEN_DIE

        # Would like to .join() here, but must use sleep() for signaling, or
        # some hogwash
        while not _ready_manifolds.empty():
            time.sleep(0.05)

    # Ditto, can't .join()
    while any(w.is_alive() for w in worker_threads):
        time.sleep(0.05)

    write_dict_to_output(output_filename, first_time = not have_written_out_already,
                         separator = csv_separator)
