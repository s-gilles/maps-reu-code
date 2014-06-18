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
THREAD_NUM = 12
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
RE_TRACE_FIELD = re.compile('.*Invariant trace field: ([-+*x0-9^]+) \[[0-9]+,([0-9]+)\] [-0-9]+ R\([-0-9]+\) = ([-+*0-9.I]+).*', re.DOTALL)

pari.set_real_precision(100)

# Each snap process is managed by one python thread.
snap_process = [ ]
for i in range(0, THREAD_NUM):
    snap_process.append(None)

# Paired with snap process creation to prevent signal propagation
def preexec_discard_signals():
    os.setpgrp()

# A thread-safe queue for distributing manifolds to python threads, that work
# may be done on them with snap
ready_manifolds = Queue.Queue()

# A thread-safe queue for distributing messages from the python threads to the
# main thread. e.g. Main thread sends SIG_MERGE to all workers, wishes to know
# when the merge actually occurs.
worker_to_main_messages = Queue.Queue()

# A list of all known hyperbolic volumes, organized by shape. The keys of this
# dictionary are polynomials defining shape fields (e.g. x^2 - x + 1), and the
# elements are lists of dictionaries. THESE dictionaries hav keys of hyperbolic
# volumes arising from that shape field, and elements of lists all known
# manifolds that exhibit that volume.
#
# e.g. full_list['x^3 - x^2 - 3'] = {'-0.43 - 1.19*I' => [L12n2018, ...],
#                                    '-0.43 + 1.19*I' => L11n371, ...],
#                                    ...
#                                   }
full_list = dict()
full_list_lock = threading.Lock()

def write_dict_to_output(output_filename = 'output.csv',  first_time = True):
    global full_list, full_list_lock
    with full_list_lock:
        if first_time:
            f = open(output_filename, 'w')
            f.write('Name,Tetrahedra,Volume,InvariantTraceField,InvariantTraceFieldDegree,Root,NumberOfComplexPlaces,Disc,DiscFactors\n')
        else:
            f = open(output_filename, 'a')

        for poly,data in sorted(full_list.items()):
            dm = re.match('x\^([0-9]+).*', poly)
            deg = '0'
            if dm is not None:
                deg = dm.group(1)
            disc = pari(poly).nfdisc()
            disc_str = str(disc)
            disc_fact_str = ''
            try:
                for p, e in disc.factor().mattranspose():
                    disc_fact_str = disc_fact_str + str(p) + '^' + str(e) + '*'
                disc_fact_str = disc_fact_str[:-1]
            except ValueError:
                disc_fact_str = disc_str

            for vol, l in sorted(data.items()):
                for rec in l:
                    #unpack tuples
                    m = rec[0]
                    ncp = rec[1]
                    root = rec[2]
                    f.write('"' + str(m) + '",')
                    f.write('"' + str(m.num_tetrahedra()) + '",')
                    f.write('"' + vol + '",')
                    f.write('"' + poly + '",')
                    f.write('"' + deg + '",')
                    f.write('"' + rec[2] + '",')
                    f.write('"' + str(rec[1]) + '",')
                    f.write('"' + disc_str + '",')
                    f.write('"' + disc_fact_str + '"\n')
        full_list = dict()
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
#     snap_process[my_unique_thread_id] = kickoff_snap()
#
# as their first action. After a suitable waiting period (say 0.1s), those
# threads corresponding to id's for which snap_process[id] is None are
# definitely hung, and may be restarted
def kickoff_snap():
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

    return cprocess

# Emulates typing string into a terminal running process. Note that since snap
# displays a prompt, entering '\n' causes snap's stdout to immediately echo the
# prompt plus string. This must be eaten, otherwise the process will hang.
def send_cmd(process, string):
    process.stdin.write((string + '\n').encode(encoding = 'UTF-8'))
    process.stdin.flush()
    return drain(process.stdout)

# Merge a dictionary to full_list. This should be called in turn by each worker
# thread once the list of manifolds to analyze has been exhausted
def merge_up_dict(local_dict):
    with full_list_lock:
        for polynomial_str, local_polydict in local_dict.items():
            fulls_polydict = full_list.setdefault(polynomial_str, dict())
            for vol, manifold_records in local_polydict.items():
                fulls_manifolds = fulls_polydict.setdefault(vol, list())
                fulls_manifolds.extend(manifold_records)

# The action that a worker thread takes.  Simply read from the Queue and
# perform computations.
def compute_shape_fields(idx):
    global SIG_FINISH, SIG_DIE, SIG_MERGE
    global snap_process
    global worker_to_main_messages
    local_dict = dict()
    fname = 'tmp_' + str(os.getpid()) + '_' + str(idx) + '.trig'
    snap_output = ''
    while True:

        # This could either be the manifold itself, or the name of a file
        # directly. It should be pretty easy to switch over once the full
        # triangulation files are saved - I'm just using temp files until they
        # get uploaded.
        manifold, sig = ready_manifolds.get()

        if sig is SIG_FINISH:
            merge_up_dict(local_dict)
            ready_manifolds.task_done()
            break
        elif sig is SIG_DIE:
            ready_manifolds.task_done()
            break
        elif sig is SIG_MERGE:
            merge_up_dict(local_dict)
            local_dict = dict()
            ready_manifolds.task_done()
            worker_to_main_messages.put((SIG_MERGE,))
            continue

        if os.path.isfile(TRIG_PATH+"/"+str(manifold)+".trig"):
            dname = TRIG_PATH+"/"+str(manifold)+".trig"
        else:
            dname = fname
            manifold.save(fname)
        while True:
            try:
                send_cmd(snap_process[idx], 'read file ' + dname)
                if snap_process[idx].poll() is not None:
                    raise IOError

                snap_output = send_cmd(snap_process[idx], 'compute invariant_trace_field')
                if snap_process[idx].poll() is not None:
                    raise IOError

                # due to O_NONBLOCK and drain(), IOError doesn't show up until
                # the call AFTER the killer
                break
            except IOError:
                print(str(manifold) + ' crashed snap in ' + str(idx) + '!')
                snap_process[idx].terminate()
                snap_process[idx] = kickoff_snap()
                continue

        while True:
            snap_output = snap_output + drain(snap_process[idx].stdout)
            if RE_INV_TRACE_FIELD_NOT_FOUND.match(snap_output):
                break
            elif RE_FUNC_REQ_GROUP.match(snap_output):
                break
            elif RE_ERROR.match(snap_output):
                print(str(manifold) + ' crashed snap in ' + str(idx) + '!')
                snap_process[idx].terminate()
                snap_process[idx] = kickoff_snap()
                print(str(idx) + ' recovered')
                break

            trace_match = RE_TRACE_FIELD.match(snap_output)
            if trace_match is not None:
                vol = str(manifold.high_precision().volume())
                polynomial = trace_match.group(1).strip()
                dm = re.match('x\^([0-9]+).*', polynomial)
                degree = 0
                ncp = int(trace_match.group(2).strip())
                root = trace_match.group(3).strip()
                if dm is not None:
                    degree = int(dm.group(1))
                if degree <= 24:
                    by_poly = local_dict.setdefault(polynomial, dict())
                    by_volume = by_poly.setdefault(vol, list())
                    # All the information to be sent back from the threads is packed in a tuple:
                    by_volume.append((manifold,ncp,root))
                break
            time.sleep(0.05)
            try:
                snap_output = snap_output + send_cmd(snap_process[idx], '')
            except IOError:
                # Force the next RE_ERROR.match to trigger
                snap_output = 'Error'

        ready_manifolds.task_done()


    if snap_process[idx] is not None:
        send_cmd(snap_process[idx], 'quit')
    snap_process[idx] = None

    return

# Wrapper around compute_shape_fields, in case any one-time startup/shutdown
# code needs to be applied.
def worker_action(idx):
    snap_process[idx] = kickoff_snap()
    compute_shape_fields(idx)
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

def beginCollection(iterator, output_filename = 'output.csv'):
    """Call this, given a batch iterator, to exhaust that batch iterator and
store the result to output_filename.  Example:

  beginCollection(BatchIterator(TorusBundleIterator), 50)

will set up the default thread state."""
    global THREAD_NUM
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
    global snap_process

    pari.set_real_precision(100)

    # Fiddle about with waiting for workers to startup
    print('Initializing...')
    worker_threads = list()
    for i in range(0, THREAD_NUM):
        new_thread = threading.Thread(group = None, target = worker_action, args = (i,))
        new_thread.daemon = True
        new_thread.start()
        while snap_process[i] is None:
            time.sleep(0.1)
        worker_threads.append(new_thread)

    have_written_out_already = False

    # To trigger these, find the PID and issue
    #
    #     kill -s SIGUSR1 $PID
    #
    # or such.
    signal.signal(signal.SIGINT, sigint_handler)
    signal.signal(signal.SIGUSR1, sigusr1_handler)
    signal.signal(signal.SIGUSR2, sigusr2_handler)

    print('Working...')

    while True:
        signal.signal(signal.SIGINT, sigint_handler)
        if main_action is ACT_DIE:
            with ready_manifolds.mutex:
                ready_manifolds.queue.clear()
            for i in range(0, THREAD_NUM):
                ready_manifolds.put((None, SIG_DIE))
            break
        elif main_action is ACT_COLLECT:
            worker_to_main_messages.queue.clear()
            for i in range(0, THREAD_NUM):
                ready_manifolds.put((None, SIG_MERGE))
            for i in range(0, THREAD_NUM):
                worker_to_main_messages.get()
            write_dict_to_output(output_filename, not have_written_out_already)
            have_written_out_already = True
            print('Wrote out current progress.')
            main_action = ACT_DISTRUBUTE_WORK
            continue
        elif main_action is ACT_COLLECT_THEN_DIE:
            for i in range(0, THREAD_NUM):
                ready_manifolds.put((None, SIG_FINISH))
            break
        elif main_action is ACT_DISTRUBUTE_WORK:
            try:
                for m in iterator.next_batch():
                    ready_manifolds.put((m, None))
            except StopIteration:
                print('Done.')
                main_action = ACT_COLLECT_THEN_DIE

        # Would like to .join() here, but must use sleep() for signaling, or
        # some hogwash
        while not ready_manifolds.empty():
            time.sleep(0.05)

    # Ditto, can't .join()
    while any(w.is_alive() for w in worker_threads):
        time.sleep(0.05)

    write_dict_to_output(output_filename, not have_written_out_already)
