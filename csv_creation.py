#!/usr/bin/python

from __future__ import print_function
from snappy import *
from cypari import *
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
THREAD_NUM = 8
CENSUS_CHUNK_SIZE = 50
SNAP_PATH='/usr/local/bin/snap'
TRIG_PATH='./triangulations/linkexteriors'

ACT_DIE = 'die now'
ACT_CENSUS = 'continually dump from census'
main_action = ACT_CENSUS

SIG_FINISH = 1
SIG_DIE = 2

RE_INV_TRACE_FIELD_NOT_FOUND = re.compile('.*Invariant trace field not found.*', re.DOTALL)
RE_FUNC_REQ_GROUP = re.compile('.*Function requires a group.*', re.DOTALL)
RE_ERROR = re.compile('.*Error.*', re.DOTALL)
RE_TRACE_FIELD = re.compile('.*Invariant trace field: ([-+*x0-9^]+) \[[0-9]+,([0-9]+)\] [-0-9]+ R\([-0-9]+\) = ([-+*0-9.I]+).*', re.DOTALL)

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

def write_dict_to_output():
    f = open('output.csv', 'w')
    f.write('Name,Tetrahedra,Volume,InvariantTraceField,InvariantTraceFieldDegree,Root,NumberOfComplexPlaces,Disc\n')
    for poly,data in sorted(full_list.items()):
        dm = re.match('x\^([0-9]+).*', poly)
        deg = '0'
        if dm is not None:
            deg = dm.group(1)
        disc = str(pari(poly).nfdisc())
        for vol, l in sorted(data.items()):
            for rec in l:
                #unpack tuples
                m = rec[0]
                ncp = rec[1]
                root = rec[2]
                f.write('"' + m.name() + '",')
                f.write('"' + str(m.num_tetrahedra()) + '",')
                f.write('"' + vol + '",')
                f.write('"' + poly + '",')
                f.write('"' + deg + '",')
                f.write('"' + rec[2] + '",')
                f.write('"' + str(rec[1]) + '",')
                f.write('"' + disc + '"\n')
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

def compute_shape_fields(idx):
    global SIG_FINISH, SIG_DIE
    local_dict = dict()
    fname = 'tmp' + str(idx) + '.trig'
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

        if os.path.isfile(TRIG_PATH+"/"+manifold.name()+".trig"):
            dname = TRIG_PATH+"/"+manifold.name()+".trig"
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
                print(manifold.name() + ' crashed snap in ' + str(idx) + '!')
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
                print(manifold.name() + ' crashed snap in ' + str(idx) + '!')
                snap_process[idx].terminate()
                snap_process[idx] = kickoff_snap()
                print(str(idx) + ' recovered')
                break

            trace_match = RE_TRACE_FIELD.match(snap_output)
            if trace_match is not None:
                vol = str(manifold.volume())
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

def worker_action(idx):
    snap_process[idx] = kickoff_snap()
    compute_shape_fields(idx)
    return

def double_sigint_handler(signum, frame):
    print()
    sys.exit(0)

def sigint_handler(signum, frame):
    global main_action
    main_action = ACT_DIE
    print('\nDying. Please wait for worker threads to terminate. (^C again REALLY kills)')
    signal.signal(signal.SIGINT, double_sigint_handler)

def sigusr1_handler(sig, frame):
    """Interrupt running process, and provide a python prompt for
    interactive debugging."""

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


if __name__ == "__main__":

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

    signal.signal(signal.SIGINT, sigint_handler)
    signal.signal(signal.SIGUSR1, sigusr1_handler)

    census_chunks_ocm = iter(LinkExteriors)

    while True:
        if main_action is ACT_DIE:
            with ready_manifolds.mutex:
                ready_manifolds.queue.clear()
            for i in range(0, THREAD_NUM):
                ready_manifolds.put((None, SIG_DIE))
            break
        elif main_action is ACT_CENSUS:
            try:
                for i in range(0,CENSUS_CHUNK_SIZE):
                    ready_manifolds.put((census_chunks_ocm.next(), None))
            except:
                print('Done. Finishing up.')
                for i in range(0, THREAD_NUM):
                    ready_manifolds.put((None, SIG_FINISH))
                break

        # Would like to .join() here, but must use sleep() for signaling, or
        # some hogwash
        while (not ready_manifolds.empty()):
            time.sleep(0.05)

    # Ditto, can't .join()
    while any(w.is_alive() for w in worker_threads):
        time.sleep(0.05)

    write_dict_to_output()
