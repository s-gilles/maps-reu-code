#!/usr/bin/python

from __future__ import print_function
from snappy import *
from cypari import *
import errno
import fcntl
import os
import Queue
import re
import subprocess
import sys
import threading
import time
import traceback

# `configuration' variables between runs/machines
THREAD_NUM = 16
SNAP_PATH='/usr/local/bin/snap'
CENSUS=LinkExteriors[:40]
TRIG_PATH='./triangulations/linkexteriors'

# Each snap process is managed by one python thread.
snap_process = [ ]
for i in range(0, THREAD_NUM):
    snap_process.append(None)

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
    f.write('Name,Tetrahedra,Volume,InvariantTraceField,InvariantTraceFieldDegree,NumberOfComplexPlaces,Disc\n')
    for poly,data in full_list.iteritems():
        dm = re.match('x\^([0-9]+).*', poly)
        deg = '0'
        if dm is not None:
            deg = dm.group(1)
        ncp = 0
        rlist = []
        for r in pari(poly).polroots():
            rlist.append(r)
        while len(rlist) > 0:
            for x in rlist[1:]:
                if x.real() == rlist[0].real() and x.imag() == -1*rlist[0].imag():
                    ncp += 1
                    rlist.remove(x)
                    break
            del rlist[0]
        disc = str(pari(poly).nfdisc())
        for vol, l in data.iteritems():
            for m in l:
                f.write('"' + m.name() + '",')
                f.write('"' + str(m.num_tetrahedra()) + '",')
                f.write('"' + vol + '",')
                f.write('"' + poly + '",')
                f.write('"' + deg + '",')
                f.write('"' + str(ncp) + '",')
                f.write('"' + disc + '"\n')
    f.close()


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
    snap_out = ''
    cprocess = subprocess.Popen([SNAP_PATH],
                                bufsize = 1,
                                shell = False,
                                stdin = subprocess.PIPE,
                                stdout = subprocess.PIPE,
                                stderr = subprocess.STDOUT)
    while not re.match('snap [0-9.]+', snap_out):
        snap_out = cprocess.stdout.readline().decode(encoding='UTF-8').rstrip()

    return cprocess

# Emulates typing string into a terminal running process. Note that since snap
# displays a prompt, entering '\n' causes snap's stdout to immediately echo the
# prompt plus string. This must be eaten, otherwise the process will hang.
def send_cmd(process, string):
    process.stdin.write((string + '\n').encode(encoding = 'UTF-8'))
    process.stdin.flush()
    process.stdout.readline()

# Merge a dictionary to full_list. This should be called in turn by each worker
# thread once the list of manifolds to analyze has been exhausted
def merge_up_dict(local_dict):
    full_list_lock.acquire()
    for polynomial_str, local_polydict in local_dict.items():
        fulls_polydict = full_list.setdefault(polynomial_str, dict())
        for vol, manifolds in local_polydict.items():
            fulls_manifolds = fulls_polydict.setdefault(vol, list())
            fulls_manifolds.extend(manifolds)
    full_list_lock.release()

def compute_shape_fields(idx):
    local_dict = dict()
    fname = 'tmp' + str(idx) + '.trig'
    while True:

        # This could either be the manifold itself, or the name of a file
        # directly. It should be pretty easy to switch over once the full
        # triangulation files are saved - I'm just using temp files until they
        # get uploaded.
        manifold, sig = ready_manifolds.get()

        if sig:
            merge_up_dict(local_dict)
            return
        if os.path.isfile(TRIG_PATH+"/"+manifold.name()+".trig"):
            dname = TRIG_PATH+"/"+manifold.name()+".trig"
        else:
            dname = fname
            manifold.save(fname)
        while True:
            try:
                send_cmd(snap_process[idx], 'read file ' + dname)
                send_cmd(snap_process[idx], 'compute shape')
                break
            except IOError:
                print(manifold.name() + ' crashed snap! [but we think we can fix it]')
                snap_process[idx].terminate()
                snap_process[idx] = kickoff_snap()
                continue

        while True:
            snap_out = snap_process[idx].stdout.readline().decode(encoding = 'UTF-8')
            if re.match('.*not found.*', snap_out):
                break
            elif re.match('.*Error.*', snap_out):
                print(manifold.name() + ' crashed snap!')
                snap_process[idx].terminate()
                snap_process[idx] = kickoff_snap()
                break

            shape_match = re.match('Shape field: ([-+*x0-9^]+) .*', snap_out)
            if shape_match is not None:
                vol = str(manifold.volume())
                polynomial = shape_match.group(1).strip()
                dm = re.match('x\^([0-9]+).*', polynomial)
                degree = 0
                if dm is not None:
                    degree = int(dm.group(1))
                if degree > 8:
                    break

                by_poly = local_dict.setdefault(polynomial, dict())
                by_volume = by_poly.setdefault(vol, list())
                by_volume.append(manifold)
                break
        ready_manifolds.task_done()

def worker_action(idx):
    snap_process[idx] = kickoff_snap()
    compute_shape_fields(idx)


if __name__ == "__main__":

    # Fiddle about with waiting for workers to startup
    print('Initializing workers')
    worker_threads = dict()
    all_started_successfully = False
    while not all_started_successfully:
        all_started_successfully = True
        for i in range(0, THREAD_NUM):
            if snap_process[i] is None:
                all_started_successfully = False
                new_thread = threading.Thread(group = None, target = worker_action, args = (i,))
                new_thread.daemon = True # Otherwise the ones that die ruin everything
                new_thread.start()
                worker_threads[i] = new_thread
        if not all_started_successfully:
            time.sleep(0.01)

    # To change this out for triangulation files, send in (pathname, False)
    print('Collecting manifold information')
    for m in CENSUS:
        ready_manifolds.put((m, False))

    # Make sure the magic join flags are there
    ready_manifolds.join()
    for i in range(0, THREAD_NUM):
        ready_manifolds.put((None, True))

    for i in range(0, THREAD_NUM):
        worker_threads[i].join()

    print('Coallating output')

    write_dict_to_output()
