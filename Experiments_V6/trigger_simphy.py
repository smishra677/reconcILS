"""SimPhy comparison using the workflow in test_reconcILS_default_40.py.

Reads g_trees.tre (one Newick tree per line), a SimPhy SQLite database, and
1/s_tree.trees from --data. Gene-tree lines must follow n_ltree/n_gtree order.
Supports one species-tree replicate. Output headers match test_reconcILS_default_40.py.
SimPhy Extra_l populates DLCILS; unavailable true NNI/hemiplasy values are blank.
Writes result, timing, memory, labeled species tree, and error files.
Uses reusable spawned workers; rerun the same command to resume checkpoints.

Example:
    python trigger_simphy.py --data /path/to/simphy_40_base
"""

import os
for _name in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS',
              'NUMEXPR_NUM_THREADS', 'BLIS_NUM_THREADS'):
    os.environ[_name] = '1'

import fcntl
import hashlib
import json
import multiprocessing as mp
import signal
import time
import argparse
import copy
import csv
import sqlite3
import sys
import traceback
import tracemalloc
from collections import Counter
from pathlib import Path
from time import perf_counter

HERE = Path(__file__).resolve().parent
# Experiments_V6 sits directly inside reconcILS/.
CODE_ROOT = HERE.parent
sys.path.insert(0, str(CODE_ROOT))
from reconcILS import *
from utils_reconcILS import *
from ete3 import PhyloTree

FIELDS = ['Process', 'Replicate', 'Gene_tree', 'Species_Tree', 'Duplication',
          'NNI', 'DLCILS', 'Loss', 'Hemiplasy', 'RHemiplasy']


def read_truth(path):
    conn = sqlite3.connect(path.resolve().as_uri() + '?mode=ro', uri=True)
    conn.row_factory = sqlite3.Row
    try:
        rows = conn.execute('''
            SELECT g.GID, g.LID, g.SID, g.n_ltree, g.n_gtree,
                   l.n_dup, l.n_loss, g.Extra_l
            FROM Gene_Trees g JOIN Locus_Trees l ON g.LID = l.LID
            ORDER BY g.SID, g.n_ltree, g.n_gtree
        ''').fetchall()
        count = conn.execute('SELECT COUNT(*) FROM Gene_Trees').fetchone()[0]
    finally:
        conn.close()
    if not rows or len(rows) != count or len({r['SID'] for r in rows}) != 1:
        raise ValueError('Expected one species replicate and a locus record for every gene tree')
    return rows


def measured(function, trace_memory=True):
    """Measure wall time and peak Python allocation, not total process memory."""
    if trace_memory:
        tracemalloc.start()
    start = perf_counter()
    try:
        result = function()
        elapsed = perf_counter() - start
        peak = tracemalloc.get_traced_memory()[1] / 1e6 if trace_memory else None
        return result, elapsed, peak
    finally:
        tracemalloc.stop()


def compare(newick, species_newick, truth, index, reader, reference, edges, trace_memory=True):
    gene = PhyloTree(newick, format=1,
                     sp_naming_function=lambda name: name.split('_')[0])
    if any(len(node.children) != 2 for node in gene.traverse() if not node.is_leaf()):
        raise ValueError('ReconcILS requires a binary gene tree')
    species = PhyloTree(species_newick, format=9,
                        sp_naming_function=lambda name: name)
    unknown = gene.get_species() - set(species.get_leaf_names())
    if unknown:
        raise ValueError(f'Gene tree contains unknown species: {sorted(unknown)}')
    topology = gene.write(format=9)
    base = dict(Replicate=index, Gene_tree=newick, Species_Tree=species_newick,
                DLCILS=0, Hemiplasy=0, RHemiplasy=0)
    rows = [dict(base, Process='True Process', Duplication=truth['n_dup'],
                 Loss=truth['n_loss'], NNI='', DLCILS=truth['Extra_l'],
                 Hemiplasy='', RHemiplasy='')]

    (reconciled, ete_events), ete_time, ete_mem = measured(lambda: gene.reconcile(species), trace_memory)
    ete_counts = Counter(reader.sp_event_ete3(reconciled, ete_events))
    rows.append(dict(base, Process='ETE3', Duplication=ete_counts['D'],
                     Loss=ete_counts['L'], NNI=0))
    ete_labels = reader.label_ete3(reconciled, {})

    sp = reader.parse(species_newick)
    tr = reader.parse(topology)
    reco = reconcils()
    sp_copy = copy.deepcopy(sp)
    sp_copy.reset()
    reco.gene_tree = copy.deepcopy(tr)
    tr.order_gene(sp)
    tr.label_internal()
    sp.label_internal()
    tr.map_gene(sp)
    reco.setCost(sp)
    #reco.I_cost=45
    sp.isRoot = tr.isRoot = sp_copy.isRoot = True
    events, rec_time, rec_mem = measured(
        lambda: reco.iterative_reconcILS(tr, sp, sp_copy, sp, []), trace_memory)
    counts = {key: [] for key in ('Process', 'Replicate', 'Duplication', 'Loss', 'NNI')}
    reader.Create_pd('reconcILS', index, events, counts)

    rec_labels = reader.label_reconcILS(copy.deepcopy(events))
    ete_labels = reader.find_replacement_for_edge(edges, ete_labels)
    rec_labels, ete_labels, _, ete_totals = reader.put_it_in_sp_tree(edges, rec_labels, ete_labels)
    # Match the consistency checks in test_reconcILS_default_40.py.
    if ete_totals['D'] != ete_counts['D'] or ete_totals['L'] != ete_counts['L']:
        raise ValueError('ETE3 branch labels disagree with its total event counts')
    labeled = copy.deepcopy(reference)
    reco.edge_to_event(labeled, copy.deepcopy(rec_labels), 0)
    reco.edge_to_event(labeled, copy.deepcopy(ete_labels), 1)
    rows.append(dict(base, Process='reconcILS', Gene_tree=reader.to_newick(reco.gene_tree),
                     Species_Tree=reader.to_newick_sp(labeled),
                     Duplication=counts['Duplication'][0], Loss=counts['Loss'][0],
                     NNI=counts['NNI'][0]))
    timing = dict(reconcILSTime_iterative=rec_time, ete_time=ete_time)
    memory = dict(reconcILS_mem=rec_mem, ete_mem=ete_mem)
    tied_nodes = 0
    stack = [sp]
    while stack:
        node = stack.pop()
        if node is not None:
            tied_nodes += node.equal
            stack.extend([node.leftChild, node.rightChild])
    equal = dict(gt=tr.to_newick(), count=tied_nodes)
    return rows, timing, memory, rec_labels, ete_labels, equal



# Reuse the same isolated worker and atomic output utilities as default40.
from test_reconcILS_default_40 import (
    worker_entry, add_events, positive_int, export_results as export_tables,
    write_text_atomic,
)


def process_tree(index, payload, species_newick, trace_memory=True):
    tree, truth = payload
    reader = readWrite.readWrite()
    reference = reader.parse(species_newick)
    reference.label_internal()
    edges = reconcils().get_edges(reference)
    rows, timing, memory, rec, ete, equal = compare(
        tree.strip(), species_newick, truth, index, reader, reference, edges, trace_memory)
    return dict(status='ok', replicate=index, rows=rows, equal=equal,
                rec_time=timing['reconcILSTime_iterative'], ete_time=timing['ete_time'],
                rec_mem=memory['reconcILS_mem'], ete_mem=memory['ete_mem'],
                edges_rec=[[list(k), rec[k]] for k in edges],
                edges_ete=[[list(k), ete[k]] for k in edges])


def export_results(output, records, suffix, aggregate):
    export_tables(output, records, '', aggregate)
    errors = ''.join(
        f"Line {r['replicate'] + 1}, GID {r.get('gid', '?')}\n{r.get('error', r['status'])}\n"
        for r in sorted(records.values(), key=lambda r: r['replicate']) if r['status'] != 'ok')
    write_text_atomic(output / 'run_error.txt', errors)
    write_text_atomic(output / 'reconcILS_labeled_sp.csv',
                      readWrite.readWrite().to_newick_sp(aggregate) + '\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--data', type=Path,
                        default=HERE / 'simphy_40_base')
    parser.add_argument('--db', type=Path, help='Database path; default: the single .db in --data')
    parser.add_argument('--species', type=Path, help='Species tree; default: DATA/1/s_tree.trees')
    parser.add_argument('--out', '--output-dir', type=Path, default=HERE / 'simphy_40_base_results')
    parser.add_argument('--limit', type=int, help='Process only the first N gene trees')
    allocated = int(os.environ.get('SLURM_CPUS_PER_TASK', '1'))
    parser.add_argument('--workers', type=positive_int, default=max(1, min(8, allocated)))
    parser.add_argument('--checkpoint-every', type=positive_int, default=100)
    parser.add_argument('--timeout', type=positive_int,
                        help='optional hard seconds per tree, including worker startup')
    parser.add_argument('--trace-memory', action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument('--retry-failed', action='store_true')
    parser.add_argument('--worker-logs', action='store_true')
    args = parser.parse_args()
    if args.limit is not None and args.limit < 1:
        parser.error('--limit must be positive')
    if args.db is None:
        databases = sorted(args.data.glob('*.db'))
        if len(databases) != 1:
            parser.error(f'Expected one .db in {args.data}; found {len(databases)}. Use --db.')
        args.db = databases[0]
    species_path = args.species or args.data / '1' / 's_tree.trees'
    truth = read_truth(args.db)
    trees = (args.data / 'g_trees.tre').read_text().splitlines()
    if len(trees) != len(truth):
        parser.error(f'{len(trees)} tree lines but {len(truth)} database records')
    if any(not tree.strip().endswith(';') for tree in trees):
        parser.error('Each gene-tree line must end with a semicolon; blank lines are not allowed')
    species = PhyloTree(species_path.read_text().strip(), format=1)
    if any(len(n.children) != 2 for n in species.traverse() if not n.is_leaf()):
        parser.error('Species tree must be binary')
    species_newick = species.write(format=9)
    total = min(args.limit or len(trees), len(trees))
    truth = [dict(row) for row in truth]
    args.out = args.out.resolve()
    args.out.mkdir(parents=True, exist_ok=True)
    state = args.out / '.simphy_state'
    state.mkdir(exist_ok=True)
    with (state / 'run.lock').open('a') as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            parser.error('another runner is using this output directory')
        run(args, trees, truth, species_newick, total, state)
    return 0


def run(args, trees, truth, species_newick, count, state):
    output = args.out
    suffix = ''
    costs = None  # worker does not override SimPhy reconciliation costs
    tasks = [(i, (trees[i], truth[i])) for i in range(count)]
    metadata = dict(version=1, species=species_newick, trace_memory=args.trace_memory,
                    input_sha256=hashlib.sha256(json.dumps([trees, truth], sort_keys=True).encode()).hexdigest())
    journal_path = state / 'results.jsonl'
    records = {}
    with journal_path.open('a+') as journal:
        journal.seek(0)
        header = journal.readline()
        if header:
            if json.loads(header) != metadata:
                raise ValueError('Inputs, species, costs, or memory-tracing setting differ from the journal; '
                                 'use another --output-dir for a fresh run.')
            while True:
                position = journal.tell()
                line = journal.readline()
                if not line:
                    break
                if not line.endswith('\n'):
                    journal.seek(position)
                    journal.truncate()
                    break
                result = json.loads(line)
                records[result['replicate']] = result
        else:
            journal.write(json.dumps(metadata) + '\n')
            journal.flush()
            os.fsync(journal.fileno())
        journal.seek(0, os.SEEK_END)
        aggregate = readWrite.readWrite().parse(species_newick)
        aggregate.label_internal()
        for result in records.values():
            add_events(aggregate, result)
        pending = iter((i, payload) for i, payload in tasks
                       if i not in records or (args.retry_failed and records[i]['status'] != 'ok'))
        print(f'{count} selected trees; {len(records)} journaled; '
              f'workers={args.workers}; trace_memory={args.trace_memory}; output={output}', flush=True)
        # Spawn safely even when imported libraries start native background
        # threads. Reuse workers so library imports are paid only once per slot.
        context = mp.get_context('spawn')
        active = {}
        idle = []
        exhausted = False
        since_export = 0
        stopped = False
        last_progress = time.monotonic()

        def stop(signum, frame):
            raise KeyboardInterrupt

        previous_term = signal.signal(signal.SIGTERM, stop)
        try:
            while active or not exhausted:
                while len(active) < args.workers and not exhausted:
                    try:
                        replicate, tree = next(pending)
                    except StopIteration:
                        exhausted = True
                        break
                    if idle:
                        process, receiver = idle.pop()
                    else:
                        receiver, sender = context.Pipe(duplex=True)
                        process = context.Process(target=worker_entry,
                                                  args=(sender, process_tree, args.trace_memory,
                                                        str(state) if args.worker_logs else None, costs))
                        process.start()
                        sender.close()
                    receiver.send((replicate, tree, species_newick))
                    active[replicate] = (process, receiver, time.monotonic())
                progressed = False
                for replicate, (process, receiver, started) in list(active.items()):
                    result = None
                    reusable = False
                    if receiver.poll():
                        try:
                            result = receiver.recv()
                            reusable = True
                        except EOFError:
                            result = dict(status='error', replicate=replicate, error='Worker exited without a result')
                    elif args.timeout is not None and time.monotonic() - started >= args.timeout:
                        process.kill()
                        result = dict(status='timeout', replicate=replicate,
                                      error=f'Exceeded {args.timeout} seconds')
                    elif not process.is_alive():
                        result = dict(status='error', replicate=replicate,
                                      error=f'Worker exited with code {process.exitcode}')
                    if result is None:
                        continue
                    if reusable:
                        idle.append((process, receiver))
                    else:
                        if process.is_alive():
                            process.kill()
                        process.join()
                        receiver.close()
                    del active[replicate]
                    result['gid'] = truth[replicate]['GID']
                    journal.write(json.dumps(result) + '\n')
                    journal.flush()
                    os.fsync(journal.fileno())
                    records[replicate] = result
                    add_events(aggregate, result)
                    since_export += 1
                    progressed = True
                    if result['status'] != 'ok':
                        print(f'Tree {replicate}: {result["status"]}: {result.get("error", "")}', flush=True)
                if since_export >= args.checkpoint_every:
                    export_results(output, records, suffix, aggregate)
                    since_export = 0
                    print(f'Checkpoint: {len(records)} recorded, {len(active)} active', flush=True)
                    last_progress = time.monotonic()
                if time.monotonic() - last_progress >= 30:
                    print(f'{len(records)} recorded; {len(active)} active', flush=True)
                    last_progress = time.monotonic()
                if not progressed and active:
                    time.sleep(0.05)
        except KeyboardInterrupt:
            stopped = True
            print('Stopping workers; exporting completed results. Rerun to resume.', flush=True)
        finally:
            for process, receiver, _ in list(active.values()) + [(p, c, 0) for p, c in idle]:
                if process.is_alive():
                    process.kill()
                process.join()
                receiver.close()
            export_results(output, records, suffix, aggregate)
            signal.signal(signal.SIGTERM, previous_term)
        counts = Counter(r['status'] for r in records.values())
        print(f'Results: {dict(counts)}', flush=True)
        if stopped:
            raise SystemExit(130)
        if any(r['status'] != 'ok' for r in records.values()):
            raise SystemExit(1)


if __name__ == '__main__':
    sys.exit(main())
