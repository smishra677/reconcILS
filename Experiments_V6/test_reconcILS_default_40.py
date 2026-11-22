"""Parallel 40-species reconciliation with legacy-compatible output tables."""
import os
# Parallelize across trees, not within imported numerical libraries.
for _name in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS',
              'NUMEXPR_NUM_THREADS', 'BLIS_NUM_THREADS'):
    os.environ[_name] = '1'

import argparse
import contextlib
import copy
import csv
import fcntl
import hashlib
import json
import math
import multiprocessing as mp
from pathlib import Path
import signal
import sys
import time
import traceback
import tracemalloc
from collections import Counter

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from reconcILS import reconcils
from utils_reconcILS import readWrite
from ete3 import PhyloTree, Tree

SPECIES = '((((M,(AF,B)),(P,(AJ,(AE,(O,D))))),((W,AH),((AN,N),((X,AG),(R,E))))),((H,(Y,AB)),((((L,V),(AC,A)),((C,Z),(F,K))),(((AM,AA),(AL,(T,(U,AI)))),((S,(G,AD)),((Q,J),(I,AK)))))));'
COLUMNS = ['Process', 'Replicate', 'Gene_tree', 'Species_Tree', 'Duplication',
           'NNI', 'DLCILS', 'Loss', 'Hemiplasy', 'RHemiplasy']
TRUE_COLUMNS = {'Duplication': 'observed_duplications', 'NNI': 'ILS',
                'DLCILS': 'ILS_DLCPAR', 'Loss': 'losses', 'Hemiplasy': 'CNH',
                'RHemiplasy': 'RKH'}


COSTS = (1.1, 1, 1)  # duplication, loss, ILS


def cost_suffix(args):
    values = (getattr(args, "duplication_cost", 1.1), getattr(args, "ils_cost", 1),
              getattr(args, "loss_cost", 1))
    return "_" + "_".join(str(int(v)) if float(v).is_integer() else str(v) for v in values)


def species_tree(newick=SPECIES):
    red = readWrite.readWrite()
    tree = red.parse_bio(newick)
    tree = tree.parse(tree.to_newick())
    tree.label_internal()
    return tree


def extract_equal(sp, events=0):
    stack = [sp]
    while stack:
        node = stack.pop()
        if node:
            events += node.equal
            stack.extend((node.leftChild, node.rightChild))
    return events


def process_tree(replicate, payload, species_newick, trace_memory=True):
    """Preserve the original calculation order, costs, and event labeling."""
    original_tree, truth = payload
    red = readWrite.readWrite()
    reco = reconcils()
    reco.D_cost, reco.L_cost, reco.I_cost = COSTS
    sp = red.parse(species_newick)
    tree = original_tree.replace('e-', '0')
    tr = red.parse(tree)
    dic = {key: [] for key in COLUMNS}
    true_row = dict(Process='True Process', Replicate=int(truth['replicate']),
                    Gene_tree=original_tree, Species_Tree=species_newick)
    true_row.update({name: int(truth[source]) for name, source in TRUE_COLUMNS.items()})
    for key in COLUMNS:
        dic[key].append(true_row[key])
    sp_small = species_tree(species_newick)
    edges = reco.get_edges(sp_small)

    started = time.perf_counter()
    if trace_memory:
        tracemalloc.start()
    genetree = PhyloTree(tr.to_newick())
    sptree = PhyloTree(species_newick)
    recon_tree, events = genetree.reconcile(sptree)
    ete_mem = tracemalloc.get_traced_memory()[1] / 10**6 if trace_memory else None
    if trace_memory:
        tracemalloc.stop()
    ete_time = time.perf_counter() - started
    dic['Gene_tree'].append(tr.to_newick())
    dic['Species_Tree'].append(species_newick)
    # sp_event_ete3 mutates loss labels; retain it before label_ete3.
    red.Create_pd_ete3('ETE3', replicate, dict(Counter(red.sp_event_ete3(recon_tree, events))), dic)
    fr = red.label_ete3(recon_tree, {})

    if trace_memory:
        tracemalloc.start()
    sp_copy = copy.deepcopy(sp)
    sp_copy.reset()
    reco.gene_tree = copy.deepcopy(tr)
    tr.order_gene(sp)
    tr.label_internal()
    sp.label_internal()
    tr.map_gene(sp)
    reco.setCost(sp)
    sp.isRoot = tr.isRoot = sp_copy.isRoot = True
    # Match the legacy timing window: preparation is traced but not timed.
    started = time.perf_counter()
    if tree == 'None':
        li = [[sp.taxa, sp.taxa, 'L', None]]
    else:
        li = reco.iterative_reconcILS(tr, sp, sp_copy, sp, [])
    rec_mem = tracemalloc.get_traced_memory()[1] / 10**6 if trace_memory else None
    if trace_memory:
        tracemalloc.stop()
    rec_time = time.perf_counter() - started
    if not li:
        return dict(status='empty', replicate=replicate,
                    error='Reconciliation returned no events')
    dic['Gene_tree'].append(red.to_newick(reco.gene_tree))
    for key in ('DLCILS', 'Hemiplasy', 'RHemiplasy'):
        dic[key].append(0)
    equal = dict(gt=tr.to_newick(), count=extract_equal(sp))
    rec = red.label_reconcILS(li)
    fr = red.find_replacement_for_edge(edges, fr)
    rec, fr, eve_rec, eve_ete = red.put_it_in_sp_tree(edges, rec, fr)
    reco.edge_to_event(sp_small, rec, 0)
    reco.edge_to_event(sp_small, fr, 1)
    dic['Species_Tree'].append(red.to_newick_sp(sp_small))
    for event, column in [('L', 'Loss'), ('D', 'Duplication')]:
        if eve_ete[event] != dic[column][-1]:
            raise ValueError(f'ETE3 {column} mismatch for replicate {replicate}: '
                             f'{eve_ete[event]} != {dic[column][-1]}')
    red.Create_pd('reconcILS', replicate, li, dic)
    if any(len(values) != 3 for values in dic.values()):
        raise ValueError(f'Incomplete result rows for replicate {replicate}')
    rows = [dict(zip(COLUMNS, values)) for values in zip(*(dic[k] for k in COLUMNS))]
    return dict(status='ok', replicate=replicate, rows=rows, equal=equal,
                rec_time=rec_time, ete_time=ete_time, rec_mem=rec_mem, ete_mem=ete_mem,
                edges_rec=[[list(k), rec[k]] for k in edges],
                edges_ete=[[list(k), fr[k]] for k in edges])


def write_csv_atomic(path, columns, rows):
    temporary = path.with_suffix(path.suffix + '.tmp')
    with temporary.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)
    os.replace(temporary, path)


def write_text_atomic(path, text):
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(text)
    os.replace(temporary, path)


def export_results(output, records, suffix, aggregate):
    def path(name):
        original = Path(name)
        return output / (original.stem + suffix + original.suffix)
    successful = sorted((r for r in records.values() if r['status'] == 'ok'),
                        key=lambda r: r['replicate'])
    write_csv_atomic(path('reconcILS_result.csv'), COLUMNS,
                     (row for result in successful for row in result['rows']))
    write_csv_atomic(path('equal_tree.csv'), ['gt', 'count'],
                     (r['equal'] for r in successful))
    write_csv_atomic(path('reconcILS_time_result.csv'), ['reconcILSTime_iterative', 'ete_time'],
                     (dict(reconcILSTime_iterative=r['rec_time'], ete_time=r['ete_time'])
                      for r in successful))
    write_csv_atomic(path('reconcILS_memory_result_with_dups_loss.csv'), ['reconcILS_mem', 'ete_mem'],
                     (dict(reconcILS_mem=r['rec_mem'], ete_mem=r['ete_mem']) for r in successful))
    write_text_atomic(path('reconcILS_labeled_sp.csv'), readWrite.readWrite().to_newick_sp(aggregate))
    errors = [f"replicate {r['replicate']}: {r['status']}: {r.get('error', '')}\n"
              for r in sorted(records.values(), key=lambda r: r['replicate']) if r['status'] != 'ok']
    write_text_atomic(path('run_error.txt'), ''.join(errors))


def add_events(aggregate, result):
    if result['status'] == 'ok':
        reco = reconcils()
        for flag, field in enumerate(('edges_rec', 'edges_ete')):
            # edge_to_event mutates dictionaries retained on the aggregate tree.
            reco.edge_to_event(aggregate, {tuple(k): dict(v) for k, v in result[field]}, flag)


def positive_int(value):
    number = int(value)
    if number < 1:
        raise argparse.ArgumentTypeError('must be at least 1')
    return number


def worker_entry(connection, function, trace_memory, log_directory, costs):
    """Reuse a clean spawned process until a task exceeds its deadline."""
    global COSTS
    COSTS = costs
    signal.signal(signal.SIGINT, signal.SIG_IGN)
    while True:
        try:
            task = connection.recv()
        except EOFError:
            break
        if task is None:
            break
        replicate, tree, species = task
        log_path = str(Path(log_directory) / f'tree_{replicate}.log') if log_directory else os.devnull
        try:
            with open(log_path, 'w') as log:
                with contextlib.redirect_stdout(log), contextlib.redirect_stderr(log):
                    result = function(replicate, tree, species, trace_memory)
        except Exception:
            result = dict(status='error', replicate=replicate, error=traceback.format_exc())
        finally:
            if tracemalloc.is_tracing():
                tracemalloc.stop()
        connection.send(result)
    connection.close()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    base = Path(__file__).resolve().parent
    parser.add_argument('--data', type=Path, default=base / 'Data_with_duplication_40_species_0.03_0.03',
                        help='folder containing trees.tre and summary.csv')
    parser.add_argument('--output-dir', type=Path, help='default: input folder (legacy layout)')
    parser.add_argument('--file_name', type=Path, default='trees.tre', help='filename')
    parser.add_argument('--species-tree', type=Path, help='Newick species tree with labels matching gene species prefixes')
    parser.add_argument('--duplication-cost', type=float, default=1.1)
    parser.add_argument('--loss-cost', type=float, default=1)
    parser.add_argument('--ils-cost', type=float, default=1)
    allocated = int(os.environ.get('SLURM_CPUS_PER_TASK', '1'))
    parser.add_argument('--workers', type=positive_int, default=max(1, min(8, allocated)))
    parser.add_argument('--limit', type=positive_int, default=1000, help='first N input rows (default: 1000)')
    parser.add_argument('--timeout', type=positive_int, default=1800,
                        help='hard seconds per tree, including worker startup and both methods')
    parser.add_argument('--checkpoint-every', type=positive_int, default=100)
    parser.add_argument('--trace-memory', action=argparse.BooleanOptionalAction, default=True,
                        help='measure memory by default; --no-trace-memory is faster but leaves blank cells')
    parser.add_argument('--retry-failed', action='store_true')
    parser.add_argument('--worker-logs', action='store_true')
    args = parser.parse_args()
    if any(not math.isfinite(v) or v <= 0 for v in
           (args.duplication_cost, args.loss_cost, args.ils_cost)):
        parser.error('Event costs must be finite and positive')
    args.species_newick = SPECIES
    if args.species_tree is not None:
        try:
            text = args.species_tree.read_text().strip()
            tree = Tree(text, format=1)
            names = tree.get_leaf_names()
            if len(names) < 2 or len(set(names)) != len(names):
                raise ValueError('at least two uniquely named species are required')
            if any('_' in name for name in names):
                raise ValueError('use species-prefix labels without underscores, matching the gene trees')
            args.species_newick = readWrite.readWrite().parse_bio(tree.write(format=9)).to_newick()
        except Exception as exc:
            parser.error(f'Invalid species tree: {exc}')
    args.data = args.data.resolve()
    file_name= args.file_name
    for name in (file_name, 'summary.csv'):
        if not (args.data / name).is_file():
            parser.error(f'Missing {args.data / name}; select the input folder with --data.')
    output = (args.output_dir or args.data).resolve()
    output.mkdir(parents=True, exist_ok=True)
    state = output / ('.default40_state' + cost_suffix(args))
    state.mkdir(exist_ok=True)
    with (state / 'run.lock').open('a') as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            parser.error('another runner is using this output directory')
        run(args, output,file_name)


def run(args, output,file_name):
    suffix = cost_suffix(args)
    species_newick = getattr(args, 'species_newick', SPECIES)
    costs = (getattr(args, 'duplication_cost', 1.1), getattr(args, 'loss_cost', 1),
             getattr(args, 'ils_cost', 1))
    gene_bytes = (args.data / file_name).read_bytes()
    summary_bytes = (args.data / 'summary.csv').read_bytes()
    trees = gene_bytes.decode().splitlines()
    summary = list(csv.DictReader(summary_bytes.decode().splitlines()))
    required = {'replicate', *TRUE_COLUMNS.values()}
    if not summary or not required.issubset(summary[0]):
        raise ValueError(f'summary.csv must contain {sorted(required)}')
    count = min(args.limit, len(summary))
    if len(trees) < count or any(not tree.strip() for tree in trees[:count]):
        raise ValueError(f'{file_name}  must contain one nonempty tree per selected summary row')
    # Positional pairing matches the original; True Process retains the summary
    # replicate label, while ETE3 and reconcILS retain zero-based row indices.
    tasks = [(i, (trees[i], summary[i])) for i in range(count)]
    metadata = dict(version=1, genes_sha256=hashlib.sha256(gene_bytes).hexdigest(),
                    summary_sha256=hashlib.sha256(summary_bytes).hexdigest(),
                    species=species_newick, costs=list(costs), trace_memory=args.trace_memory)
    state = output / ('.default40_state' + cost_suffix(args))
    state.mkdir(exist_ok=True)
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
        aggregate = species_tree(species_newick)
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
                    elif time.monotonic() - started >= args.timeout:
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
    main()
