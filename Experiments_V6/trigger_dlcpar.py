"""Parallel True Process / DLCpar / ETE3 / reconcILS comparison.

Supply a species tree with labels matching gene species prefixes. Costs are
shared by DLCpar and reconcILS; --ils-cost sets both DLCpar -C and -K.
"""
import os
for _variable in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS',
                  'NUMEXPR_NUM_THREADS', 'BLIS_NUM_THREADS'):
    os.environ[_variable] = '1'

import argparse
from collections import Counter
import contextlib
import copy
import csv
import fcntl
import hashlib
import io
import json
import math
import multiprocessing as mp
from pathlib import Path
import shutil
import signal
import subprocess
import sys
import time
import traceback
import tracemalloc

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from reconcILS import reconcils
from utils_reconcILS import readWrite
from ete3 import PhyloTree, Tree

FIELDS = ['Process', 'Replicate', 'Gene_tree', 'Species_Tree', 'Duplication',
          'NNI', 'DLCILS', 'Loss', 'Hemiplasy', 'RHemiplasy']
TRUE_FIELDS = {'Duplication': 'observed_duplications', 'NNI': 'ILS',
               'DLCILS': 'ILS_DLCPAR', 'Loss': 'losses', 'Hemiplasy': 'CNH',
               'RHemiplasy': 'RKH'}


def suffix_for(d, i, l):
    return '_' + '_'.join(str(int(x)) if float(x).is_integer() else str(x) for x in (d, i, l))


def make_species(newick):
    tree = readWrite.readWrite().parse(newick)
    tree.label_internal()
    return tree


def peak_memory(enabled):
    value = tracemalloc.get_traced_memory()[1] / 10**6 if enabled else None
    if enabled:
        tracemalloc.stop()
    return value


def process_tree(index, original, truth, config, work):
    """Retain the original four-method order and biological calculations."""
    if original.strip().rstrip(';') == 'None':
        return dict(status='skipped', replicate=index, reason='Missing gene tree')
    total = sum(int(truth[k]) for k in ('observed_duplications', 'ILS_DLCPAR', 'losses'))
    if config['max_events'] >= 0 and total > config['max_events']:
        return dict(status='skipped', replicate=index,
                    reason=f'True D + DLCILS + L = {total} exceeds {config["max_events"]}')
    red = readWrite.readWrite()
    reco = reconcils()
    reco.D_cost, reco.I_cost, reco.L_cost = config['costs']
    species_text = config['species']
    sp = red.parse(species_text)
    # Preserve the legacy gene-tree preparation convention.
    tr = red.parse(original.replace('e-', '0'))
    normalized = tr.to_newick()
    gene_path = work / f'rep_1_{index}.tre'
    dlc_gene = red.to_newick(tr)  # Retain unique copy names for DLCpar.
    gene_path.write_text(dlc_gene + '\n')
    species_path = work / 'species.tre'
    species_path.write_text(species_text + '\n')
    map_path = work / 'genes.mapsl'
    if config['species_map'] is not None:
        map_path.write_text(config['species_map'])
    else:
        species_names = set(Tree(species_text, format=1).get_leaf_names())
        leaves = Tree(dlc_gene, format=1).get_leaf_names()
        if len(set(leaves)) != len(leaves):
            raise ValueError('DLCpar requires distinct gene-copy labels')
        pairs = []
        for name in leaves:
            species_name = name.split('_', 1)[0]
            if species_name not in species_names:
                raise ValueError(f'Unknown species prefix {species_name!r} in gene {name!r}')
            pairs.append(f'{name}\t{species_name}\n')
        map_path.write_text(''.join(pairs))

    command = [config['python'], config['dlcpar']]
    d, i, l = config['costs']
    dp = command + ['dp', '-D', str(d), '-L', str(l), '-C', str(i), '-K', str(i),
                    '-s', str(species_path), '-S', str(map_path)]
    if config['seed'] is not None:
        dp += ['--seed', str(config['seed'] + index)]
    dp += [str(gene_path)]
    event_command = command + ['events', '--lct', '-s', str(species_path), '-S', str(map_path),
                              str(gene_path) + '.dlcdp.lct.tree']
    traced = config['trace_memory']
    started = time.perf_counter()
    if traced:
        tracemalloc.start()
    with (work / 'dlcpar.log').open('w') as log:
        subprocess.run(dp, cwd=work, stdout=log, stderr=log, check=True)
        with (work / 'events.tsv').open('w') as events_file:
            subprocess.run(event_command, cwd=work, stdout=events_file, stderr=log, check=True)
    dlc_mem = peak_memory(traced)
    dlc_time = time.perf_counter() - started
    with (work / 'events.tsv').open() as stream:
        table = csv.DictReader(stream, delimiter='\t')
        if not {'dup', 'loss', 'coal'}.issubset(table.fieldnames or []):
            raise ValueError('DLCpar event table is missing dup/loss/coal columns')
        event_rows = list(table)
    if not event_rows:
        raise ValueError('DLCpar returned an empty event table')
    counts = {key: sum(int(row[key]) for row in event_rows) for key in ('dup', 'loss', 'coal')}
    true = dict(Process='True Process', Replicate=int(truth['replicate']),
                Gene_tree=original, Species_Tree=species_text)
    true.update({k: int(truth[v]) for k, v in TRUE_FIELDS.items()})
    dlc = dict(Process='DLCpar', Replicate=str(index), Gene_tree=normalized, Species_Tree=species_text,
               Duplication=counts['dup'], NNI=0, DLCILS=counts['coal'], Loss=counts['loss'],
               Hemiplasy=0, RHemiplasy=0)
    dic = {field: [true[field], dlc[field]] for field in FIELDS}

    if traced:
        tracemalloc.start()
    gene_ete = PhyloTree(normalized)
    species_ete = PhyloTree(species_text)
    started = time.perf_counter()
    recon_ete, events = gene_ete.reconcile(species_ete)
    ete_mem = peak_memory(traced)
    ete_time = time.perf_counter() - started
    dic['Gene_tree'].append(normalized)
    dic['Species_Tree'].append(species_text)
    # Keep this runner's original order: label_ete3 before sp_event_ete3.
    fr = red.label_ete3(recon_ete, {})
    red.Create_pd_ete3('ETE3', index, dict(Counter(red.sp_event_ete3(recon_ete, events))), dic)

    if traced:
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
    started = time.perf_counter()
    li = reco.iterative_reconcILS(tr, sp, sp_copy, sp, [])
    rec_time = time.perf_counter() - started
    if not li:
        raise RuntimeError('reconcILS returned no event records; inspect worker.log for an internal timeout')
    reference = make_species(species_text)
    edges = reco.get_edges(reference)
    fr = red.find_replacement_for_edge(edges, fr)
    rec = red.label_reconcILS(li)
    rec, fr, _, _ = red.put_it_in_sp_tree(edges, rec, fr)
    reco.edge_to_event(reference, rec, 0)
    reco.edge_to_event(reference, fr, 1)
    rec_mem = peak_memory(traced)
    dic['Gene_tree'].append(red.to_newick(reco.gene_tree))
    dic['Species_Tree'].append(red.to_newick_sp(reference))
    for field in ('DLCILS', 'Hemiplasy', 'RHemiplasy'):
        dic[field].append(0)
    red.Create_pd('reconcILS', index, li, dic)
    rows = [dict(zip(FIELDS, values)) for values in zip(*(dic[field] for field in FIELDS))]
    if len(rows) != 4 or any(len(v) != 4 for v in dic.values()):
        raise ValueError('Incomplete four-method result')
    return dict(status='ok', replicate=index, rows=rows, rec_time=rec_time, ete_time=ete_time,
                dlcpar_time=dlc_time, rec_mem=rec_mem, ete_mem=ete_mem, dlcpar_mem=dlc_mem,
                edges_rec=[[list(k), rec[k]] for k in edges],
                edges_ete=[[list(k), fr[k]] for k in edges])


def worker(connection, function, config, state):
    # DLCpar subprocesses inherit this worker's process group, so an outer
    # deadline or interruption can kill the entire task, not leave DLCpar alive.
    os.setsid()
    signal.signal(signal.SIGINT, signal.SIG_IGN)
    while True:
        try:
            task = connection.recv()
        except EOFError:
            break
        if task is None:
            break
        index, gene, truth = task
        work = Path(state) / f'rep_{index}'
        work.mkdir(exist_ok=True)
        try:
            with (work / 'worker.log').open('w', buffering=1) as log:
                with contextlib.redirect_stdout(log), contextlib.redirect_stderr(log):
                    result = function(index, gene, truth, config, work)
        except Exception:
            result = dict(status='error', replicate=index, error=traceback.format_exc(), work=str(work))
        finally:
            if tracemalloc.is_tracing():
                tracemalloc.stop()
        connection.send(result)
    connection.close()


def close_worker(process, connection):
    try:
        os.killpg(process.pid, signal.SIGKILL)
    except ProcessLookupError:
        if process.is_alive():
            process.kill()
    process.join()
    connection.close()


def write_csv(path, fields, rows):
    temp = path.with_name(path.name + '.tmp')
    with temp.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    os.replace(temp, path)


def write_text(path, text):
    temp = path.with_name(path.name + '.tmp')
    temp.write_text(text)
    os.replace(temp, path)


def export_results(output, suffix, records, species, count):
    selected = sorted((r for idx, r in records.items() if idx < count), key=lambda r: r['replicate'])
    successful = [r for r in selected if r['status'] == 'ok']
    red, reco = readWrite.readWrite(), reconcils()
    aggregate = make_species(species)
    rows = []
    for result in successful:
        per_tree = make_species(species)
        for flag, key in enumerate(('edges_rec', 'edges_ete')):
            reco.edge_to_event(aggregate, {tuple(k): dict(v) for k, v in result[key]}, flag)
            reco.edge_to_event(per_tree, {tuple(k): dict(v) for k, v in result[key]}, flag)
        copied = [dict(row) for row in result['rows']]
        # Rebuild from saved edges so older journals also export independent labels.
        copied[-1]['Species_Tree'] = red.to_newick_sp(per_tree)
        rows.extend(copied)
    write_csv(output / f'_1_result{suffix}.csv', FIELDS, rows)
    write_csv(output / f'_1_time{suffix}.csv', ['reconcILSTime_iterative', 'ete_time', 'dlcpar_time'],
              (dict(reconcILSTime_iterative=r['rec_time'], ete_time=r['ete_time'], dlcpar_time=r['dlcpar_time'])
               for r in successful))
    write_csv(output / f'_1_mem{suffix}.csv', ['reconcILSMEM_iterative', 'ete_mem', 'dlcpar_mem'],
              (dict(reconcILSMEM_iterative=r['rec_mem'], ete_mem=r['ete_mem'], dlcpar_mem=r['dlcpar_mem'])
               for r in successful))
    write_text(output / f'_1_labeled_sp{suffix}.csv', red.to_newick_sp(aggregate))
    write_csv(output / f'status{suffix}.csv', ['Replicate', 'status', 'reason'],
              (dict(Replicate=r['replicate'], status=r['status'], reason=r.get('reason', r.get('error', '')))
               for r in selected))
    write_text(output / f'error{suffix}.txt', ''.join(
        f"replicate {r['replicate']}: {r.get('error', '')}\n" for r in selected
        if r['status'] not in ('ok', 'skipped')))


def load_journal(stream, metadata):
    stream.seek(0)
    header = stream.readline()
    records = {}
    if header:
        if json.loads(header) != metadata:
            raise ValueError('Inputs, species/map, costs, software, or settings changed; use another --output-dir')
        while True:
            offset = stream.tell()
            line = stream.readline()
            if not line:
                break
            if not line.endswith('\n'):
                stream.seek(offset)
                stream.truncate()
                break
            result = json.loads(line)
            records[result['replicate']] = result
    else:
        stream.write(json.dumps(metadata) + '\n')
        stream.flush()
        os.fsync(stream.fileno())
    stream.seek(0, os.SEEK_END)
    return records


def run_jobs(args, tasks, config, metadata, function=process_tree):
    output = args.output_dir
    suffix = suffix_for(*config['costs'])
    state = output / ('.dlcpar_state' + suffix)
    state.mkdir(parents=True, exist_ok=True)
    with (state / 'run.lock').open('a') as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            raise RuntimeError('Another runner is using this cost/output combination') from None
        with (state / 'results.jsonl').open('a+') as journal:
            records = load_journal(journal, metadata)
            pending = iter(task for task in tasks if task[0] not in records or
                           (args.retry_failed and records[task[0]]['status'] not in ('ok', 'skipped')))
            context = mp.get_context('spawn')
            active, idle = {}, []
            exhausted = False
            stopped = False
            since_export = 0
            last_progress = time.monotonic()
            print(f'{len(tasks)} selected inputs; {len(records)} journaled; workers={args.workers}; '
                  f'output={output}; suffix={suffix}', flush=True)

            def stop(signum, frame):
                raise KeyboardInterrupt

            old_term = signal.signal(signal.SIGTERM, stop)
            try:
                while active or not exhausted:
                    while len(active) < args.workers and not exhausted:
                        try:
                            task = next(pending)
                        except StopIteration:
                            exhausted = True
                            break
                        if idle:
                            process, receiver = idle.pop()
                        else:
                            receiver, sender = context.Pipe(duplex=True)
                            process = context.Process(target=worker, args=(sender, function, config, str(state)))
                            process.start()
                            sender.close()
                        receiver.send(task)
                        active[task[0]] = (process, receiver, time.monotonic())
                    progressed = False
                    for index, (process, receiver, started) in list(active.items()):
                        result = None
                        reusable = False
                        if receiver.poll():
                            try:
                                result = receiver.recv()
                                reusable = True
                            except EOFError:
                                result = dict(status='error', replicate=index, error='Worker exited without a result')
                        elif args.timeout is not None and time.monotonic() - started >= args.timeout:
                            result = dict(status='timeout', replicate=index,
                                          error=f'Exceeded {args.timeout} seconds; worker and DLCpar stopped')
                        elif not process.is_alive():
                            result = dict(status='error', replicate=index, error=f'Worker exit {process.exitcode}')
                        if result is None:
                            continue
                        if reusable:
                            idle.append((process, receiver))
                        else:
                            close_worker(process, receiver)
                        del active[index]
                        journal.write(json.dumps(result) + '\n')
                        journal.flush()
                        os.fsync(journal.fileno())
                        records[index] = result
                        since_export += 1
                        progressed = True
                        if result['status'] not in ('ok', 'skipped'):
                            print(f"Tree {index}: {result['status']}; see {state / ('rep_' + str(index))}", flush=True)
                    if since_export >= args.checkpoint_every:
                        export_results(output, suffix, records, config['species'], len(tasks))
                        since_export = 0
                    if time.monotonic() - last_progress >= 30:
                        print(f'{len(records)} recorded; {len(active)} active', flush=True)
                        last_progress = time.monotonic()
                    if not progressed and active:
                        time.sleep(.05)
            except KeyboardInterrupt:
                stopped = True
                print('Stopping workers and DLCpar subprocesses; rerun to resume.', flush=True)
            finally:
                for process, receiver, _ in list(active.values()) + [(p, c, 0) for p, c in idle]:
                    close_worker(process, receiver)
                try:
                    export_results(output, suffix, records, config['species'], len(tasks))
                finally:
                    signal.signal(signal.SIGTERM, old_term)
            counts = Counter(r['status'] for i, r in records.items() if i < len(tasks))
            print(f'Results: {dict(counts)}', flush=True)
            if stopped:
                raise SystemExit(130)
            if any(status not in ('ok', 'skipped') for status in counts):
                raise SystemExit(1)


def positive_int(text):
    value = int(text)
    if value < 1:
        raise argparse.ArgumentTypeError('must be positive')
    return value


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data', type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument('--gene-trees', '--genes', type=Path, help='default: DATA/trees.tre')
    parser.add_argument('--summary', type=Path, help='default: DATA/summary.csv (True Process observations)')
    parser.add_argument('--species-tree', type=Path, required=True)
    parser.add_argument('--species-map', type=Path, help='optional DLCpar map; otherwise derive exact gene-prefix mappings')
    parser.add_argument('--output-dir', type=Path, help='default: DATA')
    parser.add_argument('--duplication-cost', type=float, default=1.1)
    parser.add_argument('--ils-cost', type=float, default=1, help='reconcILS I cost and DLCpar -C and -K')
    parser.add_argument('--loss-cost', type=float, default=1)
    parser.add_argument('--workers', type=positive_int,
                        default=max(1, min(8, int(os.environ.get('SLURM_CPUS_PER_TASK', '1')))))
    parser.add_argument('--limit', type=positive_int, default=1000)
    parser.add_argument('--max-events', type=int, default=25,
                        help='legacy filter on true duplication + DLCILS + loss; -1 disables it')
    parser.add_argument('--checkpoint-every', type=positive_int, default=100)
    parser.add_argument('--timeout', type=positive_int, help='optional hard seconds per tree, including child processes')
    parser.add_argument('--trace-memory', action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument('--retry-failed', action='store_true')
    parser.add_argument('--seed', type=int, help='optional DLCpar seed; per-tree seed is this value plus input index')
    parser.add_argument('--dlcpar', type=Path, help='path to DLCpar script (auto-detect local installation)')
    parser.add_argument('--dlcpar-python', help='Python interpreter for DLCpar; auto-detect its Python 2 environment')
    parser.add_argument('--preflight-timeout', type=positive_int, default=300,
                        help='seconds allowed for each DLCpar startup check (default: 180)')
    args = parser.parse_args()
    if any(not math.isfinite(x) or x <= 0 for x in
           (args.duplication_cost, args.ils_cost, args.loss_cost)):
        parser.error('Event costs must be finite and positive')
    if args.max_events < -1:
        parser.error('--max-events must be -1 or nonnegative')
    args.data = args.data.resolve()
    args.output_dir = (args.output_dir or args.data).resolve()
    genes_path = (args.gene_trees or args.data / 'trees.tre').resolve()
    summary_path = (args.summary or args.data / 'summary.csv').resolve()
    try:
        species = args.species_tree.read_text().strip()
        stree = Tree(species, format=1)
        names = stree.get_leaf_names()
        if len(names) < 2 or len(names) != len(set(names)) or any('_' in name for name in names):
            raise ValueError('species labels must be distinct prefixes without gene-copy underscores')
        if any(len(node.children) != 2 for node in stree.traverse() if not node.is_leaf()):
            raise ValueError('species tree must be binary')
        text = genes_path.read_text().strip()
        if not text or not text.endswith(';'):
            raise ValueError('gene-tree file is empty or not semicolon-terminated')
        genes = [part.strip() + ';' for part in text.split(';') if part.strip()]
        summary_text = summary_path.read_text()
        summary = list(csv.DictReader(io.StringIO(summary_text)))
        if not summary or not {'replicate', *TRUE_FIELDS.values()}.issubset(summary[0]):
            raise ValueError(f'summary.csv must contain replicate and {list(TRUE_FIELDS.values())}')
        count = min(args.limit, len(genes))
        if len(summary) < count:
            raise ValueError('summary.csv has fewer rows than the selected gene trees')
        mapping = args.species_map.read_text() if args.species_map else None
    except (OSError, ValueError) as exc:
        parser.error(str(exc))
    if args.dlcpar is None:
        candidates = [parent / 'dlcpar-master/bin/dlcpar' for parent in ROOT.parents]
        executable = shutil.which('dlcpar')
        args.dlcpar = Path(executable) if executable else next((p for p in candidates if p.is_file()), None)
    if args.dlcpar is None or not args.dlcpar.is_file():
        parser.error('DLCpar was not found; supply --dlcpar /path/to/dlcpar')
    args.dlcpar = args.dlcpar.resolve()
    if args.dlcpar_python is None:
        python = args.dlcpar.parent.parent / 'py2env/bin/python'
        args.dlcpar_python = str(python) if python.is_file() else shutil.which('python2')
    python_path = shutil.which(args.dlcpar_python) if args.dlcpar_python else None
    if python_path is None:
        parser.error('DLCpar interpreter not found; supply --dlcpar-python')
    try:
        check = subprocess.run([python_path, str(args.dlcpar), '--version'],
                               capture_output=True, text=True, timeout=args.preflight_timeout, check=True)
    except subprocess.TimeoutExpired:
        parser.error(f'DLCpar version check exceeded {args.preflight_timeout} seconds; '
                     'increase --preflight-timeout for slow startup')
    except (OSError, subprocess.SubprocessError) as exc:
        parser.error(f'DLCpar preflight failed: {exc}')
    version = (check.stdout.strip() or check.stderr.strip())
    for command in ('dp', 'events'):
        try:
            subprocess.run([python_path, str(args.dlcpar), command, '--help'],
                           capture_output=True, text=True, timeout=args.preflight_timeout, check=True)
        except subprocess.TimeoutExpired:
            parser.error(f'DLCpar {command} startup check exceeded {args.preflight_timeout} seconds; '
                         'increase --preflight-timeout for slow startup')
        except subprocess.CalledProcessError as exc:
            parser.error(f'DLCpar {command} cannot start with {python_path}:\n{exc.stderr.strip()}')
        except (OSError, subprocess.SubprocessError) as exc:
            parser.error(f'DLCpar {command} preflight failed: {exc}')
    print(version, flush=True)
    config = dict(species=species, species_map=mapping,
                  costs=[args.duplication_cost, args.ils_cost, args.loss_cost],
                  python=python_path, dlcpar=str(args.dlcpar), seed=args.seed,
                  trace_memory=args.trace_memory, max_events=args.max_events)
    digest = hashlib.sha256()
    for path in [ROOT / 'reconcILS.py', *sorted((ROOT / 'utils_reconcILS').glob('*.py'))]:
        digest.update(path.read_bytes())
    metadata = dict(version=1, config=config, algorithm_sha256=digest.hexdigest(),
                    genes_sha256=hashlib.sha256(text.encode()).hexdigest(),
                    summary_sha256=hashlib.sha256(summary_text.encode()).hexdigest(),
                    dlcpar_version=version,
                    dlcpar_script_sha256=hashlib.sha256(args.dlcpar.read_bytes()).hexdigest())
    suffix = suffix_for(*config['costs'])
    destinations = [args.output_dir / f'{name}{suffix}{extension}' for name, extension in
                    [('_1_result', '.csv'), ('_1_time', '.csv'), ('_1_mem', '.csv'),
                     ('_1_labeled_sp', '.csv'), ('status', '.csv'), ('error', '.txt')]]
    inputs = {genes_path, summary_path, args.species_tree.resolve()}
    if args.species_map:
        inputs.add(args.species_map.resolve())
    if any(path.resolve() in inputs for path in destinations):
        parser.error('An output path would overwrite an input file; choose another --output-dir')
    try:
        run_jobs(args, [(i, genes[i], summary[i]) for i in range(count)], config, metadata)
    except (OSError, ValueError, RuntimeError) as exc:
        parser.exit(1, f'{exc}\n')


if __name__ == '__main__':
    main()
