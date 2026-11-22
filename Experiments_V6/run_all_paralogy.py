import os
# Parallelism is across trees. Prevent imported numerical libraries from starting
# extra native thread pools in the parent or oversubscribing each worker.
for _thread_variable in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS',
                         'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS', 'BLIS_NUM_THREADS'):
    os.environ[_thread_variable] = '1'
import copy
from ete3 import PhyloTree
import ete3
from collections import Counter
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from reconcILS import *

from utils_reconcILS import *
import tracemalloc
import pandas as pd 
import time
import re

def read_trees(i,folder):
        gene_tre= open(folder+'/rep_'+str(i)+'.tre')
        tr =gene_tre.read().strip().split('\n')
        gene_tre.close()
        return str(tr[0])

def get_edges(sp):
        return_it=[(repr(list(sorted(sp.taxa))),' to ',repr(list(sorted(sp.taxa))))]
        stack = [sp]
        while stack:
            curr=stack.pop()

            if curr.leftChild:
                return_it.append((repr(list(sorted(curr.taxa))),' to ',repr(list(sorted(curr.leftChild.taxa)))))
                stack.append(curr.leftChild)
            if curr.rightChild:
                return_it.append((repr(list(sorted(curr.taxa))),' to ',repr(list(sorted(curr.rightChild.taxa)))))
                stack.append(curr.rightChild)
        return return_it         

def edge_to_event(sp,dic):
        stack = [sp]
        sp.sp_ev_list+=[[dic[(repr(list(sorted(sp.taxa))),' to ',repr(list(sorted(sp.taxa))))],'Up']]
        while stack:
            curr=stack.pop()
            if curr.parent:
                curr.sp_ev_list+= [[dic[(repr(list(sorted(curr.parent.taxa))),' to ',repr(list(sorted(curr.taxa))))],'Up']]
            if curr.leftChild:
                stack.append(curr.leftChild)
            if curr.rightChild:
                stack.append(curr.rightChild)





def generate_species_codes(mapping,tree):
    for species, letter in mapping.items():
        tree = tree.replace(species, letter)
    return tree
species_to_letters = {
    'Carlitosyrichta': 'A',
    'Cebuscapucinus': 'B',
    'Cercocebusatys': 'C',
    'Macacafascicularis': 'D',
    'Macacanemestrina': 'E',
    'Theropithecusgelada': 'F',
    'Papioanubis': 'G',
    'Macacamulatta': 'H',
    'Mandrillusleucophaeus': 'I',
    'Chlorocebussabaeus': 'J',
    'Colobusangolensis': 'K',
    'Piliocolobustephrosceles': 'L',
    'Rhinopithecusbieti': 'M',
    'Rhinopithecusroxellana': 'N',
    'Gorillagorilla': 'O',
    'Homosapiens': 'P',
    'Panpaniscus': 'Q',
    'Pantroglodytes': 'R',
    'Pongoabelii': 'S',
    'Nomascusleucogenys': 'T',
    'Saimiriboliviensis': 'U',
    'Aotusnancymaae': 'V',
    'Callithrixjacchus': 'W',
}

list_to_prune=['Carlito_syrichta', 'Callithrix_jacchus', 'Aotus_nancymaae', 'Cebus_capucinus', 'Saimiri_boliviensis', 'Colobus_angolensis', 'Piliocolobus_tephrosceles', 'Rhinopithecus_bieti', 'Rhinopithecus_roxellana', 'Chlorocebus_sabaeus', 'Macaca_nemestrina', 'Macaca_fascicularis', 'Macaca_mulatta', 'Cercocebus_atys', 'Mandrillus_leucophaeus', 'Papio_anubis', 'Theropithecus_gelada', 'Nomascus_leucogenys', 'Pongo_abelii', 'Gorilla_gorilla', 'Pan_paniscus', 'Pan_troglodytes', 'Homo_sapiens']

def prune_trees(newick):
        red= readWrite.readWrite()
        t= ete3.Tree(newick)


       
        new_label=[]

        for i in list_to_prune:
                new_label+=[node.name for node in t.traverse() if node.name.startswith(i+'_')]


        t.prune(list(set(new_label)))
        

        newick = red.parse_bio(t.write(format=1)).to_newick()
                        
        return newick
        
def Create_pd_ete3_1(flag,i,o,dic):
        dic['Process']+=[flag]
        dic['Replicate']+=[i]
        dic['Duplication'].append(0)
        dic['Loss'].append(0)
        dic['NNI'].append(0)

        

        for i in o:
            if i in ['D','L']:
                    if i =='D':
                        dic['Duplication'][-1]=o[i]
                    elif i=='L':
                        dic['Loss'][-1]=o[i]
                    else:
                        continue

        return dic

def put_it_in_sp(df,node_to_leaves,species_edge_list):

    events={idx:{'D':0,'I':0,'L':0} for idx in species_edge_list}
    for idx, row in df.iterrows():

        if int(row['parentid'])==-1:
            if row['nodeid'].isnumeric():
                start= repr(list(sorted(node_to_leaves[str(row['nodeid'])])))
                #events[(start,' to ',start)]= {'D':row['dup'],'I':row['coal'],'L':row['loss']}
            else:
                start =repr(list({row['nodeid']}))
        
        else:
            start=repr(list(sorted(node_to_leaves[str(row['parentid'])])))
        end=''


        if row['nodeid'].isnumeric():
            end= repr(list(sorted(node_to_leaves[str(row['nodeid'])])))
        else:
            end =repr(list({row['nodeid']}))
        

        if (start,' to ',end) in events:
               events[(start,' to ',end)]= {'D':row['dup'],'I':row['coal'],'L':row['loss']}


    return events

def label_reconcILS(li):
        li.reverse()
        dic_reconcILS={}
        for i in li:
                if type(i)==list:
                        if (repr(list(sorted(i[0]))),' to ',repr(list(sorted(i[1]))))  in dic_reconcILS.keys():
                                if i[3]==True:
                                        dic_reconcILS[(repr(list(sorted(i[0]))),' to ',repr(list(sorted({i[1]}))))]+=[i[2]]
                                else:
                                        dic_reconcILS[(repr(list(sorted(i[0]))),' to ',repr(list(sorted(i[1]))))]+=[i[2]]
                        else:   
                                if i[3]==True:
                                        dic_reconcILS[(repr(list(sorted(i[0]))),' to ',repr(list(sorted({i[1]}))))]=[i[2]]
                                else:
                                        dic_reconcILS[(repr(list(sorted(i[0]))),' to ',repr(list(sorted(i[1]))))]=[i[2]]
        
        return dic_reconcILS

def label_ete3(recon_tree, fr):
    coun = 0
    coun1 = 0

    for node in recon_tree.traverse(strategy="preorder"):
        if len(node.children) >= 0:
            if hasattr(node, 'evoltype'):
                if node.evoltype in ['D', 'L']:
                    if node.evoltype == 'D':
                        coun += 1

                    node_ = sorted(node.get_species())

                    if node.evoltype == 'L' and len(node.children) > 1:
                        if not hasattr(node.children[1], 'evoltype') or not hasattr(node.children[0], 'evoltype'):
                            continue
                        else:
                            if node.up:
                                node_up = sorted(node.up.get_species())
                                if node.evoltype == 'L' and node_up == node_:
                                    continue
                                elif (repr(node_up), ' to ', repr(node_)) in fr.keys():
                                    fr[(repr(node_up), ' to ', repr(node_))] += [node.evoltype]
                                else:
                                    fr[(repr(node_up), ' to ', repr(node_))] = [node.evoltype]
                            else:
                                if (repr(node_), ' to ', repr(node_)) in fr.keys():
                                    fr[(repr(node_), ' to ', repr(node_))] += [node.evoltype]
                                else:
                                    fr[(repr(node_), ' to ', repr(node_))] = [node.evoltype]

                            if hasattr(node.children[1], 'evoltype'):
                                node.children[1].evoltype = 'S'
                            if hasattr(node.children[0], 'evoltype'):
                                node.children[0].evoltype = 'S'

                    else:
                        if node.up:
                            node_up = sorted(node.up.get_species())
                            tem_node = node.up

                            while node_ == node_up:
                                if tem_node.up:
                                    node_up = sorted(tem_node.up.get_species())
                                    tem_node = tem_node.up
                                else:
                                    break

                            if node.evoltype == 'D':
                                coun1 += 1

                            if (repr(node_up), ' to ', repr(node_)) in fr.keys():
                                fr[(repr(node_up), ' to ', repr(node_))] += [node.evoltype]
                            else:
                                fr[(repr(node_up), ' to ', repr(node_))] = [node.evoltype]
                        else:
                            if (repr(node_), ' to ', repr(node_)) in fr.keys():
                                fr[(repr(node_), ' to ', repr(node_))] += [node.evoltype]
                            else:
                                fr[(repr(node_), ' to ', repr(node_))] = [node.evoltype]
    return fr


def find_replacement_for_edge(species_edge_list,fr):
        list_to_replace ={}
        for edge in fr:
                if edge not in species_edge_list:
                        for edge_sp in species_edge_list:
                                if edge[2] in edge_sp[2]:
                                        list_to_replace[edge]=edge_sp


        for edge in list_to_replace:
                fr[list_to_replace[edge]]= fr[edge]
        
        return fr

def put_it_in_sp_tree(species_edge_list,dic_reconcILS,fr):
        for i in species_edge_list:
                rec_in_dic= {'D':0,'I':0,'L':0}
                if i in dic_reconcILS:
                        dic_in= dict(Counter(dic_reconcILS[i]))
                        for j in rec_in_dic:
                                if j in dic_in:
                                        rec_in_dic[j]=dic_in[j]
                                        eve_rec[j]+=dic_in[j]
                        dic_reconcILS[i]=rec_in_dic
                else:
                        dic_reconcILS[i]={'D':0,'I':0,'L':0}
                                                        
                                                
                ete_in_dic= {'D':0,'I':0,'L':0}

                if i in fr:
                        dic_in= dict(Counter(fr[i]))
                        for j in ete_in_dic:
                                if j in dic_in:
                                        ete_in_dic[j]=dic_in[j]
                                        eve_ete[j]+=dic_in[j]
                        fr[i]=ete_in_dic
                else:
                        fr[i]={'D':0,'I':0,'L':0}
        return dic_reconcILS,fr






def prepare_gene(tree, mode='paralogy'):
    """Keep the distinct preparation used by the original two triggers."""
    red = readWrite.readWrite()
    tr = red.parse_bio(red.to_newick(red.parse_bio(tree)))
    newick = tr.to_newick() if mode == 'sco' else red.to_newick(tr)
    if newick == ';':
        return newick
    if mode == 'sco':
        parsed = ete3.Tree(newick)
        keep = [name for name in species_to_letters if parsed.search_nodes(name=name)]
        if not keep:
            raise ValueError('SCO tree has none of the expected primate species')
        parsed.prune(list(set(keep)))
        newick = red.parse_bio(parsed.write(format=1)).to_newick()
    else:
        newick = prune_trees(newick)
    return generate_species_codes(species_to_letters, newick)


def process_tree(ili, tree, sp_str_ref, trace_memory=False, mode='paralogy'):
    """Run the original biological calculations for one isolated input tree."""
    global eve_rec, eve_ete
    reco = reconcils()
    red = readWrite.readWrite()
    tree = tree.replace('e-', '0')
    tr_str_ref = prepare_gene(tree, mode)
    sp = red.parse_bio(sp_str_ref)
    if tr_str_ref == ';':
        return dict(status='empty', replicate=ili)
    tr = red.parse(tr_str_ref)
    sp = sp.parse(sp_str_ref)
    sp_copy = copy.deepcopy(sp)

    start_time1 = time.perf_counter()
    if trace_memory:
        tracemalloc.start()
    tr.order_gene(sp)
    tr.label_internal()
    sp.label_internal()
    tr.map_gene(sp)
    reco.setCost(sp)
    sp.isRoot = tr.isRoot = sp_copy.isRoot = True
    #reco.L_cost = reco.D_cost = 2
    reco.gene_tree = copy.deepcopy(tr)
    species_edge_list = reco.get_edges(sp)
    if len(tree) < 50:
        reco.reconcILS(tr, sp, sp_copy, sp, [])
        li = red.sp_event(sp, [])
    else:
        li = reco.iterative_reconcILS(tr, sp, sp_copy, sp, [])
    rec_mem = tracemalloc.get_traced_memory()[1] / 10**6 if trace_memory else None
    if trace_memory:
        tracemalloc.stop()
    rec_time = time.perf_counter() - start_time1
    if not li:
        return dict(status='empty', replicate=ili)

    start_time2 = time.perf_counter()
    dic_reconcILS = red.label_reconcILS(li)
    if trace_memory:
        tracemalloc.start()
    genetree = PhyloTree(tr_str_ref)
    sptree = PhyloTree(sp_str_ref)
    recon_tree, events = genetree.reconcile(sptree)
    ete_mem = tracemalloc.get_traced_memory()[1] / 10**6 if trace_memory else None
    if trace_memory:
        tracemalloc.stop()
    # This also relabels loss nodes; its side effects are required by label_ete3.
    ete_counts = Counter(red.sp_event_ete3(recon_tree, events))
    if mode != 'sco':
        fr = label_ete3(recon_tree, {})
        eve_rec = {'D': 0, 'I': 0, 'L': 0}
        eve_ete = {'D': 0, 'I': 0, 'L': 0}
        fr = find_replacement_for_edge(species_edge_list, fr)
        dic_reconcILS, fr = put_it_in_sp_tree(species_edge_list, dic_reconcILS, fr)
    if mode == 'sco':
        # run_SCO uses the readWrite helpers (including removal of replaced edges).
        fr = red.label_ete3(recon_tree, {})
        fr = red.find_replacement_for_edge(species_edge_list, fr)
        dic_reconcILS, fr, eve_rec, eve_ete = red.put_it_in_sp_tree(
            species_edge_list, dic_reconcILS, fr)
        for event in ('D', 'L'):
            if eve_ete[event] != ete_counts[event]:
                raise ValueError(f'SCO ETE3 {event} count does not match branch events')
        per_tree = red.parse(sp_str_ref)
        per_tree.label_internal()
        reco.edge_to_event(per_tree, {k: dict(v) for k, v in dic_reconcILS.items()}, 0)
        reco.edge_to_event(per_tree, {k: dict(v) for k, v in fr.items()}, 1)
    rows = []
    for name, gene, events in [('reconcILS', red.to_newick(reco.gene_tree), eve_rec),
                               ('ete3', recon_tree.write(format=9), eve_ete)]:
        rows.append(dict(Process=name, Replicate=ili, Gene_tree=gene,
                         Species_Tree=red.to_newick_sp(per_tree if mode == 'sco' and name == 'reconcILS' else sp), Duplication=events['D'],
                         NNI=events['I'], Loss=events['L']))
    return dict(status='ok', replicate=ili, rows=rows, gene=tr_str_ref,
                rec_time=rec_time, post_time=time.perf_counter()-start_time2,
                rec_mem=rec_mem, ete_mem=ete_mem,
                edges_rec=[[list(k), dic_reconcILS[k]] for k in species_edge_list],
                edges_ete=[[list(k), fr[k]] for k in species_edge_list])


# The journal is authoritative; CSVs are periodically regenerated snapshots.
import argparse
import contextlib
import csv
import fcntl
import hashlib
import json
import multiprocessing as mp
import signal
import traceback

STEM = 'ZeroCol_ASTRAL_ML_ALLPARALOGS_MIN27'


def worker_entry(connection, function, trace_memory, log_directory, mode='paralogy'):
    """Reuse a clean spawned process until a task exceeds its deadline."""
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
                    result = function(replicate, tree, species, trace_memory, mode)
        except Exception:
            result = dict(status='error', replicate=replicate, error=traceback.format_exc())
        finally:
            if tracemalloc.is_tracing():
                tracemalloc.stop()
        connection.send(result)
    connection.close()


def write_csv_atomic(path, columns, rows):
    temporary = path.with_suffix(path.suffix + '.tmp')
    with temporary.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)
    os.replace(temporary, path)


def export_results(output, records, order, species_tree, stem=STEM):
    successful = sorted((r for r in records.values() if r['status'] == 'ok'),
                        key=lambda r: order[r['replicate']])
    write_csv_atomic(output / f'bio_result_{stem}_pruned67.csv',
                     ['Process', 'Replicate', 'Gene_tree', 'Species_Tree', 'Duplication', 'NNI', 'Loss'],
                     (row for r in successful for row in r['rows']))
    for suffix, field in [('pruned67', 'rec_time'), ('intro_1_write_pruned67', 'post_time')]:
        write_csv_atomic(output / f'timing_{stem}_{suffix}.csv', ['time'],
                         ({'time': r[field]} for r in successful))
    write_csv_atomic(output / f'{stem}_pruned_Gene_trees67.csv', ['Repliccate', 'Gene_trees'],
                     ({'Repliccate': r['replicate'], 'Gene_trees': r['gene']} for r in successful))
    write_csv_atomic(output / f'mem_{stem}_pruned67.csv', ['reconcILS_mem', 'ete_mem'],
                     ({'reconcILS_mem': r['rec_mem'], 'ete_mem': r['ete_mem']} for r in successful))
    write_csv_atomic(output / 'status.csv', ['Replicate', 'status', 'error'],
                     (dict(Replicate=r['replicate'], status=r['status'], error=r.get('error', ''))
                      for r in sorted(records.values(), key=lambda r: order[r['replicate']])))
    target = output / f'labeled_L_{stem}_pruned67.csv'
    temporary = target.with_suffix('.tmp')
    temporary.write_text(readWrite.readWrite().to_newick_sp(species_tree))
    os.replace(temporary, target)


def add_events(species_tree, result):
    if result['status'] == 'ok':
        reco = reconcils()
        for flag, field in enumerate(['edges_rec', 'edges_ete']):
            # edge_to_event retains and mutates the supplied dictionaries.
            reco.edge_to_event(species_tree, {tuple(k): dict(v) for k, v in result[field]}, flag)


def positive_int(value):
    number = int(value)
    if number < 1:
        raise argparse.ArgumentTypeError('must be at least 1')
    return number


def main():
    parser = argparse.ArgumentParser(description='Parallel, resumable paralogy or SCO reconciliation.')
    base = Path(__file__).resolve().parent
    parser.add_argument('--mode', choices=['auto', 'paralogy', 'sco'], default='auto',
                        help='auto selects SCO for an input filename containing SCO')
    parser.add_argument('--genes', '--gene-trees', type=Path, default=base / 'smith_et_al_2022_primates' / f'{STEM}.tre')
    parser.add_argument('--species', type=Path, default=base / 'all_species tree' / 'sp_tree_pruned.tre')
    parser.add_argument('--output-dir', type=Path, default=None)
    allocated = int(os.environ.get('SLURM_CPUS_PER_TASK', '1'))
    parser.add_argument('--workers', type=positive_int, default=max(1, min(allocated, 20)),
                        help='concurrent trees (default: allocated Slurm CPUs, capped at 20; otherwise 1)')
    parser.add_argument('--timeout', type=positive_int, default=18000000,
                        help='hard wall-clock limit per tree in seconds, including ETE3')
    parser.add_argument('--checkpoint-every', type=positive_int, default=100)
    parser.add_argument('--limit', type=positive_int, help='process only the first N trees in length order')
    parser.add_argument('--trace-memory', action='store_true', help='enable slower allocation tracing')
    parser.add_argument('--retry-failed', action='store_true', help='retry journaled errors/timeouts/empty results')
    parser.add_argument('--worker-logs', action='store_true', help='save verbose logs for each tree')
    args = parser.parse_args()
    if args.mode == 'auto':
        args.mode = 'sco' if 'SCO' in args.genes.stem.upper() else 'paralogy'
    if args.output_dir is None:
        args.output_dir = base / ('sco_efficient' if args.mode == 'sco' else 'paralogy_efficient')
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    # Hold a lock throughout execution, including exports, to prevent competing writers.
    with (output / '.run.lock').open('a') as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            parser.error('another runner is using this output directory')
        run(args, output)


def run(args, output):
    gene_bytes = args.genes.read_bytes()
    species_bytes = args.species.read_bytes()
    metadata = dict(version=1, genes_sha256=hashlib.sha256(gene_bytes).hexdigest(),
                    species_sha256=hashlib.sha256(species_bytes).hexdigest())
    if args.mode == 'sco':
        metadata.update(mode='sco', sco_version=1)
    stem = args.genes.stem if args.mode == 'sco' else STEM
    journal_path = output / 'results.jsonl'
    records = {}
    with journal_path.open('a+') as journal:
        journal.seek(0)
        header = journal.readline()
        if header:
            if json.loads(header) != metadata:
                raise ValueError('Input files differ from this journal; use another --output-dir.')
            while True:
                position = journal.tell()
                line = journal.readline()
                if not line:
                    break
                if not line.endswith('\n'):
                    # A killed writer may leave only its final record incomplete.
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
        trees = [(i, tree) for i, tree in enumerate(gene_bytes.decode().split(';')) if tree.strip()]
        trees.sort(key=lambda item: len(item[1]))
        order = {i: rank for rank, (i, _) in enumerate(trees)}
        red = readWrite.readWrite()
        species = red.parse_bio(species_bytes.decode().replace('e-', '0')).to_newick()
        species_tree = red.parse_bio(species).parse(species)
        species_tree.label_internal()
        for result in records.values():
            add_events(species_tree, result)
        selected = trees[:args.limit] if args.limit else trees
        pending = iter((i, tree) for i, tree in selected
                       if i not in records or (args.retry_failed and records[i]['status'] != 'ok'))
        print(f'{len(trees)} inputs; {len(records)} journaled; workers={args.workers}; output={output}', flush=True)
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
                                                        str(output) if args.worker_logs else None, args.mode))
                        process.start()
                        sender.close()
                    receiver.send((replicate, tree, species))
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
                    add_events(species_tree, result)
                    since_export += 1
                    progressed = True
                    if result['status'] != 'ok':
                        print(f'Tree {replicate}: {result["status"]}: {result.get("error", "")}', flush=True)
                if since_export >= args.checkpoint_every:
                    export_results(output, records, order, species_tree, stem)
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
            export_results(output, records, order, species_tree, stem)
            signal.signal(signal.SIGTERM, previous_term)
        counts = Counter(r['status'] for r in records.values())
        print(f'Results: {dict(counts)}', flush=True)
        if stopped:
            raise SystemExit(130)
        if any(r['status'] != 'ok' for r in records.values()):
            raise SystemExit(1)


if __name__ == '__main__':
    main()
