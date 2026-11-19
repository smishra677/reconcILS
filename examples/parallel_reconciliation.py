"""Run reconcILS in parallel on user-supplied gene and species trees.

Species-tree leaves use original ``SpeciesName`` labels.  Every gene-tree leaf
must be ``SpeciesName_copy``.  A deterministic two-character internal alias is
assigned to each species, copy suffixes are preserved exactly, and all
user-facing outputs are restored to the original labels.

The species tree, gene-tree file, outputs, costs, optional alias map, worker
count, timeouts, and checkpoint settings are command-line parameters.
"""
from __future__ import annotations

import argparse
import contextlib
import copy
import csv
import fcntl
import hashlib
import json
import math
import multiprocessing as mp
import os
from pathlib import Path
import signal
import sys
import time
import traceback
from collections import Counter
from typing import Any, Mapping, Sequence


PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from reconcILS.copy_labels import SpeciesCopyMap, load_supplied_mapping


THREAD_VARIABLES = (
    "OPENBLAS_NUM_THREADS",
    "OMP_NUM_THREADS",
    "MKL_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "BLIS_NUM_THREADS",
)
FIELDS = [
    "Process",
    "Replicate",
    "Gene_tree",
    "Species_Tree",
    "Duplication",
    "NNI",
    "Loss",
    "Seconds",
]


def configure_native_threads(count: int) -> None:
    """Set native-library thread counts before importing reconcILS dependencies."""
    value = str(count)
    for name in THREAD_VARIABLES:
        os.environ[name] = value


def runtime_classes():
    """Load reconcILS while supporting both an installed package and a checkout."""
    project_root = Path(__file__).resolve().parent.parent
    source_package = project_root / "reconcILS"
    if source_package.is_dir() and str(project_root) not in sys.path:
        sys.path.insert(0, str(project_root))

    import reconcILS as package

    # The legacy engine imports ``utils_reconcILS`` as a top-level package.
    package_dir = Path(package.__file__).resolve().parent
    package_path = str(package_dir)
    # Keep the project/site-packages entry ahead of the package directory. If
    # package_dir were first, a spawned worker could import reconcILS.py as a
    # top-level module and shadow the reconcILS package.
    while package_path in sys.path:
        sys.path.remove(package_path)
    sys.path.append(package_path)

    from reconcILS import reconcils
    from reconcILS.utils_reconcILS.readWrite import readWrite

    return reconcils, readWrite, package_dir


def terminal_classes():
    """Load the optional Rich presentation helpers from a checkout or installation."""
    project_root = Path(__file__).resolve().parent.parent
    source_package = project_root / "reconcILS"
    if source_package.is_dir() and str(project_root) not in sys.path:
        sys.path.insert(0, str(project_root))
    from reconcILS.rich_display import ParallelDashboard, RunDisplay

    return ParallelDashboard, RunDisplay


def positive_int(value: str) -> int:
    number = int(value)
    if number < 1:
        raise argparse.ArgumentTypeError("must be positive")
    return number


def positive_float(value: str) -> float:
    number = float(value)
    if not math.isfinite(number) or number <= 0:
        raise argparse.ArgumentTypeError("must be finite and positive")
    return number


def split_newick_records(text: str, label: str) -> list[str]:
    """Split a file containing one or more semicolon-terminated Newick trees."""
    stripped = text.strip()
    if not stripped:
        raise ValueError(f"{label} is empty")
    if not stripped.endswith(";"):
        raise ValueError(f"{label} must end with a semicolon")
    records = [part.strip() + ";" for part in stripped.split(";") if part.strip()]
    if not records:
        raise ValueError(f"{label} contains no Newick trees")
    return records


def read_species_tree(path: Path) -> tuple[str, str]:
    text = path.read_text(encoding="utf-8")
    records = split_newick_records(text, "Species-tree file")
    if len(records) != 1:
        raise ValueError("Species-tree file must contain exactly one Newick tree")
    return records[0], text


def read_gene_trees(path: Path) -> tuple[list[str], str]:
    text = path.read_text(encoding="utf-8")
    return split_newick_records(text, "Gene-tree file"), text


def prepare_species_tree(
    newick: str,
    supplied_mapping: Mapping[str, str],
    newick_format: int,
) -> tuple[str, SpeciesCopyMap]:
    """Create the shared deterministic alias map and prepare one species tree."""
    copy_map = SpeciesCopyMap.from_species_newick(
        newick,
        newick_format=newick_format,
        supplied=supplied_mapping,
    )
    return copy_map.prepare_species_newick(newick, newick_format=newick_format), copy_map


def prepare_gene_tree(newick: str, preparation: Mapping[str, Any]) -> str:
    """Preserve the supplied copy suffix while translating its species prefix."""
    copy_map = SpeciesCopyMap.from_payload(preparation["copy_map"])
    return copy_map.prepare_gene_newick(
        newick,
        newick_format=int(preparation["gene_format"]),
    )


def reconcile(
    gene_newick: str,
    species_newick: str,
    duplication: float,
    loss: float,
    ils: float,
    verbose: int,
    *,
    ils_timeout: float = 5000.0,
    include_edges: bool = False,
) -> dict[str, Any]:
    reconcils, readWrite, _ = runtime_classes()
    red = readWrite()
    reco = reconcils()
    reco.D_cost, reco.L_cost, reco.I_cost = duplication, loss, ils
    reco.time_ils = ils_timeout

    # Resolve nonbinary inputs with parse_bio, then serialize and parse again so
    # labels such as A_1 and A_2 both map to species A.
    gene_newick = red.to_newick(red.parse_bio(gene_newick))
    species_newick = red.to_newick(red.parse_bio(species_newick))
    tr, sp = red.parse(gene_newick), red.parse(species_newick)
    sp_copy = copy.deepcopy(sp)
    tr.order_gene(sp)
    tr.label_internal()
    sp.label_internal()
    tr.map_gene(sp)
    reco.setCost(sp)
    sp.isRoot = tr.isRoot = sp_copy.isRoot = True
    reco.gene_tree = copy.deepcopy(tr)
    reco.V = verbose

    started = time.perf_counter()
    events = reco.iterative_reconcILS(tr, sp, sp_copy, sp, [])
    elapsed = time.perf_counter() - started
    if not events:
        raise RuntimeError(
            "reconcILS returned no event records. The internal ILS search may have timed out, "
            "or no valid reconciliation path was produced. Inspect this tree's worker log and "
            "rerun with --verbose 1 for algorithm details."
        )
    else:
        labeled = red.label_reconcILS(events)
    counts = Counter(event for edge_events in labeled.values() for event in edge_events)
    result: dict[str, Any] = {
        "Process": "reconcILS",
        "Gene_tree": gene_newick,
        "Species_Tree": species_newick,
        "Duplication": counts["D"],
        "NNI": counts["I"],
        "Loss": counts["L"],
        "Seconds": elapsed,
    }
    if include_edges:
        reference = red.parse(species_newick)
        reference.label_internal()
        edges = reco.get_edges(reference)
        edge_events, _, _, _ = red.put_it_in_sp_tree(edges, labeled, {})
        result["_edges"] = [[list(edge), edge_events[edge]] for edge in edges]
    return result


def process_gene(
    index: int,
    gene: str,
    species: str,
    costs: Sequence[float],
    preparation: Mapping[str, Any],
    verbose: int,
) -> dict[str, Any]:
    copy_map = SpeciesCopyMap.from_payload(preparation["copy_map"])
    prepared = prepare_gene_tree(gene, preparation)
    result = reconcile(
        prepared, species, *costs, verbose,
        ils_timeout=float(preparation.get("ils_timeout", 5000.0)),
        include_edges=True,
    )
    result["Gene_tree"] = copy_map.restore_gene_newick(result["Gene_tree"])
    result["Species_Tree"] = copy_map.restore_species_newick(result["Species_Tree"])
    return dict(result, Replicate=index)


def worker_loop(
    connection,
    function,
    species: str,
    costs: Sequence[float],
    preparation: Mapping[str, Any],
    verbose: int,
    state: str,
    keep_logs: bool,
) -> None:
    """Reuse spawned workers and retain diagnostics for failed inputs."""
    signal.signal(signal.SIGINT, signal.SIG_IGN)
    while True:
        try:
            task = connection.recv()
        except EOFError:
            break
        if task is None:
            break
        index, gene = task
        log_path = Path(state) / f"tree_{index}.log"
        try:
            with log_path.open("w", buffering=1, encoding="utf-8") as log:
                with contextlib.redirect_stdout(log), contextlib.redirect_stderr(log):
                    row = function(index, gene, species, costs, preparation, verbose)
            result = dict(status="ok", replicate=index, edges=row.pop("_edges"), row=row)
        except Exception as exc:
            internal_timeout = False
            if log_path.exists():
                with log_path.open(encoding="utf-8") as stream:
                    internal_timeout = any(
                        "The function call timed out after " in line for line in stream
                    )
            result = dict(
                status="internal_timeout" if internal_timeout else "error",
                replicate=index,
                error=str(exc),
                traceback=traceback.format_exc(),
                log=str(log_path),
            )
            if internal_timeout:
                result["error"] += "; reconcILS reported an internal ILS timeout (see log)"
        if result["status"] == "ok" and not keep_logs:
            log_path.unlink(missing_ok=True)
        connection.send(result)
    connection.close()


def complete_record(result: Mapping[str, Any] | None) -> bool:
    return result is not None and result["status"] == "ok" and "edges" in result


def default_labeled_path(output: Path) -> Path:
    stem = output.stem.removesuffix("_result")
    return output.with_name(stem + "_labeled_sp.csv")


def default_alias_path(output: Path) -> Path:
    stem = output.stem.removesuffix("_result")
    return output.with_name(stem + "_species_aliases.json")


def aggregate_species(records: Sequence[Mapping[str, Any]], species_newick: str):
    reconcils, readWrite, _ = runtime_classes()
    red = readWrite()
    reco = reconcils()
    species = red.parse(species_newick)
    species.label_internal()
    edges = reco.get_edges(species)
    reco.edge_to_event(species, {edge: dict(D=0, I=0, L=0) for edge in edges}, 0)
    for result in records:
        reco.edge_to_event(
            species,
            {tuple(key): dict(value) for key, value in result["edges"]},
            0,
        )
    return species


def export_results(
    output: Path,
    labeled_output: Path,
    records: Mapping[int, Mapping[str, Any]],
    count: int,
    species_newick: str,
    copy_map_payload: Mapping[str, Any],
) -> int:
    """Write ordered results, stopping at the first incomplete record."""
    _, readWrite, _ = runtime_classes()
    output.parent.mkdir(parents=True, exist_ok=True)
    labeled_output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_name(output.name + ".tmp")
    written = 0
    included: list[Mapping[str, Any]] = []
    with temporary.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDS)
        writer.writeheader()
        for index in range(count):
            result = records.get(index)
            if not complete_record(result):
                break
            assert result is not None
            writer.writerow(result["row"])
            included.append(result)
            written += 1

    copy_map = SpeciesCopyMap.from_payload(copy_map_payload)
    tree_temporary = labeled_output.with_name(labeled_output.name + ".tmp")
    internal_labeled = readWrite().to_newick_sp(aggregate_species(included, species_newick))
    tree_temporary.write_text(
        copy_map.restore_species_newick(internal_labeled),
        encoding="utf-8",
    )
    os.replace(temporary, output)
    os.replace(tree_temporary, labeled_output)
    return written


def load_journal(journal, metadata: Mapping[str, Any]) -> dict[int, dict[str, Any]]:
    records: dict[int, dict[str, Any]] = {}
    journal.seek(0)
    header = journal.readline()
    if header:
        if json.loads(header) != metadata:
            raise ValueError(
                "Inputs, species tree, label mapping, costs, or algorithm changed since this "
                "journal; use a different --state-dir or --output for a fresh run."
            )
        while True:
            offset = journal.tell()
            line = journal.readline()
            if not line:
                break
            if not line.endswith("\n"):
                journal.seek(offset)
                journal.truncate()
                break
            result = json.loads(line)
            records[int(result["replicate"])] = result
    else:
        journal.write(json.dumps(metadata, sort_keys=True) + "\n")
        journal.flush()
        os.fsync(journal.fileno())
    journal.seek(0, os.SEEK_END)
    return records


def close_worker(process, connection) -> None:
    if process.is_alive():
        process.kill()
    process.join()
    connection.close()


def run_parallel(
    args: argparse.Namespace,
    genes: Sequence[str],
    species: str,
    preparation: Mapping[str, Any],
    metadata: Mapping[str, Any],
    copy_map_payload: Mapping[str, Any],
    function=process_gene,
) -> None:
    output = args.output.resolve()
    labeled_output = args.labeled_species_output.resolve()
    state = args.state_dir.resolve()
    state.mkdir(parents=True, exist_ok=True)
    count = len(genes)
    costs = (args.duplication_cost, args.loss_cost, args.ils_cost)
    ParallelDashboard, _ = terminal_classes()
    dashboard = ParallelDashboard(
        total=count,
        workers=args.workers,
        output=output,
        state=state,
        enabled=bool(args.verbose),
    )

    with (state / "run.lock").open("a", encoding="utf-8") as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise RuntimeError(f"Another runner is using state directory {state}") from exc

        with (state / "results.jsonl").open("a+", encoding="utf-8") as journal:
            records = load_journal(journal, metadata)
            failed = [
                record
                for index, record in records.items()
                if index < count and record["status"] != "ok"
            ]
            if failed and not args.retry_failed:
                result = min(failed, key=lambda record: record["replicate"])
                export_results(output, labeled_output, records, count, species, copy_map_payload)
                raise RuntimeError(
                    f"Gene tree {result['replicate']} previously failed: {result['error']}. "
                    "Use --retry-failed to retry it."
                )

            pending = iter(
                (index, gene)
                for index, gene in enumerate(genes)
                if not complete_record(records.get(index))
            )
            legacy = sum(
                record["status"] == "ok" and "edges" not in record
                for index, record in records.items()
                if index < count
            )

            context = mp.get_context("spawn")
            active: dict[int, tuple[Any, Any, float]] = {}
            idle: list[tuple[Any, Any]] = []
            exhausted = False
            failure: dict[str, Any] | None = None
            interrupted = False
            since_export = 0
            last_progress = time.monotonic()
            saved = sum(
                complete_record(record)
                for index, record in records.items()
                if index < count
            )
            processed = saved
            failed_this_run = 0
            dashboard.start(completed=processed, saved=saved, active=0)
            if legacy:
                dashboard.log(
                    f"Recomputing {legacy} older results without branch-event records "
                    "to build the labeled species tree.",
                    "yellow" if dashboard.rich else None,
                )

            def stop(signum, frame):
                raise KeyboardInterrupt

            previous_term = signal.signal(signal.SIGTERM, stop)
            written = saved
            try:
                while active or not exhausted:
                    while len(active) < args.workers and not exhausted:
                        try:
                            index, gene = next(pending)
                        except StopIteration:
                            exhausted = True
                            break
                        if idle:
                            process, receiver = idle.pop()
                        else:
                            receiver, sender = context.Pipe(duplex=True)
                            process = context.Process(
                                target=worker_loop,
                                args=(
                                    sender,
                                    function,
                                    species,
                                    costs,
                                    preparation,
                                    args.verbose,
                                    str(state),
                                    args.worker_logs,
                                ),
                            )
                            process.start()
                            sender.close()
                        receiver.send((index, gene))
                        active[index] = (process, receiver, time.monotonic())
                    dashboard.update(
                        completed=processed,
                        saved=saved,
                        active=len(active),
                        failed=failed_this_run,
                    )

                    progressed = False
                    for index, (process, receiver, started) in list(active.items()):
                        if failure is not None and index > failure["replicate"]:
                            close_worker(process, receiver)
                            del active[index]
                            continue

                        result = None
                        reusable = False
                        if receiver.poll():
                            try:
                                result = receiver.recv()
                                reusable = True
                            except EOFError:
                                result = dict(
                                    status="error",
                                    replicate=index,
                                    error="Worker exited without returning a result",
                                )
                        elif args.timeout is not None and time.monotonic() - started >= args.timeout:
                            result = dict(
                                status="timeout",
                                replicate=index,
                                error=f"Exceeded the requested {args.timeout:g}-second outer timeout",
                                log=str(state / f"tree_{index}.log"),
                            )
                        elif not process.is_alive():
                            result = dict(
                                status="error",
                                replicate=index,
                                error=f"Worker exited with code {process.exitcode}",
                            )
                        if result is None:
                            continue

                        if reusable:
                            idle.append((process, receiver))
                        else:
                            close_worker(process, receiver)
                        del active[index]
                        journal.write(json.dumps(result) + "\n")
                        journal.flush()
                        os.fsync(journal.fileno())
                        records[index] = result
                        since_export += 1
                        processed += 1
                        progressed = True
                        if result["status"] == "ok":
                            saved += 1
                        else:
                            failed_this_run += 1
                            if failure is None or index < failure["replicate"]:
                                failure = result
                            exhausted = True
                            log_hint = f"; log={result['log']}" if result.get("log") else ""
                            dashboard.log(
                                f"Tree {index}: {result['status']}: {result['error']}{log_hint}",
                                "red" if dashboard.rich else None,
                            )
                        dashboard.update(
                            completed=processed,
                            saved=saved,
                            active=len(active),
                            failed=failed_this_run,
                        )

                    if since_export >= args.checkpoint_every:
                        written = export_results(output, labeled_output, records, count, species, copy_map_payload)
                        dashboard.update(
                            completed=processed,
                            saved=saved,
                            active=len(active),
                            failed=failed_this_run,
                            checkpoint_rows=written,
                        )
                        dashboard.log(
                            f"Checkpoint: {written} ordered rows; {len(active)} active",
                            "cyan" if dashboard.rich else None,
                        )
                        since_export = 0
                        last_progress = time.monotonic()
                    if time.monotonic() - last_progress >= args.progress_every:
                        dashboard.update(
                            completed=processed,
                            saved=saved,
                            active=len(active),
                            failed=failed_this_run,
                        )
                        if not dashboard.rich:
                            dashboard.log(f"{saved}/{count} saved; {len(active)} active")
                        last_progress = time.monotonic()
                    if not progressed and active:
                        time.sleep(args.poll_interval)
            except KeyboardInterrupt:
                interrupted = True
                dashboard.log(
                    "Stopping workers and saving completed results; rerun to resume.",
                    "yellow" if dashboard.rich else None,
                )
            finally:
                for process, receiver, _ in list(active.values()) + [
                    (process, receiver, 0) for process, receiver in idle
                ]:
                    close_worker(process, receiver)
                try:
                    written = export_results(output, labeled_output, records, count, species, copy_map_payload)
                    dashboard.update(
                        completed=processed,
                        saved=saved,
                        active=0,
                        failed=failed_this_run,
                        checkpoint_rows=written,
                    )
                finally:
                    dashboard.close()
                    signal.signal(signal.SIGTERM, previous_term)

            if interrupted:
                dashboard.finish(
                    success=False,
                    warning=True,
                    message=(
                        "The run was interrupted after saving completed results.\n\n"
                        f"Ordered rows saved: {written}\n"
                        f"Output: {output}\n"
                        f"State directory: {state}\n"
                        "Run the same command again to resume."
                    ),
                )
                raise SystemExit(130)
            if failure is not None:
                log_hint = f"; worker log: {failure['log']}" if failure.get("log") else ""
                raise RuntimeError(
                    f"Gene tree {failure['replicate']} failed: {failure['error']}. "
                    f"Diagnostics directory: {state}{log_hint}"
                )
            dashboard.finish(
                success=True,
                message=(
                    "Parallel reconciliation completed successfully.\n\n"
                    f"Results: {output}\n"
                    f"Labeled species tree: {labeled_output}\n"
                    f"Species/copy map: {args.alias_output}\n"
                    f"Completed gene trees: {written}/{count}"
                ),
            )

def algorithm_fingerprint() -> str:
    _, _, package_dir = runtime_classes()
    digest = hashlib.sha256()
    paths = [
        package_dir / "reconcILS.py",
        package_dir / "copy_labels.py",
        *sorted((package_dir / "utils_reconcILS").glob("*.py")),
        Path(__file__).resolve(),
    ]
    for path in paths:
        digest.update(path.name.encode("utf-8"))
        digest.update(path.read_bytes())
    return digest.hexdigest()


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--gene-trees",
        type=Path,
        required=True,
        help="File of gene trees whose leaves are SpeciesName_copy",
    )
    parser.add_argument(
        "--species-tree",
        type=Path,
        required=True,
        help="File containing one species tree whose leaves are SpeciesName",
    )
    parser.add_argument("--output", type=Path, required=True, help="Result CSV path")
    parser.add_argument(
        "--labeled-species-output",
        type=Path,
        help="Labeled species-tree output; default is derived from --output",
    )
    parser.add_argument(
        "--alias-output",
        type=Path,
        help="Deterministic species/copy mapping JSON; default is derived from --output",
    )
    parser.add_argument(
        "--state-dir",
        type=Path,
        help="Checkpoint/log directory; default is a hidden directory beside --output",
    )
    parser.add_argument(
        "--mapping",
        type=Path,
        help=(
            "Optional JSON, CSV, or TSV species-to-alias mapping. CSV/TSV columns are "
            "species,alias. Missing mappings are generated automatically."
        ),
    )
    parser.add_argument("--species-format", type=int, default=1, help="ETE3 species-tree input format")
    parser.add_argument("--gene-format", type=int, default=1, help="ETE3 gene-tree input format")
    parser.add_argument("--limit", type=positive_int, help="Run only the first N gene trees")
    parser.add_argument("--duplication-cost", type=positive_float, default=1.1)
    parser.add_argument("--loss-cost", type=positive_float, default=1.0)
    parser.add_argument("--ils-cost", type=positive_float, default=1.0)
    allocated = int(os.environ.get("SLURM_CPUS_PER_TASK", "1"))
    parser.add_argument(
        "--workers",
        type=positive_int,
        default=max(1, min(8, allocated)),
        help="Worker processes (default: min(8, SLURM_CPUS_PER_TASK))",
    )
    parser.add_argument(
        "--library-threads",
        type=positive_int,
        default=1,
        help="Native numerical-library threads per worker (default: 1)",
    )
    parser.add_argument("--checkpoint-every", type=positive_int, default=100)
    parser.add_argument(
        "--progress-every",
        type=positive_float,
        default=30.0,
        help="Seconds between progress messages (default: 30)",
    )
    parser.add_argument(
        "--poll-interval",
        type=positive_float,
        default=0.05,
        help="Worker polling interval in seconds (default: 0.05)",
    )
    parser.add_argument(
        "--ils-timeout", type=positive_float, default=5000.0,
        help="Internal ILS search timeout in seconds (default: 5000)",
    )
    parser.add_argument(
        "--timeout",
        type=positive_float,
        help="Optional hard seconds per tree; default: no outer timeout",
    )
    parser.add_argument(
        "--verbose",
        type=int,
        choices=(0, 1),
        default=0,
        help="Show the DAFT-style Rich dashboard and retain verbose algorithm logs",
    )
    parser.add_argument("--retry-failed", action="store_true")
    parser.add_argument(
        "--worker-logs",
        action="store_true",
        help="Retain logs for successful trees as well as failed trees",
    )
    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    configure_native_threads(args.library_threads)
    _, RunDisplay = terminal_classes()
    display = RunDisplay(bool(args.verbose))

    for label, path in (
        ("Gene-tree file", args.gene_trees),
        ("Species-tree file", args.species_tree),
    ):
        if not path.is_file():
            parser.error(f"{label} does not exist: {path}")
    if args.mapping is not None and not args.mapping.is_file():
        parser.error(f"Mapping file does not exist: {args.mapping}")
    args.output = args.output.resolve()
    args.labeled_species_output = (
        args.labeled_species_output.resolve()
        if args.labeled_species_output is not None
        else default_labeled_path(args.output).resolve()
    )
    args.alias_output = (
        args.alias_output.resolve()
        if args.alias_output is not None
        else default_alias_path(args.output).resolve()
    )
    args.state_dir = (
        args.state_dir.resolve()
        if args.state_dir is not None
        else args.output.with_name("." + args.output.name + ".state").resolve()
    )

    # Route generated artifacts to DAFT-style folders before processing.
    # The checkpoint directory is stable across retries and never moved.
    from reconcILS.output_folders import folders
    dirs = folders(args.output)
    args.output = dirs["results"] / args.output.name
    if args.labeled_species_output.parent == dirs["results"].parent:
        args.labeled_species_output = dirs["results"] / args.labeled_species_output.name
    if args.alias_output.parent == dirs["results"].parent:
        args.alias_output = dirs["extras"] / args.alias_output.name
    if args.state_dir.parent == dirs["results"].parent:
        args.state_dir = dirs["log"] / args.state_dir.name
    for directory in dirs.values():
        directory.mkdir(parents=True, exist_ok=True)

    protected_inputs = {args.gene_trees.resolve(), args.species_tree.resolve()}
    if args.mapping is not None:
        protected_inputs.add(args.mapping.resolve())
    output_paths = (
        ("--output", args.output),
        ("--labeled-species-output", args.labeled_species_output),
        ("--alias-output", args.alias_output),
    )
    for label, path in output_paths:
        if path in protected_inputs:
            parser.error(f"{label} must differ from all input files")
    if len({path for _, path in output_paths}) != len(output_paths):
        parser.error("--output, --labeled-species-output, and --alias-output must be different")

    try:
        species_raw, species_file_text = read_species_tree(args.species_tree)
        genes, gene_file_text = read_gene_trees(args.gene_trees)
        supplied_mapping = load_supplied_mapping(args.mapping)
        species, copy_map = prepare_species_tree(
            species_raw,
            supplied_mapping,
            args.species_format,
        )
        copy_map.write_json(args.alias_output)
        copy_map_payload = copy_map.to_payload()
        preparation = {
            "copy_map": copy_map_payload,
            "gene_format": args.gene_format,
            "ils_timeout": args.ils_timeout,
        }

        # Validate every selected tree before workers start so formatting errors
        # are reported deterministically rather than after partial computation.
        selected_genes = genes[: args.limit]
        for index, gene in enumerate(selected_genes):
            try:
                copy_map.prepare_gene_newick(gene, newick_format=args.gene_format)
            except ValueError as exc:
                raise ValueError(f"Gene tree {index}: {exc}") from exc

        metadata = {
            "version": 3,
            "gene_input_sha256": hashlib.sha256(gene_file_text.encode("utf-8")).hexdigest(),
            "species_input_sha256": hashlib.sha256(species_file_text.encode("utf-8")).hexdigest(),
            "prepared_species_tree": species,
            "copy_map": copy_map_payload,
            "species_format": args.species_format,
            "gene_format": args.gene_format,
            "costs": [args.duplication_cost, args.loss_cost, args.ils_cost],
            "ils_timeout": args.ils_timeout,
            "verbose": args.verbose,
            "algorithm_sha256": algorithm_fingerprint(),
        }
        display.panel(
            "reconcILS batch run",
            "Starting parameterized parallel reconciliation.\n\n"
            f"Gene trees selected: {len(selected_genes)}\n"
            f"Workers: {args.workers}\n"
            f"Output: {args.output}",
            "green",
        )
        display.table(
            "Batch inputs and options",
            ["Item", "Value"],
            [
                ["Species-tree file", args.species_tree.resolve()],
                ["Gene-tree file", args.gene_trees.resolve()],
                ["Mapping file", args.mapping.resolve() if args.mapping else "automatic"],
                ["Species taxa", len(copy_map.original_to_internal)],
                ["Gene trees selected", len(selected_genes)],
                ["Duplication cost", args.duplication_cost],
                ["Loss cost", args.loss_cost],
                ["ILS cost", args.ils_cost],
                ["Internal ILS timeout (s)", args.ils_timeout],
                ["Outer per-tree timeout (s)", args.timeout if args.timeout is not None else "unlimited"],
                ["Workers", args.workers],
                ["Library threads per worker", args.library_threads],
                ["Result CSV", args.output],
                ["Labeled species tree", args.labeled_species_output],
                ["Species/copy map", args.alias_output],
                ["State directory", args.state_dir],
            ],
        )
        run_parallel(args, selected_genes, species, preparation, metadata, copy_map_payload)
    except (OSError, RuntimeError, ValueError, json.JSONDecodeError) as exc:
        message = (
            f"{exc}\n\n"
            f"Species-tree file: {args.species_tree}\n"
            f"Gene-tree file: {args.gene_trees}\n"
            f"State directory: {args.state_dir}"
        )
        display.panel("reconcILS batch run failed", message, "red")
        if args.verbose:
            raise SystemExit(1) from exc
        parser.exit(1, f"{exc}\n")


if __name__ == "__main__":
    main()
