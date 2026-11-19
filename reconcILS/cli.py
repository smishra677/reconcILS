"""DAFT_Test-style command line interface for the reconcILS engine."""
from __future__ import annotations

import argparse
import contextlib
import csv
import importlib
from io import StringIO
import os
from pathlib import Path
import sys
import tempfile

from .copy_labels import SpeciesCopyMap, load_supplied_mapping
from .rich_display import RunDisplay


def _existing_file(value: str) -> bool:
    """Return whether *value* names a file without rejecting long inline Newick."""
    try:
        return Path(value).is_file()
    except OSError:
        return False


def _single_tree_text(value: str, *, is_file: bool, label: str) -> str:
    text = Path(value).read_text(encoding="utf-8").strip() if is_file else value.strip()
    if not text:
        raise ValueError(f"{label} is empty")
    if not text.endswith(";") or text.count(";") != 1:
        raise ValueError(f"{label} must contain exactly one semicolon-terminated Newick tree")
    return text


def _default_alias_path(output: Path) -> Path:
    return output.with_name(output.stem + "_species_aliases.json")


def _restore_csv(path: Path, copy_map: SpeciesCopyMap, *, gene_key: str, species_key: str) -> None:
    if not path.is_file():
        return
    with path.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        fieldnames = list(reader.fieldnames or [])
        rows = list(reader)
    if not fieldnames:
        return
    for row in rows:
        if row.get(gene_key):
            row[gene_key] = copy_map.restore_gene_newick(row[gene_key])
        if row.get(species_key):
            row[species_key] = copy_map.restore_species_newick(row[species_key])
    temporary = path.with_name(path.name + ".tmp")
    with temporary.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    os.replace(temporary, path)


def _result_rows(output: str) -> list[list[str]]:
    path = Path(output)
    rows: list[list[str]] = [["Output CSV", str(path.resolve())]]
    if not path.is_file():
        return rows
    try:
        with path.open(newline="", encoding="utf-8") as stream:
            result = next(csv.DictReader(stream), None)
    except (OSError, csv.Error):
        return rows
    if result:
        for key in ("Duplication", "NNI", "Loss"):
            if key in result:
                rows.append([key, result[key]])
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(
        description="reconcILS: reconcile SpeciesName_copy gene leaves against SpeciesName leaves",
        epilog=(
            "Species-tree leaves must be SpeciesName; gene-tree leaves must be "
            "SpeciesName_copy. Use --verbose 1 for DAFT-style Rich displays."
        ),
    )
    species = parser.add_mutually_exclusive_group(required=True)
    species.add_argument("--sp", help="Species tree as Newick text (leaves: SpeciesName)")
    species.add_argument("--sp_file", help="File containing one species tree")
    parser.add_argument(
        "--gt",
        required=True,
        help="Gene tree Newick text or file (leaves: SpeciesName_copy)",
    )
    parser.add_argument("--output", required=True, help="Output CSV file name")
    parser.add_argument(
        "--mapping",
        type=Path,
        help="Optional JSON/CSV/TSV SpeciesName-to-two-character internal alias map",
    )
    parser.add_argument(
        "--alias-output",
        type=Path,
        help="Mapping report JSON; default: <output_stem>_species_aliases.json",
    )
    parser.add_argument("--species-format", type=int, default=1, help="ETE3 species-tree format")
    parser.add_argument("--gene-format", type=int, default=1, help="ETE3 gene-tree format")
    parser.add_argument("--D", type=float, help="Duplication cost")
    parser.add_argument("--L", type=float, help="Loss cost")
    parser.add_argument("--I", type=float, help="ILS cost")
    parser.add_argument("--ils-timeout", type=float, default=5000.0,
                        help="Internal ILS search timeout in seconds (default: 5000)")
    parser.add_argument(
        "--verbose",
        type=int,
        choices=(0, 1),
        default=0,
        help="Verbose/Rich display mode: 1=yes, 0=no",
    )
    args = parser.parse_args()
    import math
    if not math.isfinite(args.ils_timeout) or args.ils_timeout <= 0:
        parser.error("--ils-timeout must be a finite positive number of seconds")
    display = RunDisplay(bool(args.verbose))

    sp_is_file = args.sp_file is not None
    gt_is_file = _existing_file(args.gt)
    sp_value = args.sp_file if sp_is_file else args.sp
    assert sp_value is not None

    output = Path(args.output).resolve()
    if output.suffix.lower() != ".csv":
        parser.error("--output must end in .csv because the legacy engine derives related files from it")
    output.parent.mkdir(parents=True, exist_ok=True)
    alias_output = (
        args.alias_output.resolve() if args.alias_output is not None else _default_alias_path(output)
    )
    protected = {output, alias_output}
    input_paths = set()
    if sp_is_file:
        input_paths.add(Path(sp_value).resolve())
    if gt_is_file:
        input_paths.add(Path(args.gt).resolve())
    if args.mapping is not None:
        if not args.mapping.is_file():
            parser.error(f"Mapping file does not exist: {args.mapping}")
        input_paths.add(args.mapping.resolve())
    if protected & input_paths:
        parser.error("Output and alias-output paths must differ from all input files")
    if output == alias_output:
        parser.error("--output and --alias-output must be different files")

    species_source = str(Path(sp_value).resolve()) if sp_is_file else "inline --sp tree"
    gene_source = str(Path(args.gt).resolve()) if gt_is_file else "inline --gt tree"

    try:
        species_text = _single_tree_text(sp_value, is_file=sp_is_file, label="Species tree")
        gene_text = _single_tree_text(args.gt, is_file=gt_is_file, label="Gene tree")
        supplied = load_supplied_mapping(args.mapping)
        copy_map = SpeciesCopyMap.from_species_newick(
            species_text,
            newick_format=args.species_format,
            supplied=supplied,
        )
        prepared_species = copy_map.prepare_species_newick(
            species_text,
            newick_format=args.species_format,
        )
        prepared_gene = copy_map.prepare_gene_newick(
            gene_text,
            newick_format=args.gene_format,
        )
        copy_map.write_json(alias_output)
    except (OSError, ValueError) as exc:
        message = (
            f"reconcILS input preparation failed: {exc}\n\n"
            "Expected species leaves: SpeciesName\n"
            "Expected gene leaves: SpeciesName_copy"
        )
        display.panel("reconcILS input preparation failed", message, "red")
        parser.exit(2, message + "\n")

    from .reconcILS_validate import reconcILS_validate

    validator = reconcILS_validate()
    try:
        validator.validateme(
            species_text,
            gene_text,
            str(output),
            file_mode=False,
            verbose=bool(args.verbose),
            parser=args,
            species_source=species_source,
            gene_source=gene_source,
            copy_map=copy_map,
            species_format=args.species_format,
            gene_format=args.gene_format,
        )
    except ValueError as exc:
        if args.verbose:
            raise SystemExit(2) from exc
        parser.exit(2, f"reconcILS input validation failed: {exc}\n")

    display.panel(
        "reconcILS reconciliation",
        "Validation passed. Starting reconciliation.\n\n"
        f"Species input: {species_source}\n"
        f"Gene input: {gene_source}\n"
        f"Gene-label convention: SpeciesName_copy\n"
        f"Internal alias map: {alias_output}\n"
        f"Output: {output}",
        "blue",
    )

    with tempfile.TemporaryDirectory(prefix="reconcils_") as temporary:
        tempdir = Path(temporary)
        internal_species_path = tempdir / "species_internal.tre"
        internal_gene_path = tempdir / "gene_internal.tre"
        internal_species_path.write_text(prepared_species + "\n", encoding="utf-8")
        internal_gene_path.write_text(prepared_gene + "\n", encoding="utf-8")

        package_dir = str(Path(__file__).resolve().parent)
        if package_dir not in sys.path:
            sys.path.insert(0, package_dir)
        engine = importlib.import_module("reconcILS.reconcILS")
        legacy_args = [
            "reconcils",
            "--spTree",
            str(internal_species_path),
            "--gTree",
            str(internal_gene_path),
            "--output",
            str(output),
            "--F",
            "1",
        ]
        for name in ("D", "L", "I"):
            value = getattr(args, name)
            if value is not None:
                legacy_args.extend(["--" + name, str(value)])
        if args.verbose:
            legacy_args.extend(["--V", "1"])

        # Configure the existing algorithm instance without editing reconcILS.py.
        original_reconcils_class = engine.reconcils
        class ConfiguredReconcils(original_reconcils_class):
            def __init__(self):
                super().__init__()
                self.time_ils = args.ils_timeout
        engine.reconcils = ConfiguredReconcils
        old_argv = sys.argv
        engine_output = StringIO()
        try:
            sys.argv = legacy_args
            try:
                with contextlib.redirect_stdout(engine_output):
                    engine.main()
            except RuntimeError as exc:
                diagnostics = engine_output.getvalue().strip()
                message = (
                    "reconcILS could not complete the reconciliation.\n\n"
                    f"Reason: {exc}\n"
                    f"Species input: {species_source}\n"
                    f"Gene input: {gene_source}\n"
                    "Check that every gene leaf follows SpeciesName_copy and maps to the species tree.\n"
                    "Run again with --verbose 1 for algorithm details."
                )
                if diagnostics and args.verbose:
                    restored = copy_map.restore_gene_newick(
                        copy_map.restore_species_newick(diagnostics)
                    )
                    message += "\n\nAlgorithm output:\n" + restored
                display.panel("reconcILS reconciliation failed", message, "red")
                if args.verbose:
                    raise SystemExit(1) from exc
                parser.exit(1, message + "\n")
        finally:
            sys.argv = old_argv
            engine.reconcils = original_reconcils_class

        if args.verbose and engine_output.getvalue():
            restored = copy_map.restore_gene_newick(
                copy_map.restore_species_newick(engine_output.getvalue())
            )
            print(restored, end="" if restored.endswith("\n") else "\n")

    # The legacy engine saw only internal aliases. Restore original labels in
    # every user-facing CSV after it finishes.
    _restore_csv(output, copy_map, gene_key="Gene_tree", species_key="Species_Tree")
    legacy_log = output.with_name(output.stem + "_log.csv")
    _restore_csv(legacy_log, copy_map, gene_key="Gene_Tree", species_key="Species_Tree")

    # Organize completed artifacts like DAFT. Leave files in place on failure.
    from .output_folders import organize_single
    result_rows = _result_rows(str(output))
    try:
        directories, moved = organize_single(output, alias_output, validator)
    except (OSError, FileExistsError) as exc:
        parser.exit(1, f"Reconciliation succeeded but output organization failed: {exc}\n")
    result_output = moved.get(str(output), str(output))
    result_alias = moved.get(str(alias_output), str(alias_output))
    result_log = moved.get(str(Path(validator.log_path).resolve()), validator.log_path) if validator.log_path else None
    result_rows.extend(
        [
            ["Species/copy mapping", result_alias],
            ["Validation report", moved.get(str(Path(validator.validation_report_path).resolve()), validator.validation_report_path) if validator.validation_report_path else "not written"],
            ["Taxon/copy report", moved.get(str(Path(validator.taxon_report_path).resolve()), validator.taxon_report_path) if validator.taxon_report_path else "not written"],
            ["Run log", result_log or "not written"],
        ]
    )
    display.table("reconcILS results", ["Item", "Value"], result_rows)
    display.panel(
        "reconcILS run complete",
        "Reconciliation completed successfully.\n\n"
        "All output tree labels were restored to SpeciesName and SpeciesName_copy.\n"
        f"Results: {result_output}\n"
        f"Mapping: {result_alias}\n"
        f"Log: {result_log}",
        "green",
    )


if __name__ == "__main__":
    main()
