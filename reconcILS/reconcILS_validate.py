"""DAFT-style automatic preflight for the reconcILS command-line wrapper."""
from __future__ import annotations

from pathlib import Path
from typing import Any, Iterable, Sequence

from .copy_labels import SpeciesCopyMap
from .rich_display import RunDisplay
from .validation import validate


class reconcILS_validate:
    """Validate and report before execution; never rewrite input trees."""

    def __init__(self) -> None:
        self.species_tree_newick = None
        self.gene_tree_newick = None
        self.input_hash = None
        self.taxon_report_path = None
        self.validation_report_path = None
        self.log_path = None
        self.console = None
        self.log = None
        self.verbose = False
        self.display = RunDisplay(False)

    def log_message(self, message: str, log) -> None:
        if self.verbose:
            self.display.message(message)
        print(message, file=log, flush=True)

    @staticmethod
    def log_only(message: str, log) -> None:
        print(message, file=log, flush=True)

    def log_header(self, title: str, log) -> None:
        self.log_only("", log)
        self.log_only("=" * 80, log)
        self.log_only(title, log)
        self.log_only("-" * 80, log)

    def log_table_only(
        self,
        title: str,
        columns: Sequence[str],
        rows: Iterable[Sequence[Any]],
        log,
    ) -> None:
        self.log_header(title, log)
        self.log_only("\t".join(str(value) for value in columns), log)
        for row in rows:
            self.log_only("\t".join(str(value) for value in row), log)
        self.log_only("*" * 80, log)

    def print_panel(self, title: str, message: str, log, style: str = "green") -> None:
        self.display.panel(title, message, style)
        self.log_header(title, log)
        self.log_only(message, log)
        self.log_only("*" * 80, log)

    def print_table(
        self,
        title: str,
        columns: Sequence[str],
        rows: Iterable[Sequence[Any]],
        log,
        *,
        status_column: int | None = None,
    ) -> None:
        materialized = [list(row) for row in rows]
        self.display.table(title, columns, materialized, status_column=status_column)
        self.log_table_only(title, columns, materialized, log)

    def print_startup_report(
        self,
        *,
        parser,
        species_source: str,
        gene_source: str,
        output_prefix: str,
        log,
    ) -> None:
        log_path = Path(output_prefix + "_reconcILS_log.txt").resolve()
        self.print_panel(
            "reconcILS",
            "Starting a new reconcILS run\n\n"
            f"Output: {output_prefix}\n"
            f"Log file: {log_path}",
            log,
            "green",
        )

        if parser is not None:
            option_rows = [[key, value] for key, value in vars(parser).items()]
            self.print_table("Run options", ["Option", "Value"], option_rows, log)

        self.print_table(
            "Input sources",
            ["Item", "Value"],
            [
                ["Species tree source", species_source],
                ["Gene tree source", gene_source],
                ["Expected gene leaf", "SpeciesName_copy"],
                ["Output", output_prefix],
                ["reconcILS log file", log_path],
            ],
            log,
        )

    def print_tree_summary(self, report: dict[str, Any], log) -> None:
        rows = [
            ["Species-tree taxa", report["species_taxa_count"]],
            ["Gene-tree copies", report["gene_copy_count"]],
            ["Species represented by gene copies", len(report["gene_species_taxa"])],
            ["Species absent from gene tree", len(report["missing_species_taxa"])],
        ]
        self.print_table("Tree summary", ["Item", "Value"], rows, log)

    def validateme(
        self,
        sp,
        gene,
        output_prefix,
        *,
        file_mode: bool = False,
        verbose: bool = False,
        parser=None,
        species_source: str | None = None,
        gene_source: str | None = None,
        copy_map: SpeciesCopyMap | None = None,
        species_format: int = 1,
        gene_format: int = 1,
    ):
        prefix = str(output_prefix or "reconcILS")
        log_path = Path(prefix + "_reconcILS_log.txt")
        log_path.parent.mkdir(parents=True, exist_ok=True)
        self.verbose = bool(verbose)
        self.display = RunDisplay(self.verbose)
        self.console = self.display.console

        source_label = "file" if file_mode else "inline Newick"
        species_source = species_source or source_label
        gene_source = gene_source or source_label

        with log_path.open("w", encoding="utf-8") as log:
            self.log = log
            try:
                self.print_startup_report(
                    parser=parser,
                    species_source=species_source,
                    gene_source=gene_source,
                    output_prefix=prefix,
                    log=log,
                )

                if not sp or not gene:
                    raise ValueError("Both species-tree and gene-tree inputs are required")
                if file_mode:
                    sp = Path(sp).read_text(encoding="utf-8").strip()
                    gene = Path(gene).read_text(encoding="utf-8").strip()

                report = validate(
                    sp,
                    gene,
                    prefix,
                    copy_map=copy_map,
                    species_format=species_format,
                    gene_format=gene_format,
                )
                taxon_path = Path(prefix + "_taxon_report.tsv")
                internal = report["mapping"]["original_to_internal"]
                copies = report["copies_per_species"]
                with taxon_path.open("w", encoding="utf-8") as out:
                    out.write("species\tinternal_alias\tgene_copy_count\tin_gene_tree\n")
                    for species in report["species_taxa"]:
                        count = int(copies.get(species, 0))
                        out.write(f"{species}\t{internal[species]}\t{count}\t{int(count > 0)}\n")

                validation_rows = [
                    ["Species tree parsed", "OK", f"{report['species_taxa_count']} species"],
                    ["Gene tree parsed", "OK", f"{report['gene_copy_count']} gene copies"],
                    ["Duplicate leaves", "OK", "No duplicate species or copy labels"],
                    ["Binary rooted trees", "OK", "Both trees are strictly binary"],
                    ["Gene label format", "OK", "Every leaf is SpeciesName_copy"],
                    ["Gene species mapping", "OK", "Every copy maps to a species-tree leaf"],
                    ["Internal aliases", "OK", "Deterministic fixed two-character IDs"],
                ]
                if report["missing_species_taxa"]:
                    validation_rows.append(
                        [
                            "Species absent from gene tree",
                            "WARNING",
                            ", ".join(report["missing_species_taxa"]),
                        ]
                    )
                else:
                    validation_rows.append(["Species absent from gene tree", "OK", "None"])
                validation_rows.extend(
                    [
                        ["Taxon/copy report", "OK", str(taxon_path)],
                        ["Input fingerprint", "OK", report["input_sha256"]],
                    ]
                )

                self.print_tree_summary(report, log)
                self.print_table(
                    "Validation checks",
                    ["Check", "Status", "Details"],
                    validation_rows,
                    log,
                    status_column=1,
                )
                self.print_panel(
                    "reconcILS validation complete",
                    "Input validation completed successfully.\n\n"
                    "Public labels: SpeciesName and SpeciesName_copy\n"
                    f"Taxon/copy report: {taxon_path}\n"
                    f"Validation report: {report.get('report_file', '')}\n"
                    f"Input hash: {report['input_sha256']}",
                    log,
                    "green",
                )

                self.species_tree_newick = sp
                self.gene_tree_newick = gene
                self.input_hash = report["input_sha256"]
                self.taxon_report_path = str(taxon_path)
                self.validation_report_path = str(report.get("report_file", ""))
                self.log_path = str(log_path)
                return self
            except (ValueError, OSError, IndexError) as exc:
                message = (
                    f"Error: {exc}\n\n"
                    "Please check that:\n"
                    "- each input contains one valid Newick tree ending in ';'\n"
                    "- both trees are strictly binary\n"
                    "- species-tree leaves are unique SpeciesName labels\n"
                    "- every gene-tree leaf is SpeciesName_copy\n"
                    "- each copy suffix is non-empty and each full gene label is unique\n"
                    "- every SpeciesName prefix occurs in the species tree"
                )
                self.print_panel("reconcILS input validation failed", message, log, "red")
                raise ValueError(f"{exc} (see {log_path})") from exc
