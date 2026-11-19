"""Preflight validation for reconcILS SpeciesName/SpeciesName_copy inputs."""
from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
from pathlib import Path
import json
from typing import Any

from .copy_labels import SpeciesCopyMap, inspect_tree


def _duplicate_names(names: list[str], label: str) -> None:
    duplicate = sorted(name for name, count in Counter(names).items() if count > 1)
    if duplicate:
        raise ValueError(f"{label} has duplicate leaf labels: {', '.join(duplicate)}")


def validate(
    species: str,
    gene: str,
    output: str | None = None,
    *,
    copy_map: SpeciesCopyMap | None = None,
    species_format: int = 1,
    gene_format: int = 1,
) -> dict[str, Any]:
    """Validate one species tree and one gene tree.

    Species-tree leaves must be original species names.  Every gene-tree leaf
    must be ``SpeciesName_copy`` and its species prefix must occur in the
    species tree.  Copy suffixes are preserved exactly.
    """
    species_info = inspect_tree(species, newick_format=species_format)
    gene_info = inspect_tree(gene, newick_format=gene_format)
    species_names = list(species_info["leaf_names"])
    gene_labels = list(gene_info["leaf_names"])

    _duplicate_names(species_names, "Species tree")
    _duplicate_names(gene_labels, "Gene tree")
    if not species_info["strictly_binary"]:
        raise ValueError("Species tree must be strictly binary (including its root)")
    if not gene_info["strictly_binary"]:
        raise ValueError("Gene tree must be strictly binary (including its root)")

    mapping = copy_map or SpeciesCopyMap.from_species_names(species_names)
    if set(mapping.original_to_internal) != set(species_names):
        missing = sorted(set(species_names) - set(mapping.original_to_internal))
        extra = sorted(set(mapping.original_to_internal) - set(species_names))
        raise ValueError(
            "Species mapping does not match the supplied species tree"
            + (f"; unmapped: {', '.join(missing)}" if missing else "")
            + (f"; absent from tree: {', '.join(extra)}" if extra else "")
        )

    assignments = mapping.gene_assignments(gene_labels)
    gene_species = [assignment.species for assignment in assignments]
    copies_per_species = Counter(gene_species)
    missing_species = sorted(set(species_names) - set(gene_species))

    digest_payload = {
        "species": species.strip(),
        "gene": gene.strip(),
        "copy_separator": mapping.copy_separator,
        "mapping": mapping.to_payload(),
    }
    digest = sha256(json.dumps(digest_payload, sort_keys=True).encode("utf-8")).hexdigest()
    report: dict[str, Any] = {
        "status": "OK",
        "gene_leaf_format": "SpeciesName_copy",
        "copy_separator": mapping.copy_separator,
        "species_taxa": sorted(species_names),
        "gene_taxa": sorted(gene_labels),
        "gene_species_taxa": sorted(set(gene_species)),
        "species_taxa_count": len(species_names),
        "gene_taxa_count": len(gene_labels),
        "gene_copy_count": len(assignments),
        "copies_per_species": dict(sorted(copies_per_species.items())),
        "missing_species_taxa": missing_species,
        "input_sha256": digest,
        "mapping": mapping.to_payload(),
    }
    if output:
        path = Path(str(output) + "_validation.json")
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        report["report_file"] = str(path)
    return report


def validate_args(argv=None) -> None:
    parser = argparse.ArgumentParser(
        description="Validate reconcILS SpeciesName/SpeciesName_copy inputs without running reconciliation"
    )
    parser.add_argument("--spTree", required=True)
    parser.add_argument("--gTree", required=True)
    parser.add_argument("--F", default="0", choices=["0", "1"], help="Read trees from files")
    parser.add_argument("--output", help="Prefix for JSON validation report")
    args = parser.parse_args(argv)
    try:
        if args.F == "1":
            sp = Path(args.spTree).read_text(encoding="utf-8").strip()
            gt = Path(args.gTree).read_text(encoding="utf-8").strip()
        else:
            sp, gt = args.spTree, args.gTree
        report = validate(sp, gt, args.output)
    except (ValueError, OSError, IndexError) as exc:
        parser.exit(1, f"Validation FAILED: {exc}\n")
    print(json.dumps(report, indent=2, sort_keys=True))
