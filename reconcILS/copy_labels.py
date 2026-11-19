"""Deterministic species/copy labels for reconcILS.

Public input uses the convention::

    species tree leaf: SpeciesName
    gene tree leaf:    SpeciesName_copy

reconcILS itself receives fixed-width two-character species aliases and gene
leaves of the form ``AL_copy``.  User-facing output is restored to the original
species names.  All input renaming is performed on parsed tree leaves, never by
replacing arbitrary substrings in Newick text.
"""
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from io import StringIO
from pathlib import Path
import csv
import json
import re
from typing import Any, Iterable, Mapping, Sequence


COPY_SEPARATOR = "_"
_ALIAS_PATTERN = re.compile(r"^[A-Z0-9]{2}$")
# Unquoted Newick labels accepted by the compatibility layer.  Underscores,
# periods, dashes and accession punctuation are intentionally allowed.
_SAFE_NEWICK_TOKEN = re.compile(r"^[^\s(),:;\[\]'\"]+$")
_ALIAS_FIRST = "ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789"
_ALIAS_SECOND = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ"


class _TopologyTree:
    """Small adapter around ETE3 with a Biopython fallback for inspection/tests."""

    def __init__(self, newick: str, newick_format: int = 1) -> None:
        text = newick.strip()
        if not text:
            raise ValueError("Tree is empty")
        if not text.endswith(";"):
            raise ValueError("Tree must end with ';'")

        try:
            from ete3 import Tree as EteTree
        except ImportError:
            EteTree = None

        if EteTree is not None:
            try:
                self.kind = "ete3"
                self.tree = EteTree(text, format=newick_format)
                return
            except Exception as exc:
                raise ValueError(f"Invalid Newick tree: {exc}") from exc

        # The package requires ETE3 at runtime, but this fallback keeps the
        # mapping/validation layer usable in lightweight environments.
        try:
            from Bio import Phylo

            self.kind = "biopython"
            self.tree = Phylo.read(StringIO(text), "newick")
        except Exception as exc:
            raise ValueError(f"Invalid Newick tree: {exc}") from exc

    def leaves(self):
        if self.kind == "ete3":
            return list(self.tree.iter_leaves())
        return list(self.tree.get_terminals())

    def leaf_names(self) -> list[str]:
        return [str(leaf.name or "") for leaf in self.leaves()]

    def is_strictly_binary(self) -> bool:
        if self.kind == "ete3":
            return all(len(node.children) == 2 for node in self.tree.traverse() if not node.is_leaf())
        return all(len(node.clades) == 2 for node in self.tree.get_nonterminals())

    def write_topology(self) -> str:
        if self.kind == "ete3":
            return self.tree.write(format=9)

        def serialize(clade) -> str:
            if clade.is_terminal():
                name = str(clade.name or "")
                if not name:
                    raise ValueError("Tree has a leaf with no name")
                if not _SAFE_NEWICK_TOKEN.fullmatch(name):
                    raise ValueError(
                        f"Leaf label {name!r} requires Newick quoting, which is not supported "
                        "by the reconcILS compatibility layer"
                    )
                return name
            return "(" + ",".join(serialize(child) for child in clade.clades) + ")"

        return serialize(self.tree.root) + ";"


def alias_candidates() -> Iterable[str]:
    """Yield 1,296 fixed-width aliases, beginning A0, A1, ..., AZ, B0, ..."""
    for first in _ALIAS_FIRST:
        for second in _ALIAS_SECOND:
            yield first + second


def validate_alias(alias: str) -> None:
    if not _ALIAS_PATTERN.fullmatch(alias):
        raise ValueError(
            f"Invalid internal species alias {alias!r}; aliases must contain exactly "
            "two uppercase ASCII letters/digits and no underscore"
        )


def _validate_species_name(name: str) -> None:
    if not name:
        raise ValueError("Species-tree leaves must all have names")
    if not _SAFE_NEWICK_TOKEN.fullmatch(name):
        raise ValueError(
            f"Species name {name!r} is not a supported unquoted Newick token. "
            "Use letters, digits, underscores, periods, or dashes without whitespace "
            "or Newick delimiters."
        )


def _mapping_rows(path: Path) -> Iterable[tuple[str, str]]:
    text = path.read_text(encoding="utf-8")
    if not text.strip():
        raise ValueError(f"Mapping file is empty: {path}")

    if path.suffix.lower() == ".json":
        data = json.loads(text)
        if isinstance(data, dict) and "original_to_internal" in data:
            data = data["original_to_internal"]
        if not isinstance(data, dict):
            raise ValueError(
                "JSON mapping must be an object of species-to-alias pairs or contain "
                "an original_to_internal object"
            )
        for species, alias in data.items():
            yield str(species).strip(), str(alias).strip()
        return

    lines = [line for line in text.splitlines() if line.strip() and not line.lstrip().startswith("#")]
    if not lines:
        raise ValueError(f"Mapping file has no data rows: {path}")
    sample = "\n".join(lines[:20])
    try:
        delimiter = csv.Sniffer().sniff(sample, delimiters=",\t").delimiter
    except csv.Error:
        delimiter = "\t" if path.suffix.lower() in {".tsv", ".tab"} else ","

    rows = list(csv.reader(lines, delimiter=delimiter))
    if rows and len(rows[0]) >= 2:
        first = [cell.strip().lower() for cell in rows[0][:2]]
        if first[0] in {"species", "source", "original", "taxon"} and first[1] in {
            "alias",
            "internal",
            "target",
            "code",
            "label",
        }:
            rows = rows[1:]
    for line_number, row in enumerate(rows, 1):
        if len(row) != 2:
            raise ValueError(
                f"Mapping row {line_number} must have exactly two columns: species and alias"
            )
        yield row[0].strip(), row[1].strip()


def load_supplied_mapping(path: Path | None) -> dict[str, str]:
    """Read an optional JSON/CSV/TSV original-species to internal-alias map."""
    if path is None:
        return {}
    mapping: dict[str, str] = {}
    for species, alias in _mapping_rows(path):
        if not species or not alias:
            raise ValueError("Mapping species and alias values cannot be empty")
        if species in mapping and mapping[species] != alias:
            raise ValueError(f"Species {species!r} has more than one alias")
        validate_alias(alias)
        mapping[species] = alias
    if len(set(mapping.values())) != len(mapping):
        raise ValueError("Each supplied internal alias must identify exactly one species")
    return mapping


@dataclass(frozen=True)
class GeneCopy:
    species: str
    copy: str

    @property
    def public_label(self) -> str:
        return f"{self.species}{COPY_SEPARATOR}{self.copy}"


@dataclass(frozen=True)
class SpeciesCopyMap:
    """A deterministic, reversible map between public species and internal IDs."""

    original_to_internal: Mapping[str, str]
    copy_separator: str = COPY_SEPARATOR

    def __post_init__(self) -> None:
        if self.copy_separator != COPY_SEPARATOR:
            raise ValueError(
                "reconcILS currently requires '_' between SpeciesName and copy ID"
            )
        normalized = {str(species): str(alias) for species, alias in self.original_to_internal.items()}
        if len(normalized) < 2:
            raise ValueError("Species tree must contain at least two uniquely named species")
        for species, alias in normalized.items():
            _validate_species_name(species)
            validate_alias(alias)
        if len(set(normalized.values())) != len(normalized):
            raise ValueError("Each internal species alias must be unique")
        object.__setattr__(self, "original_to_internal", normalized)

    @property
    def internal_to_original(self) -> dict[str, str]:
        return {alias: species for species, alias in self.original_to_internal.items()}

    @classmethod
    def from_species_names(
        cls,
        species_names: Sequence[str],
        *,
        supplied: Mapping[str, str] | None = None,
    ) -> "SpeciesCopyMap":
        names = [str(name) for name in species_names]
        if len(names) < 2:
            raise ValueError("Species tree must contain at least two leaves")
        duplicates = sorted(name for name, count in Counter(names).items() if count > 1)
        if duplicates:
            raise ValueError("Species tree contains duplicate leaves: " + ", ".join(duplicates))
        for name in names:
            _validate_species_name(name)

        supplied_map = {str(key): str(value) for key, value in (supplied or {}).items()}
        unknown = sorted(set(supplied_map) - set(names))
        if unknown:
            raise ValueError(
                "Mapping contains species absent from the species tree: " + ", ".join(unknown)
            )
        for alias in supplied_map.values():
            validate_alias(alias)
        if len(set(supplied_map.values())) != len(supplied_map):
            raise ValueError("Each supplied internal alias must identify exactly one species")

        ordered = sorted(names)
        used = set(supplied_map.values())
        candidates = (candidate for candidate in alias_candidates() if candidate not in used)
        mapping: dict[str, str] = {}
        for species in ordered:
            if species in supplied_map:
                alias = supplied_map[species]
            else:
                try:
                    alias = next(candidates)
                except StopIteration as exc:
                    raise ValueError(
                        "The fixed two-character namespace supports at most 1,296 species"
                    ) from exc
                used.add(alias)
            mapping[species] = alias
        return cls(mapping)

    @classmethod
    def from_species_newick(
        cls,
        newick: str,
        *,
        newick_format: int = 1,
        supplied: Mapping[str, str] | None = None,
    ) -> "SpeciesCopyMap":
        tree = _TopologyTree(newick, newick_format)
        return cls.from_species_names(tree.leaf_names(), supplied=supplied)

    @classmethod
    def from_payload(cls, payload: Mapping[str, Any]) -> "SpeciesCopyMap":
        mapping = payload.get("original_to_internal")
        if not isinstance(mapping, Mapping):
            raise ValueError("Mapping payload is missing original_to_internal")
        separator = str(payload.get("copy_separator", COPY_SEPARATOR))
        return cls({str(key): str(value) for key, value in mapping.items()}, separator)

    def to_payload(self) -> dict[str, Any]:
        return {
            "version": 1,
            "ordering": "lexicographic by original species name",
            "alias_scheme": "fixed two-character uppercase base-36 identifiers",
            "gene_leaf_format": "SpeciesName_copy",
            "copy_separator": self.copy_separator,
            "species_count": len(self.original_to_internal),
            "original_to_internal": dict(self.original_to_internal),
            "internal_to_original": self.internal_to_original,
        }

    def write_json(self, path: Path) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(self.to_payload(), indent=2, sort_keys=True) + "\n", encoding="utf-8")

    def split_gene_label(self, label: str) -> GeneCopy:
        """Split exactly ``SpeciesName_copy``, using the longest species prefix."""
        text = str(label)
        matches = [
            species
            for species in self.original_to_internal
            if text.startswith(species + self.copy_separator)
        ]
        if not matches:
            if text in self.original_to_internal:
                raise ValueError(
                    f"Gene leaf {text!r} is missing a copy ID; expected {text}{self.copy_separator}copy"
                )
            raise ValueError(
                f"Gene leaf {text!r} does not match any species-tree leaf followed by "
                f"'{self.copy_separator}copy'"
            )
        longest = max(len(species) for species in matches)
        best = [species for species in matches if len(species) == longest]
        if len(best) != 1:
            raise ValueError(f"Gene leaf {text!r} has an ambiguous species prefix: {best}")
        species = best[0]
        copy_id = text[len(species) + len(self.copy_separator) :]
        if not copy_id:
            raise ValueError(
                f"Gene leaf {text!r} has an empty copy ID after '{self.copy_separator}'"
            )
        if not _SAFE_NEWICK_TOKEN.fullmatch(copy_id):
            raise ValueError(
                f"Copy ID {copy_id!r} in gene leaf {text!r} contains unsupported whitespace "
                "or Newick delimiters"
            )
        return GeneCopy(species, copy_id)

    def internal_gene_label(self, public_label: str) -> str:
        parsed = self.split_gene_label(public_label)
        return (
            self.original_to_internal[parsed.species]
            + self.copy_separator
            + parsed.copy
        )

    def public_gene_label(self, internal_label: str) -> str:
        alias, separator, copy_id = str(internal_label).partition(self.copy_separator)
        if not separator or not copy_id:
            raise ValueError(
                f"Internal gene label {internal_label!r} is not ALIAS{self.copy_separator}copy"
            )
        try:
            species = self.internal_to_original[alias]
        except KeyError as exc:
            raise ValueError(f"Unknown internal species alias {alias!r}") from exc
        return species + self.copy_separator + copy_id

    def prepare_species_newick(self, newick: str, *, newick_format: int = 1) -> str:
        tree = _TopologyTree(newick, newick_format)
        names = tree.leaf_names()
        if set(names) != set(self.original_to_internal) or len(names) != len(self.original_to_internal):
            missing = sorted(set(self.original_to_internal) - set(names))
            extra = sorted(set(names) - set(self.original_to_internal))
            raise ValueError(
                "Species tree does not match the mapping"
                + (f"; missing: {', '.join(missing)}" if missing else "")
                + (f"; unexpected: {', '.join(extra)}" if extra else "")
            )
        for leaf in tree.leaves():
            leaf.name = self.original_to_internal[str(leaf.name)]
        return tree.write_topology()

    def prepare_gene_newick(self, newick: str, *, newick_format: int = 1) -> str:
        tree = _TopologyTree(newick, newick_format)
        names = tree.leaf_names()
        duplicates = sorted(name for name, count in Counter(names).items() if count > 1)
        if duplicates:
            raise ValueError("Gene tree contains duplicate copy labels: " + ", ".join(duplicates))
        prepared: list[str] = []
        for leaf in tree.leaves():
            label = self.internal_gene_label(str(leaf.name))
            leaf.name = label
            prepared.append(label)
        if len(prepared) < 2:
            raise ValueError("Gene tree must contain at least two leaves")
        if len(set(prepared)) != len(prepared):
            raise ValueError("Gene-copy labels are not unique after internal mapping")
        return tree.write_topology()

    def restore_gene_newick(self, newick: str) -> str:
        """Restore internal ``AL_copy`` prefixes without touching event text."""
        reverse = self.internal_to_original
        aliases = sorted(reverse, key=len, reverse=True)
        pattern = re.compile(
            r"(?<![A-Za-z0-9_])(" + "|".join(re.escape(alias) for alias in aliases) + r")_"
        )
        return pattern.sub(lambda match: reverse[match.group(1)] + self.copy_separator, newick)

    def restore_species_newick(self, newick: str) -> str:
        """Restore internal species leaves, including event-labeled species trees."""
        reverse = self.internal_to_original
        aliases = sorted(reverse, key=len, reverse=True)
        pattern = re.compile(
            r"(?<![A-Za-z0-9_])(" + "|".join(re.escape(alias) for alias in aliases) + r")"
            r"(?=\s*[,);]|\s*$)"
        )
        return pattern.sub(lambda match: reverse[match.group(1)], newick)

    def gene_assignments(self, labels: Sequence[str]) -> list[GeneCopy]:
        return [self.split_gene_label(label) for label in labels]


def inspect_tree(newick: str, *, newick_format: int = 1) -> dict[str, Any]:
    """Return leaf names and binary status for validation/reporting."""
    tree = _TopologyTree(newick, newick_format)
    return {
        "leaf_names": tree.leaf_names(),
        "strictly_binary": tree.is_strictly_binary(),
    }
