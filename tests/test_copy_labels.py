from __future__ import annotations

import json
from pathlib import Path
import tempfile
import unittest

from reconcILS.copy_labels import SpeciesCopyMap, load_supplied_mapping
from reconcILS.validation import validate


class SpeciesCopyMapTests(unittest.TestCase):
    def test_mapping_is_independent_of_newick_child_order(self):
        first = SpeciesCopyMap.from_species_newick("((Zeta,Alpha),Beta);")
        second = SpeciesCopyMap.from_species_newick("(Beta,(Alpha,Zeta));")
        self.assertEqual(first.original_to_internal, second.original_to_internal)
        self.assertEqual(
            first.original_to_internal,
            {"Alpha": "A0", "Beta": "A1", "Zeta": "A2"},
        )

    def test_species_prefixes_and_underscores_use_longest_match(self):
        mapping = SpeciesCopyMap.from_species_newick("(A,A_B);")
        prepared = mapping.prepare_gene_newick("(A_copy1,A_B_copy1);")
        self.assertEqual(prepared, "(A0_copy1,A1_copy1);")
        self.assertEqual(mapping.restore_gene_newick(prepared), "(A_copy1,A_B_copy1);")

    def test_copy_suffix_is_preserved_exactly(self):
        species = "((Homo_sapiens,Pan_troglodytes),Mus_musculus);"
        gene = (
            "((Homo_sapiens_copy_17,Homo_sapiens_paralog_B),"
            "(Pan_troglodytes_geneX,Mus_musculus_007));"
        )
        mapping = SpeciesCopyMap.from_species_newick(species)
        prepared = mapping.prepare_gene_newick(gene)
        self.assertIn("A0_copy_17", prepared)
        self.assertIn("A0_paralog_B", prepared)
        self.assertEqual(mapping.restore_gene_newick(prepared), gene)

    def test_missing_copy_and_unknown_species_are_informative(self):
        mapping = SpeciesCopyMap.from_species_newick("(SpeciesA,SpeciesB);")
        with self.assertRaisesRegex(ValueError, "missing a copy ID"):
            mapping.prepare_gene_newick("(SpeciesA,SpeciesB_copy1);")
        with self.assertRaisesRegex(ValueError, "does not match any species-tree leaf"):
            mapping.prepare_gene_newick("(SpeciesA_copy1,Unknown_copy1);")

    def test_labeled_species_output_is_restored(self):
        mapping = SpeciesCopyMap.from_species_newick("((Alpha,Beta),Gamma);")
        labeled = "((0-1-0   A0,2-0-1   A1),0-0-0   A2);"
        self.assertEqual(
            mapping.restore_species_newick(labeled),
            "((0-1-0   Alpha,2-0-1   Beta),0-0-0   Gamma);",
        )

    def test_mapping_report_round_trip(self):
        mapping = SpeciesCopyMap.from_species_newick("((Alpha,Beta),Gamma);")
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "mapping.json"
            mapping.write_json(path)
            payload = json.loads(path.read_text())
            loaded = load_supplied_mapping(path)
            self.assertEqual(loaded, mapping.original_to_internal)
            self.assertEqual(
                SpeciesCopyMap.from_payload(payload).original_to_internal,
                mapping.original_to_internal,
            )

    def test_validation_reports_copy_counts(self):
        species = "((Alpha,Beta),Gamma);"
        gene = "((Alpha_copy1,Alpha_copy2),(Beta_copy1,Gamma_copy1));"
        report = validate(species, gene)
        self.assertEqual(report["gene_copy_count"], 4)
        self.assertEqual(report["copies_per_species"]["Alpha"], 2)
        self.assertEqual(report["missing_species_taxa"], [])


if __name__ == "__main__":
    unittest.main()
