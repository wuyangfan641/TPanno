#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
End-to-end smoke tests (Biopython/NumPy/SciPy must be installed).

Run from the project root:
    python3 test/test_pipeline.py

Coverage:
  * 00_prepare_inputs.py: cleaning, U/O -> X, short-sequence filtering, duplicate-ID error,
    and UTF-8 BOM handling
  * 03_extract_domains.py: domtblout parsing, HMM coverage filtering, domain coordinate extraction
  * bin/lib/cdhit_map.py: CD-HIT representative mapping when the representative is not the first member
  * 05_cluster_tree.py: patristic distances, tree cutting, monophyly/support parsing
  * 06_make_report.py: summary report generation and directory-name mapping
  * scripts/add_swissprot_annotation.py: BOM-prefixed FASTA and best-hit selection
"""
from __future__ import annotations

import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PY = sys.executable


def run(args, cwd, env_extra=None, check=True):
    env = os.environ.copy()
    if env_extra:
        env.update({k: str(v) for k, v in env_extra.items()})
    proc = subprocess.run([PY, *map(str, args)], cwd=cwd, env=env,
                          stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if check and proc.returncode != 0:
        raise AssertionError(f"command failed: {args}\nSTDOUT:\n{proc.stdout}\nSTDERR:\n{proc.stderr}")
    return proc


AA = "ACDEFGHIKLMNPQRSTVWY"


def seq(n=120):
    return (AA * ((n // len(AA)) + 1))[:n]


class TestPrepare(unittest.TestCase):
    def test_clean_and_filter(self):
        with tempfile.TemporaryDirectory() as td:
            td = Path(td)
            raw = td / "raw"
            raw.mkdir()
            (raw / "S1.pep.fasta").write_text(
                ">p1\n" + seq(40) + "*\n"
                ">p2\n" + seq(39) + "U\n"
                ">short\nMKTA\n"
                ">highx\n" + "X" * 60 + "\n",
                encoding="utf-8",
            )
            out = td / "prep"
            run([ROOT / "bin/00_prepare_inputs.py", "--raw-dir", raw, "--outdir", out,
                 "--write-samples", td / "samples.tsv",
                 "--min-len", 30, "--max-x", 0.5], cwd=ROOT)
            clean = out / "clean_fasta/S1.faa"
            text = clean.read_text(encoding="utf-8")
            self.assertEqual(text.count(">"), 2)
            self.assertNotIn("*", text)
            self.assertIn("X", text.split(">p2")[1])
            self.assertTrue((out / "input_qc.md").exists())
            self.assertTrue((td / "samples.tsv").exists())

    def test_utf8_bom_is_accepted(self):
        with tempfile.TemporaryDirectory() as td:
            td = Path(td)
            raw = td / "raw"
            raw.mkdir()
            payload = (">p1\n" + seq(40) + "\n>p2\n" + seq(40) + "\n").encode("utf-8")
            (raw / "S1.pep.fasta").write_bytes(b"\xef\xbb\xbf" + payload)
            out = td / "prep"
            run([ROOT / "bin/00_prepare_inputs.py", "--raw-dir", raw, "--outdir", out,
                 "--min-len", 30], cwd=ROOT)
            text = (out / "clean_fasta/S1.faa").read_text(encoding="utf-8")
            self.assertEqual(text.count(">"), 2)

    def test_uo_residues_can_be_retained(self):
        with tempfile.TemporaryDirectory() as td:
            td = Path(td)
            raw = td / "raw"
            raw.mkdir()
            (raw / "S1.pep.fasta").write_text(
                ">p1\n" + ("ACDEFGUOHIKLMNPQRSTVWY" * 2) + "\n",
                encoding="utf-8",
            )
            out = td / "prep"
            run([ROOT / "bin/00_prepare_inputs.py", "--raw-dir", raw, "--outdir", out,
                 "--no-uo-to-x", "--min-len", 30], cwd=ROOT)
            text = (out / "clean_fasta/S1.faa").read_text(encoding="utf-8")
            self.assertIn("U", text)
            self.assertIn("O", text)

    def test_duplicate_ids_fail(self):
        with tempfile.TemporaryDirectory() as td:
            td = Path(td)
            raw = td / "raw"
            raw.mkdir()
            (raw / "S1.pep.fasta").write_text(
                ">p1\n" + seq(40) + "\n>p1\n" + seq(40) + "\n", encoding="utf-8")
            proc = run([ROOT / "bin/00_prepare_inputs.py", "--raw-dir", raw,
                        "--outdir", td / "prep"], cwd=ROOT, check=False)
            self.assertNotEqual(proc.returncode, 0)
            self.assertIn("duplicate", proc.stderr)


class TestExtract(unittest.TestCase):
    def _make(self, td: Path):
        raw = td / "clean"
        raw.mkdir()
        samples = []
        for s in range(1, 5):
            fa = raw / f"S{s}.faa"
            with open(fa, "w", encoding="utf-8") as fh:
                for i in range(1, 5):
                    fh.write(f">prot_S{s}_{i}\n{seq(120)}\n")
            samples.append((f"S{s}", fa))
        with open(td / "samples.clean.tsv", "w", encoding="utf-8") as fh:
            fh.write("sample\tfasta\tn_raw\tn_clean\n")
            for s, fa in samples:
                fh.write(f"{s}\t{fa}\t4\t4\n")

        search = td / "search"
        search.mkdir()
        for s, _fa in samples:
            with open(search / f"PF_TEST__{s}.domtblout", "w", encoding="utf-8") as fh:
                fh.write("# target name accession tlen query name accession qlen E-value score bias # of c-Evalue i-Evalue score bias from to from to from to acc description\n")
                for i in range(1, 5):
                    fh.write(
                        f"prot_{s}_{i} - 120 PF_TEST.1 - 60 1e-30 120.0 0.0 1 1 "
                        f"1e-32 1e-30 110.0 0.0 5 55 30 80 28 82 0.95 test\n")
                # Low coverage: should be filtered out
                fh.write("prot_S1_1 - 120 PF_TEST.1 - 60 1e-20 50.0 0.0 2 2 "
                         "1e-22 1e-20 40.0 0.0 1 20 30 50 28 52 0.50 lowcov\n")
                # Weak E-value: should be filtered out
                fh.write("prot_S1_2 - 120 PF_TEST.1 - 60 1e-3 20.0 0.0 1 1 "
                         "1e-3 1e-3 20.0 0.0 5 55 30 80 28 82 0.90 weak\n")
        return search

    def test_extract_domains(self):
        with tempfile.TemporaryDirectory() as td:
            td = Path(td)
            search = self._make(td)
            out = td / "domains"
            run([ROOT / "bin/03_extract_domains.py", "--search-dir", search,
                 "--samples", td / "samples.clean.tsv", "--outdir", out], cwd=ROOT,
                env_extra={"MIN_SEQS": 4, "MIN_SAMPLES": 2, "HMM_COV_MIN": 0.70,
                           "HMM_EVALUE_MAX": 1e-5, "DOMAIN_MIN_LEN": 30})
            fam = out / "PF_TEST.1"
            faa = fam / "PF_TEST.1.domains.faa"
            self.assertTrue(faa.exists())
            self.assertEqual(faa.read_text().count(">"), 16)
            meta = (fam / "metadata.tsv").read_text().splitlines()
            self.assertEqual(len(meta), 17)
            self.assertIn("0.85", meta[1])
            # Extract ali 30..80, which is 51 aa long
            lines = faa.read_text().splitlines()
            self.assertEqual(len(lines[1]), 51)

    def test_bom_sample_table(self):
        with tempfile.TemporaryDirectory() as td:
            td = Path(td)
            search = self._make(td)
            table = td / "samples.clean.tsv"
            table.write_bytes(b"\xef\xbb\xbf" + table.read_bytes())
            out = td / "domains_bom"
            run([ROOT / "bin/03_extract_domains.py", "--search-dir", search,
                 "--samples", table, "--outdir", out], cwd=ROOT,
                env_extra={"MIN_SEQS": 4, "MIN_SAMPLES": 2, "HMM_COV_MIN": 0.70,
                           "HMM_EVALUE_MAX": 1e-5, "DOMAIN_MIN_LEN": 30})
            self.assertTrue((out / "PF_TEST.1/PF_TEST.1.domains.faa").exists())


class TestScreen(unittest.TestCase):
    def test_screen_and_clan_collapse(self):
        with tempfile.TemporaryDirectory() as td:
            td = Path(td)
            dom = td / "domains"
            dom.mkdir()
            (dom / "families.tsv").write_text(
                "family\tdir\tn_hits\tn_domains\tn_seqs\tn_samples\tn_missing\tstatus\n"
                "PF00001.1\tPF00001.1\t3\t3\t3\t3\t0\tok\n"
                "PF00002.1\tPF00002.1\t6\t6\t6\t3\t0\tok\n"
                "PF00003.1\tNA\t2\tNA\t2\t2\t0\tskipped: proteins<4\n",
                encoding="utf-8")
            for fam, multi in [("PF00001.1", False), ("PF00002.1", True)]:
                d = dom / fam
                d.mkdir()
                with open(d / "metadata.tsv", "w", encoding="utf-8") as fh:
                    fh.write("seq_id\tfamily\tsample\tprotein_id\tdomain_no\tn_domains\t"
                             "i_evalue\tscore\thmm_cov\thmm_from\thmm_to\tali_from\tali_to\t"
                             "env_from\tenv_to\tprotein_len\tdomain_len\tx_frac\tdescription\tsearch_file\n")
                    for i, s in enumerate(["S1", "S2", "S3"], 1):
                        fh.write(f"{s}__p1\t{fam}\t{s}\tp1\t1\t1\t1e-30\t100\t0.9\t5\t55\t"
                                 "30\t80\t28\t82\t120\t51\t0\t-\tS1.domtblout\n")
                        if multi:
                            fh.write(f"{s}__p2\t{fam}\t{s}\tp2\t1\t1\t1e-30\t100\t0.9\t5\t55\t"
                                     "30\t80\t28\t82\t120\t51\t0\t-\tS1.domtblout\n")
            (td / "clans.tsv").write_text("PF00001.1\tCL0001\tA\nPF00002.1\tCL0001\tA\n",
                                          encoding="utf-8")
            run([ROOT / "bin/03b_screen_families.py", "--domain-dir", dom,
                 "--outdir", td / "screen", "--min-samples", 3, "--min-seqs", 3,
                 "--max-families", 10, "--clan-map", td / "clans.tsv",
                 "--collapse-clans"], cwd=ROOT)
            selected = (td / "screen/family_list.txt").read_text(encoding="utf-8").split()
            self.assertEqual(selected, ["PF00001.1"])
            screen = (td / "screen/family_screen.tsv").read_text(encoding="utf-8")
            self.assertIn("clan_redundant:CL0001", screen)
            self.assertIn("03_skipped", screen)

    def test_five_column_clan_map(self):
        """The real Pfam-A.clans.tsv has accession, clan, clan name, family, description."""
        with tempfile.TemporaryDirectory() as td:
            td = Path(td)
            dom = td / "domains"
            dom.mkdir()
            (dom / "families.tsv").write_text(
                "family\tdir\tn_hits\tn_domains\tn_seqs\tn_samples\tn_missing\tstatus\n"
                "PF00001.1\tPF00001.1\t3\t3\t3\t3\t0\tok\n",
                encoding="utf-8")
            d = dom / "PF00001.1"
            d.mkdir()
            with open(d / "metadata.tsv", "w", encoding="utf-8") as fh:
                fh.write("seq_id\tfamily\tsample\tprotein_id\thmm_cov\tdomain_len\n")
                for s in ["S1", "S2", "S3"]:
                    fh.write(f"{s}__p1\tPF00001.1\t{s}\tp1\t0.9\t51\n")
            clan_text = ("# accession\tclan\tclan_name\tfamily\tdescription\n"
                         "PF00001\tCL0192\tGPCR_A\t7tm_1\t7 transmembrane receptor\n")
            (td / "clans.tsv").write_bytes(b"\xef\xbb\xbf" + clan_text.encode("utf-8"))
            run([ROOT / "bin/03b_screen_families.py", "--domain-dir", dom,
                 "--outdir", td / "screen", "--min-samples", 3, "--min-seqs", 3,
                 "--max-families", 10, "--clan-map", td / "clans.tsv",
                 "--collapse-clans"], cwd=ROOT)
            screen = (td / "screen/family_screen.tsv").read_text(encoding="utf-8")
            self.assertIn("\tCL0192\t", screen)


class TestDerep(unittest.TestCase):
    def test_representative_after_first_member(self):
        with tempfile.TemporaryDirectory() as td:
            td = Path(td)
            clstr = td / "demo.clstr"
            clstr.write_text(
                ">Cluster 0\n"
                "0\t100aa, >alpha... *\n"
                "1\t99aa, >beta... at 99%\n"
                ">Cluster 1\n"
                "0\t100aa, >gamma... at 95%\n"
                "1\t100aa, >delta... *\n",
                encoding="utf-8",
            )
            out = td / "derep.tsv"
            run([ROOT / "bin/lib/cdhit_map.py", "--clstr", clstr, "--out", out], cwd=ROOT)
            rows = [line.split("\t") for line in out.read_text(encoding="utf-8").splitlines()[1:]]
            self.assertIn(["alpha", "alpha"], rows)
            self.assertIn(["alpha", "beta"], rows)
            self.assertIn(["delta", "gamma"], rows)
            self.assertIn(["delta", "delta"], rows)


class TestCluster(unittest.TestCase):
    def test_cluster(self):
        with tempfile.TemporaryDirectory() as td:
            td = Path(td)
            tree = td / "t.nwk"
            tree.write_text("((A:0.05,B:0.05)90/100:0.1,(C:0.05,D:0.05)95/100:0.1);\n")
            meta = td / "metadata.tsv"
            with open(meta, "w", encoding="utf-8") as fh:
                fh.write("seq_id\tfamily\tsample\tprotein_id\tdomain_no\ti_evalue\n")
                for i, name in enumerate(["A", "B", "C", "D"], 1):
                    fh.write(f"{name}\tPF_TEST.1\tS{i}\tp{i}\t1\t1e-30\n")
            prefix = td / "out"
            run([ROOT / "bin/05_cluster_tree.py", "--tree", tree, "--metadata", meta,
                 "--out-prefix", prefix, "--cutoff", 0.2, "--scan", "0.1,0.2,0.5"],
                cwd=ROOT)
            rows = (td / "out.clusters.tsv").read_text().splitlines()
            self.assertEqual(len(rows), 5)
            self.assertEqual(len({r.split("\t")[1] for r in rows[1:]}), 2)
            summ = (td / "out.cluster_summary.tsv").read_text()
            self.assertIn("90", summ)
            self.assertIn("95", summ)
            self.assertIn("yes", summ)
            self.assertTrue((td / "out.cutoff_scan.tsv").exists())


class TestReport(unittest.TestCase):
    def test_report(self):
        with tempfile.TemporaryDirectory() as td:
            td = Path(td)
            (td / "03_domains").mkdir()
            (td / "03_domains/families.tsv").write_text(
                "family\tdir\tn_hits\tn_domains\tn_seqs\tn_samples\tn_missing\tstatus\n"
                "PF_TEST.1\tPF_TEST.1\t16\t16\t16\t4\t0\tok\n")
            (td / "04_phylo").mkdir()
            (td / "04_phylo/all_status.tsv").write_text(
                "family\tstatus\tn_input\tn_used\tn_aln_cols\tn_trim_cols\ttrim_fraction\tbest_model\tnotes\n"
                "PF_TEST.1\tok\t16\t16\t80\t60\t0.750\tLG+F+G4\taligner=mafft\n")
            (td / "05_clusters").mkdir()
            (td / "05_clusters/all_cluster_summary.tsv").write_text(
                "family\tcluster\tn_seqs\tn_samples\tsamples\tsample_counts\tmean_within\tmax_within\tmin_between\tmean_silhouette\tmonophyletic\tSH_aLRT\tUFBoot\trepresentative\trepresentative_sample\tpasses_min_size\n"
                "PF_TEST.1\tC01\t8\t4\tS1,S2,S3,S4\tS1:2,S2:2,S3:2,S4:2\t0.1\t0.2\t0.5\t0.7\tyes\t90\t100\tA\tS1\tyes\n")
            run([ROOT / "bin/06_make_report.py", "--results", td], cwd=ROOT)
            report = (td / "06_report/report.md").read_text(encoding="utf-8")
            self.assertIn("PF_TEST.1", report)
            self.assertTrue((td / "06_report/family_summary.tsv").exists())
            params = (td / "06_report/parameters.tsv").read_text(encoding="utf-8")
            self.assertIn("SEARCH_MODE\t", params)
            self.assertNotIn("SEARCH_MODE\tNA", params)

    def test_sanitized_directory_is_mapped_back(self):
        with tempfile.TemporaryDirectory() as td:
            td = Path(td)
            (td / "03_domains").mkdir()
            (td / "03_domains/families.tsv").write_text(
                "family\tdir\tn_hits\tn_domains\tn_seqs\tn_samples\tn_missing\tstatus\n"
                "PF_TEST.1\tPF_TEST_1\t8\t8\t8\t4\t0\tok\n")
            (td / "04_phylo").mkdir()
            (td / "04_phylo/all_status.tsv").write_text(
                "family\tstatus\tn_input\tn_used\tn_aln_cols\tn_trim_cols\ttrim_fraction\tbest_model\tnotes\n"
                "PF_TEST_1\tok\t8\t8\t80\t60\t0.750\tLG+F+G4\taligner=mafft\n")
            run([ROOT / "bin/06_make_report.py", "--results", td], cwd=ROOT)
            summary = (td / "06_report/family_summary.tsv").read_text(encoding="utf-8")
            self.assertIn("PF_TEST.1\tok\t", summary)


class TestSwissProt(unittest.TestCase):
    @unittest.skipUnless((ROOT.parent / "scripts/add_swissprot_annotation.py").exists(),
                         "top-level scripts/ directory not present")
    def test_bom_and_best_hit(self):
        with tempfile.TemporaryDirectory() as td:
            td = Path(td)
            query = td / "protein.faa"
            payload = (">p1\n" + seq(20) + "\n>p2\n" + seq(20) + "\n>p3\n" + seq(20) + "\n").encode()
            query.write_bytes(b"\xef\xbb\xbf" + payload)
            blast = td / "protein.blastp.tsv"
            blast.write_text(
                "p1\tsp|P1|A\t90.0\t20\t1e-50\t200.0\tdesc one\n"
                "p2\tsp|P2|B\t40.0\t20\t1e-10\t50.0\tdesc two\n"
                "p2\tsp|P2|C\t95.0\t20\t1e-20\t150.0\tdesc three\n",
                encoding="utf-8")
            out = td / "out.faa"
            summary = td / "out.tsv"
            run([ROOT.parent / "scripts/add_swissprot_annotation.py",
                 "--query", query, "--blast", blast, "--out", out,
                 "--summary", summary, "--mark-missing"], cwd=ROOT)
            text = out.read_text(encoding="utf-8")
            self.assertIn(">p1 swissprot=sp|P1|A", text)
            self.assertIn(">p2 swissprot=sp|P2|C", text)
            self.assertIn(">p3 swissprot=no_hit", text)
            self.assertIn("sp|P2|C", summary.read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main(verbosity=2)
