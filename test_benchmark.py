"""Check the example generator and output handling without a timing run."""
from collections import Counter
from fractions import Fraction as F
from pathlib import Path
import signal
import tempfile
import unittest

from benchmark import prepare, run, write_json
from check_precision import check
from extension_history import verify_witness
from monotone_history import CHARTS, decode


class BenchmarkTests(unittest.TestCase):
    def test_case_counts_and_ids(self):
        cases, truth = prepare()
        self.assertEqual(len(cases), 60)
        self.assertEqual(len(truth), 60)
        self.assertEqual(Counter(c["kind"] for c in cases),
                         {"nominal": 48, "ambiguous": 6, "incompatible": 6})
        self.assertEqual({c["case"] for c in cases}, set(truth))

    def test_generated_witnesses(self):
        cases, truth = prepare()
        for case in cases:
            for phases in truth[case["case"]]["histories"]:
                domain = CHARTS["negative" if phases[-1] < 0 else "positive"]
                self.assertTrue(verify_witness(phases, case["history"],
                                              case["family"], case["scale"], domain))

    def test_control_times(self):
        cases, _ = prepare()
        for case in cases:
            if case["kind"] != "nominal":
                self.assertEqual(case["history"]["times"], [F(-1, 200), F(1, 200)])

    def test_precision_check_on_generated_inputs(self):
        cases, truth = prepare()
        with tempfile.TemporaryDirectory() as folder:
            directory = Path(folder)
            write_json(directory / "inputs.json", cases)
            write_json(directory / "truth.json", truth)
            check(directory)

    @unittest.skipUnless(hasattr(signal, "setitimer"), "POSIX timer required")
    def test_does_not_overwrite_results(self):
        with tempfile.TemporaryDirectory() as folder:
            output = Path(folder) / "run"
            output.mkdir()
            sentinel = output / "keep.txt"
            sentinel.write_text("keep\n")
            with self.assertRaises(FileExistsError):
                run(output, check_only=True)
            self.assertEqual(sentinel.read_text(), "keep\n")
            self.assertEqual(list(output.iterdir()), [sentinel])

    def test_legacy_status_controls(self):
        # Manufactured controls are regression fixtures, not publication data.
        zero = [['0', '0', '0'], ['0', '0', '0']]
        controls = {
            'unique': {
                'times': ['-0.02', '0'],
                'positions': [['-0.0396', '-0.019404', '0'], ['0', '0', '0']],
                'linf_errors': ['0.0005', '0.0005'],
            },
            'ambiguous': {
                'times': ['-0.01', '0.01'], 'positions': zero,
                'linf_errors': ['0.06', '0.06'],
            },
            'incompatible': {
                'times': ['-0.01', '0.01'], 'positions': zero,
                'linf_errors': ['0.0005', '0.0005'],
            },
        }
        for name, history in controls.items():
            result = decode(history, F(1, 2))
            self.assertEqual(result["status"], name.upper())


if __name__ == "__main__":
    unittest.main()
