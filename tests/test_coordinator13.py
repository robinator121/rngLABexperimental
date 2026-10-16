"""Tests for RNGLab v0.13+ Campaign Coordinator."""
import json
from pathlib import Path
import unittest
import tempfile

from rnglab13.coordinator import CampaignRunner, CampaignConfig


class CampaignCoordinatorTests(unittest.TestCase):
    def test_campaign_runner_cycle_and_atomic_artifacts(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            run_path = Path(tmp_dir)
            cfg = CampaignConfig(run_dir=run_path, daemon_mode=False, seeds=[101])
            runner = CampaignRunner(cfg)
            runner.run()

            # Verify artifacts
            status_file = run_path / "status.json"
            summary_file = run_path / "RESEARCH_SUMMARY.md"
            self.assertTrue(status_file.exists())
            self.assertTrue(summary_file.exists())

            with open(status_file, "r") as f:
                status = json.load(f)
            self.assertFalse(status["is_running"])
            self.assertGreater(status["completed_tasks"], 0)
            self.assertIn("math01", status["active_metrics"])
            self.assertIn("math02", status["active_metrics"])
            self.assertIn("math03", status["active_metrics"])

            # Verify results directory contains task files
            results = list((run_path / "results").glob("*.json"))
            self.assertGreater(len(results), 0)

    def test_campaign_runner_respects_stop_file(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            run_path = Path(tmp_dir)
            stop_file = run_path / "STOP"
            stop_file.touch()

            cfg = CampaignConfig(run_dir=run_path, daemon_mode=True, seeds=[101])
            runner = CampaignRunner(cfg)
            self.assertTrue(runner.check_stop_condition())
            runner.run()

            with open(run_path / "status.json", "r") as f:
                status = json.load(f)
            self.assertFalse(status["is_running"])
            self.assertTrue(status["stop_requested"])
            self.assertEqual(status["completed_tasks"], 0)


if __name__ == "__main__":
    unittest.main()
