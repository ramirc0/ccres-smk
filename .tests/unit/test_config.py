"""Tests of the config template."""

from pathlib import Path

TEMPLATE = Path(__file__).parents[2] / "config" / "config.yaml.template"


def test_step_3_excludes_no_experiment():
    """Leave no filtered TF experiment out of step 3."""
    assert "tf_excluded" not in TEMPLATE.read_text()
