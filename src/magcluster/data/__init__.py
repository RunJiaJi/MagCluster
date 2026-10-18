"""Access packaged data assets (HMMs, rules, profiles)."""

from pathlib import Path
from typing import Any, Dict
import yaml

DATA_DIR = Path(__file__).resolve().parent

def get_hmm_path() -> Path:
    """Return path to packaged magnetosome.hmm file."""
    p = DATA_DIR / "hmms" / "magnetosome.hmm"
    if not p.exists():
        raise FileNotFoundError(f"Magnetosome HMM database not found at {p}")
    return p

def get_rules_path() -> Path:
    """Return path to default rules.yaml."""
    p = DATA_DIR / "rules.yaml"
    if not p.exists():
        raise FileNotFoundError(f"Rules configuration not found at {p}")
    return p

def load_rules(custom_rules_path: Path = None) -> Dict[str, Any]:
    """Load default or custom rules configuration."""
    target = custom_rules_path or get_rules_path()
    with open(target, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)

def get_profiles_summary_path() -> Path:
    """Return path to profiles_summary.tsv."""
    return DATA_DIR / "profiles_summary.tsv"
