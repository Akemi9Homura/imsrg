"""Lightweight helpers for loading the in-tree pyIMSRG extension."""

from __future__ import annotations

import importlib
from pathlib import Path
import sys
from typing import Any


def import_pyimsrg(module_dir: Path) -> Any:
    """Import pyIMSRG from one explicit build directory."""

    sys.path.insert(0, str(module_dir.resolve()))
    try:
        return importlib.import_module("pyIMSRG")
    finally:
        sys.path.pop(0)
