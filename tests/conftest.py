import sys
from pathlib import Path

# Make the synthetic-oracle generators importable in tests (they live outside the package on purpose:
# they are the correctness spec, not shipped runtime code).
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "synthetic"))
