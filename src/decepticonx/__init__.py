"""decepticonX-py-rust public API."""

from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("decepticonx-py-rust")
except PackageNotFoundError:
    __version__ = "0.1.0"

from .exceptions import (
    BackendExecutionError,
    BackendUnavailableError,
    ConsensusError,
    DecepticonXError,
    InputValidationError,
)
from .models import BranchKey, DecepticonXResult, PipelineConfig
from .pipeline import run_decepticonx

__all__ = [
    "BackendExecutionError",
    "BackendUnavailableError",
    "BranchKey",
    "ConsensusError",
    "DecepticonXError",
    "DecepticonXResult",
    "InputValidationError",
    "PipelineConfig",
    "run_decepticonx",
    "__version__",
]
