"""Public exception hierarchy."""


class DecepticonXError(RuntimeError):
    """Base exception for all pipeline failures."""


class InputValidationError(DecepticonXError, ValueError):
    """An input violates the documented data contract."""


class BackendUnavailableError(DecepticonXError, ImportError):
    """An optional algorithm backend is not installed or usable."""


class BackendExecutionError(DecepticonXError):
    """An installed backend failed while processing valid inputs."""


class ConsensusError(DecepticonXError):
    """No valid consensus can be formed from backend estimates."""
