# lotus_backend/src/infrastructure/logging.py
import logging
import sys
from pathlib import Path
from loguru import logger
from infrastructure.workspace_context import workspace_path_manager


class InterceptHandler(logging.Handler):
    """
    Redirect standard logging records to Loguru.
    """

    def emit(self, record):
        try:
            level = logger.level(record.levelname).name
        except ValueError:
            level = record.levelno

        # Find the original caller frame for correct line numbers
        frame, depth = logging.currentframe(), 2
        while frame and frame.f_code.co_filename == logging.__file__:
            frame = frame.f_back
            depth += 1

        logger.opt(depth=depth, exception=record.exc_info).log(
            level, record.getMessage()
        )


def setup_logging(debug_mode: bool = False) -> None:
    """
    Initialize application logging:
    - Console output via Loguru
    - Rotating log file in workspace
    - Intercept stdlib logging (scanpy, lotus_core, uvicorn)
    """

    # Resolve log directory (workspace-aware fallback)
    log_dir = workspace_path_manager.root / "logs" if workspace_path_manager.root else Path("./logs")
    log_dir.mkdir(exist_ok=True)
    log_file = log_dir / "lotus_app.log"

    # Remove default Loguru handlers
    logger.remove()

    # Console logger (for development)
    logger.add(
        sys.stderr,
        level="DEBUG" if debug_mode else "INFO",
        format=(
            "<green>{time:HH:mm:ss}</green> | "
            "<level>{level: <8}</level> | "
            "<cyan>{name}</cyan>:<cyan>{function}</cyan> - "
            "<level>{message}</level>"
        ),
    )

    # File logger (for diagnostics and user reports)
    logger.add(
        log_file,
        rotation="1 day",
        retention="10 days",
        compression="zip",
        level="INFO",
        encoding="utf-8",
    )

    # Intercept standard logging
    logging.basicConfig(
        handlers=[InterceptHandler()],
        level=0,
        force=True,
    )

    # Silence overly verbose libraries
    logging.getLogger("uvicorn.access").handlers = [InterceptHandler()]
    logging.getLogger("matplotlib").setLevel(logging.WARNING)

    logger.info(f"Logging initialized. Log file: {log_file}")
