# lotus_backend/src/infrastructure/logging.py
import logging
import sys
from pathlib import Path
from loguru import logger
from infrastructure.workspace_context import workspace_path_manager


class InterceptHandler(logging.Handler):
    """Redirect stdlib logging records to Loguru."""

    def emit(self, record):
        try:
            level = logger.level(record.levelname).name
        except ValueError:
            level = record.levelno

        frame, depth = logging.currentframe(), 2
        while frame and frame.f_code.co_filename == logging.__file__:
            frame = frame.f_back
            depth += 1

        logger.opt(depth=depth, exception=record.exc_info).log(
            level, record.getMessage()
        )


def _configure_log_levels(
        app_packages: list[str],
        *,
        debug_mode: bool,
) -> None:
    """
    Silence noisy third-party libs (ERROR only) and enable app logs.
    """
    # 1. Default: Silence everything (allow only ERROR/CRITICAL)
    logging.getLogger().setLevel(logging.ERROR)

    # 2. Whitelist: Enable DEBUG/INFO for our application code
    app_level = logging.DEBUG if debug_mode else logging.INFO
    for pkg in app_packages:
        logging.getLogger(pkg).setLevel(app_level)


def setup_logging(debug_mode: bool = False) -> None:
    """
    Initialize logging with specific color mapping.
    """

    # 1. Resolve log path
    try:
        log_dir = workspace_path_manager.root / "logs"
    except RuntimeError:
        log_dir = Path("./logs")

    log_dir.mkdir(exist_ok=True, parents=True)
    log_file = log_dir / "lotus_app.log"

    # 2. Reset Loguru
    logger.remove()

    # 3. Define Color Scheme
    # INFO: White, DEBUG: Cyan, WARNING: Yellow, ERROR: Red
    logger.level("INFO", color="<white>")
    logger.level("DEBUG", color="<cyan>")
    logger.level("WARNING", color="<yellow>")
    logger.level("ERROR", color="<red>")

    # 4. Console Output
    # <level> tag applies the color defined above to the content inside
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

    # 5. File Output (No colors, more retention)
    logger.add(
        log_file,
        rotation="1 day",
        retention="10 days",
        compression="zip",
        level="DEBUG",
        encoding="utf-8",
        format="{time:YYYY-MM-DD HH:mm:ss} | {level: <8} | {name}:{function} - {message}"
    )

    # 6. Intercept Standard Library Logging
    logging.basicConfig(handlers=[InterceptHandler()], level=0, force=True)

    # 7. Apply Whitelist (Silence noisy libs)
    _configure_log_levels(
        app_packages=["lotus", "lotus_core", "lotus_backend"],
        debug_mode=debug_mode,
    )
