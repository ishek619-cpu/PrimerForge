"""
PrimerForge Logger

Provides a shared logger for all PrimerForge modules.

Author: Abhishek Warad
Version: 0.1.0
"""

from pathlib import Path
import logging


def get_logger(name: str) -> logging.Logger:
    """
    Create and return a configured logger.

    Parameters
    ----------
    name : str
        Name of the module requesting the logger.

    Returns
    -------
    logging.Logger
        Configured logger instance.
    """

    # Create logs directory
    log_dir = Path("outputs/logs")
    log_dir.mkdir(parents=True, exist_ok=True)

    log_file = log_dir / "primerforge.log"

    logger = logging.getLogger(name)

    # Prevent duplicate handlers
    if logger.hasHandlers():
        return logger

    logger.setLevel(logging.INFO)

    formatter = logging.Formatter(
        "%(asctime)s | %(levelname)s | %(name)s | %(message)s"
    )

    # Console handler
    console_handler = logging.StreamHandler()
    console_handler.setFormatter(formatter)

    # File handler
    file_handler = logging.FileHandler(log_file)
    file_handler.setFormatter(formatter)

    logger.addHandler(console_handler)
    logger.addHandler(file_handler)

    return logger
