## ------------------------------------------------------------------------------------------- ##
## Logging Utilities                                                                           ##
## ------------------------------------------------------------------------------------------- ##
## @script: logging_util.py                                                                    ##
##                                                                                             ##
## @description: Configure application-wide logging with colored warnings and                  ##
##               optional handler reset.                                                       ##
##                                                                                             ##
## @author: Yousef Mustafa, Lab of Dr. William Bush.                                           ##
## ------------------------------------------------------------------------------------------- ##

"""Custom logging configuration for pathWAS."""

import logging
import sys


def configure_logging(level: int = 5, target = sys.stderr) -> None:
    """Configure root logger with optional colored warnings.

    Parameters
    ----------
    level : int, optional
        Logging level for the root logger (default 5 which corresponds to custom
        DEBUG level).
    target : IO[str], optional
        Stream target for log output, defaults to ``sys.stderr``.
    """
    logger = logging.getLogger()
    logger.setLevel(level)

    kludge = False
    if logger.hasHandlers():
        kludge = True
        logger.handlers.clear()

    RED = "\033[91m"
    RESET = "\033[0m"

    class CustomFormatter(logging.Formatter):
        def format(self, record: logging.LogRecord) -> str:
            if record.levelno == logging.WARNING:
                record.msg = f"{RED}{record.msg}{RESET}"
            return super().format(record)

    handler = logging.StreamHandler(target)
    handler.setLevel(level)
    formatter = CustomFormatter("%(levelname)s - %(message)s")
    handler.setFormatter(formatter)
    logger.addHandler(handler)

    if kludge:
        logging.info("Had to reset logging handlers!!")
