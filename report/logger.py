import os
import logging
import sys
from typing import Optional, TextIO


def configure_logger(filename: str) -> None:
    """
    Configure logger file.
    :param filename: str - filename where to log
    :return: None
    """
    # create output directory if needed
    os.makedirs(os.path.dirname(filename), exist_ok=True)

    # configure logger
    logging.basicConfig(filename=filename, level=logging.DEBUG, format='%(asctime)s %(levelname)10s: %(message)s', filemode='w')
    logging.getLogger().setLevel(logging.INFO)


def log_str(to_log: str, stdout_too: Optional[TextIO] = sys.stderr, priority: int = logging.INFO, flush: bool = False) -> None:
    """
    Write a string to log.
    :param to_log: str - string to write to log
    :param stdout_too: bool - write the string to stdout too
    :param priority: int - logging priority
    :param flush: bool - flush this logger?
    :return: None
    """
    if stdout_too is not None:
        print(to_log, file=stdout_too)

    # delete \r and make it one-liners:
    to_log = to_log.replace('\r', '')
    to_log = to_log.replace('\n', '    ')

    # finally log it!
    logging.log(priority, to_log)

    # flush it?
    if flush:
        logging.getLogger().handlers[0].flush()
