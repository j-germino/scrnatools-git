"""
sc-rna-tools package global settings manager
From scrnatools package

Created on Mon Jan 10 15:57:46 2022

@author: joe germino (joe.germino@ucsf.edu)
"""
# external package imports
import logging
import os

# scrnatools package imports
from ._utils import type_check


class Config:
    """Config manager for scrnatools"""

    def __init__(
            self,
            verbosity="warning",
            log_path="sc-rna-tools_logs",
            save_logs=False,
    ):
        self._verbosity = verbosity
        self._log_path = log_path
        self._save_logs = save_logs
        self._loggers = []  # list of loggers created by modules

    @property
    def verbosity(self) -> str:
        """The logging level, allowed are: 'debug', 'info', 'warning', 'error', 'critical'"""
        return self._verbosity

    @verbosity.setter
    def verbosity(self, verbosity: str):
        """Sets the verbosity level of the logger

        Allowed are: 'debug', 'info', 'warning', 'error', 'critical'
        """
        type_check(verbosity, "verbosity", str)
        # set the level of all the loggers
        for logger in self._loggers:
            self._set_log_level(logger, verbosity)
        self._verbosity = verbosity

    @property
    def log_path(self) -> str:
        """The path where log files are saved"""
        return self._log_path

    @log_path.setter
    def log_path(self, path: str):
        type_check(path, "path", str)
        self._log_path = path

    @property
    def save_logs(self) -> bool:
        """Determines whether logging messages are saved to a log file"""
        return self._save_logs

    @save_logs.setter
    def save_logs(self, save_logs: bool):
        type_check(save_logs, "save_logs", bool)
        # if logs are not currently being written to a file, but being updated to be written to a file, add a
        # FileHandler to each logger
        if not self.save_logs and save_logs:
            # make sure log directory exists
            self.check_log_path()
            # create the logging format and FileHandler
            formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s: %(message)s", "%Y-%m-%d %H:%M:%S")
            file_handler = logging.FileHandler(f"{self.log_path}/sc-rna-tools.log")
            file_handler.setFormatter(formatter)
            # add the FileHandler to each logger if it doesn't already have a FileHandler
            for logger in self._loggers:
                handlers = logger.handlers
                if sum([isinstance(h, logging.FileHandler) for h in handlers]) == 0:  # check if there is a FileHandler
                    logger.addHandler(file_handler)
        # if logs are currently being written to a file, but being updated to not be written to a file, remove the
        # FileHandler to each logger
        if self.save_logs and not save_logs:
            for logger in self._loggers:
                # find the FileHandler for this logger and remove it
                for h in logger.handlers:
                    if isinstance(h, logging.FileHandler):
                        logger.removeHandler(h)
        self._save_logs = save_logs

    def __str__(self) -> str:
        return f"ScoreConfig(" \
               f"verbosity: {self.verbosity}, " \
               f"log_path: {self.log_path}, " \
               f"save_logs: {self.save_logs})"

    def clear_logs(self):
        """Clears all logs saved in the logs path"""
        if os.path.exists(self._log_path):
            for f in os.listdir(self._log_path):
                # Remove all the files in the log directory
                os.remove(os.path.join(self._log_path, f))
            # remove the log directory
            os.rmdir(self._log_path)

    def check_log_path(self):
        """Checks if the log path exists, and creates it if it doesn't"""
        if not os.path.exists(self._log_path):
            os.makedirs(self._log_path)

    def create_logger(self, name: str):
        """Called by modules to create a logger for that function

        Parameters
        ----------
        name
            The name of the logger (typically the module name creating the logger)
        """
        type_check(name, 'name', str)
        # create a logger, StreamHandler, and the formatter for that logger
        logger = logging.getLogger(name)
        handler = logging.StreamHandler()
        formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s: %(message)s", "%Y-%m-%d %H:%M:%S")
        handler.setFormatter(formatter)
        handlers = logger.handlers
        if sum([isinstance(h, logging.StreamHandler) for h in handlers]) == 0:  # check if there is a StreamHandler
            logger.addHandler(handler)
        # if logs are being saved to a file as well, create a FileHandler for the new logger too
        if self.save_logs:
            self.check_log_path()
            file_handler = logging.FileHandler(f"{self.log_path}/sc-rna-tools.log")
            file_handler.setFormatter(formatter)
            if sum([isinstance(h, logging.FileHandler) for h in handlers]) == 0:  # check if there is a StreamHandler
                logger.addHandler(file_handler)
        # Set the log level to the current verbosity
        self._set_log_level(logger, self.verbosity)
        # Add the logger to the master list
        self._loggers.append(logger)
        return logger

    @staticmethod
    def _set_log_level(logger: logging.Logger, level: str):
        """Sets the log level of a modules logger

        Parameters
        ----------
        logger
            The logger to set the level of
        level
            String corresponding to the name of the logging level to set.
            Allowed: 'info', 'debug', 'warning', 'error', 'critical'
        """
        log_levels = {
            "info": logging.INFO,
            "debug": logging.DEBUG,
            "warning": logging.WARNING,
            "error": logging.ERROR,
            "critical": logging.CRITICAL
        }
        if level not in log_levels:
            raise ValueError(f"'{level}' is not a valid logging level "
                             f"('info', 'debug', 'warning', 'error', 'critical')")
        logger.setLevel(log_levels[level])


configs = Config()
