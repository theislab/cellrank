# modified from: https://github.com/theislab/scanpy/blob/master/scanpy/tests/test_logging.py
import datetime
import io
import sys

import pytest
from scanpy import Verbosity

from cellrank import logging as logg
from cellrank import settings


@pytest.fixture
def logging_state():  # noqa: PT004
    verbosity_orig = settings.verbosity
    yield
    settings.logfile = sys.stderr
    settings.verbosity = verbosity_orig


class TestLogging:
    def test_formats(self, capsys, logging_state):
        settings.logfile = sys.stderr
        settings.verbosity = Verbosity.debug
        logg.error("0")
        assert capsys.readouterr().err == "ERROR: 0\n"
        logg.warning("1")
        assert capsys.readouterr().err == "WARNING: 1\n"
        logg.info("2")
        assert capsys.readouterr().err == "2\n"
        logg.hint("3")
        assert capsys.readouterr().err == "--> 3\n"
        # TODO: this still uses scanpy's logger, but in e.g. notebooks it's fine
        # logg.debug("4")
        # assert capsys.readouterr().err == "DEBUG: 4\n"

    def test_deep(self, capsys, logging_state):
        settings.logfile = sys.stderr
        settings.verbosity = Verbosity.hint
        logg.hint("0")
        assert capsys.readouterr().err == "--> 0\n"
        logg.hint("1", deep="1!")
        assert capsys.readouterr().err == "--> 1\n"
        settings.verbosity = Verbosity.debug
        logg.hint("2")
        assert capsys.readouterr().err == "--> 2\n"
        logg.hint("3", deep="3!")
        assert capsys.readouterr().err == "--> 3: 3!\n"

    def test_logfile(self, tmp_path, logging_state):
        settings.verbosity = Verbosity.hint

        buffer = io.StringIO()
        settings.logfile = buffer
        assert settings.logfile is buffer
        assert settings.logpath is None
        logg.error("test!")
        assert buffer.getvalue() == "ERROR: test!\n"

        p = tmp_path / "test.log"
        settings.logpath = p
        assert settings.logpath == p
        assert settings.logfile.name == str(p)
        logg.hint("test2")
        logg.debug("invisible")
        assert settings.logpath.read_text() == "--> test2\n"

    def test_timing(self, monkeypatch, capsys, logging_state):
        import cellrank.logging._logging as logg

        class IncTime:
            def __init__(self):
                self.counter = 0

            def now(self, tz):
                self.counter += 1
                return datetime.datetime(2000, 1, 1, second=self.counter, microsecond=self.counter, tzinfo=tz)

            @property
            def datetime(self) -> "IncTime":
                return self

            @property
            def timezone(self):
                return datetime.timezone

        settings.logfile = sys.stderr
        t = IncTime()
        monkeypatch.setattr(logg, "datetime", t)
        settings.verbosity = Verbosity.debug

        logg.hint("1")
        assert t.counter == 1
        assert capsys.readouterr().err == "--> 1\n"
        start = logg.info("2")
        assert t.counter == 2
        assert capsys.readouterr().err == "2\n"
        logg.hint("3")
        assert t.counter == 3
        assert capsys.readouterr().err == "--> 3\n"
        logg.info("4", time=start)
        assert t.counter == 4
        assert capsys.readouterr().err == "4 (0:00:02)\n"
        logg.info("5 {time_passed}", time=start)
        assert t.counter == 5
        assert capsys.readouterr().err == "5 0:00:03\n"
