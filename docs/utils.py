import re

from enchant.tokenize import Filter


class ModnameFilter(Filter):
    """Ignore module names."""

    _pat = re.compile(r"cellrank\.+")

    def _skip(self, word: str) -> bool:
        return self._pat.match(word) is not None
