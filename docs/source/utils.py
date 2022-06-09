import re
from enchant.tokenize import Filter


class ModnameFilter(Filter):
    """Ignore module names."""

    _pat = re.compile(r"cellrank\.+")

    def _skip(self, word: str) -> bool:
        return self._pat.match(word) is not None


class SignatureFilter(Filter):
    """Ignore function signature artifacts."""

    def _skip(self, word: str) -> bool:
        # TODO(michalk8): find a better way
        return word in (
            "img[",
            "imgs[",
            "img",
            "img_key",
            "func[",
            "func",
            "combine_attrs",
            "**kwargs",
            "n_iter",
        )
