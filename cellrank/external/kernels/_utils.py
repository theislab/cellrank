from typing import Literal, Mapping, Sequence


class Markers:
    """Helper class for :class:`cellrank.external.kernels.WOT` to get proliferation/apoptosis genes."""

    _proliferation_markers: Mapping[str, Sequence[str]] = {
        "human": [],
        "mouse": [],
        "rat": [],
    }
    _apoptosis_markers: Mapping[str, Sequence[str]] = {
        "human": [],
        "mouse": [],
        "rat": [],
    }

    def proliferation_markers(
        self, organism: Literal["human", "mouse", "rat"]
    ) -> Sequence[str]:
        """Get proliferation markers for ``organism``."""
        try:
            return self._proliferation_markers[organism]
        except KeyError:
            raise NotImplementedError(
                f"Proliferation markers for `{organism}` organism are not yet implemented."
            ) from None

    def apoptosis_markers(
        self, organism: Literal["human", "mouse", "rat"]
    ) -> Sequence[str]:
        """Get apoptosis markers for ``organism``."""
        try:
            return self._apoptosis_markers[organism]
        except KeyError:
            raise NotImplementedError(
                f"Apoptosis markers for `{organism}` organism are not yet implemented."
            ) from None
