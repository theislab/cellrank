from typing import Mapping, Sequence

from typing_extensions import Literal


class MarkerGenes:
    """Helper class for :class:`cellrank.external.kernels.WOT` to get proliferation/apoptosis genes."""

    _proliferation_markers: Mapping[str, Sequence[str]] = {
        "human": (),
        "mouse": (),
    }
    _apoptosis_markers: Mapping[str, Sequence[str]] = {
        "human": (),
        "mouse": (),
    }

    @classmethod
    def proliferation_markers(
        cls, organism: Literal["human", "mouse"]
    ) -> Sequence[str]:
        """Get proliferation markers for ``organism``."""
        try:
            return cls._proliferation_markers[organism]
        except KeyError:
            raise NotImplementedError(
                f"Proliferation markers for `{organism}` are not yet implemented."
            ) from None

    @classmethod
    def apoptosis_markers(cls, organism: Literal["human", "mouse"]) -> Sequence[str]:
        """Get apoptosis markers for ``organism``."""
        try:
            return cls._apoptosis_markers[organism]
        except KeyError:
            raise NotImplementedError(
                f"Apoptosis markers for `{organism}` are not yet implemented."
            ) from None
