from cellrank.tl.kernels import ExperimentalTimeKernel


class WOTKernel(ExperimentalTimeKernel):
    """Kernel based on Waddington optimal transport [Schiebinger19]_."""

    def __init__(self, *args, **kwargs):
        # TODO: test for wot import
        super().__init__(*args, **kwargs)

    def compute_transition_matrix(self, *args, **kwargs) -> "WOTKernel":
        """TODO."""
        return self
