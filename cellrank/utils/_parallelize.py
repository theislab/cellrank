# -*- coding: utf-8 -*-
"""Module used to parallelize model fitting."""

from typing import Any, Union, Callable, Optional, Sequence
from threading import Thread
from multiprocessing import Manager

import numpy as np
import joblib as jl

from cellrank.utils._utils import _get_n_cores

_msg_shown = False


def parallelize(
    callback: Callable[[Any], Any],
    collection: Sequence[Any],
    n_jobs: Optional[int] = None,
    n_split: Optional[int] = None,
    unit: str = "",
    as_array: bool = True,
    use_ixs: bool = False,
    backend: str = "multiprocessing",
    extractor: Optional[Callable[[Any], Any]] = None,
    show_progress_bar: bool = True,
) -> Union[np.ndarray, Any]:
    """
    Parallelize function call over a collection of elements.

    Params
    ------
    callback
        Function to parallelize.
    collection
        Sequence of items which to chunkify.
    n_jobs
        Number of parallel jobs.
    n_split
        Split :paramref:`collection` into :paramref:`n_split` chunks.
        If `None`, split into :paramref:`n_jobs` chunks.
    unit
        Unit of the progress bar.
    as_array
        Whether to convert the results not :class:`numpy.ndarray`.
    use_ixs
        Whether to pass indices to the callback.
    backend
        Which backend to use for multiprocessing.
        See :class:`joblib.Parallel` for valid options.
    extractor
        Function to apply to the result after all jobs have finished.
    show_progress_bar
        Whether to show a progress bar.

    Returns
    -------
    :class:`numpy.ndarray`:w
        Result depending on :paramref:`extractor` and :paramref:`as_array`.
    """

    if show_progress_bar:
        try:
            try:
                from tqdm.notebook import tqdm
            except ImportError:
                from tqdm import tqdm_notebook as tqdm
            import ipywidgets  # noqa
        except ImportError:
            global _msg_shown
            tqdm = None

            if not _msg_shown:
                print(
                    "Unable to create progress bar. Consider installing `tqdm` as `pip install tqdm` "
                    "and `ipywidgets` as `pip install ipywidgets`.\n"
                    "Optionally, you can disable the progress bar using `show_progress_bar=False`."
                )
                _msg_shown = True
    else:
        tqdm = None

    def update(pbar, queue, n_total):
        n_finished = 0
        while n_finished < n_total:
            if queue.get() is None:
                n_finished += 1
            elif pbar is not None:
                pbar.update()

        if pbar is not None:
            pbar.close()

    def wrapper(*args, **kwargs):
        pbar = None if tqdm is None else tqdm(total=len(collection), unit=unit)

        queue = Manager().Queue()
        thread = Thread(target=update, args=(pbar, queue, len(collections)))
        thread.start()

        res = jl.Parallel(n_jobs=n_jobs, backend=backend)(
            jl.delayed(callback)(
                *((i, cs) if use_ixs else (cs,)), *args, **kwargs, queue=queue
            )
            for i, cs in enumerate(collections)
        )

        res = np.array(res) if as_array else res
        thread.join()

        return res if extractor is None else extractor(res)

    n_jobs = _get_n_cores(n_jobs, len(collection))
    if n_split is None:
        n_split = n_jobs
    collections = list(filter(len, np.array_split(collection, n_split)))

    return wrapper
