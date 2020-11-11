# -*- coding: utf-8 -*-
"""Module used to parallelize model fitting."""

from typing import Any, Union, Callable, Optional, Sequence
from threading import Thread
from multiprocessing import Manager

import joblib as jl

import numpy as np
from scipy.sparse import issparse, spmatrix

from cellrank.ul._utils import _get_n_cores


def parallelize(
    callback: Callable[[Any], Any],
    collection: Union[spmatrix, Sequence[Any]],
    n_jobs: Optional[int] = None,
    n_split: Optional[int] = None,
    unit: str = "",
    as_array: bool = True,
    use_ixs: bool = False,
    backend: str = "multiprocessing",
    extractor: Optional[Callable[[Any], Any]] = None,
    show_progress_bar: bool = True,
) -> Any:
    """
    Parallelize function call over a collection of elements.

    Parameters
    ----------
    callback
        Function to parallelize.
    collection
        Sequence of items which to chunkify or an already .
    n_jobs
        Number of parallel jobs.
    n_split
        Split ``collection`` into ``n_split`` chunks. If `None`, split into ``n_jobs`` chunks.
    unit
        Unit of the progress bar.
    as_array
        Whether to convert the results not :class:`numpy.ndarray`.
    use_ixs
        Whether to pass indices to the callback.
    backend
        Which backend to use for multiprocessing. See :class:`joblib.Parallel` for valid options.
    extractor
        Function to apply to the result after all jobs have finished.
    show_progress_bar
        Whether to show a progress bar.

    Returns
    -------
        The result depending on ``callable``, ``extractor`` and ``as_array``.
    """

    if show_progress_bar:
        from tqdm.auto import tqdm
    else:
        tqdm = None

    def update(pbar, queue, n_total):
        n_finished = 0
        while n_finished < n_total:
            try:
                res = queue.get()
            except EOFError as e:
                if not n_finished != n_total:
                    raise RuntimeError(
                        f"Finished only `{n_finished} out of `{n_total}` tasks.`"
                    ) from e
                break
            assert res in (None, (1, None), 1)  # (None, 1) means only 1 job
            if res == (1, None):
                n_finished += 1
                if pbar is not None:
                    pbar.update()
            elif res is None:
                n_finished += 1
            elif pbar is not None:
                pbar.update()

        if pbar is not None:
            pbar.close()

    def wrapper(*args, **kwargs):
        if pass_queue and show_progress_bar:
            pbar = (
                None
                if tqdm is None
                else tqdm(total=col_len, unit=unit, mininterval=0.125)
            )
            queue = Manager().Queue()
            thread = Thread(target=update, args=(pbar, queue, len(collections)))
            thread.start()
        else:
            pbar, queue, thread = None, None, None

        res = jl.Parallel(n_jobs=n_jobs, backend=backend)(
            jl.delayed(callback)(
                *((i, cs) if use_ixs else (cs,)),
                *args,
                **kwargs,
                queue=queue,
            )
            for i, cs in enumerate(collections)
        )

        res = np.array(res) if as_array else res
        if thread is not None:
            thread.join()

        return res if extractor is None else extractor(res)

    col_len = collection.shape[0] if issparse(collection) else len(collection)
    n_jobs = _get_n_cores(n_jobs, col_len)
    if n_split is None:
        n_split = n_jobs

    if issparse(collection):
        if n_split == collection.shape[0]:
            collections = [collection[[ix], :] for ix in range(collection.shape[0])]
        else:
            step = collection.shape[0] // n_split

            ixs = [
                np.arange(i * step, min((i + 1) * step, collection.shape[0]))
                for i in range(n_split)
            ]
            ixs[-1] = np.append(
                ixs[-1], np.arange(ixs[-1][-1] + 1, collection.shape[0])
            )

            collections = [collection[ix, :] for ix in filter(len, ixs)]
    else:
        collections = list(filter(len, np.array_split(collection, n_split)))

    pass_queue = not hasattr(callback, "py_func")  # we'd be inside a numba function

    return wrapper
