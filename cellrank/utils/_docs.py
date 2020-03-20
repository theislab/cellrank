# -*- coding: utf-8 -*-
from textwrap import dedent


def inject_docs(**kwargs):
    """
    Docstrings should start with "\" in the first line for proper formatting.
    """

    def decorator(obj):
        obj.__doc__ = dedent(obj.__doc__).format(**kwargs)
        return obj

    return decorator
