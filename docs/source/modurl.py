"""Based on: https://github.com/theislab/scvelo/blob/fb65bab8eb469e2b7f2793d9aca0cb9062d7e58a/docs/source/conf.py"""
from typing import Any, Tuple, Optional

import sys
from inspect import getsourcelines
from pathlib import Path, PurePosixPath


def get_obj_module(qualname: str) -> Tuple[Any, Any]:
    """Get a module/class/attribute and its original module by qualname"""
    modname = qualname
    classname, attrname = None, None
    while modname not in sys.modules:
        attrname = classname
        modname, classname = modname.rsplit(".", 1)

    # retrieve object and find original module name
    if classname:
        cls = getattr(sys.modules[modname], classname)
        modname = cls.__module__
        obj = getattr(cls, attrname) if attrname else cls
    else:
        obj = None

    return obj, sys.modules[modname]


def get_linenos(obj: Any) -> Tuple[Optional[int], Optional[int]]:
    """Get an object’s line numbers"""
    try:
        lines, start = getsourcelines(obj)
    except TypeError:
        return None, None
    else:
        return start, start + len(lines) - 1


def modurl(qualname: str) -> str:
    """Get the full GitHub URL for some object’s qualname."""
    project_dir = Path(__file__).parent.parent.parent
    obj, module = get_obj_module(qualname)
    github_url = "https://github.com/theislab/cellrank/tree/master"

    try:
        path = PurePosixPath(Path(module.__file__).resolve().relative_to(project_dir))
    except ValueError:
        return github_url

    start, end = get_linenos(obj)
    fragment = f"#L{start}-L{end}" if start and end else ""

    return f"{github_url}/{path}{fragment}"
