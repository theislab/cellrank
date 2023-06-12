from importlib import metadata

from cellrank import datasets, estimators, kernels, logging, models, pl
from cellrank._utils._lineage import Lineage
from cellrank.settings import settings

__all__ = ["datasets", "estimators", "kernels", "logging", "models", "pl", "Lineage", "settings"]

try:
    md = metadata.metadata(__name__)
    __version__ = md.get("version", "")
    __author__ = md.get("Author", "")
    __maintainer__ = md.get("Maintainer-email", "")
except ImportError:
    md = None

del metadata, md
