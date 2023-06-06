from cellrank import pl, models, kernels, logging, datasets, external, estimators
from cellrank.settings import settings
from cellrank._utils._lineage import Lineage
from importlib import metadata

try:
    md = metadata.metadata(__name__)
    __version__ = md.get("version", "")
    __author__ = md.get("Author", "")
    __maintainer__ = md.get("Maintainer-email", "")
except ImportError:
    md = None

del metadata, md
