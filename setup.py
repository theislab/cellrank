import os
from pathlib import Path

from setuptools import setup, find_packages

try:
    from cellrank import __email__, __author__, __version__, __maintainer__
except ImportError:
    __author__ = "Marius Lange, Michal Klein"
    __maintainer__ = "Marius Lange, Michal Klein"
    __email__ = "info@cellrank.org"
    __version__ = "1.2.0"


if __name__ == "__main__":
    setup(
        name="cellrank",
        use_scm_version=True,
        setup_requires=["setuptools_scm"],
        author=__author__,
        author_email=__email__,
        email=__email__,
        maintainer=__maintainer__,
        maintainer_email=__email__,
        version=__version__,
        description=Path("README.rst").read_text("utf-8").split("\n")[2],
        long_description=Path("README.rst").read_text("utf-8"),
        url="https://github.com/theislab/cellrank",
        project_urls={
            "Documentation": "https://cellrank.readthedocs.io/en/latest",
            "Source Code": "https://github.com/theislab/cellrank",
        },
        download_url="https://github.com/theislab/cellrank",
        license="BSD",
        install_requires=list(
            map(
                str.strip,
                open(os.path.abspath("requirements.txt")).read().splitlines(),
            )
        ),
        extras_require=dict(
            krylov=["pygpcca[slepc]"],
            test=[
                "pytest>=6.1.1",
                "pytest-mock>=3.1.0",
                "pytest-xdist>=2.1.0",
                "Pillow",
                "filelock",
                "mock>=4.0.2",
                "python-igraph",
                "louvain==0.6.1",
                "leidenalg==0.8.1",
                "bezier",  # curved edges for `cellrank.pl.graph`
            ],
            docs=[
                r
                for r in map(
                    str.strip,
                    open(os.path.abspath("docs/requirements.txt")).read().splitlines(),
                )
                if "requirements.txt" not in r
            ],
            dev=["pre-commit>=2.9.0"],
        ),
        zip_safe=False,
        packages=find_packages(),
        python_requires=">=3.6",
        platforms=["Linux", "MacOs", "Windows"],
        keywords=sorted(
            [
                "single-cell",
                "bio-informatics",
                "single-cell",
                "RNA velocity",
                "Markov chain",
                "GPCCA",
            ]
        ),
        classifiers=[
            "Development Status :: 5 - Production/Stable",
            "Intended Audience :: Developers",
            "Intended Audience :: Science/Research",
            "Natural Language :: English",
            "Framework :: Jupyter",
            "Operating System :: MacOS :: MacOS X",
            "Operating System :: Microsoft :: Windows",
            "Operating System :: POSIX :: Linux",
            "Typing :: Typed",
            "Programming Language :: Python :: 3",
            "Programming Language :: Python :: 3.6",
            "Programming Language :: Python :: 3.7",
            "Programming Language :: Python :: 3.8",
            "Programming Language :: Python :: 3.9",
            "Topic :: Scientific/Engineering :: Bio-Informatics",
            "Topic :: Scientific/Engineering :: Visualization",
        ],
    )
