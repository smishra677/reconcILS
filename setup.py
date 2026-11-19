from pathlib import Path
from setuptools import find_packages, setup

ROOT = Path(__file__).parent
setup(
    name="reconcils",
    version="1.1.0",
    description="Gene tree/species tree reconciliation with deterministic species-copy labels",
    long_description=(ROOT / "README.md").read_text(encoding="utf-8"),
    long_description_content_type="text/markdown",
    python_requires=">=3.9",
    packages=find_packages(include=["reconcILS", "reconcILS.*"]),
    install_requires=["pandas", "matplotlib", "igraph", "biopython", "ete3", "rich"],
    entry_points={"console_scripts": ["reconcils=reconcILS.cli:main"]},
    include_package_data=True,
    classifiers=["Intended Audience :: Science/Research", "Programming Language :: Python :: 3", "Topic :: Scientific/Engineering :: Bio-Informatics"],
)
