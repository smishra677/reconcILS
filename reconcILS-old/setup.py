from setuptools import setup, find_packages
from pathlib import Path

here = Path(__file__).parent
long_description = (here / "README.md").read_text(encoding="utf-8")

setup(
    name='reconcILS',
    version='0.1.0',
    description='A package for gene tree-species tree reconciliation',
    long_description=long_description,
    long_description_content_type='text/markdown',
    author='Sarthak Mishra',
    author_email='samishr@iu.edu',
    packages=find_packages(),
    install_requires=[
      'pandas',
      'matplotlib',
      'uuid',
      'igraph',
      
    ],
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.6',
)

