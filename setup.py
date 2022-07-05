import os
from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()
    
BUILD_ID = os.environ.get("BUILD_BUILDID", "0")

setup(
    name="Twrap",
    version="0.1" + "." + BUILD_ID,
    # Author details
    author="Sam Clark",
    author_email="sam.clark@york.ac.uk",
    description = "Wrapper package for Tfold for the folding of RNA/DNA sequences \
    and analysis of the resulting folds.",
    keywords = ["RNA", "DNA", "Folding", "Stemloops", "m&b"],
    packages=find_packages("Twrap"),
    package_dir={"": "Twrap"},
    long_description = long_description,
    url="https://github.com/PsamClark/Twrap",
    license='MIT',
    setup_requires=[
  "datetime",
  "multiprocessing",
  "os",
  "re",
  "time",
  "shutil",
  "subprocess",
  "sys",
  "pandas",
  "glob",
  "numpy",
  "Bio",
  "pathlib"
],
    tests_require=["pytest", "pytest-nunit", "pytest-cov"]
)