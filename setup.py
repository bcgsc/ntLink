import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="ntLink",
    version="1.2.0",
    author="Lauren Coombe",
    author_email="lcoombe@bcgsc.ca",
    description="Genome assembly scaffolder using long reads",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/bcgsc/ntLink",
    license="GPLv3",
    python_requires=">=3",
    install_requires=["python-igraph", "numpy"],
    scripts = ["bin/ntlink_pair.py", "bin/read_fasta.py", "bin/ntlink_utils.py"],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
    ],
)
