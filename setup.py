import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="ntLink",
    version="1.3.9",
    author="Lauren Coombe",
    author_email="lcoombe@bcgsc.ca",
    description="Genome assembly scaffolder using long reads",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/bcgsc/ntLink",
    license="GPLv3",
    python_requires=">=3",
    install_requires=["python-igraph", "numpy", "btllib"],
    scripts = ["bin/ntlink_pair.py", "bin/read_fasta.py", "bin/ntlink_utils.py",
               "bin/ntlink_filter_sequences.py", "bin/ntlink_liftover_mappings.py",
               "bin/ntlink_overlap_sequences.py", "bin/ntlink_patch_gaps.py",
               "bin/ntlink_stitch_paths.py", "bin/ntjoin_utils.py",
               "bin/ntlink_paf_output.py", "bin/path_node.py",
               "ntLink", "ntLink_rounds"],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
    ],
)
