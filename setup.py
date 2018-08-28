import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="bacpacs",
    version="0.1.0",
    author="Eran Barash",
    author_email="barashe@post.bgu.ac.il",
    description="Bacterial Pathogenicity Classification via Sparse-SVM",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/barashe/bacpacs",
    packages=setuptools.find_packages(),
    classifiers=(
        "Programming Language :: Python :: 2",
        "License :: OSI Approved :: MIT License",
        "Operating System :: POSIX :: Linux",
    ),
    install_requires=[
        "numpy",
        "pandas",
        "scikit-learn",
        "biopython",

    ],
)