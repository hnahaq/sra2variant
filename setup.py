from setuptools import setup, find_packages

from sra2variant import __version__

setup(
    name="sra2variant",
    version=__version__,
    description="Convert SRA file to variant",
    packages=find_packages(),
    python_requires='>=3.7',
    install_requires=["pyvcf"],
    url="https://github.com/wuaipinglab/ncov_sequencing_variant",
    author="Chengyang Ji",
    author_email="chengyang.ji12@alumni.xjtlu.edu.cn",
    entry_points={
        "console_scripts": [
            "sra2variant-WGS-PE = sra2variant.WGS_PE:main",
            "fastq2variant-WGS-PE = sra2variant.WGS_PE:main_fastq",
            "sra2variant-ARTIC-PE = sra2variant.ARTIC_PE:main",
            "fastq2variant-ARTIC-PE = sra2variant.ARTIC_PE:main_fastq",
        ]
    },
    zip_safe=False,
    license="MIT",
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3 :: Only",
        "Operating System :: Unix",
        "Operating System :: MacOS"
    ]
)
