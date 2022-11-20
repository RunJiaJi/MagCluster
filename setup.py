from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()
setup(
    name='magcluster',
    version='0.2.11',
    description="Magnetosome Gene Clusters Analyzer",
    url="https://github.com/RunJiaJi/magcluster",
    author='Runjia Ji',
    author_email='jirunjia@gmail.com',
    # py_modules=["magcluster", 'args', 'capture_args', 'maga', 'magm', 'magsc', 'main'],
    # package_dir={'': 'src'},
    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3.0",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "License :: OSI Approved :: GNU General Public License v2 or later (GPLv2+)",
        "Operating System :: POSIX :: Linux",
    ],
    long_description=long_description,
    long_description_content_type="text/markdown",
    install_requires=[
        "clinker",
        'pandas',
    ],
    entry_points={
        'console_scripts':[
            'magcluster = magcluster:main',
        ]
    },
    package_data={
        'magcluster':['data/*.fasta'],
    }
)