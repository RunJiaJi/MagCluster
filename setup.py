from setuptools import setup

with open("README.md", "r") as fh:
    long_description = fh.read()
setup(
    name='magcluster',
    version='0.0.2',
    description="Magnetosome gene cluster annotation, screening and mapping tool",
    url="https://github.com/RunJiaJi/magcluster",
    author='Runjia Ji',
    author_email='jirunjia@gmail.com',
    py_modules=["magcluster"],
    package_dir={'': 'src'},
    classifiers=[
        "Programming Language :: Python :: 3.0",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "License :: OSI Approved :: GNU General Public License v2 or later (GPLv2+)",
        "Operating System :: OS Independent",
    ],
    long_description=long_description,
    long_description_content_type="text/markdown",
    entry_points={
        'console_scripts':[
            'magcluster = magcluster:main'
        ]
    }
)