from setuptools import setup, find_packages

setup(
    name='regroup',
    version='0.0.1',
    author='Jack B. Greisman',
    author_email='greisman@g.harvard.edu',
    packages=find_packages(),
    description='Determine new space groups for analyzing pump-probe crystallography experiments',
    install_requires=[
        "pandas",
        "numpy",
    ],
    entry_points={
        'console_scripts': [
            'regroup=regroup.regroup:main',
        ]
    }
)
