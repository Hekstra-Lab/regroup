from setuptools import setup, find_packages

setup(
    name='regroup',
    version='0.1.0',
    author='Harrison K. Wang, Jack B. Greisman',
    author_email='hwang1@g.harvard.edu',
    packages=find_packages(),
    description='Space group conversion upon directional perturbation',
    install_requires=[
        "pandas",
        "numpy",
        "reciprocalspaceship",
        "gemmi",
    ],
    entry_points={
        'console_scripts': [
            'regroup=regroup.regroup:main',
            'regroup.low_sym=regroup.low_sym:main'
        ]
    }
)
