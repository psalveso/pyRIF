from setuptools import setup


__version__ ='1.0.2'


with open("README.md", "r") as readme_file:
    readme = readme_file.read()


setup(
    name = 'pyRIF',
    version = __version__,
    author = 'Patrick Salveson',
    author_email = 'patricksalveson@gmail.com',
    description = 'Using Rotamer Interaction Fields from RIFGen/Dock in python',
    packages = ['pyRIF'],
    package_dir={'pyRIF' : 'pyRIF'},
    install_requires = ['numpy', 'xbin', 'getpy', 'pynerf', 'h5py'],
    zip_safe = False,
    long_description=readme,
    long_description_content_type='text/markdown',
    url='https://github.com/psalveso/pyRIF',
    classifiers=[
        'Programming Language :: Python :: 3',
        'Topic :: Scientific/Engineering',
        'Topic :: Software Development',
        'License :: OSI Approved :: MIT License',
        'Operating System :: Unix'
    ],
)
