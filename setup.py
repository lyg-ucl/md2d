from setuptools import setup, find_packages
import pathlib

here = pathlib.Path(__file__).parent.resolve()

long_description = (here / 'README.md').read_text(encoding='utf-8')

setup(
    name='md2d',  # Required
    version='1.2.0',  # Required

    description='A python module for accurate determination of diffusion coefficient from molecular dynamics',
    long_description=long_description,
    long_description_content_type='text/markdown', 

    url='https://github.com/lyg-ucl/md2d',
    license='GNU General Public License v3 (GPLv3)',
    author='Yunguo Li', 
    author_email='liyungo@ustc.edu.cn', 

    classifiers=[ 
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        'Development Status :: 3 - Alpha',
    ],

    keywords='diffusion coefficient, molecular dynamics', 

    py_modules=["md2d","einstein","vis","plot","md2vasp"],
    #package_dir={'': 'src'}, 
    packages=find_packages(),  # Required

    python_requires='>=3.6, <4',

    install_requires=['numpy','scipy','matplotlib','mdanalysis'], 

    entry_points={
        'console_scripts': [
            'md2d=md2d:md2dexe',
            'md2vasp=md2d:md2v',
        ],
    },

)
