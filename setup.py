from setuptools import setup, find_packages

setup(
    name='cdem3D_viewer',
    version='0.1.0',
    packages=find_packages(),
    include_package_data=False,
    install_requires=[
        'numpy',
        'pandas',
        'matplotlib',
        'pyvista',
    ],
    author='Nestor Cardozo',
    description='A 3D viewer for cdem3d simulation data',
    classifiers=[
        'Programming Language :: Python :: 3',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.7',
)