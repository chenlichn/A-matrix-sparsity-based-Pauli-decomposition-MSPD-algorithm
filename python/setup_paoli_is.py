from setuptools import setup, Extension
import numpy as np

paoli_is_ext = Extension(
    name='paoli_is',
    sources=['paoli_is.c'],
    include_dirs=[np.get_include()],
    define_macros=[
        ('NPY_NO_DEPRECATED_API', 'NPY_1_7_API_VERSION'),
        ('_CRT_SECURE_NO_WARNINGS', '1'),
        ('WIN64', '1')
    ],
    extra_compile_args=[
        '/O2', '/std:c11', '/W2', '/utf-8'
    ],
    extra_link_args=['/STACK:10000000']
)

setup(
    name='paoli_is',
    version='1.0',
    ext_modules=[paoli_is_ext]
)