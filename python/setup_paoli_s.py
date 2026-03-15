from setuptools import setup, Extension
import numpy as np

# 定义扩展模块
paoli_s_ext = Extension(
    name='paoli_s',  # 模块名，对应import paoli_s
    sources=['paoli_s.c'],  # 源文件
    include_dirs=[np.get_include()],  # 包含numpy头文件
    define_macros=[
        ('NPY_NO_DEPRECATED_API', 'NPY_1_7_API_VERSION'),  # 禁用numpy废弃API
        ('_CRT_SECURE_NO_WARNINGS', '1'),  # 禁用MSVC安全警告
        ('WIN64', '1')  # 标记64位Windows
    ],
    extra_compile_args=[
        '/O2',          # 最高优化级别
        '/std:c11',     # C11标准
        '/W2',          # 中等警告级别
        '/utf-8',       # UTF-8编码
        '/arch:AVX2'    # 启用AVX2指令集加速
    ],
    extra_link_args=[
        '/STACK:10000000'  # 增大栈空间，避免内存溢出
    ]
)

# 执行编译
setup(
    name='paoli_s',
    version='1.0',
    description='对称矩阵泡利分解C扩展模块',
    ext_modules=[paoli_s_ext]
)