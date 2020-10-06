from setuptools import setup


setup(
    name='cldfbench_myphoible',
    py_modules=['cldfbench_myphoible'],
    include_package_data=True,
    zip_safe=False,
    entry_points={
        'cldfbench.dataset': [
            'myphoible=cldfbench_myphoible:Dataset',
        ]
    },
    install_requires=[
        'cldfbench',
        'pyglottolog',
        'pyclts',
    ],
    extras_require={
        'test': [
            'pytest-cldf',
        ],
    },
)
