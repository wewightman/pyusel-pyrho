from setuptools import Extension, setup

# Compile the C extension and put it in the "rho" folder
rho = Extension(
    name="pyrho.rho.__rho__",
    include_dirs=["pyrho/rho"],
    depends=["pyrho/rho/rho.h"],
    sources=["pyrho/rho/rho.c"]
)

# run setup tools
setup(
    name='pyusel-pyrho',
    description="C-Backed ultrasound coherence toolkit",
    author_email="wew12@duke.edu",
    packages=['pyrho', 'pyrho.rho', 'pyrho.rho.__rho__'],
    package_dir={
        'pyrho':'pyrho', 
        'pyrho.rho':'pyrho/rho',
        'pyrho.rho.__rho__':'pyrho/rho'
    },
    install_requires = [
        "numpy",
        "pyusel-cinpy @ https://github.com/wewightman/pyusel-cinpy/archive/main.tar.gz",
    ],
    license="MIT",
    ext_modules=[rho],
    version="0.0.0"
)