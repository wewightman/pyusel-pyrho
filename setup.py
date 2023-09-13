from setuptools import Extension, setup

# Compile the C extension and put it in the "rho" folder
rho = Extension(
    name="pyrho.rho.__rho__",
    include_dirs=["pyrho/rho"],
    depends=["pyrho/rho/rho.h", "pyrho/rho/cubic.h"],
    sources=["pyrho/rho/rho.c", "pyrho/rho/cubic.c"]
)

# Compile the C extension and put it in the "rho" folder
trig = Extension(
    name="pyrho.trig.__trig__",
    include_dirs=["pyrho/trig"],
    depends=["pyrho/trig/trig.h"],
    sources=["pyrho/trig/trig.c"]
)

# run setup tools
setup(
    name='pyusel-pyrho',
    description="C-Backed ultrasound coherence toolkit",
    author_email="wew12@duke.edu",
    packages=[
        'pyrho', 
        'pyrho.rho', 
        'pyrho.rho.__rho__', 
        'pyrho.trig', 
        'pyrho.trig.__trig__', 
        'pyrho.processors'
    ],
    package_dir={
        'pyrho':'pyrho', 
        'pyrho.rho':'pyrho/rho',
        'pyrho.rho.__rho__':'pyrho/rho',
        'pyrho.trig':'pyrho/trig', 
        'pyrho.trig.__trig__':'pyrho/trig',
        'pyrho.processors':'pyrho/processors',
    },
    install_requires = [
        "numpy",
        "pyusel-cinpy @ https://github.com/wewightman/pyusel-cinpy/archive/main.tar.gz",
        "pyusel-interp @ https://github.com/wewightman/pyusel-interp/archive/main.tar.gz",
        "pyusel-types @ https://github.com/wewightman/pyusel-types/archive/main.tar.gz",
    ],
    license="MIT",
    ext_modules=[rho, trig],
    version="0.0.0"
)