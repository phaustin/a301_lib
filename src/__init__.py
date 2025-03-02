from importlib.metadata import version, PackageNotFoundError

try:
    print(f"{in inint {version('a301_lib')=}")
    __version__ = version("a301_lib")
except PackageNotFoundError:
    # package is not installed
    pass
