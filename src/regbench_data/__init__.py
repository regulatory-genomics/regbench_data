import importlib.metadata
from importlib.resources import files
import pooch

__version__ = importlib.metadata.version("regbench_data")

POOCH = pooch.create(
    path=pooch.os_cache("regbench_data"),
    registry=None,
    base_url='',
)
POOCH.load_registry(files('regbench_data.data').joinpath('registry.txt'))