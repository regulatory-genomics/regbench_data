from pathlib import Path
import pooch
from dataclasses import dataclass

@dataclass
class OsfObject:
    """ Represents an object stored on the Open Science Framework (OSF).

    Attributes
    ----------
    id : str
        The unique identifier for the OSF object.
    name : str
        The name of the file to be downloaded.
    hash : str | None
        The SHA256 hash of the file for integrity checking. If None, no hash check is
        performed.
    """
    id: str
    name: str
    hash: str | None = None

    def __post_init__(self):
        self.url = f"https://osf.io/download/{self.id}"

    def fetch(self, cache_dir: Path | None = None) -> Path:
        """ Fetches the OSF object and returns the local path to the downloaded file.

        Parameters
        ----------
        cache_dir: Path | None, optional
            The location of the cache folder on disk. This is where the file will be saved.
            If None, will save to a folder in the default cache location for your
            operating system (see pooch.os_cache).
        """
        return pooch.retrieve(
            self.url,
            known_hash=self.hash,
            fname=self.name,
            path=cache_dir,
            progressbar=True,
        )