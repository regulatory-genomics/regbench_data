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

    def fetch(self, path: Path | None = None) -> Path:
        """ Fetches the OSF object and returns the local path to the downloaded file.

        Parameters
        ----------
        path : Path | None, optional
            The directory where the file should be downloaded. If None, the default
            download directory is used.
        """
        return pooch.retrieve(
            self.url,
            known_hash=self.hash,
            fname=self.name,
            path=path,
            progressbar=True,
        )