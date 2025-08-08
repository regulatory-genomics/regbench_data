from pathlib import Path
from .utils import OsfObject

def get_CAGE(name: str) -> list[Path]:
    registy = {
        "K562" [
            OsfObject(
                id="54278",
                name="K562_+.w5z",
                hash="sha256:301e52eb63aff6ec442d7d81fcd13f2a3ee19def735ec5e5662b3a3c396cf000",
            ),
            OsfObject(
                id="mh6ka",
                name="K562_-.w5z",
                hash="sha256:00185eeb469bae9acde78113013e547cfd7951aee5d635fa87c70e35c4116356",
            ),
        ],
    }

    if name not in registy:
        raise ValueError(f"Unknown CAGE dataset: {name}")
    return registy[name].fetch()