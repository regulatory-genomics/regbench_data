from pathlib import Path

from regbench_data import POOCH

CAGE_DATA = {
    'K562': {
        'plus': 'CAGE_K562_+.w5z',
        'minus': 'CAGE_K562_-.w5z',
    }
}

def list_cage() -> list[str]:
    """Lists all available datasets."""
    return list(CAGE_DATA.keys())

def retrieve_cage(
    id: str | list[str],
) -> dict[str, tuple[Path, Path]]:
    """Retrieves all datasets.

    Parameters
    ----------
    id : str | list[str]
        The ID or list of IDs of the datasets to retrieve.

    Returns
    -------
    dict[str, Dataset]
        A dictionary mapping dataset IDs to Dataset objects.
    """

    if isinstance(id, str):
        id = [id]

    datasets = {}
    for dataset_id in id:
        if dataset_id not in CAGE_DATA:
            raise ValueError(f"Dataset ID {dataset_id} not found. Available datasets: {list()}")
        dataset = CAGE_DATA[dataset_id]
        datasets[dataset_id] = (
            Path(POOCH.fetch(dataset['plus'], progressbar=True)),
            Path(POOCH.fetch(dataset['minus'], progressbar=True))
        )
    return datasets

RNA_DATA = {
    'adipose_Subcutaneous': {
        'plus': 'total_RNA_seq_subcutaneous_adipose_tissue_+.w5z',
        'minus': 'total_RNA_seq_subcutaneous_adipose_tissue_-.w5z',
    }
}

def list_rna() -> list[str]:
    """Lists all available datasets."""
    return list(RNA_DATA.keys())

def retrieve_rna(
    id: str | list[str],
) -> dict[str, tuple[Path, Path]]:
    """Retrieves all datasets.

    Parameters
    ----------
    id : str | list[str]
        The ID or list of IDs of the datasets to retrieve.

    Returns
    -------
    dict[str, Dataset]
        A dictionary mapping dataset IDs to Dataset objects.
    """

    if isinstance(id, str):
        id = [id]

    datasets = {}
    for dataset_id in id:
        if dataset_id not in RNA_DATA:
            raise ValueError(f"Dataset ID {dataset_id} not found. Available datasets: {list()}")
        dataset = RNA_DATA[dataset_id]
        datasets[dataset_id] = (
            Path(POOCH.fetch(dataset['plus'], progressbar=True)),
            Path(POOCH.fetch(dataset['minus'], progressbar=True))
        )
    return datasets