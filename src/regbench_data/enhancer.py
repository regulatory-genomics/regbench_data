from pathlib import Path
import polars as pl
import yaml
from dataclasses import dataclass

from regbench_data import POOCH

ENHANCER_DATA = {
    'Gasperini2019': {
        'metadata_file': 'Gasperini_2019_Cell_metadata.yaml',
        'data': [
            {'name': 'Gasperini2019.tsv.gz', 'data_file': 'Gasperini_2019_Cell_Gasperini2019.tsv.gz'}
        ]
    },
    'Nasser2021': {
        'metadata_file': 'Nasser_2021_Nature_metadata.yaml',
        'data': [
            {'name': 'Nasser2021.tsv.gz', 'data_file': 'Nasser_2021_Nature_Nasser2021.tsv.gz'}
        ]
    },
    'Schraivogel2020': {
        'metadata_file': 'Schraivogel_2020_NatMethods_metadata.yaml',
        'data': [
            {'name': 'Schraivogel2020.tsv.gz', 'data_file': 'Schraivogel_2020_NatMethods_Schraivogel2020.tsv.gz'}
        ] 
    }
}

def list_enhancer() -> list[str]:
    """Lists all available datasets."""
    return list(ENHANCER_DATA.keys())

def retrieve_enhancer(
    id: str | list[str] | None = None,
    p_value: float | None = None
) -> dict[str, 'Dataset']:
    """Retrieves all datasets.

    Parameters
    ----------
    id : str | list[str] | None, optional
        The ID or list of IDs of the datasets to retrieve. If None, retrieves all datasets.
    p_value : float | None, optional
        If provided, modify the 'label' column in each ScreeningResult based on the adjusted p-value.
        Otherwise, the 'label' column is left unchanged.

    Returns
    -------
    dict[str, Dataset]
        A dictionary mapping dataset IDs to Dataset objects.
    """

    if id is None:
        id = list_enhancer()
    elif isinstance(id, str):
        id = [id]

    datasets = {}
    for dataset_id in id:
        if dataset_id not in ENHANCER_DATA:
            raise ValueError(f"Dataset ID {dataset_id} not found. Available datasets: {list_enhancer()}")
        dataset = ENHANCER_DATA[dataset_id]
        metadata = POOCH.fetch(dataset['metadata_file'])
        for d in dataset['data']:
            data = [(d['name'], POOCH.fetch(d['data_file'], progressbar=True))]
            datasets[dataset_id] = Dataset.load(metadata, data, p_value=p_value)
    return datasets

@dataclass
class ScreeningResult:
    """Represents a screening result from a dataset.

    Attributes
    ----------
    sample_term_id : str
        The term ID of the sample.
    sample_name : str
        The name of the sample.
    assembly : str
        The genome assembly used for the sample.
    result : pl.DataFrame
        The results of the screening, loaded as a Polars DataFrame.
        It contains the following columns:
            - chrom: str
            - chrom_start: int
            - chrom_end: int
            - gene_symbol: str
            - gene_chrom: str
            - gene_TSS: int
            - label: str (categorical, e.g., '0' for non-significant and '1' for significant)
            - effect_size: float | None
            - adjusted_p_value: float | None
    """

    sample_term_id: str
    sample_name: str
    assembly: str
    result: pl.DataFrame

    def __len__(self):
        return len(self.result)

    def label_counts(self) -> dict[str, int]:
        return self.result["label"].value_counts().sort("label")
    
    def distance_to_tss(self) -> pl.Series:
        """Calculates the distance from the enhancer center to the gene TSS."""
        return ((self.result["chrom_start"] + self.result["chrom_end"]) / 2 - self.result["gene_TSS"]).abs()

    def load_data(
        self, csv: Path, p_value: float | None = None
    ) -> "ScreeningResult":
        """Loads the screening results from an OSF object.

        Parameters
        ----------
        object : OsfObject
            The OSF object containing the screening results.
        p_value : float | None, optional
            If provided, modify the 'label' column based on the adjusted p-value.
            Otherwise, the 'label' column is left unchanged.
        """
        df = pl.read_csv(
            csv,
            separator="\t",
            schema_overrides={
                "chrom": pl.String,
                "chrom_start": pl.UInt64,
                "chrom_end": pl.UInt64,
                "gene_symbol": pl.String,
                "gene_chrom": pl.String,
                "gene_TSS": pl.UInt64,
                "label": pl.Int64,
            },
            null_values={"effect_size": "NA", "adjusted_p_value": "NA"},
        )
        if p_value is not None:
            df = df.with_columns(
                pl.when(pl.col("adjusted_p_value") <= p_value)
                .then(1)
                .otherwise(0)
                .cast(pl.Int64)
                .alias("label")
            )
        self.result = df


def concatenate(
    results: list[ScreeningResult] | list[tuple[str, ScreeningResult]],
) -> ScreeningResult:
    """Concatenates multiple ScreeningResult objects into a single one.

    The sample_term_id, sample_name and assembly must be the same for all results.

    Parameters
    ----------
    results
        The list of ScreeningResult objects to concatenate.
        Each object can optionally be associated with a string key. The keys are added
        as a new column 'source' in the resulting DataFrame.

    Returns
    -------
    ScreeningResult
        A new ScreeningResult object containing the concatenated results.
    """
    assert len(results) > 0, "No results to concatenate"

    if isinstance(results[0], tuple):
        for key, r in results:
            r.result = r.result.with_columns(pl.lit(key).alias("source"))
        results = [r for _, r in results]

    sample_term_id = results[0].sample_term_id
    sample_name = results[0].sample_name
    assembly = results[0].assembly
    for r in results:
        if r.sample_term_id != sample_term_id:
            raise ValueError("All results must have the same sample_term_id")
        if r.sample_name != sample_name:
            raise ValueError("All results must have the same sample_name")
        if r.assembly != assembly:
            raise ValueError("All results must have the same assembly")
    concatenated_df = pl.concat([r.result for r in results], how="vertical")
    return ScreeningResult(sample_term_id, sample_name, assembly, concatenated_df)


@dataclass
class Dataset:
    """Represents a dataset containing screening results.

    Attributes
    ----------
    id : str
        The unique identifier for the dataset.
    results : list[ScreeningResult]
        A list of screening results contained in the dataset.
    """

    id: str
    results: list[ScreeningResult]

    @classmethod
    def load(
        cls,
        metadata: Path,
        data: list[tuple[str, Path]],
        p_value: float | None = None,
    ) -> "Dataset":
        with open(metadata, "r") as f:
            metadata = yaml.safe_load(f)

        files = {x["file"]: x for x in metadata["data"]}
        results = []
        for filename, d in data:
            info = files.get(filename)
            if info is None:
                raise ValueError(f"File {filename} not found in metadata")
            res = ScreeningResult(
                info["sample_term_id"],
                info["sample_name"],
                info["assembly"],
                None,
            )
            res.load_data(d, p_value=p_value)
            results.append(res)

        return Dataset(metadata["id"], results)

if __name__ == "__main__":
    datasets = retrieve_enhancer()
    for id, dataset in datasets.items():
        print(f"Dataset ID: {id}, Number of Results: {len(dataset.results)}")
        for result in dataset.results:
            print(f"  Sample: {result.sample_name}, Labels: {result.label_counts()}")
