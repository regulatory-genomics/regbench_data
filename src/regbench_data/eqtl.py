from pathlib import Path
import polars as pl
import yaml
from dataclasses import dataclass

from regbench_data import POOCH

EQTL_DATA = {
    'adipose_Subcutaneous': 'GTEx_v10_SuSiE_eQTL_Adipose_Subcutaneous.v10.eQTLs.SuSiE_summary.parquet',
}

def list_eqtl() -> list[str]:
    """Lists all available datasets."""
    return list(EQTL_DATA.keys())

def retrieve_eqtl(
    id: str | list[str],
) -> dict[str, pl.DataFrame]:
    """Retrieves all datasets.

    Parameters
    ----------
    ids : str | list[str]
        The ID or list of IDs of the datasets to retrieve.

    Returns
    -------
    dict[str, pl.DataFrame]
        A dictionary mapping dataset IDs to Polars DataFrames. The DataFrames contain the following columns:
            - gene_id: GENCODE/Ensembl gene ID or RNAcentral URS ID
            - phenotype_id: Phenotype ID, e.g., intron coordinates and cluster combined with gene ID for sQTLs
            - gene_name: GENCODE gene name
            - biotype: gene or transcript classification (protein coding, lncRNA, etc.)
            - variant_id: variant ID in the format {chr}_{pos_first_ref_base}_{ref_seq}_{alt_seq}_b38
            - pip: posterior inclusion probability (PIP)
            - af: allele frequency of the ALT allele (in-sample)
            - cs_id: credible set ID (number)
            - cs_size: credible set size
            - afc: allelic fold change (aFC) of the lead variant (highest PIP) in the credible set
            - afc_se: standard error of the aFC of the lead variant (highest PIP) in the credible set
    """

    if isinstance(id, str):
        id = [id]

    datasets = {}
    for dataset_id in id:
        if dataset_id not in EQTL_DATA:
            raise ValueError(f"Dataset ID {dataset_id} not found. Available datasets: {list_eqtl()}")
        data_file = POOCH.fetch(EQTL_DATA[dataset_id], progressbar=True)
        datasets[dataset_id] = pl.read_parquet(data_file)
    return datasets

if __name__ == "__main__":
    datasets = list_eqtl()
    print(datasets)
    datasets = retrieve_eqtl(datasets[0])
    print(datasets)