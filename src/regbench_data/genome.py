from pathlib import Path
import pooch
from .utils import OsfObject

def fetch_genome_fasta(name: str) -> Path:
    registy = {
        "GRCh38": (
            'gencode_v41_GRCh38.fa.gz',
            'https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_41/GRCh38.primary_assembly.genome.fa.gz',
            'sha256:4fac949d7021cbe11117ddab8ec1960004df423d672446cadfbc8cca8007e228',
        ),
    }
    if name not in registy:
        raise ValueError(f"Unknown genome: {name}")
    
    file_name, url, hash_value = registy[name]
    return pooch.retrieve(
        url,
        known_hash=hash_value,
        fname=file_name,
        path=None,
        progressbar=True,
        processor=pooch.Decompress(method = "gzip"),
    )

def fetch_genome_annotation(name: str) -> Path:
    registy = {
        "GRCh38": (
            'gencode_v41_GRCh38.gff3.gz',
            "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_41/gencode.v41.basic.annotation.gff3.gz",
            'sha256:b82a655bdb736ca0e463a8f5d00242bedf10fa88ce9d651a017f135c7c4e9285',
        ),
    }
    if name not in registy:
        raise ValueError(f"Unknown genome: {name}")
    
    file_name, url, hash_value = registy[name]
    return pooch.retrieve(
        url,
        known_hash=hash_value,
        fname=file_name,
        path=None,
        progressbar=True,
    )