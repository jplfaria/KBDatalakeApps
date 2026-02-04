from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class PathsPangenome:

    # project root (e.g., where you run from)
    root: Path = Path("./pangenome/").resolve()

    # ---- set in __post_init__ ----
    genome_dir: Path = None
    assembly_dir: Path = None

    ref_master_faa_protein_fitness: Path = None
    ref_master_faa_protein_phenotype: Path = None

    # generate by genome prep step
    genome_prep_clade_data: Path = None

    out_pangenome_library: Path = None
    out_pangenome_skani_database: Path = None
    out_master_faa_pangenome_members: Path = None
    out_master_faa: Path = None
    out_mmseqs_dir: Path = None

    def __post_init__(self):
        # local dirs
        object.__setattr__(self, "genome_dir", self.root / "genome")
        object.__setattr__(self, "assembly_dir", self.root / "assembly")

        object.__setattr__(self,
                           "out_master_faa_pangenome_members",
                           self.root / "master_faa_pangenome_members.faa")
        object.__setattr__(self, "out_master_faa", self.root / "master_faa.faa")
        object.__setattr__(self, "out_mmseqs_dir", self.root / "mmseqs2")

        object.__setattr__(
            self, "out_pangenome_library",
            self.root / "library_members.txt")
        object.__setattr__(
            self, "out_pangenome_skani_database",
            self.root / "skani_fast_pangenome_members")

        reference_root = Path("/data/reference_data")
        object.__setattr__(self, "ref_master_faa_protein_fitness", reference_root / "proteins_fitness_genomes.faa")
        object.__setattr__(self, "ref_master_faa_protein_phenotype", reference_root / "proteins_phenotype_genomes.faa")

        object.__setattr__(self, "genome_prep_clade_data", self.root / "../user_to_clade.json")

    def ensure(self) -> "PathsPangenome":
        for p in [self.ref_master_faa_protein_fitness,
                  self.ref_master_faa_protein_phenotype,
                  self.genome_prep_clade_data]:
            if not p.exists():
                raise FileNotFoundError(f"Missing reference data: {p}")
        # create required local dirs
        for p in [
            self.genome_dir,
            self.assembly_dir,
        ]:
            p.mkdir(parents=True, exist_ok=True)
        return self
