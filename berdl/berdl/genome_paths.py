from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class GenomePaths:
    # project root (e.g., where you run from)
    root: Path = Path(".").resolve()
    reference_root: Path = Path("/data/reference_data")

    # ---- set in __post_init__ ----
    genome_dir: Path = None
    assembly_dir: Path = None
    pangenome_dir: Path = None
    ani_dir: Path = None
    skani_library_dir: Path = None

    library_user_genomes: Path = None
    json_user_to_clade = None
    ani_kepangenomes_out: Path = None
    ani_fitness_out: Path = None
    ani_phenotypes_out: Path = None

    skani_fast_kepangenome_db: Path = None
    skani_fast_fitness_db: Path = None
    skani_fast_phenotypes_db: Path = None

    def __post_init__(self):
        # local dirs
        object.__setattr__(self, "genome_dir", self.root / "genome")
        object.__setattr__(self, "assembly_dir", self.root / "assembly")
        object.__setattr__(self, "ani_dir", self.root / "ani")
        object.__setattr__(self, "skani_library_dir", self.root / "library")
        object.__setattr__(self, "pangenome_dir", self.root / "pangenome")

        # output files
        object.__setattr__(
            self, "library_user_genomes",
            self.skani_library_dir / "user_genome.txt"
        )
        object.__setattr__(
            self, "ani_kepangenomes_out",
            self.ani_dir / "kepangenomes_fast.out"
        )
        object.__setattr__(
            self, "ani_fitness_out",
            self.ani_dir / "fitness_fast.out"
        )
        object.__setattr__(
            self, "ani_phenotypes_out",
            self.ani_dir / "phenotypes_fast.out"
        )
        object.__setattr__(
            self, "json_user_to_clade",
            self.pangenome_dir / "user_to_clade.json"
        )

        # reference data (derived from configurable root)
        skani_root = self.reference_root / "skani_library"
        object.__setattr__(
            self, "skani_fast_kepangenome_db",
            skani_root / "fast_kepangenome"
        )
        object.__setattr__(
            self, "skani_fast_fitness_db",
            skani_root / "fast_fitness"
        )
        object.__setattr__(
            self, "skani_fast_phenotypes_db",
            skani_root / "fast_phenotypes"
        )

    def ensure(self) -> "GenomePaths":
        for p in [self.reference_root, self.skani_fast_kepangenome_db, self.skani_fast_fitness_db, self.skani_fast_phenotypes_db]:
            if not p.exists():
                raise FileNotFoundError(f"Missing reference data: {p}")
        # create required local dirs
        for p in [
            self.genome_dir,
            self.assembly_dir,
            self.ani_dir,
            self.skani_library_dir,
            self.pangenome_dir
        ]:
            p.mkdir(parents=True, exist_ok=True)
        return self
