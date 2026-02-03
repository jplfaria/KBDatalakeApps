from abc import ABC, abstractmethod
from typing import Iterable, Protocol, Optional, Any, Union


class QueryPangenomeABC(ABC):

    @abstractmethod
    def get_clade_gene_clusters(self, clade_id: str) -> Any:
        """Return all gene clusters for a given clade."""
        raise NotImplementedError

    @abstractmethod
    def get_clade_by_member(self, member_id: str) -> Any:
        """Return the clade for a given genome/member (or clade metadata)."""
        raise NotImplementedError

    @abstractmethod
    def get_cluster_members(self, cluster_id: str) -> Any:
        """Return rows mapping genes to the given cluster."""
        raise NotImplementedError

    @abstractmethod
    def get_clusters_members(self, cluster_ids: Iterable[str]) -> Any:
        """Return rows mapping genes to any of the given clusters."""
        raise NotImplementedError

    @abstractmethod
    def get_clade_members(self, clade_id: str) -> Any:
        """Return genomes/members belonging to the given clade."""
        raise NotImplementedError

    @abstractmethod
    def get_member_representative(self, member_id: str) -> str:
        """Return the GTDB species clade id for the given genome/member."""
        raise NotImplementedError

    @abstractmethod
    def get_member_ani_matrix(self, member_id: str) -> Any:
        """Return ANI rows where member_id is genome1_id or genome2_id."""
        raise NotImplementedError

    @abstractmethod
    def get_clade_metadata(self, clade_id: str) -> Any:
        """Return metadata rows for a given clade."""
        raise NotImplementedError
