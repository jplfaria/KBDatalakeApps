import hashlib


def _hash_string(s):
    return hashlib.sha256(s.encode("utf-8")).hexdigest()


class HashSeq(str):

    def __new__(cls, v, strip_ending_star=True):
        # print('validate!!', v)
        if strip_ending_star:
            v = HashSeq.strip_ending_star(v)
        instance = super().__new__(cls, v.upper())
        return instance

    @property
    def hash_value(self):
        h = _hash_string(self)
        return h

    @staticmethod
    def strip_ending_star(s):
        if s.endswith('*'):
            return s[:-1]
        return s


class HashSeqList(list):

    def append(self, o, /):
        if type(o) is str:
            super().append(HashSeq(o))
        elif type(o) is HashSeq:
            super().append(o)
        else:
            raise ValueError('bad type')

    @property
    def hash_value(self):
        h_list = [x.hash_value for x in self]
        hash_seq = "_".join(sorted(h_list))
        return _hash_string(hash_seq)


class ProteinSequence(HashSeq):

    _STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")
    _EXTENDED_ONLY = {"U", "O"}  # selenocysteine, pyrrolysine
    _AMBIGUOUS_AA = set("BJZX*")

    def __new__(cls, sequence, strip_ending_star=True):
        obj = super().__new__(cls, sequence, strip_ending_star=strip_ending_star)
        return obj

    @property
    def sequence(self) -> str:
        return str(self)

    def is_standard(self) -> bool:
        """True if the sequence contains only standard amino acids."""
        return bool(self) and set(self) <= self._STANDARD_AA

    def is_extended(self) -> bool:
        """True if all residues are standard or extended (U, O)."""
        allowed = self._STANDARD_AA | self._EXTENDED_ONLY
        return bool(self) and set(self) <= allowed

    def is_ambiguous(self) -> bool:
        """True if the sequence contains any ambiguous residue codes."""
        return any(res in self._AMBIGUOUS_AA for res in self)

    def is_valid(self) -> bool:
        """True if all residues are recognized AA codes (standard, extended, or ambiguous)."""
        allowed = self._STANDARD_AA | self._EXTENDED_ONLY | self._AMBIGUOUS_AA
        return bool(self) and set(self) <= allowed

    def z_compress(self):
        import zlib
        return zlib.compress(str(self).encode("utf-8"))
