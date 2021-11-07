import subprocess
from pathlib import Path


class Kmer():
    def __init__(self, query, mer, elements, out="temp") -> None:
        self.query = query
        self.mer = mer
        self.elements = elements
        out = Path(out)
        Path.mkdir(out, parents=True, exist_ok=True)
        self.out = str(out.joinpath("kmer_count.jf"))

    def get_kmer(self):
        subprocess.run(["jellyfish", "count", "-m", self.mer, "-s",
                       self.elements, "-t", "10", "-C", self.query, "-o", self.out])

    def get_histogram(self):
        subprocess.run(["jellyfish", "histo", self.out])

    def get_histogram(self):
        subprocess.run(["jellyfish", "histo", self.out])
