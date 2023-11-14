from filemanager import FileManager
from numpy import ndarray
import numpy as np

class ChispitApp:
    def __init__(self) -> None:
        self.processed_files:FileManager
        self.RSCU_per_gene:list[list[dict[str, float]]]
        self.RSCU_per_aa:list[list[dict[str, float]]]
        self.sequences_lenght:list[list[int]]
        
    def set_processed_files(self, processed_files:FileManager) -> None:
        self.processed_files = processed_files
        
    def get_processed_files(self) -> FileManager:
        return self.processed_files
    
    def set_RSCU_per_gene(self, rscu_per_gene:list[list[dict[str, float]]]) -> None:
        self.RSCU_per_gene = rscu_per_gene
        
    def get_RSCU_per_gene(self) -> list[list[dict[str, float]]]:
        return self.RSCU_per_gene
    
    def set_RSCU_per_aa(self, RSCU_per_aa:list[list[dict[str, float]]]):
        self.RSCU_per_aa = RSCU_per_aa
    
    def get_RSCU_per_aa(self) -> list[list[dict[str, float]]]:
        return self.RSCU_per_aa
    
    def set_sequences_lenght(self, sequences_lenght:list[list[int]]) -> None:
        self.sequences_lenght = sequences_lenght
        
    def get_sequences_lenght(self) -> list[list[int]]:
        return self.sequences_lenght
    
    def calculate_RSCU_per_gene(self) -> None:
        files:FileManager = self.get_processed_files()
        RSCU_per_gene:list[list[dict[str, float]]] = [[gene.get_RSCU() for gene in genome] for genome in files.get_converted_genomes()]
        RSCU_per_aa:list[list[dict[str, list[float]]]] = [[gene.get_RSCU_per_aa() for gene in genome] for genome in files.get_converted_genomes()]
        self.set_RSCU_per_gene(RSCU_per_gene)
        self.set_RSCU_per_aa(RSCU_per_aa)
    
    def process_files(self, file_path:str, genetic_code:dict[str, str]) -> None:
        files:FileManager = FileManager()
        files.set_file_path(file_path)
        files.prepare_files_to_read()
        files.get_sequences_from_multiple_files()
        files.convert_to_genemodel()
        [[gene.process_gene(genetic_code) for gene in genome] for genome in files.get_converted_genomes()]
        self.set_sequences_lenght([[gene.get_sequence_lenght()      for gene in genes ] for genes  in files.get_converted_genomes()])
        self.set_processed_files(files)
        self.calculate_RSCU_per_gene()
    