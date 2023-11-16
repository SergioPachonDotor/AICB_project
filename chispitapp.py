from filemanager import FileManager
from codegenerator import CodeGenerator
from numpy import ndarray
import numpy as np
from tqdm import tqdm
import sys
from geneticcode import genetic_code

class ChispitApp:
    def __init__(self) -> None:
        self.processed_files:FileManager
        self.RSCU_per_gene:list[list[dict[str, float]]]
        self.RSCU_per_aa:list[list[dict[str, float]]]
        self.sequences_lenght:list[list[int]]
        self.codes_to_sequences_lenght_array:dict[str, list[int]]
        self.generated_codes_per_genome:list[list[str]]
        self.result_matrix:list[list[float]]
        
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
    
    def set_mapped_codes_to_sequences_lenght_array(self, mapped_codes:list[dict[str, list[int]]]) -> None:
        self.codes_to_sequences_lenght_array = mapped_codes
        
    def get_mapped_codes_to_sequences_lenght_array(self) -> list[dict[str, list[int]]]:
        return self.codes_to_sequences_lenght_array
    
    def set_generated_codes_per_genome(self, generated_codes_per_genome:list[list[str]]) -> None:
        self.generated_codes_per_genome = generated_codes_per_genome
        
    def get_generated_codes_per_genome(self) -> list[list[str]]:
        return self.generated_codes_per_genome
    
    def set_result_matrix(self, result_matrix:list[list[float]]) -> None:
        self.result_matrix = result_matrix
        
    def get_result_matrix(self) -> list[list[float]]:
        return self.result_matrix
    
    def calculate_RSCU_per_gene(self) -> None:
        files:FileManager = self.get_processed_files()
        RSCU_per_gene:list[list[dict[str, float]]] = [[gene.get_RSCU() for gene in genome] for genome in files.get_converted_genomes()]
        RSCU_per_aa:list[list[dict[str, list[float]]]] = [[gene.get_RSCU_per_aa() for gene in genome] for genome in files.get_converted_genomes()]
        self.set_RSCU_per_gene(RSCU_per_gene)
        self.set_RSCU_per_aa(RSCU_per_aa)
        
    def map_codes_to_sequences_lenght(self):
        code_generator:CodeGenerator = CodeGenerator()
        sequences_lenght:list[list[int]] = self.get_sequences_lenght()
        RSCU_per_aa:list[list[dict[str, float]]] = self.get_RSCU_per_aa()
        generated_codes:list[list[str]] = code_generator.calculate_RSCU_scores_per_genome(RSCU_per_aa)
        genome_maps:list[dict[str, list[int]]] = []
        for genome_index, genome_codes in enumerate(generated_codes):
            codes_to_sequences_lenght_map:dict[str, list[int]] = {}
            for index, gene_code in enumerate(genome_codes):
                
                if codes_to_sequences_lenght_map.get(gene_code, 'unknown') == 'unknown':
                    codes_to_sequences_lenght_map[gene_code] = []
                
                if codes_to_sequences_lenght_map.get(gene_code, 'unknown') != 'unknown':
                    codes_to_sequences_lenght_map[gene_code].append(sequences_lenght[genome_index][index])
                    
            genome_maps.append(codes_to_sequences_lenght_map)

        self.set_mapped_codes_to_sequences_lenght_array(genome_maps)
        self.set_generated_codes_per_genome(generated_codes)
        
        
    def compare_pair_of_genomes(self, genome1:dict[str, list[int]], genome2:dict[str, list[int]]) -> float:
        shared_genomes:int= 0
    
    
        for k1,v1 in genome1.items():
            for k2, v2 in genome2.items():
                
                if k1 == k2:
                    contiguous_sum:int = (np.sum(v1) + np.sum(v2))
                    shared_genomes += contiguous_sum
                    
        return shared_genomes
        
    def compare_genomes(self) -> None:
        genomes:list[dict[str, list[int]]] = self.get_mapped_codes_to_sequences_lenght_array()
        num_genomes = len(genomes)
        result_matrix = [[0] * num_genomes for _ in range(num_genomes)]

        for i in tqdm(range(num_genomes)):
            for j in range(i + 1, num_genomes):
                if i != j:
                    total_genome_sum:int= 0
                    
                    for v1 in genomes[i].values():
                        total_genome_sum += np.sum(v1)
                        
                    for v2 in genomes[j].values():
                        total_genome_sum += np.sum(v2)
                        
                    shared_genomes = self.compare_pair_of_genomes(genomes[i], genomes[j])
                    result_matrix[i][j] = round(shared_genomes/total_genome_sum, 6)
                    
        for i in tqdm(range(num_genomes)):
            for j in range(i, num_genomes):
                result_matrix[i][j] = 1
                
        for i in tqdm(range(num_genomes)):
            for j in range(i):
                result_matrix[i][j] = result_matrix[j][i]
        
        self.set_result_matrix(result_matrix)
        
    def save_matrix(self, filename:str='output') -> None:
        col_names:list[str] = self.get_processed_files().get_genome_names()
        row_names:list[str] = self.get_processed_files().get_genome_names()
        matrix:list[list[float]] = self.get_result_matrix()
        with open(f'{filename}.txt', 'w') as file:
            file.write("\t" + "\t".join(col_names) + "\n")
            for i, row in enumerate(matrix):
                file.write(row_names[i] + "\t" + "\t".join(map(str, row)) + "\n")
        
    def process_files(self, file_path:str, genetic_code:dict[str, str], output_file_name:str) -> None:
        files:FileManager = FileManager()
        files.set_file_path(file_path)
        files.prepare_files_to_read()
        files.get_sequences_from_multiple_files()
        print('\nPreparing files to be processed ...')
        files.convert_to_genemodel()
        print('\nFiles processed successfully ...')
        [[gene.process_gene(genetic_code) for gene in genome]  for genome in files.get_converted_genomes()]
        self.set_sequences_lenght([[gene.get_sequence_lenght() for gene in genes ] for genes  in files.get_converted_genomes()])
        self.set_processed_files(files)
        self.calculate_RSCU_per_gene()
        self.map_codes_to_sequences_lenght()
        print('\nComparing genomes ...')
        self.compare_genomes()
        self.save_matrix(output_file_name)
    
    
if __name__ == '__main__':
    
    if len(sys.argv) != 3:
        print("USage: python chispisApp.py <data_path> <output_name>")
        
    else:
        data_path:str = str(sys.argv[1])
        output_name:str = str(sys.argv[2])
        chispis = ChispitApp()
        chispis.process_files(data_path, genetic_code, output_name)
    