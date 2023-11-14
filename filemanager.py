import os
from genemodel import GeneModel

class FileManager:
    def __init__(self) -> None:
        self.file_path:str
        self.filenames:list[str]
        self.sequences_array:list[dict[str, str]]
        self.converted_genomes:list[list[GeneModel]]

    def get_file_path(self) -> str:
        return self.file_path

    def set_file_path(self, file_path:str) -> None:
        self.file_path = file_path
    
    def set_converted_genomes(self, converted_genomes:list[list[GeneModel]]) -> None:
        self.converted_genomes = converted_genomes
    
    def get_converted_genomes(self) -> list[list[GeneModel]]:
        return self.converted_genomes
    
    def get_filenames(self) -> list[str]:
        return self.filenames
        
    def set_filenames(self, filenames:list[str]) -> None:
        self.filenames = filenames
        
    def set_sequences_array(self, sequences_array:list[dict[str, str]]) -> None:
        self.sequences_array = sequences_array
        
    def get_sequences_array(self) -> list[dict[str, str]]:
        return self.sequences_array
        
    def prepare_files_to_read(self) -> None:
        file_path:str = self.get_file_path()
        file_to_read:list[str] = os.listdir(file_path)
        self.set_filenames(file_to_read)

    def get_genes_from_file(self, file_path:str) -> dict[str, str]:
        genes_from_fasta:dict[str, str] = {}
        with open(file_path, 'r') as multifasta:
            sequence_id:str = None
            for line in multifasta:
                line:str = line.strip()
                if line.startswith('>'):
                    sequence_id = line[1:]
                    genes_from_fasta[sequence_id] = ''
                else:
                    genes_from_fasta[sequence_id] += line
        return genes_from_fasta
    
    def get_sequences_from_multiple_files(self) -> None:
        files_to_read:list[str] = self.get_filenames()
        files_path:str = self.get_file_path()
        complete_path_array:list[str] = [os.path.join(files_path,file) for file in files_to_read]
        sequences:list[dict[str, str]] = [self.get_genes_from_file(file_path) for file_path in complete_path_array]
        self.set_sequences_array(sequences)
        self.set_filenames(complete_path_array)

    def convert_to_genemodel(self) -> None:
        genomes:list[dict[str, str]] = self.get_sequences_array()
        converted_genomes:list[list[GeneModel]] = []
        for genome in genomes:
            genome_objects:list[list[GeneModel]] = []
            for id, sequence in genome.items():
                newGenomeModel:GeneModel = GeneModel()
                newGenomeModel.set_gene_id(id)
                newGenomeModel.set_complete_seq(sequence)
                newGenomeModel.set_sequence_lenght()
                genome_objects.append(newGenomeModel)
            converted_genomes.append(genome_objects)
        self.set_converted_genomes(converted_genomes)

if __name__ == '__main__':
    files = FileManager()
    files.set_file_path('./data/')
    files.prepare_files_to_read()
    files.get_sequences_from_multiple_files()
    print(list(files.get_sequences_array()[0].keys())[0])
    files.convert_to_genemodel()
    print(files.get_converted_genomes()[0][1].get_gene_id())
    
    
