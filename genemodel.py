from collections import Counter

class GeneModel:
    def __init__(self) -> None:
        self.gene_id:str
        self.sequence:str
        self.sequence_length:int
        self.codons:list[str]
        self.codon_frequencies:Counter
        self.genetic_code:dict[str, str]
        self.RSCU:dict[str, float]
        self.unknown_codons:list[str]
        self.aminoacid_sequence:str
        self.codons_per_aa:dict[str, list[str]]
        self.RSCU_per_aa:dict[str, list[float]]
    
    def set_genetic_code(self, genetic_code:dict[str, str]) -> None:
        self.genetic_code = genetic_code
    
    def get_genetic_code(self) -> dict[str, str]:
        return self.genetic_code
    
    def set_gene_id(self, gene_id: list[str]) -> None:
        self.gene_id = gene_id
        
    def get_gene_id(self) -> list[str]:
        return self.gene_id
    
    def set_sequence_lenght(self) -> None:
        sequence:str = self.get_complete_seq()
        self.sequence_length = len(sequence)
        
    def get_sequence_lenght(self) -> int:
        return self.sequence_length
    
    def set_complete_seq(self, sequence: str) -> None:
        self.sequence = sequence
        
    def get_complete_seq(self) -> str:
        return self.sequence
    
    def set_codons(self, codons:list[str]) -> None:
        self.codons = codons
        
    def get_codons(self) -> list[str]:
        return self.codons
    
    def set_codon_frequencies(self, codon_frequencies:Counter) -> None:
        self.codon_frequencies = codon_frequencies
        
    def get_codon_frequencies(self) -> Counter:
        return self.codon_frequencies
    
    def set_RSCU(self, RSCU:dict[str, float]) -> None:
        self.RSCU = RSCU

    def get_RSCU(self) -> dict[str, float]:
        return self.RSCU
    
    def set_unknown_codons(self, unknown_codons:list[str]) -> None:
        self.unknown_codons = unknown_codons
    
    def get_unknown_codons(self) -> list[str]:
        return self.unknown_codons
    
    def set_aminoacid_sequence(self, aminoacid_sequence:str) -> None:
        self.aminoacid_sequence = aminoacid_sequence
        
    def get_aminoacid_sequence(self) -> str:
        return self.aminoacid_sequence
    
    def set_codons_per_aa(self, codons_per_aa:dict[str, list[str]]) -> None:
        self.codons_per_aa = codons_per_aa
        
    def get_codons_per_aa(self) -> dict[str, list[str]]:
        return self.codons_per_aa
    
    def set_RSCU_per_aa(self, codons_per_aa:dict[str, list[float]]) -> None:
        self.codons_per_aa = codons_per_aa

    def get_RSCU_per_aa(self) -> dict[str, list[float]]:
        return self.codons_per_aa

    def calculate_codons(self) -> None:
        DNA_sequence:str = self.get_complete_seq()
        codons:list[str] = [DNA_sequence[i:i+3] for i in range(0, len(DNA_sequence), 3)]
        self.set_codons(codons)
    
    def get_complementary_sequence(self) -> str:
        complement_dict:dict[str, str] = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
        DNA_sequence:str = self.get_complete_seq()
        reverse_complement:str = ''.join(complement_dict[base] for base in reversed(DNA_sequence))
        return reverse_complement
    
    def translate_gene(self):
        filtered_genetic_code:dict[str, str] = self.get_genetic_code()
        codons:list[str] = self.get_codons()
        amino_acids:list[str] = []
        
        for codon in codons:
            amino_acid = filtered_genetic_code.get(codon, 'X')
            amino_acids.append(amino_acid)
            
        amino_acid_sequence:str = ''.join(amino_acids)
        self.set_aminoacid_sequence(amino_acid_sequence)

    def calculate_codon_frequencies(self) -> None:
        codons:list[str] = self.get_codons()
        frequencies = Counter(codons)
        self.set_codon_frequencies(frequencies)
    
    def calculate_RSCU(self) -> None:
        codon_frequencies:Counter = self.get_codon_frequencies()
        genetic_code:dict[str, str] = self.get_genetic_code()
        unknown_codons:list[str] = []
        rscu:dict[str, float] = {}
        
        for codon, observed_frequency in codon_frequencies.items():
            aa = genetic_code.get(codon, 'Unknown')
            
            if aa != 'Unknown' and aa != '*' and aa != 'M':
                synonymous_codons = [observed_codon for observed_codon in codon_frequencies if genetic_code.get(observed_codon, 'Unknown') == aa]
                n_i = len(synonymous_codons)
                rscu[codon] = round(observed_frequency / (1 / n_i * sum(codon_frequencies[syn_codon] for syn_codon in synonymous_codons)), 4)
            
            else:
                unknown_codons.append(codon)
                
        self.set_unknown_codons(unknown_codons)                
        self.set_RSCU(rscu)
        self.obtain_RSCU_per_aa()
        
    def obtain_RSCU_per_aa(self) -> None:
        rscu_per_aa:dict[str, list[float]] = {}
        codons_per_aa:dict[str, list[float]] = self.get_codons_per_aa()
        calculated_RSCU:dict[str, float] = self.get_RSCU()
        
        for aa, codons in codons_per_aa.items():
            rscu_per_aa[aa] = [calculated_RSCU[codon] if codon in calculated_RSCU else 0.0 for codon in codons]
            
        self.set_RSCU_per_aa(rscu_per_aa)

    def process_genetic_code(self, genetic_code:dict[str, str]) -> None:
        filtered_genetic_code:dict[str, str] = {codon : aa for codon, aa in genetic_code.items() if aa != '*' and aa != 'M'}
        codons_per_aminoacid:dict[str, list] = {}
        
        for codon, aa in filtered_genetic_code.items():
            
            if aa not in codons_per_aminoacid.keys():
                codons_per_aminoacid[aa]:list[str] = []
                codons_per_aminoacid[aa].append(codon)
                
            elif aa in codons_per_aminoacid.keys():
                codons_per_aminoacid[aa].append(codon)
        
        self.set_genetic_code(filtered_genetic_code)
        self.set_codons_per_aa(codons_per_aminoacid)
        
    def process_gene(self, genetic_code:dict[str, str]) -> None:
        self.process_genetic_code(genetic_code)
        self.calculate_codons()
        self.calculate_codon_frequencies()
        self.calculate_RSCU()
        self.get_RSCU_per_aa()