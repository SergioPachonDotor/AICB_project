class CodeGenerator:
    def __init__(self) -> None:
        self.model:dict = {
                        'A':'',
                        'C':'',
                        'D':'',
                        'E':'',
                        'F':'',
                        'G':'',
                        'H':'',
                        'I':'',
                        'K':'',
                        'L':'',
                        'M':'',
                        'N':'',
                        'P':'',
                        'Q':'',
                        'R':'',
                        'S':'',
                        'T':'',
                        'V':'',
                        'W':'',
                        'Y':'',
                        }

    def set_model(self, model:dict[str, str]) -> None:
        self.model = model
    
    def get_model(self) -> dict[str, str]:
        return self.model
    
    def create_code_per_gene(self, calculated_results:dict[str, int]) -> str:
        code_model:dict[str, str] = self.get_model()
        model_aminoacids:list[str] = list(code_model.keys())
        
        for aminoacid, score in calculated_results.items():
            if aminoacid not in model_aminoacids:
                code_model[aminoacid] = '0'
            elif aminoacid in model_aminoacids:
                code_model[aminoacid] = f'{score}'
        generated_code:str = ''.join(list(code_model.values()))
        return generated_code

    def calculate_aa_score(self, arr) -> int:
        max_value = arr[0]
        max_index = 0
        for i in range(1, len(arr)):
            if arr[i] > max_value:
                max_value = arr[i]
                max_index = i
        if max_value <= 0:
            return 0
        elif max_value < 1.2:
            return 1
        elif max_value >= 1.2:
            return max_index + 2

    def calculate_aa_scores_array(self, gene_RSCU:dict[str, list[float]]) -> dict[str, int]:
        RSCU_values:list[list[float]] = list(gene_RSCU.values())
        RSCU_aa_keys:list[str] = list(gene_RSCU.keys())
        output_scores:list[float] = [self.calculate_aa_score(score_arr) for score_arr in RSCU_values]
        return dict(zip(RSCU_aa_keys, output_scores))
    
    # Use This Function
    def calculate_RSCU_scores_per_genome(self, genomes:list[list[dict[str, float]]]) -> list[list[str]]:
        RSCU_scores_per_genome:list[list[str]] = []
        for genome in genomes:
            RSCU_in_genome:list[str] = []
            for gene in genome:
                gene_score:dict[str, int] = self.calculate_aa_scores_array(gene)
                gene_score_string:str = self.create_code_per_gene(gene_score)
                RSCU_in_genome.append(gene_score_string)
            RSCU_scores_per_genome.append(RSCU_in_genome)    
        return RSCU_scores_per_genome
