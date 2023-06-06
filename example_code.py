import Bio.Entrez as Entrez
import time
from src.blastn_search import (
    read_config,
    blastn_and_parse,
)

# Read config and set Entrez details
email, api_key = read_config("/home/kjeong/.config/entrez/entrez_keys")
Entrez.email = email
Entrez.api_key = api_key

# Initialize the file path and the number of hit sequences to be saved
fasta_file_path = "./blastn_trial.fas"  # your 35 sequences fasta file
num_hit_seqs = 5  # number of hit sequences to save

# Call the function
start_time = time.time()
blastn_and_parse(fasta_file_path, num_hit_seqs)
end_time = time.time()
duration = end_time - start_time
print(f"Total {duration:.2f} seconds elapsed.")
print("Search and save completed.")
