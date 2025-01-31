# coordinates_mapper

# T2G and G2T coordinates mapper

## About

This script maps **transcript coordinates to genomic coordinates** (T2G) and **genomic coordinates to transcript coordinates** (G2T) using CIGAR strings. 
It also supports **reverse strand mapping** for transcripts on the `-` strand.


## Design/Features

- **Parallel Processing** – Uses multi-threading for faster query processing
- **Memory Efficient** – Handles large files with **Pandas chunking** and **garbage collection**
- **Graceful Termination** – Supports **Ctrl+C (`SIGINT`) & `SIGTERM` signals** for cleanup
- **Logging & Debugging** – Structured logging for better debugging
- **Retry Mechanism** – Automatic retries for robustness
- **Scalability** – Designed to handle millions of records efficiently


## Flow of execution

- Script starts execution with if __name__ == "__main__" and reading the transcript and query file inputs
- The TranscriptMapper class reads the transcript file, parses CIGAR strings and supports **both forward (`+`) and reverse (`-`) strands**.
- The QueryProcessor class reads the query file and determines whether it's a **T2G (Transcript to Genome)** or **G2T (Genome to Transcript)** query, retrieves and maps coordinates, saves output to an output file
- Results are saved in text file and the program exits gracefully

## Setup

**1. Install Dependencies:**

bash:

```bash
pip install pandas numpy pytest
```

**2. Clone git repository:**

bash:

```bash
git clone https://github.com/mulini/coordinates_mapper.git
cd coordinates_mapper
```

## Usage

**1. File formats:**

**Inputs:**

Transcript file:

**Added strand for a plus or minus strand**

```txt
Transcript	Chromosome	Genomic_Start	CIGAR	Strand
TR1	CHR1	3	8M7D6M2I2M11D7M	+
TR2	CHR2	10	5M2D3M	-
```

Query file:

**Added type T2G or G2T for either transcript to genome or genome to transcript mapping**

```txt
Type	Transcript	Transcript_Coord	Chromosome	Genome_Coord
T2G	TR1	2		
T2G	TR2	3		
G2T		CHR1	5
G2T		CHR2	8
```

**2. Running the script:**


bash:

```bash
python coord_mapper.py --transcript_file sample_transcripts.txt --query_file sample_queries.txt --output_file results.txt
```


- --transcript_file transcripts.txt - Path to the transcript mapping file
- --query_file queries.txt - Path to the query file
- --output_file results.txt - Path to save the output results

**3. Output file:**

```txt
Type	Transcript	Transcript_Coord	Chromosome	Genome_Coord
T2G	TR1	2	CHR1	5
T2G	TR2	3	CHR2	8
G2T	TR1	2	CHR1	5
G2T	TR2	3	CHR2	8
```


**4. Running in debug mode:**

bash:

```bash
python coord_mapper.py --transcript_file sample_transcripts.txt --query_file sample_queries.txt --output_file results.txt 2>&1 | tee debug.log
```

**5. Running with parallel processing:**

bash:

```bash
export NUM_THREADS=8
python coord_mapper.py --transcript_file sample_transcripts.txt --query_file sample_queries.txt --output_file results.txt
```


## Testing with unit and integration tests

**1. Install pytest:**

bash:

```bash
pip install pytest
```

**2. Run all test:**

bash:

```bash
pytest -v test_coord_mapper.py
```

Tests included:

Unit tests:

- test_parse_cigar() - Validate CIGAR parsing
- test_build_transcript_to_genome_map() - tests forward strand transcript to genome mapping is accurate
- test_build_reverse_transcript_map()	 - tests reverse strand transcript to genome mapping is accurate
- test_get_genomic_coordinates() - Checks genomic coordinates are correctly returned
- test_get_transcript_coordinates() - Checks transcript coordinates are correctly returned

Integration tests:

- test_query_processor() - Tests entire workflow from query to output
- test_empty_transcript_file() - System behavior when transcript file is empty
- test_empty_query_file() - System behavior when query file is empty
