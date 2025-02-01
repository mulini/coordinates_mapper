import os
import re
import sys
import pandas as pd
import logging
import argparse

# Configure logging for debugging
logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")


class TranscriptMapper:
    """
    Maps transcript coordinates to genomic coordinates and vice versa.
    Handles CIGAR strings, insertions, deletions, and reverse strand mapping.
    """

    def __init__(self, transcript_file):
        self.transcript_file = transcript_file
        self.transcript_mapping = {}  # Stores transcript → genome mappings
        self.genome_to_transcript_mapping = {}  # Stores genome → transcript mappings
        self.load_transcripts()

    def load_transcripts(self):
        """Loads transcript mappings from a file and processes CIGAR strings."""
        try:
            logging.info(f"Reading transcript data from {self.transcript_file}...")

            if not os.path.exists(self.transcript_file) or os.stat(self.transcript_file).st_size == 0:
                raise FileNotFoundError(f"Transcript file '{self.transcript_file}' is missing or empty.")

            df = pd.read_csv(self.transcript_file, sep="\t", dtype=str).fillna("")
            logging.info(f"Loaded {df.shape[0]} transcript entries.")

            for _, row in df.iterrows():
                try:
                    transcript_name, chromosome, genomic_start, cigar, strand = row
                    genomic_start = int(genomic_start)

                    if not cigar.strip():
                        logging.warning(f"Skipping {transcript_name}: Missing CIGAR string.")
                        continue

                    # Choose mapping function based on strand
                    if strand == "-":
                        mapping = self.build_reverse_transcript_map(genomic_start, cigar)
                    else:
                        mapping = self.build_transcript_to_genome_map(genomic_start, cigar)

                    self.transcript_mapping[transcript_name] = (chromosome, mapping)

                    # Reverse mapping: genome → transcript
                    for transcript_pos, genome_pos in mapping.items():
                        if isinstance(genome_pos, int):  # Ignore insertions
                            self.genome_to_transcript_mapping[(chromosome, genome_pos)] = (transcript_name, transcript_pos)

                except ValueError as ve:
                    logging.warning(f"Skipping invalid row: {row} -> {ve}")

        except Exception as e:
            logging.error(f"Error loading transcript file: {e}")
            sys.exit(1)

    @staticmethod
    def parse_cigar(cigar):
        """Parses a CIGAR string into a list of operations."""
        parsed_cigar = re.findall(r'(\d+)([MID])', cigar)
        if not parsed_cigar:
            raise ValueError(f"Invalid CIGAR format: '{cigar}'")
        return [(op, int(length)) for length, op in parsed_cigar]

    def build_transcript_to_genome_map(self, genomic_start, cigar):
        """
        Creates a mapping of transcript coordinates to genomic positions (forward strand).
        """
        transcript_pos = 0  # 0-based index (fixing off-by-one error)
        genomic_pos = genomic_start
        transcript_to_genome = {}

        for op, length in self.parse_cigar(cigar):
            if op == 'M':  # Matches (aligns directly)
                for _ in range(length):
                    transcript_to_genome[transcript_pos] = genomic_pos
                    transcript_pos += 1
                    genomic_pos += 1
            elif op == 'D':  # Deletion in transcript (skip genomic positions)
                genomic_pos += length
            elif op == 'I':  # Insertion in transcript (extra bases)
                for _ in range(length):
                    transcript_to_genome[transcript_pos] = f"Insertion before {genomic_pos}"
                    transcript_pos += 1

        return transcript_to_genome

    def build_reverse_transcript_map(self, genomic_start, cigar):
        """
        Creates a mapping for transcripts on the reverse strand.
        """
        transcript_pos = 0  # 0-based index
        genomic_pos = genomic_start
        transcript_to_genome = {}

        parsed_cigar = self.parse_cigar(cigar)

        for op, length in parsed_cigar:
            if op == 'M':  # Matches (aligns directly)
                for _ in range(length):
                    transcript_to_genome[transcript_pos] = genomic_pos
                    transcript_pos += 1
                    genomic_pos -= 1  # Reverse direction
            elif op == 'D':  # Deletion in transcript (skip genomic positions)
                genomic_pos -= length
            elif op == 'I':  # Insertion in transcript (extra bases)
                for _ in range(length):
                    transcript_to_genome[transcript_pos] = f"Insertion after {genomic_pos}"
                    transcript_pos += 1

        return transcript_to_genome

    def get_genomic_coordinates(self, transcript_name, transcript_coord):
        """Finds the genomic coordinate for a given transcript position."""
        if transcript_name in self.transcript_mapping:
            chromosome, mapping = self.transcript_mapping[transcript_name]
            return chromosome, mapping.get(transcript_coord, "Not Found")
        return "Not Found", "Not Found"

    def get_transcript_coordinates(self, chromosome, genome_coord):
        """Finds the transcript coordinate for a given genomic position."""
        return self.genome_to_transcript_mapping.get((chromosome, genome_coord), ("Not Found", "Not Found"))


class QueryProcessor:
    """
    Handles query processing for transcript-to-genome (T2G) and genome-to-transcript (G2T) mappings.
    """

    def __init__(self, mapper, query_file, output_file):
        self.mapper = mapper
        self.query_file = query_file
        self.output_file = output_file

    def process_queries(self):
        """Processes T2G and G2T queries."""
        if not os.path.exists(self.query_file) or os.stat(self.query_file).st_size == 0:
            raise ValueError(f"Query file '{self.query_file}' contains no valid rows.")

        df = pd.read_csv(self.query_file, sep="\t")

        output_data = []

        for _, row in df.iterrows():
            try:
                query_type = row["Type"]

                if query_type == "T2G":
                    transcript_name = row["Transcript"]
                    transcript_coord = int(row["Transcript_Coord"])
                    chromosome, genome_coord = self.mapper.get_genomic_coordinates(transcript_name, transcript_coord)

                elif query_type == "G2T":
                    chromosome = row["Chromosome"]
                    genome_coord = int(row["Genome_Coord"])
                    transcript_name, transcript_coord = self.mapper.get_transcript_coordinates(chromosome, genome_coord)

                output_data.append([query_type, transcript_name, transcript_coord, chromosome, genome_coord])

            except Exception as e:
                logging.warning(f"Skipping query row {row} -> {e}")

        self.save_output(output_data)

    def save_output(self, output_data):
        """Saves the processed queries to the output file."""
        if not output_data:
            logging.warning("No valid results. Check input files.")
            return

        df = pd.DataFrame(output_data, columns=["Type", "Transcript", "Transcript_Coord", "Chromosome", "Genome_Coord"])
        df.to_csv(self.output_file, sep='\t', index=False)
        logging.info(f"Results saved to {self.output_file}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Map transcript coordinates to genomic coordinates.")
    parser.add_argument("--transcript_file", required=True)
    parser.add_argument("--query_file", required=True)
    parser.add_argument("--output_file", required=True)
    args = parser.parse_args()

    transcript_mapper = TranscriptMapper(args.transcript_file)
    query_processor = QueryProcessor(transcript_mapper, args.query_file, args.output_file)
    query_processor.process_queries()

    logging.info("Processing completed successfully.")
