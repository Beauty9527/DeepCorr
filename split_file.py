from filesplit.split import Split
import argparse

split = Split("mpileup_genome.pileup","./split_mpileup")
split.bylinecount(linecount=10000000)

def build_parser():
    """Setup option parsing for sample."""
    parser = argparse.ArgumentParser(description='Split the mpileup_genome.pileup by line.')
    parser.add_argument('--mpileup', type=str, help='The complete mpileup file')
    parser.add_argument('--output-folder', type=str, help='Path to output split_mpileup file.')

    args = parser.parse_args()

    return args

if __name__ == '__main__':
    parsed_args = build_parser()
    split = Split(parsed_args.mpileup, parsed_args.output_folder)
    split.bylinecount(linecount=10000000)
