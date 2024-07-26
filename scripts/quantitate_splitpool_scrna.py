#!/usr/bin/env python3
import gzip
import argparse
import sys
from pathlib import Path
import pysam

# This script is from the GitHub repository https://github.com/s-andrews/splitpoolquantitation 

def main():
    global options
    options = get_options()
    genes = parse_gtf(options.gtf_file)

    counts = count_overlaps(options.bam_file,genes)

    if options.output == "matrix":
        write_matrix_output(counts,options.outfile)

    elif options.output == "sparse":
        write_sparse_output(counts,genes,options.outfile)

def write_sparse_output(counts, allgenes, outdir):
    if not options.quiet:
        print("Writing counts to sparse folder",outdir, file=sys.stderr)

    # Create the folder if it doesn't exist
    outpath = Path(outdir)
    outpath.mkdir(parents=True, exist_ok=True)

    # Collect the full list of genes
    genes = set()
    item_count = 0
    for cell in counts.keys():
        for gene in counts[cell].keys():
            item_count += 1
            genes.add(gene)

    genes = sorted(list(genes))

    # Get the mapping of genes to ids
    genes_to_ids = {}
    for chr in allgenes.keys():
        for gene in allgenes[chr]:
            genes_to_ids[gene["name"]] = gene["id"]

    # Write out the feature list
    with gzip.open(outpath / "features.tsv.gz", "wt", encoding="utf-8") as out:
        for gene in genes:
            line = [genes_to_ids[gene], gene, "Gene Expression"]
            print("\t".join(line), file=out)

    # Collect all of the cell ids
    cell_ids = sorted(list(counts.keys()))

    # Write out the cell list
    with gzip.open(outpath / "barcodes.tsv.gz", "wt", encoding="utf-8") as out:
        for cell_id in cell_ids:
            print(cell_id, file=out)

    # Write out the counts as sparse values
    with gzip.open(outpath / "matrix.mtx.gz", "wt", encoding="utf-8") as out:
        print("%%MatrixMarket matrix coordinate integer general", file=out)
        print("%metadata_json: {\"software_version\": \"quantitate_splitpool_scrna\", \"format_version\": 2}", file=out)

        # Print totals of rows(genes), columns(cells), entries(measures)
        line = [len(genes),len(cell_ids),item_count]
        print("\t".join([str(x) for x in line]), file=out)

        for i,cell in enumerate(cell_ids):
            for j,gene in enumerate(genes):
                if gene in counts[cell]:
                    line = [j+1,i+1,len(counts[cell][gene])]
                    print("\t".join([str(x) for x in line]),file=out)

def write_matrix_output(counts, outfile):
    if not options.quiet:
        print("Writing counts to matrix",outfile, file=sys.stderr)

    # Collect the full list of genes
    genes = set()
    for cell in counts.keys():
        for gene in counts[cell].keys():
            genes.add(gene)

    genes = sorted(list(genes))

    with open(outfile,"wt") as out:

        # Print the header first
        cell_ids = sorted(list(counts.keys()))

        line=["Gene"]
        line.extend(cell_ids)
        print("\t".join(line), file=out)


        # Print the data
        for gene in genes:
            line = [gene]

            for cell in cell_ids:
           # bug? cell order are not match to printed 
           # for cell in counts.keys():
                if gene in counts[cell]:
                    line.append(str(len(counts[cell][gene])))
                else:
                    line.append("0")

            print("\t".join(line), file=out)





def count_overlaps(file,genes):
    if not options.quiet:
        print("Counting "+file)

    # We save and restore the verbosity so we can suppress an otherwise
    # annoying warning that we don't have an index (which we don't expect)
    save = pysam.set_verbosity(0)
    bamfile = pysam.AlignmentFile(file,"rb")
    pysam.set_verbosity(save)

    # We're going to keep track of what we've seen in this structure
    # we index by cell identity then record the Position+UMIs we've seen for
    # each gene so we can deduplicate

    counts = {}

    current_chromosome = ""
    chr_genes = []
    start_index = 0

    cell_count = 0

    for x,read in enumerate(bamfile.fetch(until_eof=True)):
        if x%100000 == 0:
            print("Searched",x,"reads",file=sys.stderr)

        # TESTING
        #if x==1000000:
        #    break

        chromosome = read.reference_name
        #check chromosome is not None, if is None, give a value, otherwise will cause AttributeError: 'NoneType' object has no attribute 'lower' 
        if chromosome is not None:
            chromosome = chromosome
        else:
            print("chromosome error")
            chromosome = "c"

        if chromosome.lower().startswith("chr"):
            chromosome = chromosome[3:]

        if chromosome != current_chromosome:
            if not options.quiet:
                print("Searching chromosome",chromosome)
            current_chromosome = chromosome
            if chromosome in genes:
                chr_genes = genes[chromosome]
            else:
                chr_genes = []

            start_index = 0

        for i in range(start_index,len(chr_genes)):
            # See if this is a feature we never need to look at again
            if i==start_index and read.reference_start > chr_genes[i]["end"]:
                start_index += 1
                continue
        
            # Can we skip based on strand
            if read.is_reverse != chr_genes[i]["is_reverse"]:
                continue

            # Can we stop looking
            if chr_genes[i]["start"] > read.reference_end:
                break

            # Do we overlap
            if read.reference_start < chr_genes[i]["end"] and read.reference_end > chr_genes[i]["start"]:
                # We need to check whether an alignment segment actually overlaps, rather
                # than just an intron
                found_block = False
                for rs,re in read.blocks:
                    if rs < chr_genes[i]["end"] and re > chr_genes[i]["start"]:
                        found_block = True
                        break

                if not found_block:
                    # None of our blocks actually overlap this gene
                    continue

                # We need to extract the UMI and cell identity from the read
                sections = read.query_name.split(":")
                umi = sections[-1]

                # We add the read start position to the umi for completeness
                if read.is_reverse:
                    umi += str(read.reference_end)
                else:
                    umi += str(read.reference_start)

                # We take the number of barcode sections from the options
                cell_id = ":".join(sections[(-1 - options.barcodes):-1])
                
                if not cell_id in counts:
                    cell_count += 1
                    counts[cell_id] = {}

                if not chr_genes[i]["name"] in counts[cell_id]:
                    counts[cell_id][chr_genes[i]["name"]] = set()

                if umi in counts[cell_id][chr_genes[i]["name"]]:
                    continue

                counts[cell_id][chr_genes[i]["name"]].add(umi)               
 
    if not options.quiet:
        print("Found",cell_count,"cells", file=sys.stderr)
    return counts

def parse_gtf(file):
    if not options.quiet:
        print(f"Parsing genes from {file}", file=sys.stderr)

    genes = {}

    gene_count = 0
    failed_count = 0

    infh = None

    if file.lower().endswith(".gz"):
        infh = gzip.open(file,"rt", encoding="utf8")

    else:
        infh = open(file,"rt",encoding="utf-8")

    for line in infh:
        if line.startswith("#"):
            continue

        sections = line.strip().split("\t")

        if sections[2] != "gene":
            continue

        chromosome = sections[0]
        if chromosome.lower().startswith("chr"):
            chromosome = chromosome[3:]

        start = int(sections[3])
        end = int(sections[4])
        is_reverse = sections[6] == "-"
        name = None
        id = None

        ##### TESTING ONLY #####
        #if chromosome != "1":
        #    break

        annotations = sections[8].split(";")
        for annot in annotations:
            if annot.strip().startswith("gene_id"):
                id = annot.strip().split(" ",2)[1].replace("\"","").strip()

            if annot.strip().startswith("gene_name"):
                name = annot.strip().split(" ",2)[1].replace("\"","").strip()
                
            if name is not None and id is not None:
                break

        if name is None:
            failed_count += 1
            continue

        if id is None:
            raise Exception("Couldn't find a gene id in '"+line+"'")

        else:
            gene_count += 1

        if not chromosome in genes:
            genes[chromosome] = []

        genes[chromosome].append({"start":start, "end":end, "is_reverse": is_reverse, "name":name, "id":id})

    infh.close()

    # It looks like the genes aren't listed in order in the gtf file so
    # we need to fix that.

    for n in genes.keys():
        g = genes[n]
        g.sort(key=lambda x: x["start"])


    if gene_count == 0:
        print(f"No genes found.  Rejected {failed_count} genes for having no gene_name attribute", file=sys.stderr)

    else:
        if not options.quiet:
            print(f"Found {gene_count} genes and rejected {failed_count} genes with no name", file=sys.stderr)

    return genes


            


def get_options():

    parser = argparse.ArgumentParser(description="Quantitate single cell spatial transcriptomic data")

    parser.add_argument("gtf_file", type=str, help="The GTF file of features")
    parser.add_argument("bam_file", type=str, help="The BAM file to quantitate")
    parser.add_argument("outfile", type=str, help="The output file (matrix) or folder (sparse) to write to")

    parser.add_argument("--quiet",action="store_true", default=False, help="Suppress progress messages")
    parser.add_argument("--barcodes",type=int, default=4, help="Number of embedded barcodes in cell id")
    parser.add_argument("--output", type=str, default="matrix", help="The format of the output - matrix(default) or sparse")

    return(parser.parse_args())


if __name__ == "__main__":
    main()
