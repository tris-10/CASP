import argparse
import sys

def main():
    parser = argparse.ArgumentParser(description="This script parses the stats file generated by methphaser and"
                                                 "reports more detailed information for each reported phase block pair")
    parser.add_argument("input_stats",nargs="*")
    parser.add_argument("output_stats")
    parser.add_argument("--min_reads",default=10,type=int)
    parser.add_argument("--min_percentage",default=0.5,type=float)
    args = parser.parse_args()

    with open(args.output_stats, "w") as opf:
        for f in args.input_stats:
            with open(f, "r") as  ipf:
                header = ipf.readline()
                if not header.startswith(",snp_phased_block_1"):
                    continue
                ipf.readline()
                opf.write("Block1\tBlock2\tDistance\tCall\tTotal\tSame\tFlip\tDiff\tCpG1\tCpG2\tBalance\n")
                for line in ipf:
                    items = line.strip().split(",")
                    start1 = int(items[1][2:])
                    end1 = int(items[2][:-2])
                    start2 = int(items[3][2:])

                    distance = start2-end1
                    call = items[9] if items[9] != "not same" else "flip"
                    same = items[10]
                    flip = items[11]
                    cpg1 = items[14]
                    cpg2 = items[15]
                    total = "NA"
                    diff = "NA"
                    bias = "NA"
                    if same != "":
                        total = int(same) + int(flip)
                        same_p = int(same) / total
                        flip_p = int(flip) / total
                        diff = max([flip_p, same_p]) - min([flip_p, same_p])
                        bias = min([int(cpg1), int(cpg2)]) / (int(cpg1) + int(cpg2))
                        if total <= args.min_reads:
                            call = "low_overlap"
                        elif diff <= args.min_percentage:
                            call = "low_difference"
                        opf.write(
                            "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}%\t{8}\t{9}\t{10}%\n".format(start1, start2, distance, call,
                                                    total,same, flip, round(diff * 100, 2), cpg1, cpg2, round(bias * 100, 2)))
                    else:
                        call = "no_overlap"
                        opf.write(
                            "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\n".format(start1, start2, distance, call,
                                                    total,same, flip, diff, cpg1, cpg2, bias))



if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("User interrupted, exiting")
        sys.exit(1)
