# =========================================================================================
# This file is used to pre-proccess the bcontr.tmp file and find simple simplifications.
# So far this includes collecting lines together which are the same and introducing
# a factor.
# =========================================================================================
import argparse     # Parse arguments

# Parse arguments from gecco
parser = argparse.ArgumentParser(
                description="""Simplify ITF algo code, but collecting lines together which
                               are the same and proceed one another""",
                formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-i','--input',default=None,help='ITF binary contraction file')
parser.add_argument('-o','--output',default=None,help='ITF binary contraction file')
args = parser.parse_args()

if args.input is None:
    print("Error in python ITF processor: Must provide input file")
    exit(1)
if args.output is None:
    print("Error in python ITF processor: Must provide output file")
    exit(1)

inp = args.input
out = args.output

# Open bcontr.tmp file to read from
f1 = open(inp,"r")
# Open bcontr2.tmp file to read from
f2 = open(out,"w")

# Each line will be saved a compared to the previous one
old_line = ""
factor = 1.0

# Read each line of bcontr.tmp
for line in f1:

    if ("Error" in line):
        print("Error in translating GeCCo code to ITF, check bcontr.tmp file for more details", file=f2)
        quit()

    # If the current line is the same as the precceding line, increment a factor
    # + do not print out the line
    if (line == old_line):
        factor = factor + 1.0
        continue

    # line != old_line, so print out precceding with new factor
    if (old_line != "" or old_line != "BEGIN" or old_line != "END"):
        if (factor > 1.0):
            words = old_line.split()
            # TODO: check if factor already there...
            # if simplify to 1.0, dont print
            if ("*" in words[2]):
                wordsf = words[2].split("*",1)
                factor = factor * float(wordsf[0])
                if (factor == 1.0):
                    # No need to print factor
                    wordsf[0] = ""
                else:
                    wordsf[0] = str(factor) + "*"

                words[2] = "".join(wordsf)
            else:
                # No factor already, so simply add factor to the string
                words[2] = str(factor) + "*" + words[2]

            words[0] = " " + words[0]
            old_line = " ".join(words)
            print(old_line, file=f2)
            factor = 1.0
        else:
            # Lines are different, so no factor needed
            print(old_line, end="", flush=True, file=f2)

    old_line = line

f1.close()
f2.close()
