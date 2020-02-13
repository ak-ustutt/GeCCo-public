# =========================================================================================
#
# =========================================================================================
import argparse     # Parse arguments
import os

def readlines_reverse(filename):
    with open(filename) as qfile:
        qfile.seek(0, os.SEEK_END)
        position = qfile.tell()
        line = ''
        while position >= 0:
            qfile.seek(position)
            next_char = qfile.read(1)
            if next_char == "\n":
                yield line[::-1]
                line = ''
            else:
                line += next_char
            position -= 1
        yield line[::-1]


parser = argparse.ArgumentParser(
                description="""Remove extra lines.""",
                formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('file',nargs='+',help='input file')
args = parser.parse_args()
files = vars(args)['file']

out = 'bcontr3.tmp'

temp_out = 'tmp.tmp'
f2 = open(temp_out,"w+")

begin = False
end = False

inters =[]

for line in readlines_reverse(files[0]):
    if line == '': continue

    if 'END' in line:
        print(line, file=f2)
        continue

    if 'BEGIN' in line:
        print(line, file=f2)
        continue

    words = line.split()
    if 'STIN' not in words[0]:
        begin = True

        if end:
            inters = []
            end = False

        for w in words:
            if 'STIN' in w:
                tensor = w.split('[')[0]
                inters.append(tensor)
        print(line, file=f2)

    elif 'STIN' in words[0]:
        end = True


        #print(inters)
        combined_inters = '\t'.join(inters)
        if words[0].split('[')[0].replace('.','') in combined_inters:
            print(line, file=f2)

            # Intermediates may depend on other intermediates
            for w in words:
                if 'STIN' in w:
                    tensor = w.split('[')[0]
                    inters.append(tensor)
            #print('Needed ', line, words[0].split('[')[0], combined_inters)

f2.close()

f2 = open(out,"w+")
# Read though temp file and reverse order
for line in readlines_reverse(temp_out):
    if line == '': continue

    print(line, file=f2)

f2.close()
os.remove('tmp.tmp')
