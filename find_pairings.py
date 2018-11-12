import argparse
import os
import difflib

def createParser():
    parser = argparse.ArgumentParser(description='Find file pairs in directory')
    parser.add_argument('-I', '--input', type=str, nargs=1, 
                        help='Directory to find pairs. If left blank, current directory is used',
                        dest='input')
    parser.add_argument('-f', '--filter', type=str, nargs='*',
                        help='Limit which filetypes are searched. Will search only for file types specified',
                        dest='filter')

    args = parser.parse_args()

    return args

def parseArgs():
    args = createParser()

    if args.input is None:
        workingDirectory = os.getcwd()
    else:
        workingDirectory = ''.join(args.input)
    
    if args.filter is None:
        fileTypes = []
    else:
        fileTypes = args.filter

    arguments = {'workingDirectory' : workingDirectory,
                    'fileTypes' : fileTypes}

    return arguments

def findPairs(directory, fileTypes, N=1):
    # This only works if all paired files either end with _1 or _2 before the file extension
    files = [f for f in os.listdir(directory) if os.path.isfile(os.path.join(directory, f))]
    if len(fileTypes) > 0:
        for f in files:
            if not any(f.endswith(fileType) for fileType in fileTypes):
                files.remove(f)
        
    pairs = []
    for f in files:
        if '_1.' in f:
            match = f[0:f.rfind('_1.')] + "_2" + f[f.find('.'):]
            if match in files:
                pairs.append((os.path.join(directory, f), os.path.join(directory, match), os.path.join(directory, f[0:f.find('_1.')])))
    
    return pairs

def printPairs(pairs):
    for p in pairs:
        print('{}\t{}\t{}'.format(p[0], p[1],p[2]))

def main():
    args = parseArgs()
    pairs = findPairs(args['workingDirectory'], args['fileTypes'])
    printPairs(pairs)

if __name__ == '__main__':
    main()