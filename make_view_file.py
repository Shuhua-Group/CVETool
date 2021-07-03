from sys import argv

infile=argv[1]
outputfile=argv[2]

def get_len(file):
    with open(file,'r') as f:
        n=0
        length='len'
        for line in f:
            if n==0:
                x=line.split()
                length=len(x)
            n=n+1
    return length
    
def x(file,out,a):
    with open(file,'r') as f:
        outfile=open(out,'w')
        n=0
        genotype=[]
        for line in f:
            x=line.split()
            if n==0:
                for i in range(0,a):
                    svnum=x[i].split('SV')[1].split('_',1)[0]
                    snvnum=x[i].split('SNV')[1].split('_',1)[0]
                    samplenum=x[i].split('_')[2]
                    genotype.append('SV('+svnum+')/SNV('+snvnum+'):'+samplenum)
            else:
                for w in range(0,a):
                    if x[w]!='NA':
                        outfile.write(genotype[w]+'\t'+x[w]+'\n')
            n=n+1
    return 0
    
len=get_len(infile)
result=x(infile,outputfile,len)