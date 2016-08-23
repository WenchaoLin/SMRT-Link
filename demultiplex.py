from pbcore.io import IndexedBamReader,FastaReader,FastqWriter
import sys

def main(parser):

    args = parser.parse_args()

    def makeFqName(bcPair):
        return '{}/{}--{}.fastq'.format(args.outDir,*[bcNames[i] for i in bcPair])

    bcNames     = {i:rec.name for i,rec in enumerate(FastaReader(args.barcodeFasta))}
    bcNames[-1] = 'NoBC'
    bam         = IndexedBamReader(args.ccsBAM)
    
    for bcPair in set(zip(bam.bcForward,bam.bcReverse)):
         with FastqWriter(makeFqName(bcPair)) as writer:
            for rec in bam[(bam.bcForward==bcPair[0])&(bam.bcReverse==bcPair[1])]:
                header = rec.readName 
                if not args.noBcQual: 
                    header += ' bq=%i'%rec.bcQual 
                writer.writeRecord(header,rec.read(aligned=False),rec.peer.query_qualities)
    
if __name__ == '__main__':

    import argparse
    import os

    parser = argparse.ArgumentParser(prog='demultiplex.py', description='Generate fastq files per barcode from barcoded BAM (with bc/bq tags)')
    parser.add_argument('ccsBAM', metavar='ccsBAM', type=str,
                    help='pacbio barcoded BAM file with .pbi index')
    parser.add_argument('barcodeFasta', metavar='barcodeFasta', type=str,
                    help='fasta of barcode sequences used to make bam (should match BAM header)')
    parser.add_argument('--outDir', dest='outDir', type=str, default=os.getcwd(),
                    help='output directory. default cwd' )
    parser.add_argument('--noBcQual', dest='noBcQual', action='store_true', default=False,
                    help='do not include barcde quality score in header. default include quailty')

    main(parser)

