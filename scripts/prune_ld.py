from optparse import OptionParser
import os 
import gzip 


parser = OptionParser()
parser.add_option("-f", "--gwas_file", dest="gwas_file")
parser.add_option("-w", "--window", dest="window",default=5000)
parser.add_option("--bp_head", dest="bp_head",default='BP')

(options, args) = parser.parse_args()
gwas_file = options.gwas_file 
BP_head = options.bp_head
LD_WINDOW = int(options.window) 

bp_old = 0 
chr_old = 0

pruned_lines = [] 

with open(gwas_file, 'rb') as fp:

   # read in line from sumstats file 
   line = fp.readline()

   # get index from header 
   header = line.split() 
   SNP_index = header.index("SNP") 
   BP_index = header.index(BP_head)
   BETAS_index = header.index("BETA_STD") 
   CHR_index = header.index("CHR")

   while line: 
      params = line.split()
      bp_new = (params[BP_index])
      chr_new = (params[CHR_index])

      # write header 
      if bp_new == BP_head:
         pruned_lines.append(line)

      # if not the header 
      if bp_new != BP_head:
         bp_new = int(bp_new)
         chr_new = int(chr_new) 

         # if base pairs are a window apart or a different chromosome 
         if abs(int(bp_new) - int(bp_old)) >= LD_WINDOW or chr_old is not chr_new:

            bp_old = bp_new 
            chr_old = chr_new 
            
            # write whole line of snp info
            pruned_lines.append(line)

      line = fp.readline()

fp.close() 
 
# write file (overwrites intermediated gwas file)
f = open(gwas_file,'w')

for line in pruned_lines:
   f.write(line)

f.close()





