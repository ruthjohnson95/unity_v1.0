from unity_v1 import *
import pandas as pd 
import os
from optparse import OptionParser

# global logging
logging.basicConfig(format='%(asctime)s - %(message)s', datefmt='%d-%b-%y %H:%M:%S', level=logging.INFO)

def main():
    parser = OptionParser()
    parser.add_option("--s", "--seed", dest="seed", default="7")
    parser.add_option("--H", "--H", dest="H")
    parser.add_option("--N", "--N", dest="N", default=1000)
    parser.add_option("--id", "--id", dest="id", default="unique_id")
    parser.add_option("--its", "--ITS", dest="its", default=250)
    parser.add_option("--gwas_file", dest="gwas_file")
    parser.add_option("--outdir", dest="outdir")

    (options, args) = parser.parse_args()

    # set seed
    seed = int(options.seed)
    np.random.seed(seed)

    # get experiment params
    H = float(options.H)
    N = int(options.N)
    ITS = int(options.its)
    id = options.id
    gwas_file = options.gwas_file
    outdir = options.outdir

    logging.info("Using gwas file: %s" % gwas_file) 

    gwas_df = pd.read_csv(gwas_file,sep=' ')
    z = gwas_df['BETA_STD']

    p_init = smart_start(z, N)

    p_est, p_var, p_list = run_MCMC(p_init, H, N, z, ITS)

    logfile = os.path.join(outdir, id + '.' + str(seed)+".log")

    f = open(logfile, 'w')
    print_func("Estimate p: %.4f" % p_est, f)
    print_func("SD p: %.4g" % math.sqrt(p_var), f)

    f.close()


if __name__ == "__main__":
    main()
