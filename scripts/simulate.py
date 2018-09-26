#!/usr/bin/env python

from optparse import OptionParser
import numpy as np
import scipy.stats as st
import math
import sys
import os
import logging
import pandas as pd
import re

# setup global logging
logging.basicConfig(format='%(asctime)s - %(message)s', datefmt='%d-%b-%y %H:%M:%S', level=logging.INFO)


# simulates gwas effect sizes assuming an infintiesimal variance (h/M) with LD
def simulate_ivar_LD(p_sim, sigma_g_m, N, M, V):
    c = st.bernoulli.rvs(p=p_sim, size=M)
    true_p = (np.sum(c)/float(M))

    sd = math.sqrt(sigma_g_m)
    gamma = st.norm.rvs(loc=0, scale=sd, size=M)
    beta = np.multiply(gamma, c)
    h_sim = sigma_g_m * M
    sigma_e = (1-h_sim)/float(N)

    mu = np.matmul(V, beta)
    env_cov = np.multiply(np.eye(M), sigma_e)
    cov = np.matmul(V, env_cov)
    z = st.multivariate_normal.rvs(mean=mu, cov=cov)

    return z, c, gamma, true_p


# simulate gwas effects with infinitesimal variance (h/M) WITHOUT LD
def simulate_ivar_gw_noLD(p_sim, h_sim, N, M, sim_name, outdir):
    sigma_g_m = h_sim/float(M)

    c = st.bernoulli.rvs(p=p_sim, size=M)
    true_p = (np.sum(c)/float(M))

    sd = math.sqrt(sigma_g_m)
    gamma = st.norm.rvs(loc=0, scale=sd, size=M)
    beta = np.multiply(gamma, c)
    h_sim = sigma_g_m * M
    sigma_e = (1-h_sim)/float(N)

    mu = beta
    env_sd = np.sqrt(sigma_e)

    z = st.norm.rvs(loc=mu, scale=env_sd, size=M)

    logging.info("True prop causals: %.4f" % true_p)

    outfile = sim_name+'.gwas'
    full_outfile = os.path.join(outdir, outfile)

    z_df = pd.DataFrame(data=z, columns=['BETA_STD'])
    z_df.to_csv(full_outfile, sep=' ', index=False)
    logging.info("Saving simulated gwas to: %s" % outfile)

    return


# simulates GWAS across entire genome
def simulate_ivar_gw_LD(p_sim, h_sim, N, ld_list_file, rsid_list_file, ld_dir, bim_dir, outdir, M=None):

    total_causals_gw = 0
    M_gw = 0

	# count how many SNPs if not provided
    if M is None:
        with open(ld_list_file, 'r') as ld_list:
            for ld_file in ld_list:
                ld_file = os.path.basename(ld_file)
                ld_file = ld_file.rstrip()
                ld_file_b = os.path.join(ld_dir, ld_file)
                ld_b = np.loadtxt(ld_file_b)
                M_b = ld_b.shape[0]

                # add to running total of SNPs
                M_gw += M_b
                logging.info("Calculating total number of SNPs...")

    else: # user provided M
		M_gw = M

    with open(ld_list_file, 'r') as ld_list:
        with open(rsid_list_file, 'r') as rsid_list:
            for ld_file, rsid_file in zip(ld_list, rsid_list):
                ld_file = os.path.basename(ld_file)
                ld_file_b = os.path.join(ld_dir, ld_file)
                ld_file_b= ld_file_b.rstrip()

                sigma_g_m = h_sim/float(M_gw)
                ld_b = np.loadtxt(ld_file_b)

                logging.info("Simulating effect sizes using ld matrix: %s" % os.path.basename(ld_file_b))

                M_b = ld_b.shape[0]
                z_b, c_b, gamma_b, true_p = simulate_ivar_LD(p_sim, sigma_g_m, N, M_b, ld_b)

    			# keep running total of total causals
                total_causals_gw += np.sum(c_b)

                # open bim file containing BP/CHR info
                rsid_file = os.path.basename(rsid_file)
                rsid_file = rsid_file.rstrip() 
                bim_prefix = rsid_file.split('.')
#                CHR = re.sub("\D", "", chr_prefix[0])
                full_bim_file =  os.path.join(bim_dir, bim_prefix[0]+'.'+bim_prefix[1]+'.'+bim_prefix[2]+'.bim') 

                logging.info("Using bim file for SNP/BP info: %s" % full_bim_file) 

                df = pd.read_csv(full_bim_file, header=None, sep='\t')

                df.columns = ['CHR', 'SNP', 'POS', 'BP', 'Allele1', 'Allele2']
                df['BETA_STD'] = z_b

                locus_fname = bim_prefix[0]+'.'+bim_prefix[1]+'.'+bim_prefix[2]+'.gwas'

                locus_full_fname = os.path.join(outdir, locus_fname)

                df.to_csv(locus_full_fname, sep=' ', index=False)

                logging.info("Saving locus to: %s" % locus_fname)

    true_p = total_causals_gw/float(M_gw)
    logging.info("True prop causals: %.4f" % true_p)

    return


# header statement
def print_header(sim_name, p_sim, h_sim, N, outdir, rsid_list, ld_list, bim_dir, ld_dir, M):
    logging.info("Simulation name: %s" % sim_name)
    logging.info("Simulating with params -- p: %.4f, H2: %.4f, N: %d" % (p_sim, h_sim, N))

    if None in [rsid_list, ld_list, bim_dir, ld_dir]: # missing LD param
        logging.info("Did not find information for LD...simulating without LD")
        logging.info("Going to simulate %d SNPs" % M)


    else:
        logging.info("Simulating with LD found in dir: %s" % ld_dir)

    logging.info("Outputing simulated gwas to: %s" % outdir)

    return


def main():
    parser = OptionParser()
    parser.add_option("--sim_name", dest="sim_name", default="test_pipeline")
    parser.add_option("--h_sim", dest="h_sim", default=0.20)
    parser.add_option("--p_sim", dest="p_sim", default=0.05)
    parser.add_option("--N", dest="N", default=100000)
    parser.add_option("--rsid_list", dest="rsid_list")
    parser.add_option("--ld_list", dest="ld_list")
    parser.add_option("--outdir", dest="outdir")
    parser.add_option("--bim_dir", dest="bim_dir")
    parser.add_option("--ld_dir", dest="ld_dir")
    parser.add_option("--seed", dest="seed", default=100)
    parser.add_option("--M", dest="M")

    (options, args) = parser.parse_args()
    sim_name = options.sim_name
    h_sim = float(options.h_sim)
    p_sim = float(options.p_sim)
    N = int(options.N)
    seed = int(options.seed)
    outdir = options.outdir
    ld_list = options.ld_list
    rsid_list = options.rsid_list
    bim_dir = options.bim_dir
    ld_dir = options.ld_dir
    M = options.M
    if M != 'NA':
        M = int(M)

    # set the seed
    np.random.seed(seed)

    # print the header for user
    print_header(sim_name, p_sim, h_sim, N, outdir, rsid_list, ld_list, bim_dir, ld_dir, M)

    if 'NA' in [rsid_list, ld_list, bim_dir, ld_dir]:
        # no LD
        simulate_ivar_gw_noLD(p_sim, h_sim, N, M, sim_name, outdir)
    else:
        # LD
        simulate_ivar_gw_LD(p_sim, h_sim, N, ld_list, rsid_list, ld_dir, bim_dir, outdir)

    logging.info("FINISHED simulating")
    logging.info("Simulations can be found at: %s" % outdir)


if __name__== "__main__":
  main()
