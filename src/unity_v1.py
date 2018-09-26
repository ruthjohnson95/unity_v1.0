from auxilary import *


def log_q_pdf(p_star, p_old):
    B = 10 # noise constant
    alpha1 = beta_lam + p_old*B
    alpha2 = beta_lam + (1-p_old)*B
    log_q = st.beta.logpdf(x=p_star, a=alpha1, b=alpha2)
    return log_q


def q_rand(p_old):
    B = 10 # noise constant
    alpha1 = beta_lam + p_old*B
    alpha2 = beta_lam + (1-p_old)*B

    p_star = st.beta.rvs(a=alpha1, b=alpha2)

    return p_star


def accept_prob(z, p_star, p_old, H_gw, N):

    log_q_star = log_q_pdf(p_star, p_old)
    log_q_old = log_q_pdf(p_old, p_star)
    log_p_star = log_p_pdf_noLD(p_star, H_gw, z, N)
    log_p_old = log_p_pdf_noLD(p_old, H_gw, z, N)

    # DOUBLE CHECK THIS!!!
    r = (log_p_star - log_p_old) + (log_q_old - log_q_star)

    if r < EXP_MAX:
        R = math.exp(r)
    else:
        R = 100

    accept = min(1, R)

    return accept


def log_p_pdf_noLD(p, H_gw, z, N):

    M_gw = len(z)
    sigma_g = H_gw / float(M_gw)
    sigma_e = (1-H_gw)/float(N)

    mu = 0
    sd_1 = np.sqrt(sigma_g + sigma_e)
    sd_0 = np.sqrt(sigma_e)

    log_pdf_0 = st.norm.logpdf(x=z, loc=mu, scale=sd_0) # non-causal
    log_pdf_1 = st.norm.logpdf(x=z, loc=mu, scale=sd_1) # causal

    d_0 = np.add(log_pdf_0, np.log(1-p))
    d_1 = np.add(log_pdf_1, np.log(p))

    snp_pdfs = logsumexp_vector([d_0, d_1])

    log_pdf = np.sum(snp_pdfs)

    log_p = st.beta.logpdf(x=p, a=beta_lam,b=beta_lam)

    log_p_pdf = log_pdf + log_p

    return log_p_pdf


def smart_start(z, N):
    # concert betas to zscore
    zscores = np.multiply(z, np.sqrt(N))
    causal_inds_list = np.where(zscores >= 3.0)
    causal_inds = causal_inds_list[0]
    M = len(z)
    p_init = len(causal_inds)/float(M)
    return p_init


def run_MCMC(p_init, H, N, z, ITS):
    BURNIN = ITS/BURN

    # calculating acceptance probabilities
    ACCEPT_P = 0
    p_list = []

    # initialize
    p_old = p_init

    for i in range(0, ITS):

        # sample p
        p_star = q_rand(p_old)
        accept_p = accept_prob(z, p_star, p_old, H, N)

        # decide to accept or reject
        u = st.uniform.rvs(size=1)
        if u < accept_p: # ACCEPT
            p = p_star
            ACCEPT_P += 1
        else: # REJECT
            p = p_old

        # update with sample
        p_old = p

        # debugging
        if i % 10 == 0:
            logging.info("Iteration %d, p(t): %.4g" % (i, p_old))

        # save the values
        if i >= BURN:
            p_list.append(p_old)

        # end iteration loop

    # return posterior mean
    p_est = np.mean(p_list)
    p_std = np.std(p_list)

    return p_est, p_std, p_list
