r"""
This module encodes a model for finding "edges" in photometric data, a common
step in the photometric calibration of images and catalogs.  In the following,
we describe the photometric model for the data with edges.

We assume we are measuring photometry in some collection of bands; thus we
measure a set of magnitudes for each object, 

.. math::
    \vec{A}_\mathrm{obs} = \left( A_0, A_1, \ldots, A_{N-1} \right)^T.

There is an associated "color" space, defined by 

.. math::
    \vec{C}_\mathrm{obs} \equiv \left( A_1 - A_0, A_2 - A_1, \ldots , A_{N-1} - A_{N-2} \right)^T
    
(note that the "vectors" :math:`\vec{A}` and :math:`\vec{C}` have different
dimension!). There are associated "true" magnitudes and colors for each object
:math:`\vec{A}` and :math:`\vec{C}` (these are "latent"---not directly
observed---due to noise in the instrument). 

We are looking for an "edge" in the photometry: there is a population of objects
with a "brightest" object, such that the population has support only for 

.. math::
    \vec{w} \cdot \vec{A} \geq e 
    
with 

.. math::
    \vec{w} = \left( 1 - c_0, c_0 - c_1,\ldots, c_{N-3}-c_{N-2}, c_{N-2}\right)^T 
    
or, equivalently, 

.. math::
    A_0 + \vec{c} \cdot \vec{C} \geq e. 
    
The vector :math:`\vec{c}` plays the role of a "color correction" to the
limiting magnitude in the 0th photometry band (which is singled out relative to
the other bands in the limiting magnitude :math:`e`).  

In fact, the code is a bit more sophisticated than this, and works with
"centered" colors, so that 

.. math::
    A_0 + \vec{c} \cdot \left( \vec{C} - \vec{C}_\mathrm{center} \right) \geq e_\mathrm{center}

is the edge.  The model outputs both :math:`e` and :math:`e_\mathrm{center}`.

We imagine that the population with the edge follows a multivariate Gaussian
with mean :math:`\vec{\mu}` and covariance :math:`\mathbf{\Sigma}` up to the
edge: 

.. math::
    p\left( \vec{A} \mid \vec{\mu}, \mathbf{\Sigma}, \vec{c}, e \right) \propto 
        \begin{cases} 
            N\left( \vec{A} \mid \vec{\mu}, \mathbf{\Sigma} \right) & \vec{w} \cdot \vec{A} \geq e \\ 
            0 & \mathrm{otherwise} 
        \end{cases}

We further imagine there is a second population of contaminants that do not
belong to the population with the edge, but nevertheless appear in our
photometric catalog.  For these we also assign a multivariate Gaussian
population, now with no edge or cut on the photometry 

.. math::
    p\left( \vec{A} \mid \vec{\mu}', \mathbf{\Sigma}' \right) \propto 
        N\left( \vec{A} \mid \vec{\mu}', \mathbf{\Sigma}' \right). 
        
The full catalog is a mixture of the two populations, with contamination
fraction :math:`f_\mathrm{bg}`: 

.. math::
    p\left( \vec{A} \mid \vec{\mu}, \vec{\mu}', \mathbf{\Sigma}, \mathbf{\Sigma}', \vec{c}, e, f_\mathrm{bg} \right) \propto 
        \left( 1 - f_\mathrm{bg} \right) p\left(\vec{A} \mid \vec{\mu}, \mathbf{\Sigma}, \vec{c}, e \right) + 
        f_\mathrm{bg} p\left( \vec{A} \mid \vec{\mu}', \mathbf{\Sigma}' \right)

If we assume Gaussian uncertainties on the photometry, 

.. math::
    p\left(\vec{A}_\mathrm{obs} \mid \vec{A} \right) \propto 
        N\left( \vec{A}_\mathrm{obs} \mid \vec{A}, \mathbf{\Sigma}_\mathrm{obs} \right),
        
which (for now, though it is easy to generalize) we will take to have
independent noise in each band, 

.. math::
    \mathbf{\Sigma}_\mathrm{obs} = \mathrm{diag} \, \left( \vec{\sigma}_\mathrm{obs}^2 \right), 
    
then we can analytically marginalize over the latent true photometric
magnitudes, :math:`\vec{A}`.  

For the foreground popluation with an edge, we obtain

.. math::
    p\left(\vec{A}_\mathrm{obs} \mid \vec{\mu}, \mathbf{\Sigma}, \vec{c}, e \right) = 
        \int \mathrm{d} \vec{A} \, p\left( \vec{A}_\mathrm{obs} \mid \vec{A} \right) p\left(
\vec{A} \mid \vec{\mu}, \mathbf{\Sigma}, \vec{c} \right) = 
        N\left(\vec{A}_\mathrm{obs} \mid \vec{\mu}, \mathbf{\Sigma} +
        \mathbf{\Sigma}_\mathrm{obs} \right) \alpha 
        
where the normalization :math:`\alpha` is 

.. math::
    \alpha = \frac{1 + \mathrm{erf}\,\left(\frac{\left( e_\mathrm{obs} - e\right) \sigma_e^2 + \left(\mu_e - e \right) \sigma_{e,\mathrm{obs}}^2}{\sqrt{2} \sigma_e \sigma_{e,\mathrm{obs}} \sqrt{\sigma_e^2 + \sigma_{e,\mathrm{obs}}^2}} \right)}
        {1 + \mathrm{erf}\left( \frac{\mu_e - e}{\sqrt{2} \sigma_e} \right)},

where 

.. math::
    \mu_e & = \vec{w} \cdot \vec{\mu}, \\
    e_\mathrm{obs} & = \vec{w} \cdot \vec{A}_\mathrm{obs}, \\
    \sigma_e^2 & = \vec{w}^T \mathbf{\Sigma} \vec{w}, \\
    \sigma_{e,\mathrm{obs}^2} & = \vec{w}^T \mathbf{\Sigma}_\mathrm{obs}\vec{w} = \sum_i w_i^2 \sigma_{\mathrm{obs},i}^2 
    
(the last equality holding only for independent photometric uncertainties).

For the background population, we obtain

.. math::
    p\left( \vec{A}_\mathrm{obs} \mid \vec{\mu}', \mathbf{\Sigma}'\right) = 
        \int \mathrm{d} \vec{A} \, p\left(\vec{A}_\mathrm{obs} \mid \vec{A} \right) p\left( \vec{A} \mid \vec{\mu}', \mathbf{\Sigma}' \right) = N\left( \vec{A}_\mathrm{obs} \mid \vec{\mu}', \mathbf{\Sigma}' + \mathbf{\Sigma}_\mathrm{obs} \right).

Putting the mixture together, the marginal likelihood for an observation
:math:`\vec{A}_\mathrm{obs}` under our population model is 

.. math::
    p\left(\vec{A}_\mathrm{obs} \mid \vec{\mu}, \vec{\mu}', \mathbf{\Sigma}, \mathbf{\Sigma}', \vec{c}, e, f_\mathrm{bg} \right) \propto 
        \left( 1 - f_\mathrm{bg} \right) N\left( \vec{A}_\mathrm{obs} \mid \vec{\mu}, \mathbf{\Sigma} + \mathbf{\Sigma}_\mathrm{obs} \right) \alpha + 
        f_\mathrm{bg} N\left(\vec{A}_\mathrm{obs} \mid \vec{\mu}', \mathbf{\Sigma}' + \mathbf{\Sigma}_\mathrm{obs} \right)

There is one final trick for numerical stability of the gradient.  When
evaluating the likelihood function, we are going to have to compute
:math:`\alpha`, which will involve computation of functions like

.. math::
    f(x) = \log\left( 1 + \mathrm{erf}\,\left(x\right) \right). 
    
For some objects, there will be very small values of :math:`x`, because they
will be overwhelmingly unlikely to be in the foreground component.  As :math:`x
\to -\infty`, :math:`\mathrm{erf}\,\left( x\right) \to -1`, and the above
computation becomes subject to a lot of roundoff error. Additionally, long
before :math:`f(x)` is not representable the gradient computed naively will have
a :math:`0/0` singularity. To avoid these failures, we switch to a series
expansion and custom derivative of the above function whenever :math:`x < 3`: 

.. math::
    f(x) = \begin{cases} 
        \log\left( 1 + \mathrm{erf}\,\left(x\right) \right) & x \geq -3 \\ 
        -x^2 - \log\left( - \sqrt{\pi} x \right) - \frac{1}{2 x^2} & x < -3
    \end{cases}
    
whence 

.. math::
    f'(x) = \begin{cases} 
        \frac{2}{\sqrt{\pi}} \frac{e^{-x^2}}{1 + \mathrm{erf}\,\left( x \right)} & x \geq -3 \\
        -2 x - \frac{1}{x} + \frac{1}{x^3} & x < -3 
    \end{cases} 
"""

__all__ = ['edge_model', 'jax_prng_key']

from jax import custom_jvp
import jax
import jax.numpy as jnp
import jax.scipy.special as jss
import numpy as np
import numpyro
import numpyro.distributions as dist

def jax_prng_key(seed=None):
    if seed is None:
        seed = np.random.randint(1<<32)
    return jax.random.PRNGKey(seed)

@custom_jvp
def log1p_erf(x):
    x = jnp.array(x)
    return jnp.where(x < -4.0, -x*x - jnp.log(-jnp.sqrt(np.pi)*x) + 1/(x*x)*(-0.5 + 1/(x*x)*(5.0/8.0 - 37.0/(24.0*x*x))), jnp.log1p(jss.erf(x)))

@log1p_erf.defjvp
def log1p_erf_jvp(primals, tangents):
    x, = primals
    dx, = tangents

    ans = jnp.where(x < -4.0, -x*x - jnp.log(-jnp.sqrt(np.pi)*x) + 1/(x*x)*(-0.5 + 1/(x*x)*(5.0/8.0 - 37.0/(24.0*x*x))), jnp.log1p(jss.erf(x)))
    ans_dot = jnp.where(x < -4.0, -2*x + 1/x*(-1 + 1/(x*x)*(1 + 1/(x*x)*(-5.0/2.0 + 37.0/(4.0*x*x)))), 2/jnp.sqrt(np.pi)*jnp.exp(-x*x)/(1 + jss.erf(x)))
    return ans, ans_dot*dx

def log_edge_normalization_factor(e, mu_e, sigma_e, e_obs, sigma_e_obs):
    sigma_e2 = sigma_e*sigma_e
    sigma_e_obs2 = sigma_e_obs*sigma_e_obs
    
    obs_arg = ((e_obs - e)*sigma_e2 + (mu_e - e)*sigma_e_obs2)/(jnp.sqrt(2)*sigma_e*sigma_e_obs*jnp.sqrt(sigma_e2 + sigma_e_obs2))
    pop_arg = (mu_e-e)/(jnp.sqrt(2)*sigma_e)

    log_numer = log1p_erf(obs_arg)
    log_denom = log1p_erf(pop_arg)

    return log_numer - log_denom

def mean_sample(postfix, mean_vec, scale_vec):
    N = mean_vec.shape[0]
    mu_unit = numpyro.sample('mu_unit_' + postfix, dist.Normal(loc=0, scale=1), sample_shape=(N,))
    mu = numpyro.deterministic('mu_' + postfix, mean_vec + scale_vec*mu_unit)

    return mu, mu_unit

def covariance_sample(postfix, scale_vec, eta=1):
    N = scale_vec.shape[0]

    scale_unit = numpyro.sample('scale_unit_' + postfix, dist.HalfNormal(scale=1), sample_shape=(N,))
    scale = numpyro.deterministic('scale_' + postfix, scale_unit * scale_vec)
    corr_cholesky = numpyro.sample('corr_cholesky_' + postfix, dist.LKJCholesky(N, eta))
    cov_cholesky = numpyro.deterministic('cov_cholesky_' + postfix, scale[:,None]*corr_cholesky)
    cov = numpyro.deterministic('cov_' + postfix, jnp.matmul(cov_cholesky, cov_cholesky.T))

    return cov, cov_cholesky, corr_cholesky, scale, scale_unit

def edge_model(Aobs, cov_obs, e_center_mu=0.0, e_center_sigma=1.0, c_mu=None, c_sigma=None, c_center=None, mu_bg=None, cov_bg=None, f_bg=None, nu_lkj=1):
    Aobs = np.array(Aobs)
    cov_obs = np.array(cov_obs)

    nobs, nband = Aobs.shape
    assert cov_obs.shape == (nobs, nband, nband), 'size mismatch between `Aobs` and `cov_obs`'

    A_mu = np.mean(Aobs, axis=0)
    sigma_A = np.std(Aobs, axis=0)

    if f_bg is None:
        f_bg = numpyro.sample('f_bg', dist.Uniform())

    if c_mu is None:
        c_mu = np.zeros(nband-1)
    if c_sigma is None:
        c_sigma = np.zeros(nband-1)
    if c_center is None:
        c_center = np.zeros(nband-1)

    c_unit = numpyro.sample('c_unit', dist.Normal(loc=0, scale=1), sample_shape=(nband-1,))
    c = numpyro.deterministic('c', c_mu + c_sigma*c_unit)
    if nband == 2:
        w = numpyro.deterministic('w', jnp.array([1-c[0], c[0]]))
    else:
        w = numpyro.deterministic('w', jnp.concatenate((jnp.array([1-c[0]]), -jnp.diff(c), jnp.array([c[-1]]))))

    e_unit = numpyro.sample('e_unit', dist.Normal(loc=0, scale=1))
    e_centered = numpyro.deterministic('e_centered', e_center_mu + e_center_sigma*e_unit)
    e = numpyro.deterministic('e', e_centered + jnp.dot(c, c_center))

    mu_fg, _ = mean_sample('fg', A_mu, sigma_A)
    cov_fg, _, _, _, _ = covariance_sample('fg', sigma_A, nu_lkj)

    if mu_bg is None and cov_bg is None:
        mu_bg, _ = mean_sample('bg', A_mu, sigma_A)
        cov_bg, _, _, _, _ = covariance_sample('bg', sigma_A, nu_lkj)
    elif mu_bg is None or cov_bg is None:
        raise ValueError('either both `mu_bg` and `cov_bg` must be `None` or neither can be `None`')

    mu_e = jnp.dot(w, mu_fg)
    e_obs = jnp.dot(Aobs, w)
    sigma_e2 = jnp.sum(w[:,None]*cov_fg*w[None,:])
    sigma_e = jnp.sqrt(sigma_e2)
    sigma_e_obs2 = jnp.sum(w[None,:,None]*w[None,None,:]*cov_obs, axis=(1, 2))
    sigma_e_obs = jnp.sqrt(sigma_e_obs2)

    log_alpha = numpyro.deterministic('log_alpha', log_edge_normalization_factor(e, mu_e, sigma_e, e_obs, sigma_e_obs))

    cov_fg_obs = cov_fg[None,:,:]+cov_obs
    cov_bg_obs = cov_bg[None,:,:]+cov_obs

    log_f_fg = jnp.log1p(-f_bg)
    log_f_bg = jnp.log(f_bg)

    logp_fg = numpyro.deterministic('logp_fg', log_f_fg + dist.MultivariateNormal(loc=mu_fg[None,:], covariance_matrix=cov_fg_obs).log_prob(Aobs) + log_alpha)
    logp_bg = numpyro.deterministic('logp_bg', log_f_bg + dist.MultivariateNormal(loc=mu_bg[None,:], covariance_matrix=cov_bg_obs).log_prob(Aobs))
    logp_total = jnp.logaddexp(logp_fg, logp_bg)
    
    log_fg_prob = numpyro.deterministic('log_fg_prob', logp_fg - logp_total)

    numpyro.factor("likelihood", jnp.sum(logp_total))
