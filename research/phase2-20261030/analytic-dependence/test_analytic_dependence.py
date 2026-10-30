"""Analytic identities and independent small-matrix checks, not coverage tests."""
import math
import numpy as np
import analytic_dependence as a


def test_exact_finite_n_variance_and_third_cumulant_against_dense_gaussian_form():
    p=a.parameters()
    for n in (1,2,9,32):
        for phi in (0.,.4,.997):
            r=phi**np.abs(np.arange(n)[:,None]-np.arange(n)[None,:])
            centroid_cov=r/p['beta']
            vc=.5*np.trace(centroid_cov@centroid_cov)/n**2
            third=np.trace(centroid_cov@centroid_cov@centroid_cov)/n**3+p['shape_third_cumulant']/n**2
            true=vc+p['shape_variance']/n
            assert math.isclose(a.exact_mc_variance(n,phi,p),true,rel_tol=1e-13)
            assert math.isclose(a.exact_mean_skewness(n,phi,p),third/true**1.5,rel_tol=1e-13)


def test_expected_bm_variance_against_dense_batch_covariance():
    p=a.parameters();n,b,phi=48,8,.97;blocks=n//b
    d=np.abs(np.arange(n)[:,None]-np.arange(n)[None,:])
    gamma=p['centroid_component_variance']*phi**(2*d)+np.eye(n)*p['shape_variance']
    averaging=np.zeros((blocks,n))
    for i in range(blocks):averaging[i,i*b:(i+1)*b]=1/b
    batch_cov=averaging@gamma@averaging.T
    centering=np.eye(blocks)-np.ones((blocks,blocks))/blocks
    expected=np.trace(centering@batch_cov)/(blocks*(blocks-1))
    assert math.isclose(a.exact_expected_bm_mcvariance(n,b,phi,p),expected,rel_tol=1e-13)


def test_independent_mode_moments_match_dense_harmonic_covariance():
    p=a.parameters();m=p['M'];beta=p['beta'];d=np.roll(np.eye(m),-1,axis=0)-np.eye(m)
    covariance=np.linalg.inv(beta/m*np.eye(m)+m/beta*d.T@d)
    shape=np.eye(m)/m-np.ones((m,m))/m**2
    cv=np.eye(m)/m-np.ones((m,m))/(2*m*m)
    assert math.isclose(p['shape_mean'],np.trace(shape@covariance),rel_tol=1e-13)
    assert math.isclose(p['shape_variance'],2*np.trace(shape@covariance@shape@covariance),rel_tol=1e-13)
    assert math.isclose(p['CV_variance'],2*np.trace(cv@covariance@cv@covariance),rel_tol=1e-13)
    assert math.isclose(a.exact_mc_variance(8192,0,p),p['CV_variance']/8192,rel_tol=1e-13)


def test_rng_roles_and_replicate_seeds_disjoint_and_generator_loop_matches_scalar():
    seeds={a.seed(c,i,role) for c in ('phi0p0','phi0p997','phi0p999') for i in range(20) for role in ('initial','centroid','shape')}
    assert len(seeds)==180
    p=a.parameters();n,phi=13,.997;c=.2
    compiled=a.simulate(np.random.default_rng(51),np.random.default_rng(52),c,n,phi,2.,p['shape_weights'])
    cg=np.random.default_rng(51);sg=np.random.default_rng(52);expected=[]
    for _ in range(n):
        c=phi*c+math.sqrt((1-phi*phi)/2)*cg.standard_normal()
        shape=sum(w*sg.standard_normal()**2 for w in p['shape_weights'])
        expected.append(.25+.5*c*c+shape)
    np.testing.assert_allclose(compiled,expected,rtol=0,atol=1e-15)
