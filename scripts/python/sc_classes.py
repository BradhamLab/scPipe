import pomegranate as pmg
import numpy as np

class SignalNoiseModel(object):    

    def __init__(self, X):
        """
        Fit data to a mixture model to distinguish signal from technical noise.

        Fit data to a mixture model where technical noise of scRNAseq profiles
        is modelled using a Poisson distribution, and true expression is
        modelled as a Gaussian Distribution. 
        
        Parameters
        ----------
        X : numpy.array
            Single-cell expression profile across cells. Values are assumed to
            be log-transformed.
        
        Attributes
        ----------
        gmm: pomegranate.GeneralMixtureModel 
            Mixture of a Poisson and Normal Distribution to de-convolve noise
            and signal.
        range: list
            Two element list with range of possible values. Minimum is always
            zero, and max is set to the maximum observed value + spread.
            #TODO, find this by increasing x by dx until cdf(x) - 1 < eps

        Methods
        -------
        pdf : Calculate the probability of values within an array.
        cdf : Calculate cumulative density probabilities of values in an array.
        threshold: Find the first value such that P(Noise) < P(Signal)
        """
        #  use count data for mixtures
        counts, bins = np.histogram(X, bins=30)

        # estimate center as maximum non-zero count
        mode_idx = np.where(counts == np.max(counts[1:]))[0][0]
        # estimate 2SD as center - end value --> find values in normal dist.
        normal_spread = bins[-1] - bins[mode_idx]

        # find minimum value heuristically expected to be in signal
        noise_indices = np.where(bins < bins[mode_idx] - normal_spread)[0]
        if len(noise_indices) > 0:
            normal_min = bins[noise_indices[-1]]
        else:
            # no values below expected threshold, set to first non-zero value
            normal_min = bins[1]

        # estimate Normal distribution from likely normal samples
        signal = pmg.NormalDistribution.from_samples(X[X >= normal_min])
        # estimate Poisson from likely non-normal samples
        pois_samples = X[X < normal_min]
        percent_zeros = sum(X == 0) / len(X) 
        pois_lambda = percent_zeros
        if len(pois_samples) != 0:
            pois_lambda = max(np.mean(pois_samples), percent_zeros)
        noise = pmg.PoissonDistribution(pois_lambda)
        # instantiate and fit mixture to data
        gmm = pmg.GeneralMixtureModel([noise, signal])
        self.gmm = gmm.fit(X)
        self.range = [0, np.max(X) + normal_spread]

    def pdf(self, X):
        return self.gmm.probability(X)

    def cdf(self, X):
        """
        Calculate cumulative probabilities for values within an array.
        
        Parameters
        ----------
        X : np.array
            Array of values defined in the mixture domain.
        
        Returns
        -------
        np.array
            Array of cumulative densities for each provided value.
        """
        values = np.zeros_like(X)
        for i, x in enumerate(X):
            space = np.arange(0, x + 0.01, 0.01)
            values[i] = np.sum(self.pdf(space)*0.01)
        return values

    def threshold(self):
        space = np.arange(self.range[0], self.range[1], 0.01).reshape(-1, 1)
        p_x = self.gmm.predict_proba(space)
        idx = 0
        while p_x[idx][0] > p_x[idx][1] and idx < p_x.shape[0]:
            idx += 1
        return space[idx][0]