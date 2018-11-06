import collections

import pomegranate as pmg
import numpy as np
import sklearn

import scanpy.api as sc
from tqdm import tqdm

import sc_utils
import sc_expression
import sc_clustering

cluster_spaces = {'controls': {'neighbor_params': {'metric': ['euclidean',
                                                              'cosine',
                                                              'correlation'],
                                                   'n_neighbors': [3, 5, 7] + 
                                                    list(range(10, 50 + 1, 5)),
                                                   'use_rep': ['X', None]},
                                'on_off_params': [{'on_off': False},
                                                  {'on_off': True,
                                                   'params':
                                                        {'method': 'mixture',
                                                         'test_fits': True}},
                                                  {'on_off': True,
                                                   'params':
                                                        {'method': 'mixture',
                                                         'test_fits': False}},
                                                  {'on_off': True,
                                                   'params':
                                                        {'method': 'otsu'}}],
                                 'hvg_params': {'method': ['gini',
                                                           'dispersion'],
                                                'percentile':
                                                    np.arange(0.0, 0.5, 0.05),
                                                'ignore_zeros': [True, False]}},

                  'clusterA': {'neighbor_params': {'metric': ['euclidean',
                                                              'cosine',
                                                              'correlation'],
                                                   'n_neighbors': [3, 5, 7] + 
                                                       list(range(10, 30, 5)),
                                                   'use_rep': ['X', None]},
                               'on_off_params': [{'on_off': False},
                                                 {'on_off': True,
                                                  'params':
                                                       {'method': 'mixture',
                                                        'test_fits': True}},
                                                 {'on_off': True,
                                                  'params':
                                                       {'method': 'mixture',
                                                        'test_fits': False}},
                                                 {'on_off': True,
                                                  'params':
                                                       {'method': 'otsu'}}],
                                'hvg_params': {'method':
                                                   ['gini', 'dispersion'],
                                               'percentile':
                                                   np.arange(0.0, 0.5, 0.05),
                                               'ignore_zeros':
                                                   [True, False]}}}

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

class Individual(object):

    def __init__(self, chromosome=None):
        self.chromosome = chromosome
        self.fitness = None

    def __repr__(self):
        msg = "Individual with the following genotype:\n"
        for key, item in self.chromosome.items():
            msg += "  {}: {}\n".format(key, item)
        msg += "  Fitness: {}\n".format(self.fitness)
        return msg
    
    def cross(self, partner):
        # ensure individuals are comparable
        if self.chromosome.keys() != partner.chromosome.keys():
            raise ValueError("Incomparable indivuals. Genes between" +
                             "individuals do not match.")
        # randomly choose breakpoint for genetic cross
        genes = list(self.chromosome.keys())
        break_point = np.random.choice(range(1, len(genes) - 1))

        # create chromosome for offsprings
        child1 = {x: None for x in genes}
        child2 = child1.copy()

        # cross genes between parents and children
        for x in genes[0:break_point]:
            child1[x] = self.chromosome[x]
            child2[x] = partner.chromosome[x]
        for y in genes[break_point:]:
            child1[y] = partner.chromosome[y]
            child2[y] = self.chromosome[y]

        return Individual(child1), Individual(child2)


class ClusterFitness(object):

    def __init__(self, data):
        # remove genes with all zereos, just in case
        sc.pp.filter_genes(data, min_cells=3)

        # pre-compute stuff
        self.ranked_genes = {'dispersion': None, 'gini':None}
        # rank hvg genes
        print('Pre-computing high variable genes via dispersion...')
        self.ranked_genes['dispersion'] = sc_utils.variable_genes(data,
                                                            method='dispersion',
                                                            percentile=0)
        print('Pre-computing high variable genes via gini metric...')
        self.ranked_genes['gini'] = sc_utils.variable_genes(data, method='gini',
                                                       percentile=0)

        print('Pre-computing principle components for raw data..')
        components = min(50, *data.shape)
        sc.pp.pca(data, n_comps=components)
        # pre-compute data
        self.data = {'raw': data, 'mixture.on_off': None, 'mixture.fit': None,
                     'otsu': None}
        print("Pre-computing on-off via mixture models...")
        self.data['mixture.on_off'] = sc_expression.set_on_off(data,
                                                               method='mixture', 
                                                               test_fits=False)
        print("Pre-computing on-off via mixture models and assessing fit...")
        self.data['mixture.fit'] = sc_expression.set_on_off(data,
                                                            method='mixture',
                                                            test_fits=True)
        print("Pre-computing on-off via otsu's method...")
        self.data['otsu'] = sc_expression.set_on_off(data, method='otsu')
        
    def score(self, individual):
        cluster_input = ClusterFitness.__parse_chromosome(individual.chromosome)
        # get expression data based off of input parameters
        if cluster_input['on_off_params']['on_off']:
            params = cluster_input['on_off_params']['params']
            if params['method'] == 'mixture':
                if params['test_fits']:
                    data = self.data['mixture.fit']
                else:
                    data = self.data['mixture.on_off']
            else:
                data = self.data[params['method']]
        else:
            data  = self.data['raw']
        
        # get genes
        genes =  self.ranked_genes[cluster_input['hvg_params']['method']]
        percentile = cluster_input['hvg_params']['percentile']
        idx = int(percentile*len(genes))
        
        # pass clustering parameters with variable genes specified
        hvg_df = sc_clustering.cluster_cells(data,
                               neighbor_kwargs=cluster_input['neighbor_params'],
                               genes=genes[idx:])

        individual.fitness = sklearn.metrics.silhouette_score(
                                                 self.data['raw'].obsm['X_pca'],
                                                 hvg_df.obs['louvain'])

    @staticmethod
    def __parse_chromosome_key(keys, value):
        if len(keys) > 1:
            new_value = {keys[-1]: value}
            return ClusterFitness.__parse_chromosome_key(keys[:-1], new_value)
        else:
            return {keys[0]: value}
    
    @staticmethod
    def __merge_dict(d1, d2):
        new_dict = {}
        for key in d1.keys():
            if key in d2.keys():
                if isinstance(d1[key], collections.MutableMapping) and\
                isinstance(d2[key], collections.MutableMapping):
                    new_dict[key] = {**d1[key], **d2[key]}
            else:
                new_dict[key] = d1[key]
        for key in d2.keys():
            if key not in d1.keys():
                new_dict[key] = d2[key]
        return new_dict

    @staticmethod
    def __parse_chromosome(chromosome):
        parameter_kwargs = {}
        for key in chromosome.keys():
            if '.' in key:
                kwargs = ClusterFitness.__parse_chromosome_key(key.split('.'),
                                                             chromosome[key])
                parameter_kwargs = ClusterFitness.__merge_dict(parameter_kwargs,
                                                               kwargs)
            else:
                parameter_kwargs[key] = chromosome[key]

        return parameter_kwargs
    

class GeneticAlgorithm(object):

    def __init__(self, parameter_space, pop_size=100, mutation_rate=0.03,
                 generations=500, verbose=True):
        """[summary]
        
        Parameters
        ----------
        parameter_space : dict
            A possibly nested dictionary containing parameter values to be
            evaluated. Keys pointing to non-dictionary entries should point to
            lists containing all possible test values. All keys should be
            strings with no "."s, as these are used when flattening input.

            Example:
                {'x': {'range': [0, 1, 2],
                       'domain': [3, 4, 5]},
                 'y': [-1, -2, -3]}
        pop_size : int, optional
            Total population size. The default is 100.
        mutation_rate : float, optional
            Probability to induce random mutation after a cross. The default is
            0.03, which will randomly mutate a single gene in a child in 3% of
            the crosses. 
        generations : int, optional
            Total number of generation to spawn. The default is 500.
        
        Attributes
        ----------
            genomic_space : dict
                Dictionary where each key represents a gene. Values are lists
                of possible values (alleles), each gene can exhibit. Dictionary
                is essentially a flattened version of `parameter_space`. 

                Example
                    {'x.range': [0, 1, 2],
                     'x.domain': [3, 4, 5],
                     'y': [-1, -2, -3]}

            population : list
                List of current individuals.

        Methods
        -------
            breed:
            mutate:
            optimize:
            random_individual:
        """

        self.__set_genomic_space(parameter_space)
        self.__initial_population(pop_size)
        self.generations = generations
        self.mutation_rate = mutation_rate
        self.verbose = verbose

    @staticmethod
    def __evaluate_key(key):
        if not isinstance(key, str):
            raise ValueError("Keys must be strings for easy concatentation.")
        if '.' in key:
            raise ValueError("Cannot have `.` in dictionary keys.")

    @staticmethod
    def flatten(d, parent_key='', sep='.'):
        """[summary]
        
        Parameters
        ----------
        d : [type]
            [description]
        parent_key : str, optional
            [description] (the default is '', which [default_description])
        sep : str, optional
            [description] (the default is '.', which [default_description])
        
        Returns
        -------
        [type]
            [description]

        References
        ----------

        Taken shamelessly from here:
            https://stackoverflow.com/questions/6027558/flatten-nested-python-dictionaries-compressing-keys
        """

        items = []
        for k, v in d.items():
            GeneticAlgorithm.__evaluate_key(k)
            new_key = parent_key + sep + k if parent_key else k
            if isinstance(v, collections.MutableMapping):
                items.extend(GeneticAlgorithm.flatten(v, new_key,
                                                      sep=sep).items())
            else:
                items.append((new_key, v))
        return dict(items)

    def __set_genomic_space(self, parameter_space):
        self.genomic_space = self.flatten(parameter_space, sep='.')

    def __initial_population(self, n=1000):
        """Randomally initialize a population of set size."""
        self.population = [None]*n
        for i in range(n):
            self.population[i] = self.random_individual()
    
    @staticmethod
    def __score2prob(scores):
        """
        Convert fitness scores to probabilities using the softmax function.
        
        Parameters
        ----------
        scores : numpy.array
            Fitness scores for each individual in the current population
        
        Raises
        ------
        ValueError
            Raised if `scores` is not castable to a numpy.array.
        
        Returns
        -------
        numpy.array
            Probability to select each individual within the population 
        """

        try:
            scores = np.array(scores)
        except:
            raise ValueError("`scores` must be castable to numpy array.")
        ex_scores = np.exp(scores - np.max(scores))
        return ex_scores / ex_scores.sum()

    def random_individual(self):
        """
        Create a random individual.

        Returns
        -------
        Individual
            Individual with randomally selected values for
            `Individual.chromosome`. 
        """

        chromosome = {key: np.random.choice(self.genomic_space[key])\
                      for key in self.genomic_space}
        return Individual(chromosome)

    def best_performer(self):
        scores = [x.fitness for x in self.population]
        return self.population[np.argsort(scores)[0]]

    def breed(self, fitness_function):
        for i in tqdm(range(self.generations)):
            new_population = []
            # calculate fitness for current population
            for each in self.population:
                fitness_function.score(each)
            # convert fitness scores to probabilities for pairing selection
            p_select = self.__score2prob([x.fitness for x in self.population])
            # select pairs
            pairs = np.random.choice(self.population,
                                     (int(len(self.population) / 2), 2),
                                     replace=True, p=p_select)
            # create children by crossing parents
            for parent1, parent2 in pairs:
                child1, child2 = parent1.cross(parent2)
                r1, r2 = np.random.random(2)
                if r1 <= self.mutation_rate:
                    child1 = self.mutate(child1)
                if r2 <= self.mutation_rate:
                    child2 = self.mutate(child2)
                new_population += [child1, child2]
            self.population = new_population

        # score final population
        for each in self.population:
            fitness_function.score(each)

    def mutate(self, individual):
        # choose random 'gene' to change
        key = np.random.choice(list(self.genomic_space.keys()))
        # mutate current gene value to a new value. 
        value = np.random.choice([x for x in self.genomic_space[key]\
                                  if individual.chromosome[key] != x])
        individual.chromosome[key] = value
        return individual
