import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize

from outils import *


class gp_class:

    def init_hypers(self, case='canonical'):
        self.x = None
        self.y = None
        if case == 'random':
            self.sigma = np.random.rand()*2
            self.gamma = np.random.rand()
            self.mu = np.random.rand()*0.1
            self.sigma_n = np.random.rand()*0.1
        elif case == 'canonical':
            self.sigma = 10
            self.gamma = 1/2
            self.mu = 0.1
            self.sigma_n = 0.1
            
    def show_hypers(self):
        print(f'gamma: {self.gamma}, i.e., lengthscale = {np.sqrt(1/(2*self.gamma))}')
        print(f'sigma: {self.sigma}')
        print(f'sigma_n: {self.sigma_n}')
        print(f'mu: {self.mu}')
        
        
    def sample(self, how_many=1):
        samples = np.random.multivariate_normal(self.mean, self.cov, size=how_many)
        self.samples = samples.T
        return self.samples

        
    def plot_samples(self, linestyle='-', v_axis_lims=None):
        if v_axis_lims == None:
            v_axis_lims = np.max(np.abs(self.samples))
        plt.figure(figsize=(9, 4))
        error_bars = 2 * self.sigma
        plt.fill_between(self.time, - error_bars, error_bars,
                         color='blue', alpha=0.1, label='95% barres d\'erreur')
        plt.plot(self.time, np.zeros_like(
            self.time), alpha=0.7, label='moyenne')
        if self.samples.shape[1] == 1:
            plt.plot(self.time, self.samples, linestyle, c='r', alpha=1)
        else:
            plt.plot(self.time, self.samples, linestyle, alpha=1)
        plt.title('observations du GP')
        plt.xlabel('temps')
        plt.legend(loc=1)
        plt.xlim([min(self.time), max(self.time)])
        plt.ylim([-v_axis_lims, v_axis_lims])
        plt.tight_layout()
        
        
    def load(self, x, y):
        self.Nx = len(x)
        self.x = x
        self.y = y
        
        
        
    def plot_data(self):
            
        plt.figure(figsize=(9, 3))
        plt.plot(self.x, self.y, '.r', markersize=8, label='données')
        plt.xlabel('temps')
        plt.legend(loc=1)
        plt.xlim([min(self.x), max(self.x)])
        plt.tight_layout()
        
        
    def plot_posterior(self, n_samples=0, v_axis_lims=None):

        plt.figure(figsize=(9, 3))
        plt.plot(self.time, self.mean, 'b', label='posterieur')

        plt.plot(self.x, self.y, '.r', markersize=8, label='données')
        error_bars = 2 * np.sqrt((np.diag(self.cov)))
        plt.fill_between(self.time, self.mean - error_bars, self.mean +
                         error_bars, color='blue', alpha=0.1, label='95% barres d\'erreur')
        if n_samples > 0:
            self.compute_posterior(where=self.time)
            self.sample(how_many=n_samples)
            plt.plot(self.time, self.samples, alpha=0.7)
        plt.title('Posterieur')
        plt.xlabel('temps')
        plt.legend(loc=1, ncol=3)
        plt.xlim([min(self.x), max(self.x)])
        plt.ylim([-35, 35])
        plt.tight_layout()
