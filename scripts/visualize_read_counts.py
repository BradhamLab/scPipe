"""
Script to visualize and organize read counts after quality control with fastp.

@author: Dakota Hawkins
@date: July 27, 2018
"""

# system imports
import os

# file io imports 
import json

# patten matching
import re

# numerical libraries
import numpy as np
from scipy.interpolate import interp1d
import statsmodels.api as sm
import pandas as pd

# plotting libraries
import seaborn as sns
from matplotlib import pyplot as plt
import matplotlib.patches as mpatches

class SummarizeFastpReads(object):
    """
    Class to summarize output from read quality control using `fastp`.
    """
    
    def __init__(self, output_dir, sample_regex, bad_threshold, ugly_threshold):
        """
        Class to summarize output from read quality control using `fastp`

        Creates a SnakeMake html report with three charts:
            
            1. A violin plot showing read count distributions between treatment
               groups after `fastp` quality control.
            2. A CDF plot showing percentage of cells that fall into 'good',
               'bad', and 'ugly' groups depending on provided thresholds of
               read counts.
            3. A stacked bar plot showing proportion of cells from each
               treatment placed in 'good', 'bad', and 'ugly' groups.
        
        A .csv file is also created with summarizing read count, read length,
        treatment, and quality placement.

        Args:
            output_dir (string): path to output directory containing `fastp`
                output. Directory is assumed to use the following format:

                    <output_dir>/<sample_id/>

                Where a file, `fastp.json` exists within each <sample_id>
                subdirectory.
            sample_regex (list, string): list of regex patterns to extract all
                possible treatments from a sample name.
            bad_threshold (int): maximum number of reads for a cell to be
                considered as 'bad'.
            ugly_threshold (int): maximum number of reads for a cell to be
                considered as 'ugly'.
        """
        self.output_dir = output_dir
        self.sample_regex = [re.compile(x) for x in sample_regex]
        self.bad_threshold = bad_threshold
        self.ugly_threshold = ugly_threshold
        self.read_df = self.generate_dataframe()


    def generate_dataframe(self):
        """
        Create a dataframe of read data from `fastp` json files.

        Return:
            (pd.DataFrame)
            Dataframe with 7 columns:
                pre.filter.total: total number of reads for a given cell pre
                    quality control.
                pre.filter.total.log10: log10 of total number of reads for a
                    given cell pre quality control.
                pre.filter.length: average read length of both read ends for
                    each cell pre quality control.
                post.filter.total: total number of reads for given cell post
                    quality control.
                post.filter.total.log10: log10 of total number of reads for a
                    given cell post quality control.
                post.filter.length: average read length of both read ends for
                    each cell post quality control.
                treatment: experimental treatment applied to each cell. 
        """
        all_samples = {}
        for path, dirs, files in os.walk(self.output_dir):
            for each in files:
                if each == 'fastp.json':
                    sample_name = os.path.basename(path)
                    each = os.path.join(path, each)
                    all_samples = {**all_samples,
                                   **self.get_read_data(each, sample_name)}
        df = pd.DataFrame(all_samples).T
        df['pre.filter.total.log10'] = [np.log10(x) for x in\
                                        df['pre.filter.total']]
        df['post.filter.total.log10'] = [np.log10(x) for x in \
                                         df['post.filter.total']]
        return df

    def treatment_from_sample_name(self, sample_name):
        """
        Retrieve treatment from a sample name.

        Args:
            sample_name (string): sample name/id.
        Output:
            (string): treatment applied to `sample`. 
        """

        for pattern in self.sample_regex:
            pattern_match = re.search(pattern, sample_name)
            if pattern_match is not None:
                return pattern_match.group()
        
        return None

    def get_read_data(self, fastp_json, sample_name):
        """
        Get sample read data from a `fastp` json file.

        Get sample read data from a `fastp` json file. This function is used to
        create a dataframe of read information for every sample. See 
        `generate_dataframe()` for more information.

        Args:
            fastp_json (string): json file created by `fastp`.
            sample_name (string): name of the current sample
        Return:
            (dict): dictionary with keys, values:
                1. pre.filter.total | total number of reads pre qc filtering.
                2. pre.filter.length | average read length pre qc filtering. 
                3. post.filter.total | total number of reads post qc filtering.
                4. post.filter.length | average read length post qc filtering.
                5. treatment | treatment applied to cell.
        """
        sample_dict = {sample_name: {}}
        read_dict = {}
        with open(fastp_json, 'r') as f:
            read_dict = json.load(f)

        pre_filter = read_dict['summary']['before_filtering']
        r_tot = pre_filter['total_reads']
        r_len = np.mean([pre_filter['read1_mean_length'],
                        pre_filter['read2_mean_length']])
        
        post_filter = read_dict['summary']['after_filtering']
        fr_tot = post_filter['total_reads']
        fr_len = np.mean([post_filter['read1_mean_length'],
                        post_filter['read2_mean_length']])

        treatment = self.treatment_from_sample_name(sample_name)
        key_values = zip(['pre.filter.total', 'pre.filter.length',
                          'post.filter.total', 'post.filter.length',
                          'treatment'],
                         [r_tot, r_len, fr_tot, fr_len, treatment])
        sample_dict[sample_name] = {x:y for x, y in key_values}
        return sample_dict

    def create_violin_plot(self):
        """
        Visualize post qc log10 read distribution between treatment types.
        """ 
        sns.set(style='whitegrid')
        n_values = self.read_df['treatment'].value_counts()
        colors = sns.color_palette()
        treatments = sorted(n_values.index.values)
        patches = []
        for i, each in enumerate(treatments):
            color = colors[i % len(treatments)]
            label = "{} (n={})".format(each, n_values[treatments[i]])
            patches.append(mpatches.Patch(color=color, label=label))
        ax = sns.violinplot(data=self.read_df, x='treatment',
                            y='post.filter.total.log10')
        ax.set(xlabel='Treatment', ylabel='$\log_{10}($# of reads$)$',
            title='Read Distributions per Treatment')
        plt.legend(handles=patches)
        
        return ax

    @staticmethod
    def _interpolate_values(function, start, end, include_end=False):
        """Perform linear interpolation."""
        new_x = np.linspace(start, end, endpoint=include_end)
        new_y = function(new_x)
        return(new_x, new_y)

    def create_cdf_plot(self):
        """
        Visualize percentage of cells falling into "bad", "ugly" and "good"
        categories. 
        """
        palette = sns.color_palette('bright')
        red_yellow_green = [palette[2], palette[4], palette[1]]
        # red_yellow_green = ['red', 'yellow', 'green']
        reads = self.read_df['post.filter.total.log10']
        ecdf = sm.distributions.ECDF(reads)
        x = np.linspace(0, max(reads), len(reads))
        y = ecdf(x)

        # reads below "bad" threshold
        bad = np.where(x < self.bad_threshold)[0]
        percent_bad = y[bad[-1]]

        # reads below "ugly" threshold, but above "bad"
        ugly = np.where((x >= self.bad_threshold) &
                        (x < self.ugly_threshold))[0]
        percent_ugly = y[ugly[-1]]

        # reads above "ugly" threshold
        good = np.where(x >= self.ugly_threshold)[0]

        estimate_cdf = interp1d(x, y)
        fig, ax = plt.subplots(1,1)
        for i, each in enumerate([bad, ugly, good]):
            plot_x = x[each]
            plot_y = y[each]
            # bad read coverage, estimate from data to bad threshold
            if i == 0:
                end_x, end_y = self._interpolate_values(estimate_cdf,
                                                        np.max(plot_x),
                                                        self.bad_threshold)
                plot_x = np.hstack((plot_x, end_x))
                plot_y = np.hstack((plot_y, end_y))

            # ugly read coverage, estimate from bad and to ugly thresholds
            elif i == 1:
                start_x, start_y = self._interpolate_values(estimate_cdf,
                                                            self.bad_threshold,
                                                            np.min(plot_x))
                end_x, end_y = self._interpolate_values(estimate_cdf,
                                                        np.max(plot_x),
                                                        self.ugly_threshold)
                plot_x = np.hstack((start_x, plot_x, end_x))
                plot_y = np.hstack((start_y, plot_y, end_y))

            # good reads, estimate from ugly to data 
            else:
                start_x, start_y = self._interpolate_values(estimate_cdf,
                                                            self.ugly_threshold,
                                                            np.min(plot_x))
                plot_x = np.hstack((start_x, plot_x))
                plot_y = np.hstack((start_y, plot_y))

            # draw lines to bad and ugly thresholds 
            if i < 2:
                xmax = [self.bad_threshold, self.ugly_threshold][i] / x[-1]
                ax.axhline(y=np.max(plot_y), xmin=0, xmax=xmax,
                           color=red_yellow_green[i], linestyle='--')
                                                    
            plt.plot(plot_x, plot_y, color=red_yellow_green[i])
            ax.fill_between(plot_x, plot_y, facecolor=red_yellow_green[i],
                            interpolate=True, alpha=0.25)

        y_w_thresholds = sorted(np.hstack((np.arange(0, 1.2, 0.2), 
                                        np.array([percent_bad, percent_ugly]))))
        labels = ['{:.1f}%'.format(x*100) for x in y_w_thresholds]
        bad_loc = np.where(y_w_thresholds)[0][0]
        if y_w_thresholds[bad_loc + 1] != percent_ugly and\
        (percent_ugly - percent_bad) < 0.25:
            labels[bad_loc + 1] = ''
        
        plt.yticks(y_w_thresholds, labels)
        plt.xlim((0, np.max(reads)))
        plt.ylim((0, 1))
        ax.set(xlabel='$\log_{10}($# of reads$)$', ylabel='Percent of Cells',
               title='Cumulative Distribution of Reads per Cell')
        
        return ax