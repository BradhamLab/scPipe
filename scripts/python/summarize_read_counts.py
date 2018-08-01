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

# snakemake import
from snakemake.utils import report

class SummarizeFastpReads(object):
    """
    Class to summarize output from read quality control using `fastp`.
    """
    
    def __init__(self, fastp_dir, output_dir, treatment_regex, bad_threshold,
                 ugly_threshold):
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
            fastp_dir (string): path to `fastp` output directory containing
                `fastp` output. Directory is assumed to use the following
                format:

                    <fastp_dir>/<sample_id/>

                Where a file, `fastp.json` exists within each <sample_id>
                subdirectory.
            output_dir (string): directory where output should be written to. 
            treatment_regex (list, string): list of regex patterns to extract
                all possible treatments from a sample name.
            bad_threshold (int): maximum number of reads for a cell to be
                considered as 'bad'.
            ugly_threshold (int): maximum number of reads for a cell to be
                considered as 'ugly'.
        """
        self.fastp_dir = fastp_dir
        self.output_dir = output_dir
        self.treatment_regex = [re.compile(x) for x in treatment_regex]
        self.bad_threshold = bad_threshold
        self.log10_bad = np.log10(bad_threshold)
        self.ugly_threshold = ugly_threshold
        self.log10_ugly = np.log10(ugly_threshold)
        self.read_df = self.generate_dataframe()
        self.plot_df = self.read_df.replace([np.inf, -np.inf], np.nan)
        self.plot_df.dropna(inplace=True)

    def snakemake_report(self):
        """
        Generate a SnakeMake report.

        Generates an html report containing diagnostic plots for read coverage
        quality. All images are saved to a `plots` subdirecty in the specified
        output directory from `self.output_dir`. Finally, a `.csv` file
        containing read summaries is written to `self.output_dir`. 
        """

        # check if output dir exists
        if not os.path.exists(self.output_dir):
            os.mkdir(self.output_dir)

        # check if plot dir exists
        if not os.path.exists('plots'):
            os.mkdir('plots')

        # output paths
        csv_path = os.path.join(self.output_dir, 'read_summary.csv')

        # report loc
        html_path = os.path.join(self.output_dir, 'report.html')

        # plot loc
        plot_dir = os.path.join(self.output_dir, 'plots')
        if not os.path.exists(plot_dir):
            os.makedirs(plot_dir)

        violin_path = os.path.join(plot_dir, 'violin.png')
        cdf_path = os.path.join(plot_dir, 'cdf.png')
        stacked_bar_path = os.path.join(plot_dir, 'stacked_barplot.png')

        # write csv file
        self.read_df.to_csv(csv_path)


        # clear any pyplot figures just in case
        plt.cla()

        # create plots
        self.create_violin_plot()
        plt.savefig(violin_path)
        plt.cla()

        self.create_cdf_plot()
        plt.savefig(cdf_path)
        plt.cla()

        self.create_stack_barplot()
        plt.savefig(stacked_bar_path)
        plt.cla()
  
        rst_markup = """
        ============================================
        Read Summary Following Fastp Quality Control
        ============================================

        Read Counts per Treatment
        =========================
        .. image:: plots/violin.png

        Read Coverage Quality
        =====================
        .. image:: plots/cdf.png

        Coverage Quality per Treatment
        ==============================
        .. image:: plots/stacked_barplot.png
        """.format(violin_path, cdf_path, stacked_bar_path)
        report(rst_markup, html_path,
               metadata="Author: Dakota Hawkins (dyh0110@bu.edu)")

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
        for path, dirs, files in os.walk(self.fastp_dir):
            for each in files:
                if each == 'fastp.json':
                    sample_name = os.path.basename(path)
                    each = os.path.join(path, each)
                    all_samples = {**all_samples,
                                   **self.get_read_data(each, sample_name)}
        df = pd.DataFrame(all_samples).T
        df.to_csv('test.csv')
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

        for pattern in self.treatment_regex:
            pattern_match = re.search(pattern, sample_name)
            if pattern_match is not None:
                return pattern_match.group()
        
        return None

    def _get_good_bad_ugly(self, read_count):
        """
        Get good-bad-ugly group based on read count.
        """
        if read_count < self.bad_threshold:
            return 'Bad'
        elif read_count < self.ugly_threshold:
            return 'Ugly'
        return 'Good'

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
                6. quality | quality of read coverage defined by threshold
                    values. Can be either 'Good', 'Bad' or 'Ugly'.
                7. pre.filter.total.log10 | log10 of total reads pre filtering.
                8. post.filter.total.log10 | log10 of total reads post filtering.
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
        gbu_group = self._get_good_bad_ugly(fr_tot)

        # avoid NaNs/infinites for plotting

        treatment = self.treatment_from_sample_name(sample_name)
        key_values = zip(['pre.filter.total', 'pre.filter.length',
                          'post.filter.total', 'post.filter.length',
                          'treatment', 'quality'],
                         [r_tot, r_len, fr_tot, fr_len, treatment, gbu_group])
        sample_dict[sample_name] = {x:y for x, y in key_values}
        return sample_dict

    def create_violin_plot(self):
        """
        Visualize post qc log10 read distribution between treatment types.
        """ 
        sns.set(style='whitegrid')
        n_values = self.plot_df['treatment'].value_counts()
        colors = sns.color_palette()
        treatments = n_values.index.values
        patches = []
        for i, each in enumerate(treatments):
            color = colors[i % len(treatments)]
            label = "{} (n={})".format(each, n_values[treatments[i]])
            patches.append(mpatches.Patch(color=color, label=label))
        
        ax = sns.violinplot(data=self.plot_df, x='treatment',
                            y='post.filter.total.log10')
        ax.set(xlabel='Treatment', ylabel='$\log_{10}($# of reads$)$',
               title='Read Distributions per Treatment')
        plt.legend(handles=patches, loc='upper center',
                   bbox_to_anchor=(0.5, -0.1), fancybox=False,
                   ncol=len(treatments), fontsize='x-small',
                   frameon=False, shadow=False)
        plt.gcf().subplots_adjust(bottom=0.15)
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
        red_yellow_green = ['#E06666', '#FFDC77', '#93C47D']
        reads = self.plot_df['post.filter.total.log10']
        ecdf = sm.distributions.ECDF(reads)
        x = np.linspace(0, max(reads), len(reads))
        y = ecdf(x)

        # reads below "bad" threshold
        bad = np.where(x < self.log10_bad)[0]
        percent_bad = 0
        if len(bad) > 0:
            percent_bad = y[bad[-1]]

        # reads below "ugly" threshold, but above "bad"
        ugly = np.where((x >= self.log10_bad) &
                        (x < self.log10_ugly))[0]
        percent_ugly = 0
        if len(ugly) > 0: 
            percent_ugly = y[ugly[-1]]

        # reads above "ugly" threshold
        good = np.where(x >= self.log10_ugly)[0]

        estimate_cdf = interp1d(x, y)
        fig, ax = plt.subplots(1,1)
        for i, each in enumerate([bad, ugly, good]):
            plot_x = x[each]
            plot_y = y[each]

            # check to see if there's anything to plot:
            if len(plot_x) > 0:

                # bad read coverage, estimate from data to bad threshold
                if i == 0:
                    end_x, end_y = self._interpolate_values(estimate_cdf,
                                                            np.max(plot_x),
                                                            self.log10_bad)
                    plot_x = np.hstack((plot_x, end_x))
                    plot_y = np.hstack((plot_y, end_y))

                # ugly read coverage, estimate from bad and to ugly thresholds
                elif i == 1:
                    start_x, start_y = self._interpolate_values(estimate_cdf,
                                                                self.log10_bad,
                                                                np.min(plot_x))
                    end_x, end_y = self._interpolate_values(estimate_cdf,
                                                            np.max(plot_x),
                                                            self.log10_ugly)
                    plot_x = np.hstack((start_x, plot_x, end_x))
                    plot_y = np.hstack((start_y, plot_y, end_y))

                # good reads, estimate from ugly to data 
                else:
                    start_x, start_y = self._interpolate_values(estimate_cdf,
                                                                self.log10_ugly,
                                                                np.min(plot_x))
                    plot_x = np.hstack((start_x, plot_x))
                    plot_y = np.hstack((start_y, plot_y))

                # draw lines to bad and ugly thresholds 
                if i < 2:
                    xmax = [self.log10_bad, self.log10_ugly][i] / x[-1]
                    ax.axhline(y=np.max(plot_y), xmin=0, xmax=xmax,
                            color=red_yellow_green[i], linestyle='--')
                                                        
                plt.plot(plot_x, plot_y, color=red_yellow_green[i])
                ax.fill_between(plot_x, plot_y, facecolor=red_yellow_green[i],
                                interpolate=True, alpha=0.5)

        # format y-axis tick marsk to include threshold percents.
        y_w_thresholds = sorted(np.hstack((np.arange(0, 1.2, 0.2), 
                                        np.array([percent_bad, percent_ugly]))))
        labels = ['{:.1f}%'.format(x*100) for x in y_w_thresholds]
        bad_loc = np.where(y_w_thresholds)[0][0]
        if y_w_thresholds[bad_loc + 1] != percent_ugly and\
        (percent_ugly - percent_bad) < 0.25:
            labels[bad_loc + 1] = ''

        # create legend
        cumulative_percentages = [0, percent_bad, percent_ugly, 1]
        patches = []
        for i, each in enumerate(['Bad', 'Ugly', 'Good']):
            color = red_yellow_green[i]
            label = "{} (%={:.1f})".format(each, (cumulative_percentages[i + 1]
                                           - cumulative_percentages[i])
                                           * 100)
            patches.append(mpatches.Patch(color=color, label=label))
        
        plt.yticks(y_w_thresholds, labels)
        plt.xlim((0, np.max(reads)))
        plt.ylim((0, 1))
        ax.set(xlabel='$\log_{10}($# of reads$)$', ylabel='Percent of Cells',
               title='Cumulative Distribution of Reads per Cell')
        return ax

    def create_stack_barplot(self):
        """
        Visualize percentage of cells falling into 'Good', 'Bad', and 'Ugly'
        catergories for each treatment group.
        """
        color_dict = {x:y for x, y in zip(['Bad', 'Ugly', 'Good'],\
                                          ['#E06666', '#FFDC77', '#93C47D'])}
        by_treatment = self.plot_df.groupby('treatment')
        quality_by_treatment = by_treatment['quality'].value_counts().unstack()
        percentages = quality_by_treatment.apply(lambda x: x / sum(x), axis=1)
        
        treatments = percentages.index.values
        qualities = quality_by_treatment.columns.values
        width = 0.35
        ind = np.arange(len(treatments))
        plots = []
        fig, ax = plt.subplots(1,1)
        
        for i, each in enumerate(qualities):
            if i == 0:
                plots.append(ax.bar(ind, percentages[each],
                                    color=color_dict[each], width=width))
            else:
                cum_prob = sum([percentages[x] for x in qualities[0:i]])
                plots.append(ax.bar(ind, percentages[each],
                                     bottom=cum_prob,
                                     color=color_dict[each], width=width))

        plt.ylabel('Percentage')
        plt.title('Quality Percentages by Treatment')
        plt.xticks(ind, treatments)
        plt.yticks(np.arange(0, 1.2, 0.2),
                   ['{:.0f}%'.format(x*100) for x in np.arange(0, 1.2, 0.2)])
        plt.legend([x[0] for x in plots], qualities,
                   loc='upper center', bbox_to_anchor=(0.5, -0.05),
                   fancybox=False, ncol=len(treatments),
                   frameon=False, shadow=False)
        plt.ylim(0, 1)
        plt.gcf().subplots_adjust(bottom=0.15)
        return ax

if __name__ == "__main__":
    snakemake_exists = True
    try:
        snakemake
    except NameError:
        snakemake_exists = False
    
    if snakemake_exists:
        read_summary = SummarizeFastpReads(snakemake.params['fastp'],
                                           snakemake.params['outdir'],
                                           snakemake.params['regex'],
                                           snakemake.params['bad'],
                                           snakemake.params['ugly'])
        read_summary.snakemake_report()