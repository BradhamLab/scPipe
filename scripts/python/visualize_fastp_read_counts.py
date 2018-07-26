import os
import json

import re

import numpy as np
from scipy.interpolate import interp1d
import statsmodels.api as sm
import pandas as pd

import seaborn as sns
from matplotlib import pyplot as plt
import matplotlib.patches as mpatches

def generate_dataframe(output_dir):
    all_samples = {}
    for path, dirs, files in os.walk(output_dir):
        for each in files:
            if each == 'fastp.json':
                sample_name = os.path.basename(path)
                each = os.path.join(path, each)
                all_samples = {**all_samples, **get_read_data(each, sample_name)}
    df = pd.DataFrame(all_samples).T
    df['pre.filter.total.log10'] = [np.log10(x) for x in df['pre.filter.total']]
    df['post.filter.total.log10'] = [np.log10(x) for x in df['post.filter.total']]

def treatment_from_sample_name(sample_name):
    asw = re.compile('^ASW')
    chlorate = re.compile('^Chlorate')
    dmso = re.compile('^DMSO')
    mk886 = re.compile('^MK886')

    for pattern in [asw, chlorate, dmso, mk886]:
        pattern_match = re.search(pattern, sample_name)
        if pattern_match is not None:
            return pattern_match.group()
    
    return None


def get_read_data(fastp_json, sample_name):
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

    treatment = treatment_from_sample_name(sample_name)
    key_values = zip(['pre.filter.total', 'pre.filter.length',
                      'post.filter.total', 'post.filter.length',
                      'treatment'],
                      [r_tot, r_len, fr_tot, fr_len, treatment])
    sample_dict[sample_name] = {x:y for x, y in key_values}
    return sample_dict

def create_violin_plot(read_df): 
    sns.set(style='whitegrid')
    n_values = read_df['treatment'].value_counts()
    colors = sns.color_palette()
    treatments = sorted(n_values.index.values)
    patches = []
    for i, each in enumerate(treatments):
        color = colors[i % len(treatments)]
        label = "{} (n={})".format(each, n_values[treatments[i]])
        patches.append(mpatches.Patch(color=color, label=label))
    ax = sns.violinplot(data=read_df, x='treatment', y='post.filter.total.log10')
    ax.set(xlabel='Treatment', ylabel='$\log_{10}($# of reads$)$',
           title='Read Distributions per Treatment')
    plt.legend(handles=patches)
    
    return ax

def interpolate_values(function, start, end, include_end=False):
    new_x = np.linspace(start, end, endpoint=include_end)
    new_y = function(new_x)
    return(new_x, new_y)

def create_cdf_plot(read_df, bad_threshold, ugly_threshold):
    pastels = sns.color_palette('bright')
    red_yellow_green = [pastels[2], pastels[4], pastels[1]]
    # red_yellow_green = ['red', 'yellow', 'green']
    reads = read_df['post.filter.total.log10']
    ecdf = sm.distributions.ECDF(reads)
    x = np.linspace(min(reads), max(reads), len(reads))
    y = ecdf(x)

    # reads below "bad" threshold
    bad = np.where(x < bad_threshold)[0]
    percent_bad = y[bad[-1]]

    # reads below "ugly" threshold, but above "bad"
    ugly = np.where((x >= bad_threshold) & (x < ugly_threshold))[0]
    percent_ugly = y[ugly[-1]]

    # reads above "ugly" threshold
    good = np.where(x >= ugly_threshold)[0]

    estimate_cdf = interp1d(x, y)
    fig, ax = plt.subplots(1,1)
    for i, each in enumerate([bad, ugly, good]):
        plot_x = x[each]
        plot_y = y[each]
        # bad read coverage, estimate from data to bad threshold
        if i == 0:
            end_x, end_y = interpolate_values(estimate_cdf,
                                              np.max(plot_x), bad_threshold)
            plot_x = np.hstack((plot_x, end_x))
            plot_y = np.hstack((plot_y, end_y))

        # ugly read coverage, estimate from bad and to ugly thresholds
        elif i == 1:
            start_x, start_y = interpolate_values(estimate_cdf,
                                                  bad_threshold,
                                                  np.min(plot_x))
            end_x, end_y = interpolate_values(estimate_cdf,
                                              np.max(plot_x), ugly_threshold)
            plot_x = np.hstack((start_x, plot_x, end_x))
            plot_y = np.hstack((start_y, plot_y, end_y))

        # good reads, estimate from ugly to data 
        else:
            start_x, start_y = interpolate_values(estimate_cdf,
                                                  ugly_threshold,
                                                  np.min(plot_x))
            plot_x = np.hstack((start_x, plot_x))
            plot_y = np.hstack((start_y, plot_y))

        # draw lines to bad and ugly thresholds 
        if i < 2:
            xmax = [bad_threshold, ugly_threshold][i] / x[-1]
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
           title='CDF of $\log_{10}($# of reads$)$')
    plt.show()