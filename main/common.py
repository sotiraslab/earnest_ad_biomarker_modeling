#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 12 15:02:45 2024

@author: tom.earnest
"""

import argparse
import os
import pickle
import warnings

import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
import pandas as pd

def atn_subscripts(a='bin', t='bin', n='bin', upper=True):
    
    def proc(x, upper=upper):
        if isinstance(x, int):
            val = {0: None, 1: 'bin', 2: 'cat', 3: 'con', 4: 'svm'}[x]
        else:
            val = x
            
        val = val.upper() if val is not None and upper else val
        return val
        
    a = proc(a)
    t = proc(t)
    n = proc(n)
    
    alab = 'A$_{\\rm %s}$' % a if a is not None else ' - '
    tlab = 'T$_{\\rm %s}$' % t if t is not None else ' - '
    nlab = 'N$_{\\rm %s}$' % n if n is not None else ' - '
    included = [a is not None, t is not None, n is not None]
    
    if sum(included) == 1:
        idx = included.index(True)
        full = [alab, tlab, nlab][idx]
    else:
        full = alab + '/' + tlab + '/' + nlab
    return full

def contains_results(experiment):
    return os.path.exists(os.path.join(get_outputs_path(), experiment, 'results.csv'))

def get_outputs_path():
    this_dir = os.path.dirname(os.path.abspath(__file__))
    outputs_dir = os.path.join(this_dir, 'outputs')
    return outputs_dir

def load_results(experiment, result, check_short=True):
    outputs_dir = get_outputs_path()
    path = os.path.join(outputs_dir, experiment, result)
    if not os.path.exists(path) and not check_short:
        raise RuntimeError(f'"{result}" for "{experiment}" are missing.')

    if not os.path.exists(path):
        short_folder = experiment + '_short'
        short_path = os.path.join(outputs_dir, short_folder, result)

        if os.path.exists(short_path):
            warnings.warn(f'The full results are missing for "{experiment}", but the short results are present.  Using the latter!',
                          RuntimeWarning)
            path = short_path
        else:
            raise RuntimeError(f'"{result}" for "{experiment}" are missing.')

    ext = os.path.splitext(path)[1].lower().strip('.')
    if ext == 'pickle':
        with open(path, 'rb') as f:
            output = pickle.load(f)
    elif ext == 'csv':
        output = pd.read_csv(path)
    else:
        raise RuntimeError(f'Cannot load result with extension "{ext}"')

    return output

def locate_outfolder(experiment):
    outputs = os.listdir(get_outputs_path())
    long_outputs = [s for s in outputs if '_short' not in s]
    search1 = [s for s in long_outputs if experiment in s and contains_results(s)]
    search2 = [s for s in outputs if experiment in s and contains_results(s)]
    if search1:
        return os.path.join(get_outputs_path(), search1[0])
    elif search2:
        warnings.warn(f'Only found short results for "{experiment}".  Using those...')
        return os.path.join(get_outputs_path(), search2[0])
    else:
        raise ValueError(f'Cannot find any results for "{experiment}".')

def parse():
    parser = argparse.ArgumentParser()
    parser.add_argument('--short', action='store_true')
    return parser.parse_args()

def set_font_properties(arial_font=None):

    if arial_font is None:
        plt.rcParams.update({'font.family': 'arial'})
    else:
        try:
            font_prop = fm.FontProperties(fname=arial_font)
            plt.rcParams.update({
                'font.family': font_prop.get_name()})
        except:
            pass
    plt.rcParams.update({'font.size': 14})
    
def set_labels_baseline_exp(fig):
    ax = fig.axes[0]
    f = atn_subscripts
    labs = [
        'Baseline',
        f(1, 0, 0),
        f(0, 1, 0),
        f(0, 0, 1),
        f(1, 1, 1),
        f(2, 0, 0),
        f(0, 2, 0),
        f(0, 0, 2),
        f(2, 2, 2),
        f(3, 0, 0),
        f(0, 3, 0),
        f(0, 0, 3),
        f(3, 3, 3),
        f(4, 0, 0),
        f(0, 4, 0),
        f(0, 0, 4),
        f(4, 4, 4),]
    ax.set_yticklabels(labs)
    
def set_labels_binary_exp(fig):
    ax = fig.axes[0]
    f = atn_subscripts
    labs = [
        f(1, 1, 1),
        f(2, 1, 1),
        f(1, 2, 1),
        f(1, 1, 2),
        f(2, 2, 2),
        f(3, 1, 1),
        f(1, 3, 1),
        f(1, 1, 3),
        f(3, 3, 3),
        f(4, 0, 0),
        f(0, 4, 0),
        f(0, 0, 4),
        f(4, 4, 4),]
    ax.set_yticklabels(labs)
    
def set_labels_combo_vs_csf(fig):
    ax = fig.axes[0]
    f = atn_subscripts
    labs = [
        'Baseline',
        f(1, 1, 1) + ' [Imaging]',
        f(1, 1, 1) + ' [CSF]',
        f(2, 2, 2) + ' [Imaging]',
        f(2, 2, 2) + ' [CSF]',
        f(3, 3, 3) + ' [Imaging]',
        f(3, 3, 3) + ' [CSF]']
    ax.set_yticklabels(labs)
    
def setup_output(call_file, short=False):
    outputs_dir = get_outputs_path()
    foldername = os.path.splitext(os.path.basename(call_file))[0]
    if short:
        foldername += '_short'
    requested_dir = os.path.join(outputs_dir, foldername)
    if not os.path.isdir(requested_dir):
        os.mkdir(requested_dir)
    return requested_dir