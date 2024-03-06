#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  6 09:12:35 2024

@author: tom.earnest
"""

import os

from experiments import experiment_test_all_atn_predictors

def main(rerun=False, replot=True):
    
    # setup output
    experiment_name = os.path.splitext(os.path.basename(__file__))[0]
    this_file = os.path.abspath(__file__)
    output_folder = os.path.abspath(os.path.join(this_file, '..', '..', 'outputs', experiment_name))
    
    if not os.path.isdir(output_folder):
        os.mkdir(output_folder)
    
    
if __name__ == '__main__':
    main()