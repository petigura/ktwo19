#!/usr/bin/env python
from argparse import ArgumentParser
import glob
import os
from matplotlib.pylab import *
import cPickle as pickle
from collections import OrderedDict


import ktwo19.io
#import ktwo19.tables
#import ktwo19.values
#import ktwo19.plotting._corner
#import ktwo19.plotting.omc

def main():
    psr = ArgumentParser()
    subpsr = psr.add_subparsers(title="subcommands", dest='subcommand')
    psr_parent = ArgumentParser(add_help=False)

    psr2 = subpsr.add_parser('photodyn-prepare', parents=[psr_parent], )
    psr2.set_defaults(func=photdyn_prepare)

    psr2 = subpsr.add_parser('create-val', parents=[psr_parent], )
    psr2.add_argument('name',type=str)
    psr2.set_defaults(func=create_val)

    psr2 = subpsr.add_parser('create-plot', parents=[psr_parent], )
    psr2.add_argument('name',type=str)
    psr2.set_defaults(func=create_plot)

    psr2 = subpsr.add_parser('create-table', parents=[psr_parent], )
    psr2.add_argument('name',type=str)
    psr2.set_defaults(func=create_table)

    psr2 = subpsr.add_parser('update-paper', parents=[psr_parent])
    psr2.set_defaults(func=update_paper)

    args = psr.parse_args()
    args.func(args)

def photdyn_prepare(args):
    df = ktwo19.io.load_table('ktwo-everest',cache=1)
    fn = 'analysis/photodyn/photometry-ktwo.tsv'
    df.to_csv(fn,sep='\t',index=None)
    print "created {}".format(fn)

    df = ktwo19.io.load_table('rv',cache=2)
    fn = 'analysis/photodyn/rv.tsv'
    df.to_csv(fn,sep='\t',index=None)
    print "created {}".format(fn)

    df = ktwo19.io.load_table('rv-trend-removed',cache=2)
    fn = 'analysis/photodyn/rv-trend-removed.tsv'
    df.to_csv(fn,sep='\t',index=None)
    print "created {}".format(fn)

def create_table(args):
    w = Workflow()
    w.create_file('table', args.name ) 

def create_plot(args):
    w = Workflow()
    w.create_file('plot', args.name ) 

def create_val(args):
    w = Workflow()
    w.create_file('val',args.name) 

def update_paper(args):
    w = Workflow()
    w.update_paper() 

class Workflow(object):
    def __init__(self):
        # plots
        d = OrderedDict()
        d['dumbbell'] = ktwo24.plotting.dumbbell.fig_dumbbell
        d['lithwick-corner4'] = ktwo24.plotting._corner.fig_lithwick4
        d['ttvfast-corner6'] = ktwo24.plotting._corner.fig_ttvfast6
        d['lithwick-samples'] = ktwo24.plotting.omc.fig_lithwick_samples
        d['ttvfast-samples'] = ktwo24.plotting.omc.fig_ttvfast_samples
        self.plot_dict = d

        d = OrderedDict()
        d['rv'] = ktwo24.tables.tab_rv
        d['rv-stub'] = lambda : ktwo24.tables.tab_rv()[:10]
        d['transit-times-predict'] = ktwo24.tables.tab_transit_times_predict
        d['transit-times-predict-stub1'] = \
                lambda : ktwo24.tables.tab_transit_times_predict()[:5]
        d['transit-times-predict-stub2'] = \
                lambda : ktwo24.tables.tab_transit_times_predict()[-5:]
        self.table_dict = d

        d = OrderedDict()
        d['stat'] = ktwo24.values.val_stat
        d['fit'] = ktwo24.values.val_fit
        d['lithwick'] = ktwo24.values.val_lithwick
        d['ttvfast'] = ktwo24.values.val_ttvfast
        self.val_dict = d

        d = OrderedDict()
        d['table'] = self.table_dict
        d['plot'] = self.plot_dict
        d['val'] = self.val_dict
        self.all_dict = d

    def key2fn(self, key, kind):
        if kind=='plot':
            return 'fig_'+key+'.pdf'
        if kind=='table':
            return 'tab_'+key+'.tex'
        if kind=='val':
            return 'val_'+key+'.tex'
            
    def create_file(self, kind, name):
        i = 0

        for key, func in self.all_dict[kind].iteritems():
            if kind=='plot':
                if name=='all':
                    func()
                elif name==key:
                    func()
                else:
                    continue

                fn = self.key2fn(key, 'plot')
                gcf().savefig(fn)

            elif kind=='table':
                if name=='all':
                    lines = func()
                elif name==key:
                    lines = func()
                else:
                    continue
                    
                # Remove last \\
                fn = self.key2fn(key, 'table')
                with open(fn,'w') as f:
                    f.writelines("\n".join(lines))

            elif kind=='val':
                fn = self.key2fn(key, 'val')
                if name=='all':
                    lines = func()
                elif name==key:
                    lines = func()
                else:
                    continue

                lines1 = [
                    "\\newcommand{\%s}[1]{%%" % key,
                    "\IfEqCase{#1}{%",
                ]

                lines2 = [
                    "}[\PackageError{tree}{Undefined option to tree: #1}{}]%",
                    "}%"
                ]
                lines = lines1 + lines + lines2

                with open(fn,'w') as f:
                    f.writelines("\n".join(lines))

            i+=1

        if i==0:
            assert False, name + " not a valid key"

    def update_paper(self):
        for kind, d in self.all_dict.iteritems():
            for key, val in d.iteritems():
                fn = self.key2fn(key, kind)
                cmd = 'cp {} paper/'.format(fn)
                print cmd
                os.system(cmd)

if __name__=="__main__":
    main()

