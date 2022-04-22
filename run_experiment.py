# core python modules:
import os
import sys
from functools import partial
from contextlib import contextmanager
# 3rd party modules:
import multiprocessing as mp
import shutil as sh
import glob
import pickle
# local modules:
import fourdenvar
import experiment_setup as es
import observations
import run_jules as rjda


@contextmanager
def poolcontext(*args, **kwargs):
    """
    Function to control the parallel run of other functions
    :param args:
    :param kwargs:
    :return:
    """
    pool = mp.Pool(*args, **kwargs)
    yield pool
    pool.terminate()


def ens_member_run(ens_number_xi, seed_val=0, params=None, ens_dir_out=None):
    """
    Function to run a prior or posterior ensemble member
    :param ens_number_xi: tuple of ensemble member number (int) and corresponding parameter vector (arr)
    :param seed_val: seed value used for any perturbations within experiment (int)
    :param params: parameters to update in ensemble member run (lst)
    :param ens_dir_out: output directory for this ensemble (str)
    :return: string confirming ensemble member has run (str)
    """
    ens_number = ens_number_xi[0]
    xi = ens_number_xi[1]
    nml_dir = os.path.join(es.work_directory, "output_seed%d_xb_ens%d" % (seed_val, ens_number))
    if not os.path.exists(nml_dir):
        os.makedirs(nml_dir)
    for file in glob.glob(es.nml_directory + '/*.nml'):
        sh.copy(file, nml_dir)
    try:
        rj = rjda.RunJulesDa(params=params, values=xi, nml_dir=nml_dir)
        rj.run_jules_dic(output_name='ens'+str(ens_number), out_dir=ens_dir_out)
        dump_file_pattern = os.path.join(ens_dir_out, "ens%d.dump*" % ens_number)
        dumps = glob.glob(dump_file_pattern)
        for f in dumps:
            os.remove(f)
        sh.rmtree(nml_dir)
    except ValueError:
        sh.rmtree(nml_dir)
        print 'Something went wrong at: ' + str(ens_number)
    return 'ensemble member '+ str(ens_number) + ' run!'


def ens_run(x_ens, seed_val=0, params=None, ens_dir_out=None):
    """
    Perform a parallel run of JULES models given an ensemble of paramter vectors
    :param x_ens: ensemble of paramter vectors (arr)
    :param seed_val: seed value used for any perturbations in the experiment (int)
    :param xa: switch if this is a prior or posterior ensemble run (bool)
    :param ens_dir_out: directory name for ensemble output (str)
    :return: string confirming if the ensemble has been run (str)
    """
    print 'Running ensemble'

    if not os.path.exists(ens_dir_out):
        os.makedirs(ens_dir_out)

    mp.freeze_support()
    with poolcontext(processes=es.num_processes) as pool:
        res = pool.map(partial(ens_member_run, seed_val=seed_val, params=params, ens_dir_out=ens_dir_out), enumerate(x_ens))
    pool.close()
    pool.join()
    return 'Ensemble has been run'


if __name__ == "__main__":

    seed_val = int(sys.argv[1])

    # Output for a particular run is put into this directory structure:
    # 
    # OUT_DIR/work/
    # OUT_DIR/seed/background/xb.daily.nc
    # OUT_DIR/seed/ensemble_xb/ensNN.daily.nc
    # OUT_DIR/seed/ensemble_xa/ensNN.daily.nc
    # OUT_DIR/seed/xa_ensMM.pkl
    # OUT_DIR/seed/plot/

    out_dir_run = os.path.join(es.output_directory, "%d" % seed_val)

    out_dir_background = os.path.join(out_dir_run, "background")
    out_dir_xb = os.path.join(out_dir_run, "ensemble_xb")
    out_dir_xa = os.path.join(out_dir_run, "ensemble_xa")
    out_dir_plot = os.path.join(out_dir_run, "plot")

    file_posterior = os.path.join(out_dir_run, "xa_ens%d.pkl" % es.ensemble_size)

    # Override the loader functions from experiment_setup now that we know the
    # seed value.
    es.jules_hxb = partial(observations.extract_jules_hxb, out_dir=out_dir_background)
    es.jules_hxb_ens = partial(observations.extract_jules_hxb_ens, out_dir=out_dir_xb)

    # instantiate JULES data assimilation class
    jda = fourdenvar.FourDEnVar(seed_val=seed_val)
    params = jda.p_keys

    # if 'run_xb' is in system arguments then run JULES with prior parameters
    if 'run_xb' in sys.argv:

        if not os.path.exists(out_dir_background):
            os.makedirs(out_dir_background)

        nml_dir = os.path.join(es.work_directory, "output_seed%d_xb" % seed_val)
        if not os.path.exists(nml_dir):
            os.makedirs(nml_dir)
        for file in glob.glob(es.nml_directory + '/*.nml'):
            sh.copy(file, nml_dir)

        rj = rjda.RunJulesDa(params=params, values=jda.xb, nml_dir=nml_dir)
        rj.run_jules_dic(output_name='xb', out_dir=out_dir_background)

        sh.rmtree(nml_dir)
    # remove any old output in folders
    old_outs = glob.glob(os.path.join(out_dir_xb, '*.nc'))
    for f in old_outs:
        os.remove(f)
    # run prior ensemble
    ens_run(jda.xbs, seed_val=jda.seed_val, params=params, ens_dir_out=out_dir_xb)
    # if 'run_xa' is in system arguments then run posterior ensemble
    if 'run_xa' in sys.argv:
        jda = fourdenvar.FourDEnVar(assim=True, seed_val=seed_val)
        params = jda.p_keys
        # find posterior estimate and posterior ensemble
        xa = jda.find_min_ens_inc()
        xa_ens = jda.a_ens(xa[1])
        # pickle posterior parameter ensemble array
        with open(file_posterior, "wb") as f:
            pickle.dump(xa_ens, f)
        # remove any old output in folders
        old_outs = glob.glob(os.path.join(out_dir_xa, '*.nc'))
        for f in old_outs:
            os.remove(f)
        # run posterior ensemble
        ens_run(xa_ens, seed_val=jda.seed_val, params=params, ens_dir_out=out_dir_xa)
    if 'plot' in sys.argv:
        es.save_plots(file_posterior, out_dir_xa, out_dir_xb, out_dir_plot)
    print 'Experiment has been run'
