import math
import functools
import itertools
import sys
import typing
from copy import copy

from scipy.stats import binom
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patheffects as patheffects
import plotly.graph_objects as go

from src.annotation import Annotation
from src.report.html_templates import float_to_str

# type ConfidenceType = tuple[float, float, float | str, float, float, float, float]

def combine_distribs(deletes, inserts):
    """
    Combine insert and delete models/distributions
    :param deletes: ndarray - delete distribution
    :param inserts: ndarray - insert distribution
    :return: ndarray - combined array of the same length
    """
    # how much to fill?
    to_fill = sum(deletes == 0.0) + 1
    while to_fill < len(inserts) and inserts[to_fill] > 0.0001:
        to_fill += 1

    # create the end array
    end_distr = np.zeros_like(deletes, dtype=float)

    # fill it!
    for i, a in enumerate(inserts[:to_fill]):
        end_distr[i:] += (deletes * a)[:len(deletes) - i]

    return end_distr


def const_rate(n, p1=0.0, p2=1.0, p3=1.0):
    """
    Constant rate function.
    :param n: int - allele number (unused)
    :param p1: float - constant parameter
    :param p2: float - linear parameter (unused)
    :param p3: float - additional parameter (unused)
    :return: float - p1
    """
    return p1


def linear_rate(n, p1=0.0, p2=1.0, p3=1.0):
    """
    Linear rate function.
    :param n: int - allele number
    :param p1: float - constant parameter
    :param p2: float - linear parameter
    :param p3: float - additional parameter (unused)
    :return: float - p1 + p2 * n
    """
    return p1 + p2 * n


def n2_rate(n, p1=0.0, p2=1.0, p3=1.0):
    """
    Quadratic rate function.
    :param n: int - allele number
    :param p1: float - constant parameter
    :param p2: float - linear parameter
    :param p3: float - quadratic parameter
    :return: float - p1 + p2 * n + p3 * n * n
    """
    return p1 + p2 * n + p3 * n * n


def exp_rate(n, p1=0.0, p2=1.0, p3=1.0):
    """
    Exponential rate function.
    :param n: int - allele number
    :param p1: float - constant parameter
    :param p2: float - linear parameter
    :param p3: float - exponential parameter
    :return: float - p1 + p2 * e^(p3 * n)
    """
    return p1 + p2 * math.exp(p3 * n)


def clip(value, minimal, maximal):
    """
    Clips value to range <minimal, maximal>
    :param value: ? - value
    :param minimal: ? - minimal value
    :param maximal: ? - maximal value
    :return: ? - clipped value
    """
    return min(max(minimal, value), maximal)


def model_full(rng, model_params, n, rate_func=linear_rate):
    """
    Create binomial model for both deletes and inserts of STRs
    :param rng: int - max_range of distribution
    :param model_params: 4-tuple - parameters for inserts and deletes
    :param n: int - target allele number
    :param rate_func: function - rate function for deletes
    :return: ndarray - combined distribution
    """
    p1, p2, p3, q = model_params
    deletes = binom.pmf(np.arange(rng), n, clip(1 - rate_func(n, p1, p2, p3), 0.0, 1.0))
    inserts = binom.pmf(np.arange(rng), n, q)
    return combine_distribs(deletes, inserts)


def model_template(rng, model_params, rate_func=linear_rate):
    """
    Partial function for model creation.
    :param rng: int - max_range of distribution
    :param model_params: 4-tuple - parameters for inserts and deletes
    :param rate_func: function - rate function for deletes
    :return: partial function with only 1 parameter - n - target allele number
    """
    return functools.partial(model_full, rng, model_params, rate_func=rate_func)


def generate_models(min_rep: int, max_rep: int, multiple_backgrounds: bool = True) -> typing.Iterator[typing.Tuple[int | str, int | str]]:
    """
    Generate all pairs of alleles (models for generation of reads).
    :param min_rep: int - minimal number of repetitions
    :param max_rep: int - maximal number of repetitions
    :param multiple_backgrounds: bool - whether to generate all background states
    :return: generator of allele pairs (numbers or 'E' or 'B')
    """
    for model_index1 in range(min_rep, max_rep):
        for model_index2 in range(model_index1, max_rep):
            yield model_index1, model_index2
        yield model_index1, 'E'
        if multiple_backgrounds:
            yield 'B', model_index1

    yield 'B', 'B'
    yield 'E', 'E'


def generate_models_one_allele(min_rep: int, max_rep: int) -> typing.Iterator[typing.Tuple[int | str, int | str]]:
    """
    Generate all pairs of alleles (models for generation of reads).
    :param min_rep: int - minimal number of repetitions
    :param max_rep: int - maximal number of repetitions
    :return: generator of allele pairs (numbers or 'E' or 'B'), 'X' for non-existing allele
    """
    for model_index1 in range(min_rep, max_rep):
        yield model_index1, 'X'

    yield 'B', 'X'
    yield 'E', 'X'


class Inference:
    """ Class for inference of alleles. """

    MIN_REPETITIONS = 1

    # default parameters for inference (miSeq default)
    DEFAULT_MODEL_PARAMS = (0.00716322, 0.000105087, 0.0210812, 0.0001648)
    DEFAULT_FIT_FUNCTION = 'linear'

    def __init__(self, read_distribution, params_file, str_rep=3, minl_primer1=5, minl_primer2=5, minl_str=5, p_bckg_closed=None,
                 p_bckg_open=None, p_expanded=None):
        """
        Initialization of the Inference class + setup of all models and their probabilities.
        :param read_distribution: ndarray(int) - read distribution
        :param params_file: str - filename of parameters, None for defaults
        :param str_rep: int - length of the STR
        :param minl_primer1: int - minimal length of the left primer
        :param minl_primer2: int - minimal length of the right primer
        :param minl_str: int - minimal length of the STR
        :param p_bckg_closed: float - probability of the background model for closed observation
        :param p_bckg_open: float - probability of the background model for open observation
        :param p_expanded: float - probability of the expanded model (if None it is equal to other models)
        """
        # assign variables
        self.str_rep = str_rep
        self.minl_primer1 = minl_primer1
        self.minl_primer2 = minl_primer2
        self.minl_str = minl_str
        self.read_distribution = read_distribution
        self.sum_reads_log = np.log(np.sum(read_distribution))
        self.sum_reads = np.sum(read_distribution)
        self.params_file = params_file
        self.p_expanded = p_expanded
        self.p_bckg_closed = p_bckg_closed
        self.p_bckg_open = p_bckg_open

    def construct_models(self, min_rep, max_rep, e_model):
        """
        Construct all models needed for current inference.
        :param min_rep: int - minimal allele to model
        :param max_rep: int - maximal allele to model
        :param e_model: int - model for expanded alleles
        :return: None
        """
        # extract params
        model_params, rate_func_str = self.read_params(self.params_file)
        str_to_func = {'linear': linear_rate, 'const': const_rate, 'exponential': exp_rate, 'square': n2_rate}
        rate_func = const_rate
        if rate_func_str in str_to_func.keys():
            rate_func = str_to_func[rate_func_str]

        # save min_rep and max_rep
        self.min_rep = min_rep
        self.max_rep = max_rep  # non-inclusive
        self.max_with_e = e_model + 1  # non-inclusive

        # get models
        mt = model_template(self.max_with_e, model_params, rate_func)
        self.background_model = np.concatenate([np.zeros(self.min_rep, dtype=float),
                                                np.ones(self.max_with_e - self.min_rep, dtype=float) / float(self.max_with_e - self.min_rep)])
        self.expanded_model = mt(self.max_with_e - 1)
        self.allele_models = {i: mt(i) for i in range(min_rep, max_rep)}
        self.models = {'E': self.expanded_model, 'B': self.background_model}
        self.models.update(self.allele_models)

        # get model likelihoods
        open_to_closed = 10.0

        l_others = 1.0
        l_bckg_open = 0.01
        l_exp = 1.01

        l_bckg_model_open = 1.0

        if self.p_expanded is None:
            self.p_expanded = l_exp
        if self.p_bckg_open is None and self.p_bckg_closed is None:
            self.p_bckg_open = l_bckg_open
            self.p_bckg_closed = self.p_bckg_open / open_to_closed
        if self.p_bckg_closed is None:
            self.p_bckg_closed = self.p_bckg_open / open_to_closed
        if self.p_bckg_open is None:
            self.p_bckg_open = self.p_bckg_closed * open_to_closed

        self.model_probabilities = {'E': self.p_expanded, 'B': l_bckg_model_open}
        self.model_probabilities.update({i: l_others for i in self.allele_models.keys()})

    def read_params(self, params_file):
        """
        Reads all parameters written with write_params(print_all=True)
        :param params_file: str - filename to read parameters from, if None, load default params
        :return: 4-tuple, 2-tuple, function - parameters for model, read count drop, and error function for model distributions
        """
        if params_file is None:
            return self.DEFAULT_MODEL_PARAMS, self.DEFAULT_FIT_FUNCTION

        # read 2nd and last line of the file
        with open(params_file) as f:
            lines = f.readlines()
            fit_function = lines[1].strip().split()[1]
            split = list(map(float, lines[-1].strip().split()))

        if len(split) < 4:
            print('ERROR: parameters were not read successfully, using defaults!', file=sys.stderr)
            return self.DEFAULT_MODEL_PARAMS, self.DEFAULT_FIT_FUNCTION

        # extract parameters from last line of file
        model_params = tuple(split[0:4])

        return model_params, fit_function

    def likelihood_rl(self, rl):
        """
        Likelihood of a read with this length.
        :param rl: int - read length
        :return: float - likelihood of a read this long
        """
        # print('rl', self.read_distribution[rl] / float(self.sum_reads))
        return self.read_distribution[rl] / float(self.sum_reads)

    @staticmethod
    def likelihood_model(model, g):
        """
        Likelihood of a generated allele al from a model of
        :param model: ndarray - model that we evaluate
        :param g: int - observed read count
        :return: float - likelihood of a read coming from this model
        """
        return model[g]

    @staticmethod
    def likelihood_intersection(model_i, model_j, g):
        return min(model_i[g], model_j[g])

    def likelihood_coverage(self, true_length, rl, closed=True):
        """
        Likelihood of generating a read with this length and this allele.
        :param true_length: int - true number of repetitions of an STR
        :param rl: int - read length
        :param closed: bool - if the read is closed - i.e. both primers are there
        :return: float - likelihood of a read being generated with this attributes
        """
        whole_inside_str = max(0, true_length * self.str_rep + self.minl_primer1 + self.minl_primer2 - rl + 1)
        # closed_overlapping = max(0, rl - self.minl_primer1 - self.minl_primer2 - true_length * self.str_rep + 1)
        open_overlapping = max(0, rl + true_length * self.str_rep - 2 * self.minl_str + 1)

        assert open_overlapping > whole_inside_str, '%d open %d whole inside %d %d %d' % (
            open_overlapping, whole_inside_str, true_length, rl, self.minl_str)

        return 1.0 / float(open_overlapping - whole_inside_str)

    def likelihood_read_allele(self, model, observed, rl, closed=True):
        """
        Likelihood of generation of read with observed allele count and rl.
        :param model: ndarray - model for the allele
        :param observed: int - observed allele count
        :param rl: int - read length
        :param closed: bool - if the read is closed - i.e. both primers are there
        :return:
        """
        if closed:
            return self.likelihood_rl(rl) * self.likelihood_model(model, observed) * self.likelihood_coverage(observed, rl, True)
        else:
            number_of_options = 0
            partial_likelihood = 0
            for true_length in itertools.chain(range(observed, self.max_rep), [self.max_with_e - 1]):
                partial_likelihood += self.likelihood_model(model, true_length) * self.likelihood_coverage(true_length, rl, False)
                number_of_options += 1

            return self.likelihood_rl(rl) * partial_likelihood / float(number_of_options)

    def likelihood_read_intersection(self, model_i, model_j, observed, rl, closed=True):
        """
        Likelihood of generation of read with observed allele count and rl.
        :param model_i: ndarray - model for the left allele
        :param model_j: ndarray - model for the right allele
        :param observed: int - observed allele count
        :param rl: int - read length
        :param closed: bool - if the read is closed - i.e. both primers are there
        :return:
        """
        if closed:
            return self.likelihood_rl(rl) * self.likelihood_intersection(model_i, model_j, observed) * self.likelihood_coverage(observed, rl, True)
        else:
            number_of_options = 0
            partial_likelihood = 0
            for true_length in itertools.chain(range(observed, self.max_rep), [self.max_with_e - 1]):
                partial_likelihood += self.likelihood_intersection(model_i, model_j, true_length) * self.likelihood_coverage(true_length, rl, False)
                number_of_options += 1

            return self.likelihood_rl(rl) * partial_likelihood / float(number_of_options)

    @functools.lru_cache()
    def likelihood_read(self, observed: int, rl: int, model_index1: int, model_index2: int = None, closed: bool = True) -> float:
        """
        Compute likelihood of generation of a read from either of those models.
        :param observed: int - observed allele count
        :param rl: int - read length
        :param model_index1: char/int - model index for left allele
        :param model_index2: char/int - model index for right allele or None if mono-allelic
        :param closed: bool - if the read is closed - i.e. both primers are there
        :return: float - likelihood of this read generation
        """
        # TODO: tuto podla mna nemoze byt len tak +, chyba tam korelacia modelov, ale v ramci zjednodusenia asi ok
        allele1_likelihood = self.model_probabilities[model_index1] * self.likelihood_read_allele(self.models[model_index1], observed, rl, closed)
        allele2_likelihood = 0.0 if model_index2 is None else self.model_probabilities[model_index2] * self.likelihood_read_allele(
            self.models[model_index2], observed, rl, closed)
        p_bckg = self.p_bckg_closed if closed else self.p_bckg_open
        bckgrnd_likelihood = p_bckg * self.likelihood_read_allele(self.models['B'], observed, rl, closed)

        # alleles_intersection = min(model_prob_j, model_prob_i) * self.likelihood_read_intersection(model_i, model_j, observed, rl, closed)
        # if alleles_intersection > 0.0:
        #    print('%g %g %g %s %s %d' % (alleles_intersection, allele2_likelihood, allele1_likelihood,
        #    str(model_index1), str(model_index2), observed))

        assert not np.isnan(allele2_likelihood)
        assert not np.isnan(allele1_likelihood)
        assert not np.isnan(bckgrnd_likelihood)
        # assert alleles_intersection <= max(allele1_likelihood, allele2_likelihood), '%g %g %g %s %s %d' % (
        #    alleles_intersection, allele2_likelihood, allele1_likelihood, str(model_index1), str(model_index2), observed)

        # print('read_%s' % (str(closed)), observed, 'all1_lh', allele1_likelihood, 'all2_lh', allele2_likelihood)

        return allele1_likelihood + allele2_likelihood + bckgrnd_likelihood  # - alleles_intersection

    def infer(self, annotations: list[Annotation], filt_annotations: list[Annotation], index_rep: int,
              verbose: bool = True, monoallelic: bool = False) -> dict[tuple[int | str, int | str]: float]:
        """
        Does all the inference, computes for which 2 combination of alleles are these annotations and parameters the best.
        argmax_{G1, G2} P(G1, G2 | AL, COV, RL) ~ P(AL, COV, RL | G1, G2) * P(G1, G2) = prod_{read_i} P(al_i, cov_i, rl_i | G1, G2) * P(G1, G2)
         = independent G1 G2 = prod_{read_i} P(al_i, cov_i, rl_i | G1) * P(al_i, cov_i, rl_i | G2) * P(G1) * P(G2) {here G1, G2 is from possible
         alleles, background, and expanded, priors are from params}
         P(al_i, cov_i, rl_i | G1) - 2 options: 1. closed evidence (al_i = X), we know X; 2. open evidence (al_i >= X), cl_i == True if i is closed
         1.: P(al_i, cov_i, rl_i, cl_i | G1) = P(rl_i from read distrib.) * p(allele is al_i | G1) * P(read generated closed evidence | rl_i, al_i)
         2.: P(rl_i is from r.distr.) * P(allele is >= al_i | G1) * P(read generated open evidence | rl_i, al_i)
        :param annotations: list(Annotation) - closed annotated reads (both primers set)
        :param filt_annotations: list(Annotation) - open annotated reads (only one primer set)
        :param index_rep: int - index of a repetition
        :param verbose: bool - print more stuff?
        :param monoallelic: bool - do we have a mono-allelic motif (i.e. chrX/chrY and male sample?)
        :return: dict(tuple(int, int):float) - directory of model indices to their likelihood
        """
        # generate closed observed and read_length arrays
        observed_annots = np.array([ann.module_repetitions[index_rep] for ann in annotations])
        rl_annots = np.array([len(ann.read_seq) for ann in annotations])
        closed_annots = np.ones_like(observed_annots, dtype=bool)

        # generate open observed and read_length arrays
        observed_fa = np.array([ann.module_repetitions[index_rep] for ann in filt_annotations])
        rl_fa = np.array([len(ann.read_seq) for ann in filt_annotations])
        closed_fa = np.zeros_like(observed_fa, dtype=bool)

        # join them and keep the information if they are open or closed
        observed_arr = np.concatenate((observed_annots, observed_fa)).astype(int)
        rl_arr = np.concatenate((rl_annots, rl_fa)).astype(int)
        closed_arr = np.concatenate([closed_annots, closed_fa]).astype(bool)

        # generate the boundaries:
        overhead = 3
        if len(observed_annots) == 0:
            max_rep = max(observed_fa) + overhead  # non-inclusive
            min_rep = max(self.MIN_REPETITIONS, max(observed_fa) - overhead)  # inclusive
        else:
            max_rep = max(observed_annots) + overhead + 1  # non-inclusive
            min_rep = max(self.MIN_REPETITIONS, min(observed_annots) - overhead)  # inclusive

        # expanded allele
        e_allele = max_rep
        if len(observed_fa) > 0:
            e_allele = max(max_rep, max(observed_fa) + 1)

        # generate all the models
        self.construct_models(min_rep, max_rep, e_allele)

        # go through every model and evaluate:
        evaluated_models = {}
        for m1, m2 in generate_models_one_allele(min_rep, max_rep) if monoallelic else generate_models(min_rep, max_rep, multiple_backgrounds=True):
            evaluated_models[(m1, m2)] = 0.0
            # go through every read
            for obs, rl, closed in zip(observed_arr, rl_arr, closed_arr):
                lh = self.likelihood_read(obs, rl, m1, None if m2 == 'X' else m2, closed=closed)
                # TODO weighted sum according to the closeness/openness of reads?
                evaluated_models[(m1, m2)] += np.log(lh)

        return evaluated_models

    def print_pcolor(self, lh_dict: dict[tuple[int | str, int | str]: float], display_file: str | None,
                     name: str, lognorm: bool = True) -> tuple[np.array, tuple[int, int]]:
        """
        Get maximum likelihood option and alternatively print it to image file.
        :param lh_dict: dict(tuple(int, int):float) - directory of model indices to their likelihood
        :param display_file: str|None - filename for pcolor image output
        :param name: str - name to use in title
        :param lognorm: bool - use loglog scale in displaying likelihood array
        :return: tuple(int, int) - option with the highest likelihood
        """
        # convert to a numpy array:
        lh_array = np.zeros((self.max_rep, self.max_rep + 1))
        for (k1, k2), v in lh_dict.items():
            if k2 == 'X':  # if we have mono-allelic
                k2 = k1
            if k2 == 'B' or k1 == 'E' or (isinstance(k1, int) and isinstance(k2, int) and k2 < k1):  # B is the smallest, E is the largest!
                k1, k2 = k2, k1
            if k1 == 'B':
                k1 = 0
            if k2 == 'B':
                k2 = 0
            if k1 == 'E':  # only if k2 is 'E' too.
                k1 = 0
            if k2 == 'E':
                k2 = self.max_rep
            lh_array[k1, k2] = v

        # get minimal and maximal likelihood
        ind_good = (lh_array < 0.0) & (lh_array > -1e10) & (lh_array != np.nan)
        if len(lh_array[ind_good]) == 0:
            return lh_array, (0, 0)
        lh_array[~ind_good] = np.NINF
        z_min, z_max = min(lh_array[ind_good]), max(lh_array[ind_good])

        max_str = len(lh_array)

        # generate image file if specified:
        if display_file is not None:
            plt.figure()

            if lognorm:
                lh_view = -np.log(-lh_array)
                z_min = -np.log(-z_min)
                z_max = -np.log(-z_max)
            else:
                lh_view = lh_array.copy()

            # background (B, i) - copy it below min_rep
            lh_view[self.min_rep - 1, :] = lh_view[0, :]

            lh_copy = lh_view.copy()
            lh_copy[-1, self.min_rep] = lh_copy[0, 0]
            lh_copy[-1, self.min_rep + 1] = lh_copy[0, self.max_rep]

            # background (B,B)
            bg_size = max(2, (len(lh_view) - self.min_rep) // 6)
            if len(lh_view) - self.min_rep <= 6:
                bg_size = 1
            lh_view[-bg_size:, self.min_rep:self.min_rep + bg_size] = lh_view[0, 0]
            # expanded (E,E)
            lh_view[-bg_size:, self.min_rep + bg_size:self.min_rep + 2 * bg_size] = lh_view[0, self.max_rep]

            # plotting
            title = '%s likelihood of options (%s)' % ('Loglog' if lognorm else 'Log', name)
            plt.title(title)
            plt.xlabel('2nd allele')
            plt.ylabel('1st allele')
            start_ticks = 5
            step_ticks = 5
            plt.xticks(
                np.concatenate([np.array(range(start_ticks - self.min_rep, max_str - self.min_rep, step_ticks)), [max_str - self.min_rep]]) + 0.5,
                list(range(start_ticks, max_str, step_ticks)) + ['E(>%d)' % (self.max_with_e - 2)])
            plt.yticks(np.concatenate([np.array(range(start_ticks - self.min_rep + 1, max_str - self.min_rep + 1, step_ticks)), [0]]) + 0.5,
                       list(range(start_ticks, max_str, step_ticks)) + ['B'])
            palette = copy(plt.cm.jet)
            palette.set_under('gray', 1.0)
            plt.pcolor(lh_view[self.min_rep - 1:, self.min_rep:], cmap=palette, vmin=z_min, vmax=z_max)
            plt.colorbar()

            # draw dividing line(s):
            plt.plot([max_str - self.min_rep, max_str - self.min_rep], [0, max_str - self.min_rep + 1], 'k', linewidth=3)
            plt.plot([0, max_str - self.min_rep + 1], [1, 1], 'k', linewidth=3)

            # text background:
            plt.text(float(bg_size) / 2.0, max_str - self.min_rep + 1 - float(bg_size) / 2.0, 'BG', size=20, horizontalalignment='center',
                     verticalalignment='center', path_effects=[patheffects.withStroke(linewidth=2.5, foreground="w")])
            # text expanded
            plt.text(bg_size + float(bg_size) / 2.0, max_str - self.min_rep + 1 - float(bg_size) / 2.0, 'Exp', size=20, horizontalalignment='center',
                     verticalalignment='center', path_effects=[patheffects.withStroke(linewidth=2.5, foreground="w")])

            # save
            plt.savefig(display_file + '.pdf')
            plt.savefig(display_file + '.png')
            plt.close()

            # ----- PLOTLY HISTOGRAM -----
            text = [['' for _ in range(max_str - self.min_rep + 1)] for _ in range(max_str - self.min_rep + 1)]
            text[-1][0] = 'B'
            text[-1][1] = 'E'

            hovertext = [[f'{j}/{i}' for i in list(range(self.min_rep, max_str)) + ['E']] for j in ['B'] + list(range(self.min_rep, max_str))]
            hovertext[0][-1] = 'E/E'
            hovertext[-1][0] = 'B'
            hovertext[-1][1] = 'E'

            fig = go.Figure()
            fig.add_trace(go.Heatmap(z=lh_copy[self.min_rep - 1:, self.min_rep:],
                                     text=text, name='', hovertext=hovertext,
                                     showscale=True, colorscale='Jet'))
            fig.add_vline(x=max_str - self.min_rep - 0.5, line_width=5, line_color='black', opacity=1)
            fig.add_hline(y=0.5, line_width=5, line_color='black', opacity=1)

            fig.update_traces(texttemplate='%{text}', textfont_size=15,
                              hovertemplate='<b>%{{hovertext}} - {log} likelihood:\t%{{z}}</b>'.format(log='Loglog' if lognorm else 'Log'))
            fig.update_layout(width=500, height=450,
                              template='simple_white',
                              yaxis_fixedrange=True, xaxis_fixedrange=True,
                              title=title)
            fig.update_yaxes(title_text='1st allele', tickmode='array',
                             tickvals=np.concatenate([np.array(range(start_ticks - self.min_rep + 1, max_str - self.min_rep + 1, step_ticks)), [0]]),
                             ticktext=list(range(start_ticks, max_str, step_ticks)) + ['B'])
            fig.update_xaxes(title_text='2nd allele', tickmode='array',
                             tickvals=np.concatenate(
                                 [np.array(range(start_ticks - self.min_rep, max_str - self.min_rep, step_ticks)), [max_str - self.min_rep]]),
                             ticktext=list(range(start_ticks, max_str, step_ticks)) + ['E(>%d)' % (self.max_with_e - 2)])

            with open(display_file + '.json', 'w') as f:
                f.write(fig.to_json())

            # fig.write_image(display_file + '_plotly.png')

        # output best option
        best = sorted(np.unravel_index(np.argmax(lh_array), lh_array.shape))
        return lh_array, (int(best[0]), int(best[1]))

    def convert_to_sym(self, best: tuple[int, int], monoallelic: bool) -> tuple[int | str, int | str]:
        """
        Convert numeric alleles to their symbolic representations.
        :param best: (int, int) - numeric representation of alleles
        :param monoallelic: bool - if this is monoallelic version
        :return: (int|str, int|str) - symbolic representation of alleles
        """
        # convert it to symbols
        if best[0] == 0 and best[1] == self.max_rep:
            best_sym = ('E', 'E')
        else:
            best_sym = tuple(map(lambda x: 'E' if x == self.max_rep else 'B' if x == 0 else x, best))

        # if mono-allelic return 'X' as second allele symbol
        if monoallelic:
            best_sym = (best_sym[0], 'X')

        return best_sym

    def get_confidence(self, lh_array: np.ndarray, predicted: tuple[int, int], monoallelic: bool = False) -> tuple[float, float, float | str, float, float, float, float]:
        """
        Get confidence of a prediction.
        :param lh_array: 2D-ndarray - log likelihoods of the prediction
        :param predicted: tuple(int, int) - predicted alleles
        :param monoallelic: bool - do we have a mono-allelic motif (i.e. chrX/chrY and male sample?)
        :return: tuple[float, float, float | str, float, float, float, float] - prediction confidence of all, first, and second allele(s), background and expanded states
        """
        # get confidence
        lh_corr_array = lh_array - np.max(lh_array)
        lh_sum = np.sum(np.exp(lh_corr_array))
        confidence = np.exp(lh_corr_array[predicted[0], predicted[1]]) / lh_sum
        if predicted[0] == predicted[1]:  # same alleles - we compute the probability per allele
            confidence1 = np.sum(np.exp(lh_corr_array[predicted[0], :])) / lh_sum
            confidence2 = np.sum(np.exp(lh_corr_array[:, predicted[1]])) / lh_sum
        elif predicted[1] == lh_corr_array.shape[0]:  # expanded allele - expanded is only on one side of the array
            confidence1 = (np.sum(np.exp(lh_corr_array[predicted[0], :])) + np.sum(np.exp(lh_corr_array[:, predicted[0]])) - np.exp(
                lh_corr_array[predicted[0], predicted[0]])) / lh_sum
            confidence2 = np.sum(np.exp(lh_corr_array[:, predicted[1]])) / lh_sum
        else:  # normal behavior - different alleles , no expanded, compute all likelihoods of the alleles
            confidence1 = (np.sum(np.exp(lh_corr_array[predicted[0], :])) + np.sum(np.exp(lh_corr_array[:, predicted[0]])) - np.exp(
                lh_corr_array[predicted[0], predicted[0]])) / lh_sum
            confidence2 = (np.sum(np.exp(lh_corr_array[:, predicted[1]])) + np.sum(np.exp(lh_corr_array[predicted[1], :])) - np.exp(
                lh_corr_array[predicted[1], predicted[1]])) / lh_sum

        confidence_back = np.exp(lh_corr_array[0, 0]) / lh_sum
        confidence_back_all = np.sum(np.exp(lh_corr_array[0, :])) / lh_sum
        confidence_exp = np.exp(lh_corr_array[0, self.max_rep]) / lh_sum
        confidence_exp_all = np.sum(np.exp(lh_corr_array[:, self.max_rep])) / lh_sum

        if monoallelic:
            confidence2 = '---'

        return confidence, confidence1, confidence2, confidence_back, confidence_back_all, confidence_exp, confidence_exp_all

    @staticmethod
    def write_output(file_desc: str | typing.TextIO, predicted: tuple[int | str, int | str], conf: tuple[float, float, float | str, float, float, float, float], name: str | int):
        """
        Write result of one prediction.
        :param file_desc: file descriptor - where to write to
        :param predicted: tuple(int/char, int/char) - predicted alleles
        :param conf: tuple[float, float, float | str, float, float, float, float] - confidence of prediction (whole, 1st allele, 2nd allele, background and expanded alleles)
        :param name: str/int - name/number of the sample
        :return: None
        """
        def write_output_fd(f, predicted, conf, name):
            print(f'Predicted alleles for {name}: (confidence = {float_to_str(conf[1], percents=True)})', file=f)
            print(f'\t{str(predicted[0]):3s} (confidence = {float_to_str(conf[1], percents=True)})', file=f)
            print(f'\t{str(predicted[1]):3s} (confidence = {float_to_str(conf[2], percents=True)})', file=f)
            print(f'B   B   {float_to_str(conf[3], percents=True)}', file=f)
            print(f'all B   {float_to_str(conf[4], percents=True)}', file=f)
            print(f'B   E   {float_to_str(conf[5], percents=True)}', file=f)
            print(f'all E   {float_to_str(conf[6], percents=True)}', file=f)

        if type(file_desc) is str:
            with open(file_desc, 'w') as f:
                write_output_fd(f, predicted, conf, name)
        else:
            write_output_fd(file_desc, predicted, conf, name)

    def genotype(self, annotations: list[Annotation], filt_annotations: list[Annotation], index_rep: int, file_pcolor: str | None,
                 file_output: str | None, name: str, monoallelic: bool = False) -> tuple[tuple[str | int, str | int], tuple[float, float, float | str, float, float, float, float]]:
        """
        Genotype based on all annotations - infer likelihoods, print pcolor and write output
        :param annotations: list(Annotation) - good (blue) annotations
        :param filt_annotations: list(Annotation) - (grey) annotations with one primer
        :param index_rep: i
        nt - index of a repetition
        :param file_pcolor: str - file prefix for a pcolor image
        :param file_output: str - file for genotyping output
        :param name: str - name of the sample
        :param monoallelic: bool - do we have a mono-allelic motif (i.e. chrX/chrY and male sample?)
        :return: tuple - predicted symbols and confidences
        """
        # if we do not have any good annotations, then quit
        if len(annotations) == 0 and len(filt_annotations) == 0:
            return ('B', 'B'), (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

        # infer likelihoods
        lh_dict = self.infer(annotations, filt_annotations, index_rep, verbose=False, monoallelic=monoallelic)

        # print pcolor image
        lh_array, predicted = self.print_pcolor(lh_dict, file_pcolor, name)

        # adjust for no spanning reads (should output Background)
        if len(annotations) == 0:
            predicted = (0, 0)

        # convert numbers to symbols
        predicted_sym = self.convert_to_sym(predicted, monoallelic)

        # get confidence of our prediction
        confidence = self.get_confidence(lh_array, predicted, monoallelic)

        # write output
        if file_output is not None:
            self.write_output(file_output, predicted_sym, confidence, name)

        # return predicted and confidence
        return predicted_sym, confidence
