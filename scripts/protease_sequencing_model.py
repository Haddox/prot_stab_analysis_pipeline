# coding: utf-8



import logging
logger = logging.getLogger(__name__)

import copy

import traitlets

import numpy
import numpy as np

import theano.tensor as T
import theano

import pymc3
theano.config.compute_test_value = "ignore"

from scipy.special import expit
import scipy.stats.distributions

import scipy.stats
import scipy.optimize

from pymc3.distributions import Continuous
from pymc3.distributions.continuous import get_tau_sd
from pymc3.distributions.continuous import bound

from utility import resolve_subclass, SubclassName

class FlatNormal(Continuous):

	def __init__(self, mu=0.0, w=0.0, tau=None, sd=None, *args, **kwargs):
		super(FlatNormal, self).__init__(*args, **kwargs)
		self.mean = self.median = self.mode = self.mu = mu
		self.tau, self.sd = get_tau_sd(tau=tau, sd=sd)
		self.variance = 1. / self.tau
		self.w = w

	def logp(self, value):
		tau = self.tau
		sd = self.sd
		mu = self.mu
		w = self.w
		return bound(
			(-tau * ( T.maximum(
				T.minimum(value - mu+w, 0)**2.0,
				T.maximum(value - mu-w, 0)**2.0
			)) + T.log(tau / numpy.pi / 2.)) / 2.,
			tau > 0,
			sd > 0
		)

def unorm(v):
	return v / v.sum()

class SelectionResponseFn(object):
	@property
	def population_params(self):
		return []

	# Given ec50 and protease concentration, what fraction
	# of cells are going to survive?
	def selection_mass(self, **kwargs):
		# Split two ways here
		# 	Either we're going to do this full theano
		# 	or we're doing this full numpy
		if any(isinstance(v, (T.TensorVariable, T.TensorConstant)) for v in list(kwargs.values())):
			kwargs = { k : T.as_tensor_variable(v) for k, v in list(kwargs.items()) }
			return self.selection_mass_impl(num=T, **kwargs)
		else:
			return self.selection_mass_impl(num=numpy, **kwargs)

class LogisticResponse(SelectionResponseFn):
	def selection_mass_impl(self, num, sel_level, sel_k, sel_ec50):
		sel_xs = sel_k * (sel_level - sel_ec50)
		return 1 - 1 / (1 + num.exp(-sel_xs))

class NormalSpaceErfResponse(SelectionResponseFn):
	@property
	def population_params(self):
		return ["conc_factor"]

	def selection_mass_impl(self, num, sel_level, sel_k, sel_ec50, conc_factor, div_parent_is_expressor):
		if ( div_parent_is_expressor ):
			# runs only for round1 (i.e. where the parent is the 0 protease sort)
			return self.selection_mass_from_naive(num, sel_level, sel_k, sel_ec50, conc_factor) / \
				   self.selection_mass_from_naive(num, -100, sel_k, 0, conc_factor)
		else:
			# runs for rounds 2 and 3 (i.e. where the parent is a previously sorted pool that was grown up)
			return self.selection_mass_from_naive(num, sel_level, sel_k, sel_ec50, conc_factor)

	# Given ec50 and protease concentration, what fraction
	# of cells are going to survive?
	def selection_mass_impl_from_naive(self, num, sel_level, sel_k, sel_ec50, conc_factor):
		sel_xs = sel_k * (conc_factor ** (sel_level - sel_ec50) - 1.0)   # Bracketed term of Eq. 10
																		 # Note that [E]/EC_50 has been replaced by
		if num == numpy:                                                 # conc_factor ** (sel_level - sel_ec50)
			erf = scipy.special.erf                                      # because enzyme concentrations (sel_level)
		else:                                                            # are passed to this function on a log scale, with the base
			erf = T.erf                                                  # defined by conc_factor.


		return (1.0 - erf(sel_xs)) / 2.0                                # The full Eq. 10


class KdResponse(SelectionResponseFn):
	@property
	def population_params(self):
		return []

	def selection_mass_impl(self, num, conc, cs, mu, sigma, kd, viable_frac, csmu):

		if num == numpy:           
			erf = scipy.special.erf
			ln = np.log
		else:
			erf = T.erf
			ln = T.log

		return  (1/2 - 1/2 * erf( ( csmu + ln( ( kd + conc ) / conc ) ) / ( 1.4142 * sigma ) ))


class NormalSpaceLogisticResponse(SelectionResponseFn):
	@property
	def population_params(self):
		return ["conc_factor"]

	def selection_mass_impl(self, num, sel_level, sel_k, sel_ec50, conc_factor):
		sel_xs = sel_k * (conc_factor ** (sel_level - sel_ec50) - 1.0)
		return 1 - 1 / (1 + num.exp(-sel_xs * 2.45))

import sys
def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


# This class gets the following traits passed in when initalized:
    # response_fn = ("NormalSpaceErfResponse",),  # Given ec50 and protease concentration, what fraction of cells are going to survive?
    # min_selection_mass = [5e-7],
    # min_selection_rate = [0.0001],
    # outlier_detection_opt_cycles = [3],
    # sel_k = [0.8]
class FractionalSelectionModel(traitlets.HasTraits):
	response_fn = SubclassName(SelectionResponseFn)

	@property
	def response_impl(self):
	   return resolve_subclass(SelectionResponseFn, self.response_fn)()

	#sel_k = traitlets.Dict(
	#    traits = dict(
	#        __class__ = SubclassName(Continuous)
	#    ),
	#    default_value = dict(
	#        __class__="FlatNormal",
	#        mu=1.63, #4.0, #1.63,
	#        sd=0.00002, #0.0001, #0.00002,
	#        w=1.02 #2.5 #1.02,
	#    )
	#)
	#
	#@property
	#def sel_k_class(self):
	#    return resolve_subclass(Continuous, self.sel_k["__class__"])
	#
	#@property
	#def sel_k_kwargs(self):
	#    kwargs = dict(self.sel_k)
	#    kwargs.pop("__class__")
	#    return kwargs

	# sel_k = traitlets.Float()#min_selection_rate

	sigma = traitlets.Float()
	cs = traitlets.Float()


	min_selection_rate = traitlets.Union([
		traitlets.Bool(),
		traitlets.Float()
	])

	min_selection_mass = traitlets.Union([
		traitlets.Enum(["global", "per_selection"]),
		traitlets.Float()
	])

	homogenous_k = traitlets.Bool(default_value=True)

	outlier_detection_opt_cycles = traitlets.Integer(default_value=1)

	def __init__(self, **kwargs):
		# Override 'super' error-handling logic in HasTraits base __init__
		# __init__ swallows errors from unused kwargs until v4.3
		for key in kwargs:
			if not self.has_trait(key):
				raise TypeError("__init__() got an unexpected keyword argument '%s'" % key )
		traitlets.HasTraits.__init__(self, **kwargs)

	@staticmethod
	def lognorm_params(sd, mode):
		return dict(
			tau = sd**-2.,
			mu = numpy.log(mode) + sd ** 2,
		)

	# How many parents does this sort round have?
	#
	# pop_specs is just self.population_data
	# start_key comes from self.population_data.keys()
	@staticmethod
	def parent_depth(start_key, pop_specs):
		seen_keys = set()
		cur_depth = 0
		key = start_key

		while pop_specs[key]["parent"] is not None:
			seen_keys.add(key)
			parent = pop_specs[key]["parent"]

			if parent in seen_keys:
				raise ValueError(
					"Cycle in population parents. start_key: %s cycle_key: %s" %
					(start_key, pop_specs[key]["parent"]))
			if parent not in pop_specs:
				raise ValueError(
					"Invalid parent specified: %s pop_specs: %s" %
					(parent, list(pop_specs.keys())))
			cur_depth += 1
			key = parent

		return cur_depth

	def generate_data(self, pop_specs, sel_k, sel_ec50, init_pop):
		populations = {}

		for pkey in sorted(list(pop_specs.keys()), key=lambda pkey: self.parent_depth(pkey, pop_specs)):
			p = pop_specs[pkey]
			if p["parent"] is None:
				start_pop = init_pop
			else:
				start_pop = populations[p["parent"]]["P_sel"]

			start_dist = unorm(start_pop)
			if p["selection_level"] is not None:
				src_dist = start_dist * self.response_impl.selection_mass(
					sel_level = p["selection_level"], sel_k = sel_k, sel_ec50 = sel_ec50,
					div_parent_is_expressor = self.population_data[pdat['parent']]["selection_level"] is None,
					**{param : p[param] for param in self.response_impl.population_params}
				)
				fraction_selected = src_dist.sum() / start_dist.sum()
			else:
				src_dist = start_dist
				fraction_selected = 1

			selected = numpy.random.multinomial(p["P_sel"], unorm(src_dist)).astype(float)

			populations[pkey] = {}
			populations[pkey].update(p)
			populations[pkey].update({
				"P_sel" : selected,
				"Frac_sel_pop" : fraction_selected
			})

		return populations


	# First, just add the variable to the pymc3 model with the given distribution
	# But then, if the dist has a .transform (which I think they all do...)
	#  ??????
	def add_fit_param(self, name, dist):
		var = self.model.Var(name, dist, data=None)

		if dist.transform:
			forward_trans = dist.transform.forward(var)  # what does thes even do???
			self.to_trans[name] = (
				"%s_%s__" % (name, dist.transform.name),
				lambda v: forward_trans.eval({var : v})
			)


		self.fit_params[name] = var

		return var

	# selection_level == concentration
	# population_data looks like this:
	#    {
    #       0: 
    #       1: {parent: 0, concentration:1, num_selected:7600000, Frac_sel_pop: 0.97, conc_factor: 3, 
    #                      name: [], seq_counts: [], P_sel: []}
    #       2:
    #      ...
    #    }
	def build_model(self, population_data):

		# Check for differerences between population_data.items()
		#  and self.response_impl.population_params and throw warnings if there are differences
		for k, p in list(population_data.items()):
			unused_keys = set(p.keys()).difference(
				["P_sel", "Frac_sel_pop", "concentration", "parent"] +
				list(self.response_impl.population_params)
			)
			if unused_keys:
				logger.warning("Unused keys in population_data[%r] : %s", k, unused_keys)

		# Make sure that there are same number of designs in each sort round
		num_members = set( len(p["P_sel"]) for p in list(population_data.values()) )
		assert len(num_members) == 1, "Different observed population memberships: %s" % num_members
		self.num_members = num_members.pop()

		# Remake population_data but with "selection_level" as key and "P_sel" as value. counts0 gets excluded
		# selected_observations looks like this:
		# {
		#	1: [],
		#   2: [],
		#   ...
		# }
		selected_observations = {
			v["concentration"] : v["P_sel"] for v in list(population_data.values()) if v["concentration"] is not None
		}

		# First, assign all ec50 values to selection_level 0

		# Next, assign ec50 values to 1 - (last_round_with_counts)
		# for sl in selected_observations:
		# 	start_ec50[ ((sl - 1) > start_ec50) & (selected_observations[sl] > 0) ] = sl - 1

		self.model = pymc3.Model()

		self.to_trans = {}
		self.fit_params = {}

		self.model_populations = {}
		self.population_data = population_data

		# sort population_data.keys() by number of parents. i.e. naive, sort1, sort1, sort1, sort2, sort2, sort3, sort3
		pops_by_depth = sorted(
			list(population_data.keys()),
			key=lambda pkey: self.parent_depth(pkey, population_data))

		# Literally gets set to [1, 2, 3, 4, 5, 6] # but specifically not round 0
		self.modeled_populations = [ p for p in pops_by_depth if population_data[p]["parent"] is not None ]

		with self.model:
			#sel_k = self.add_fit_param(
			#    "sel_k",
			#    self.sel_k_class.dist(**self.sel_k_kwargs))

			#sel_k = self.add_fit_param(
			#            "sel_k",
			#            FlatNormal.dist(**default_selk_dict[self.sel_k_dict]))

			# sel_k=self.sel_k


			# sel_values = protease concentrations [0, 1, 2, 3, 4, 5, 6]
			sel_values = set(
				float(p["concentration"])
				for p in list(self.population_data.values())
				if p["concentration"] is not None
			)
			# sel_mag = 6
			# sel_mag = max(sel_values) - min(sel_values)
			# self.sel_range = {
			#  lower: -3,
			#  upper: 9
			# }
			self.sel_range = dict(lower=min(sel_values) / 50, upper=max(sel_values) * 10)
			logger.info("Inferred kd range: %s", self.sel_range)
			print(self.sel_range['lower'])

			start_kd = numpy.full_like(list(selected_observations.values())[0], self.sel_range['upper'] / 2 )

			# == sel_ec50 ==
			# Create a pymc3 random variable from -3 to 9 with guess = start_ec50
			kd = self.add_fit_param(
				"kd",
				pymc3.Uniform.dist(
					shape=self.num_members,
					testval=start_kd,
					**self.sel_range
				)
			)

			viable_frac = self.add_fit_param(
				"viable_frac",
				pymc3.Uniform.dist(
					shape=self.num_members,
					testval=0.15,
					lower=0.05,
					upper=0.50
				)
			)
			# kd = self.add_fit_param(
			# 	"kd",
			# 	pymc3.Lognormal.dist(
			# 		shape=self.num_members,
			# 		testval=start_kd,
			# 		mu=(np.log(max(sel_values)) + np.log(min(sel_values))) / 2,
			# 		sd=(np.log(max(sel_values)) - np.log(min(sel_values))) / 4)
			# 	)


			# == sel_ec50 ==
			# Create a pymc3 random variable from -3 to 9 with guess = start_ec50
			mu = self.add_fit_param(
				"mu",
				pymc3.Uniform.dist(
					# shape=self.num_members,
					testval=9.5,
					# sd=0.3,
					# mu=9.5,
					lower=8.5,
					upper=10.5)
				)

			csmu = self.add_fit_param(
				"csmu",
				pymc3.Uniform.dist(
					# shape=self.num_members,
					testval=0.1,
					# sd=0.3,
					# mu=9.5,
					lower=0.01,
					upper=1)
				)

			# mu = self.add_fit_param(
			# 	"mu",
			# 	pymc3.Normal.dist(
			# 		shape=self.num_members,
			# 		testval=9.5,
			# 		sd=0.3,
			# 		mu=9.5
			# 		)
			# 	)


			# sigma = self.add_fit_param(
			# 			"sigma",
			# 			pymc3.HalfNormal.dist(sd=1, testval=1))

			sigma = self.add_fit_param(
				"sigma",
				pymc3.Uniform.dist(
					# shape=self.num_members,
					testval=1,
					# sd=0.3,
					# mu=9.5,
					lower=0.8,
					upper=10)
				)

			cs = self.add_fit_param(
				"cs",
				pymc3.Uniform.dist(
					# shape=self.num_members,
					testval=5,
					# sd=0.3,
					# mu=9.5,
					lower=1,
					upper=9)
				)

			# cs = self.add_fit_param(
			# 			"cs",
			# 			pymc3.HalfNormal.dist(sd=5, testval=5))


			# add a random rate called min_selection_rate at which each cell can randomly get accross
			if self.min_selection_rate:
				if isinstance(self.min_selection_rate, bool):
					logger.info("Adding adaptive min_selection_rate.")
					# min_selection_rate = self.add_fit_param(
					# 	"min_selection_rate",
					# 	pymc3.HalfNormal.dist(sd=.0002, testval=.0001))

					min_selection_rate = self.add_fit_param(
						"min_selection_rate",
						pymc3.Uniform.dist(
							# shape=self.num_members,
							testval=0.0001,
							# sd=0.3,
							# mu=9.5,
							lower=0,
							upper=0.01)
						)
				else:
					logger.info("Adding const min_selection_rate: %.03f" % self.min_selection_rate)
					min_selection_rate = float(self.min_selection_rate)
			else:
				min_selection_rate = 0.0

			# idk what this does but it's set to a constant 5e-7 and takes the float() path
			if self.min_selection_mass:
				if self.min_selection_mass == "global":
					logger.info("Adding global adaptive min_selection_mass.")
					min_selection_mass = self.add_fit_param(
						"min_selection_mass",
						pymc3.HalfNormal.dist(sd=1e-3, testval=1e-12))
				elif self.min_selection_mass == "per_selection":
					logger.info("Adding per-selection adaptive min_selection_mass.")
					min_selection_mass = self.add_fit_param(
						"min_selection_mass",
						pymc3.HalfNormal.dist(shape=len(self.modeled_populations), sd=1e-3, testval=1e-12))
				else:
					logger.info("Adding const min_selection_mass: %.03f" % self.min_selection_mass)
					min_selection_mass = float(self.min_selection_mass)
			else:
				min_selection_mass = 0.0

			# pidx = [0, 1, 2, 3, 4, 5], pkey = [1, 2, 3, 4, 5, 6]
			for pidx, pkey in enumerate(self.modeled_populations): # Sum over all rounds as in Eq. 15

				# get dict of info for this selection step
				pdat = population_data[pkey]

				# == constant 5e-7
				p_min_selection_mass = (
					min_selection_mass[pidx]
						if self.min_selection_mass == "per_selection" else
					min_selection_mass
				)

				# start_pop = [counts in parent pool]
				start_pop = population_data[pdat["parent"]]["P_sel"].astype(float)
				# P_in = [frac of total parent pool]
				P_in = unorm(start_pop)

				# this is always true
				if pdat["concentration"] is not None:

					# Frac_sel is a pymc3 deterministic variable using the NormalSpaceErfResponse
					#  equation, and is based on sel_ec50
					# Indicates fraction of cells expected to survive given the ec50 and concentration
					Frac_sel = self.response_impl.selection_mass(                                  # Eq. 10, where "sel_k" is K_sel and
						conc = pdat["concentration"], mu = mu, sigma = sigma, cs = cs, kd = kd, viable_frac = viable_frac, csmu = csmu,
						**{param : pdat[param] for param in self.response_impl.population_params}  # selection_level is the enzyme "concentration" on a
					)                                                                              # log scale, with base specified in the input.

					# Frac_sel but with the random carryover
					Frac_sel_star = min_selection_rate + (Frac_sel * (1 - min_selection_rate))     # Eq. 11, with min_selection_rate as "a"

					# [Fraction of each design in the pool we expect to see, unnormalized]
					P_cleave_prenormalized = P_in * Frac_sel_star    # The numerator of Eq. 12, normalized shortly

					# Fraction of expressing cells that are expected survive in pymc3 variables
					#    P_in.sum() == 1
					Frac_sel_pop = P_cleave_prenormalized.sum() / P_in.sum() # Precalculate this for Eq. 14

				# never runs. Implies that we have a parent but didn't have protease
				else:
					P_cleave = P_in ## this gets overwritten
					Frac_sel_pop = 1.0

				#multinomial formula returns nan if any p == 0, file bug?
				# Add epsilon to selection prob to avoid nan-results when p==0

				# P_cleave is fraction of each design in the output pool we expect to see
				# clip P_cleave_prenormalized to 1
				#  bcov doesn't undestand lowerbound here
				P_cleave = unorm(                                                                     # Finish calculating P_cleave per Eq. 12 by normalizingNormalize P_cleave and ensure all probabilities
					T.clip(P_cleave_prenormalized, (p_min_selection_mass + 1e-9) * Frac_sel_pop, 1))  # it to sum to 1. A lower threshold is added to all probabilities to
				                                                                                      # prevent crashes. This low threshold is rarely relevant because
				# pop_mask = [ indices of designs that were present in parent pool ]
				pop_mask = numpy.flatnonzero(start_pop > 0)

				# n_sel = total number of cells collected for things that were in the parent pool
				n_sel = pdat["P_sel"][pop_mask].sum()

				# selected_%i = [ expected number of cells for this design collected given fraction we expect to see and total number of cells ]
				#   IMPORTANT NOTE: This is where the experimental counts is added with observed=pdat["P_sel"]
				selected = pymc3.distributions.Multinomial(     # Eq. 13
					name = "selected_%s" % pkey,                #
					n=n_sel ,                                   #
					p=P_cleave[pop_mask],
					observed=pdat["P_sel"][pop_mask]
				)

				if pdat.get("Frac_sel_pop", None) is not None:
					# n_sel = total cells collected that made it through the gate
					n_sel = pdat["P_sel"].sum()
					# n_assay = number of cells that were expressing proteins even including noise cells
					n_assay = numpy.floor(float(n_sel) / pdat["Frac_sel_pop"]) #/ pdat["parent_expression"])

					# Make a variable for total number of cells that made it through the gate
					#   given that we have the number that were expressing and the fraction
					#   that were collected relative to expressing
					
					# bcov modified this to be a selection from the parent pool rather than from the 
					#    sister-no-protease pool. The reason for this is that in very strong pools,
					#    Frac_sel_pop gets close to 1 and this Binomial has a lot of trouble.
					#    It's better to allow for the larger pool size of the true parent to make 
					#    the log-likihood slope gentler near p=1
					total_selected = pymc3.distributions.Binomial(    #
						name = "total_selected_%s" % pkey,            # Eq. 14
						n = n_assay,                                  #
						p = Frac_sel_pop,#*pdat['parent_expression'],   #
						observed = n_sel)                             #
				else:
					total_selected = pdat["P_sel"].sum()

				self.model_populations[pkey] = {
					"selection_mass" : self._function(P_cleave * Frac_sel_pop), # ???
					"P_cleave" : self._function(P_cleave),         # fraction of each design we expect to see
					"Frac_sel_pop" : self._function(Frac_sel_pop), # fraction of expressing cells that survived protease
					"P_sel" : self._function(selected),        # total cell of each design through the gate
					"n_sel" : self._function(total_selected),  # total cells through the gate
				}

		# make sure all fit params have been wrapped (why?)
		self.fit_params = { k : self._function(v) for k, v in list(self.fit_params.items()) }
		self.logp = self._function(self.model.logpt)

		return self

	def get_frac_sel_pop_llh(self, params):

		llhs = numpy.zeros(len(self.modeled_populations))

		for pidx, pkey in enumerate(self.modeled_populations):
			pdat = self.population_data[pkey]

			n_sel = pdat["P_sel"].sum()
			n_assay = numpy.floor(float(n_sel) / pdat["Frac_sel_pop"] )
			llh = scipy.stats.binom.logpmf(
				n_sel,
				n=n_assay,
				p=self.model_populations[pkey]["Frac_sel_pop"](params)
			)

			llhs[pidx] = llh
		return llhs.sum()


	def get_ec50_llh(self, params):
		return self.ec50_logp_traces(params, params['kd'][...,None], False).sum()



	def optimize_params(self, start = None):
		logger.info("optimize_params: %i members", self.num_members)
		if start is not None:
			start = self.to_transformed(start)
			for k in self.model.test_point:
				if k not in start:
					start[k] = self.model.test_point[k]
		# pymc3.sample()
		MAP = pymc3.find_MAP(start=start, model=self.model, method="L-BFGS-B") #fmin=scipy.optimize.fmin_l_bfgs_b)

		return { k : v(MAP) for k, v in list(self.fit_params.items()) }

	# exactly the same as opt_ec50_cred_outliers except faster
	def opt_ec50_cred_outliers2(self, src_params):
		logger.info("scan_ec50_outliers: %i members", self.num_members)
		params = copy.deepcopy(src_params)

		xs = self.get_ec50_trace_range()
		ec50_logp_traces = numpy.nan_to_num(self.ec50_logp_traces(params, xs))

		# most likely ec50 index
		max_points = numpy.argmax(ec50_logp_traces, axis=-1)

		# current ec50 index
		currents = numpy.searchsorted(xs, src_params["kd"], "left")

		# If we're off by 1, call it an outlier
		is_outlier = (max_points < currents-1) | ( max_points > currents )
		num_outlier = is_outlier.sum()

		# replace all the outlier points with the ec50 values of the max_point
		params["kd"][is_outlier] = xs[max_points[is_outlier]]

		params["kd"] = params["kd"].clip(xs[0], xs[-1])


		print("min_sel_rate: %.4f -- csmu: %.1f -- sigma: %.2f "%(params['min_selection_rate'],
					params['csmu'], params['sigma']) )

		if ( num_outlier / self.num_members > 0.2 ):
			print("resetting parameters")
			params['min_selection_rate'] = 0.001
			params['mu'] = 9.5
			params['sigma'] = 1
			params['viable_frac'][:] = 0.15
			params['cs'] = 5
			params['csmu'] = 0.1
		else:
			params['min_selection_rate'] = np.clip(params['min_selection_rate'], 0.0000000001, 0.0095)
			params['mu'] = np.clip(params['mu'], 8.7, 10.3)
			params['sigma'] = np.clip(params['sigma'], 0.9, 9)
			params['cs'] = np.clip(params['cs'], 2, 8)
			params['csmu'] = np.clip(params['csmu'], 0.05, 0.8)

		logger.info(
			"Modified %.3f outliers. (%i/%i)",
			 num_outlier / self.num_members, num_outlier, self.num_members)

		return params


	def opt_ec50_cred_outliers(self, src_params):
		logger.info("scan_ec50_outliers: %i members", self.num_members)
		params = copy.deepcopy(src_params)

		# bcov -- massive speed boost. Also probably helps with stability since it doesn't get modified as we update outliers
		precomputed_Frac_sel_pop = {
			pkey : self.model_populations[pkey]["Frac_sel_pop"](params)
			for pkey in self.modeled_populations
		}
		print(precomputed_Frac_sel_pop)

		num_outlier = 0

		for i in range(self.num_members):
			if i % 1000 == 0:
				logger.info("scan_ec50_outliers: %i / %i  outlier count: %s", i, self.num_members, num_outlier)

			cred_summary = self.estimate_ec50_cred(params, i, precomputed_Frac_sel_pop=precomputed_Frac_sel_pop)

			# figure out which X we correspond to
			current = numpy.searchsorted(cred_summary["xs"], cred_summary["sel_ec50"], "left")
			#rb = numpy.searchsorted(cred_summary["xs"], cred_summary["sel_ec50"], "right")

			# what is the probability of the most likely bin
			m_pmf = cred_summary["pmf"].argmax()

			if ( i == 31227 ):
				print(cred_summary["xs"][m_pmf])

			# Basically, are we in the most probably bin or right next to it. Otherwise call it an outlier and assign most probable
			if m_pmf < current - 1 or m_pmf > current:
				if ( i == 31227 ):
					print(cred_summary["xs"][m_pmf])
				num_outlier += 1
				params["sel_ec50"][i] = cred_summary["xs"][m_pmf]

		logger.info(
			"Modified %.3f outliers. (%i/%i)",
			 num_outlier / self.num_members, num_outlier, self.num_members)

		return params

	def find_MAP(self, start = None):
		params = self.optimize_params(start)

		for _ in range(self.outlier_detection_opt_cycles):

			resampled = self.opt_ec50_cred_outliers2(params)
			params = self.optimize_params(resampled)


		params = self.opt_ec50_cred_outliers2(params)

		return params

	def get_ec50_trace_range(self):
		# return numpy.linspace(self.sel_range['lower']+1,self.sel_range['upper']-1, (self.sel_range['upper'] - self.sel_range['lower'] - 2)*10 + 1
		return np.exp(numpy.linspace(np.log(self.sel_range['lower']+1),np.log(self.sel_range['upper']-1), 
					(np.log(self.sel_range['upper']) - np.log(self.sel_range['lower']) -2 )*10 + 1))

	# calculate all ec50 traces at the same time
	def ec50_logp_traces(self, base_params, ec50_range, subtract_max=True ):
		llh_by_ec50_gen = numpy.zeros((len(base_params['kd']), ec50_range.shape[-1], len(self.model_populations)))

		if ( len(ec50_range.shape) > 1 ):
			assert( ec50_range.shape[0] == len(base_params['kd']) )

		if self.min_selection_rate:
			if self.min_selection_rate == True:
				min_selection_rate = base_params["min_selection_rate"]
			else:
				min_selection_rate = self.min_selection_rate
		else:
			min_selection_rate = 0

		if self.min_selection_mass:
			if isinstance(self.min_selection_mass, str):
				min_selection_mass = base_params["min_selection_mass"]
			else:
				min_selection_mass = self.min_selection_mass
		else:
			min_selection_mass = 0

		# for loop over 6 rounds of selection
		for pidx, pkey in enumerate(self.modeled_populations):
			pdat = self.population_data[pkey]

			p_min_selection_mass = (
				min_selection_mass[pidx]
					if self.min_selection_mass == "per_selection" else
				min_selection_mass
			)

			# this design's fraction in parent population
			parent_pop_fraction = unorm(self.population_data[pdat['parent']]["P_sel"])

			in_parent_pool = parent_pop_fraction > 0

			## calculate selection results for full ec50 range
			## base_params['sel_k'] * (pdat['conc_factor'] ** (pdat['selection_level'] - ec50_range) - 1.0 ))

			# sel_k = 0.8
			# if "sel_k" in base_params:
			# 	sel_k = base_params["sel_k"]
			# else:
			# 	sel_k = self.sel_k

			selected_fraction = self.response_impl.selection_mass(                                  # Eq. 10, where "sel_k" is K_sel and
				conc = pdat["concentration"], mu = base_params["mu"], sigma = base_params["sigma"], cs = base_params["cs"], kd = ec50_range,
				viable_frac = base_params['viable_frac'][...,None], csmu = base_params['csmu'],
				**{param : pdat[param] for param in self.response_impl.population_params}  # selection_level is the enzyme "concentration" on a
					) 

			# fraction surviving protease for the whole ec50 sweep
			# selected_fraction = self.response_impl.selection_mass(
			# 	sel_level = pdat["selection_level"], sel_k = sel_k, sel_ec50 = ec50_range,
			# 	div_parent_is_expressor = self.population_data[pdat['parent']]["selection_level"] is None,
			# 	**{param : pdat[param] for param in self.response_impl.population_params}
			# )

			selected_fraction = min_selection_rate + (selected_fraction * (1 - min_selection_rate))

			# fraction we expect to see in this population. selected_fraction / Frac_sel_pop is enrichment
			sel_pop_fraction = parent_pop_fraction[...,None] * (selected_fraction / self.model_populations[pkey]["Frac_sel_pop"](base_params) )

			# log-likelihood of seeing the actual counts given total counts and our expected fraction in population
			sample_llhs = scipy.stats.binom.logpmf(
				pdat["P_sel"][...,None],  # counts of this design in this selection step
				n=pdat["P_sel"].sum(),    # total counts in this selection step
				p=numpy.clip(sel_pop_fraction, p_min_selection_mass + 1e-9, 1.0) 
			)

			llh_by_ec50_gen[in_parent_pool,:,pidx] = sample_llhs[in_parent_pool]

		llh_by_ec50 = llh_by_ec50_gen.sum(axis=-1)

		if ( subtract_max ):
			return llh_by_ec50 - numpy.nanmax(llh_by_ec50, axis=-1)[...,None]
		else:
			return llh_by_ec50

	# compute log-probability of each ec50 for a design given a sweep of ec50s
	# bcov changed include_global_terms to false because it appears to hurt more than help
	def ec50_logp_trace(self, base_params, sample_i, ec50_range, include_global_terms=False, precomputed_Frac_sel_pop=None):
		llh_by_ec50_gen = numpy.zeros((len(ec50_range), len(self.model_populations)))

		if self.min_selection_rate:
			if self.min_selection_rate == True:
				min_selection_rate = base_params["min_selection_rate"]
			else:
				min_selection_rate = self.min_selection_rate
		else:
			min_selection_rate = 0

		if self.min_selection_mass:
			if isinstance(self.min_selection_mass, str):
				min_selection_mass = base_params["min_selection_mass"]
			else:
				min_selection_mass = self.min_selection_mass
		else:
			min_selection_mass = 0

		# for loop over 6 rounds of selection
		for pidx, pkey in enumerate(self.modeled_populations):
			pdat = self.population_data[pkey]

			p_min_selection_mass = (
				min_selection_mass[pidx]
					if self.min_selection_mass == "per_selection" else
				min_selection_mass
			)

			# this design's fraction in parent population
			parent_pop_fraction = unorm(self.population_data[pdat['parent']]["P_sel"])[sample_i]

			if parent_pop_fraction == 0:
				continue

			## calculate selection results for full ec50 range
			## base_params['sel_k'] * (pdat['conc_factor'] ** (pdat['selection_level'] - ec50_range) - 1.0 ))

			# sel_k = 0.8
			if "sel_k" in base_params:
				sel_k = base_params["sel_k"]
			else:
				sel_k = self.sel_k

			# fraction surviving protease for the whole ec50 sweep
			selected_fraction = self.response_impl.selection_mass(
				sel_level = pdat["selection_level"], sel_k = sel_k, sel_ec50 = ec50_range,
				div_parent_is_expressor = self.population_data[pdat['parent']]["selection_level"] is None,
				**{param : pdat[param] for param in self.response_impl.population_params}
			)

			selected_fraction = min_selection_rate + (selected_fraction * (1 - min_selection_rate))

			# fraction we expect to see in this population. selected_fraction / Frac_sel_pop is enrichment
			sel_pop_fraction = parent_pop_fraction * (selected_fraction / precomputed_Frac_sel_pop[pkey]) # self.model_populations[pkey]["Frac_sel_pop"](base_params)

			# log-likelihood of seeing the actual counts given total counts and our expected fraction in population
			sample_llhs = scipy.stats.binom.logpmf(
				pdat["P_sel"][sample_i],  # counts of this design in this selection step
				n=pdat["P_sel"].sum(),    # total counts in this selection step
				p=numpy.clip(sel_pop_fraction, p_min_selection_mass + 1e-9, 1.0) 
			)

			# Redo calculations, but include effects on Frac_sel_pop that this selection of ec50 would bring
			# bcov believes that this does more harm than good as it allows for wacky situations to seem ok
			if include_global_terms and pdat.get("Frac_sel_pop") is not None:

				# calculate fraction that survive protease using current model sel_ec50
				prev_selected_fraction = self.response_impl.selection_mass(
					sel_level = pdat["selection_level"], sel_k = sel_k, sel_ec50 = base_params['sel_ec50'][sample_i],
					div_parent_is_expressor = self.population_data[pdat['parent']]["selection_level"] is None,
					**{param : pdat[param] for param in self.response_impl.population_params}
				)

				# selected mass based on current model
				prev_selected_mass = parent_pop_fraction * prev_selected_fraction

				# selected masses based on parameter sweep
				selected_mass = parent_pop_fraction * selected_fraction

				selected_count = pdat["P_sel"].sum()

				# expected counts in parent pool?
				source_count = numpy.floor(float(selected_count) / pdat["Frac_sel_pop"])

				# Frac_sel_pop adjusted for how much setting this ec50 would change it by 
				modified_global_selection_fractions = (
					precomputed_Frac_sel_pop[pkey] #self.model_populations[pkey]["Frac_sel_pop"](base_params)
					+ selected_mass - prev_selected_mass
				)

				sample_llhs += scipy.stats.binom.logpmf(
					selected_count,
					n=source_count,
					p=modified_global_selection_fractions
				)


			llh_by_ec50_gen[:,pidx] = sample_llhs

		llh_by_ec50 = llh_by_ec50_gen.sum(axis=1)

		return llh_by_ec50 - numpy.nanmax(llh_by_ec50)

	# calculate all ec50 cred intervals at the same time
	def estimate_ec50_creds(self, base_params, cred_spans = [.68, .95], super_span=-15):
		"""Estimate EC50 credible interval for a single ec50 parameter via model probability."""
		#xs = numpy.arange(self.sel_range["lower"]+0.1, self.sel_range["upper"]-0.1, .1)
		xs = self.get_ec50_trace_range()

		# get the logp value for each ec50 in xs
		logps = numpy.nan_to_num(self.ec50_logp_traces(base_params, xs))

		# compute probability mass function
		pmfs = numpy.exp(logps) / numpy.sum(numpy.exp(logps), axis=-1)[...,None]
		cdfs = numpy.cumsum(pmfs, axis=-1)

		# There's no way to vectorize searcsorted so this will have to do
		cred_summaries = []
		for ec50_i in range(self.num_members):

			logp = logps[ec50_i]
			cdf = cdfs[ec50_i]
			pmf = pmfs[ec50_i]

			# get credibility intervals by just looking at where the cdf crosses 0.05 and 0.95
			cred_intervals = {}
			for cred_i in cred_spans:
				cdf_b = (1 - cred_i) / 2
				l_b = xs[numpy.searchsorted(cdf, cdf_b, side="left")]
				u_b = xs[numpy.searchsorted(cdf, 1 - cdf_b, side="right")]
				cred_intervals[cred_i] = (l_b, u_b)

			super_low = xs[np.searchsorted(np.log(cdf), super_span, side="left").clip(0, len(xs)-1)]
			super_high = xs[np.searchsorted(-np.log(1-cdf), -super_span, side="right").clip(0, len(xs)-1)]

			kd = base_params["kd"][ec50_i]
			if ( super_low == xs[0] ):
				super_low = np.exp( np.log(kd) - ( np.log(super_high) - np.log(kd) ) )
			if ( super_high == xs[-1] ):
				super_high = np.exp( np.log(kd) + ( np.log(kd) - np.log(super_low) ) )

			cred_summaries.append(dict(
				xs = xs,
				pmf = pmf,
				cdf = cdf,
				logp = logp,
				kd = base_params["kd"][ec50_i],
				cred_intervals = cred_intervals,
				super_span = (super_low, super_high)
			))
		return cred_summaries

	def estimate_ec50_cred(self, base_params, ec50_i, cred_spans = [.68, .95], precomputed_Frac_sel_pop=None):
		"""Estimate EC50 credible interval for a single ec50 parameter via model probability."""
		#xs = numpy.arange(self.sel_range["lower"]+0.1, self.sel_range["upper"]-0.1, .1)
		xs = self.get_ec50_trace_range()

		# get the logp value for each ec50 in xs
		logp = numpy.nan_to_num(self.ec50_logp_trace(base_params, ec50_i, xs, False, precomputed_Frac_sel_pop=precomputed_Frac_sel_pop))

		# compute probability mass function
		pmf = numpy.exp(logp) / numpy.sum(numpy.exp(logp))
		cdf = numpy.cumsum(pmf)

		# get credibility intervals by just looking at where the cdf crosses 0.05 and 0.95
		cred_intervals = {}
		for cred_i in cred_spans:
			cdf_b = (1 - cred_i) / 2
			l_b = xs[numpy.searchsorted(cdf, cdf_b, side="left")]
			u_b = xs[numpy.searchsorted(cdf, 1 - cdf_b, side="right")]
			cred_intervals[cred_i] = (l_b, u_b)

		return dict(
			xs = xs,
			pmf = pmf,
			cdf = cdf,
			logp = logp,
			sel_ec50 = base_params["sel_ec50"][ec50_i],
			cred_intervals = cred_intervals
		)

	@staticmethod
	def plot_cred_summary(ec50_cred, ax=None):
		if ax is None:
			from matplotlib import pylab
			ax = pylab.gca()

		ax.plot( ec50_cred["xs"], ec50_cred["pmf"], label="pmf" )
		ax.plot( ec50_cred["xs"], ec50_cred["cdf"], label="cdf" )
		ax.axvline(ec50_cred["sel_ec50"], alpha=.5, label="sel_e50: %.2f" % ec50_cred["sel_ec50"])

		for ci, (cl, cu) in list(ec50_cred["cred_intervals"].items()):
			ax.axvspan(cl, cu, color="red", alpha=.2, label="%.2f cred" % ci)

	def model_selection_summary(self, params):
		def normed_pop(v):
			return v / v.sum()

		return {
			pkey : {
				"P_sel"  : self.population_data[pkey]["P_sel"],
				"selected_fraction" : normed_pop(self.population_data[pkey]["P_sel"].astype(float)),
				"P_cleave" : normed_pop(self.model_populations[pkey]["P_cleave"](params))
			}
			for pkey in self.model_populations
		}


	# make a really nice graph for a single design
	def plot_kdfit_summary(model, i, fit):
		import scipy.stats
		from matplotlib import pylab
		import matplotlib.pyplot as plt


		sel_sum = model.model_selection_summary(fit)

		sel_levels = {
			k : p["concentration"] if p["concentration"] else 0
			for k, p in list(model.population_data.items())}

		sel_fracs = {
			k : p["P_sel"][i] / p["P_sel"].sum()
			for k, p in list(model.population_data.items())}

		pylab.xticks(
			list(sel_levels.values()), list(sel_levels.keys()))
		# pylab.xlim((-1, 7))

		porder = [
			k for k, p in
			sorted(model.population_data.items(), key=lambda x: 9e9 if x[1]["concentration"] is None else x[1]["concentration"])]

		pylab.plot(
			[sel_levels[k] for k in porder],
			[sel_fracs[k] for k in porder],
			"-o",
			color="black", label="observed")

		lbl = False
		for k in sel_sum:
			n = sel_sum[k]["P_sel"].sum()
			p = sel_sum[k]["selected_fraction"][i]
			sel_level = model.population_data[k]["concentration"]
			counts=sel_sum[k]["P_sel"][i]
			plt.text(sel_levels[k] + 0.2, sel_fracs[k], '%.0f' % counts)

			if p<=0:
				continue

			bn = scipy.stats.binom(n=n, p=p)

			parkey = model.population_data[k]["parent"]
			pylab.plot(
				[sel_levels[parkey], sel_levels[k]],
				[sel_fracs[parkey], float(bn.ppf(.5)) / n],
				"--", color="red", alpha=.25
			)




			for ci in (.68, .95, .99):
				pylab.plot(
					[sel_level] * 2, bn.ppf([ci, 1-ci]) / n,
					linewidth=10, color="red", alpha=.25,
					label="predicted" if not lbl else None
				)
				lbl=True

		pylab.xscale("log")

		expected_f = []
		xs = []
		for pkey in porder:
			pdat = model.population_data[pkey]
			if ( pdat['parent'] is None ):
				continue
			xs.append(pdat['concentration'])

			parent_pop_fraction = unorm(model.population_data[pdat['parent']]["P_sel"])[i]
			selected_fraction = model.response_impl.selection_mass(                                  # Eq. 10, where "sel_k" is K_sel and
						conc = pdat["concentration"], mu = fit['mu'], sigma = fit['sigma'], cs = fit['cs'], kd = fit['kd'][i],
						viable_frac = fit['viable_frac'][i], csmu = fit['csmu'],
						**{param : pdat[param] for param in model.response_impl.population_params}  # selection_level is the enzyme "concentration" on a
					) 
			# selected_fraction = model.response_impl.selection_mass(
			# 	sel_level = pdat["selection_level"], sel_k = model.sel_k, sel_ec50 = fit['sel_ec50'][i],
			# 	div_parent_is_expressor = model.population_data[pdat['parent']]["selection_level"] is None,
			# 	**{param : pdat[param] for param in model.response_impl.population_params}
			# )
			selected_fraction = fit["min_selection_rate"] + (selected_fraction * (1 - fit["min_selection_rate"]))
			sel_pop_fraction = parent_pop_fraction * (selected_fraction / model.model_populations[pkey]["Frac_sel_pop"](fit) )

			expected_f.append(sel_pop_fraction)

		pylab.plot(
			xs,
			expected_f,
			"-o",
			color="blue", label="model")
		pylab.legend(fontsize="large", loc="best")

		pylab.twinx()
		xs = numpy.linspace(-2, 8)
		sel_ec50 = fit["kd"][i]
		# min_sel = fit["min_selection_rate"][i]
		# mu = fit['mu'][i]
		via_frac = fit['viable_frac'][i]
		# sel_k = fit["sel_k"][i] if len(fit["sel_k"]) > 1 else fit["sel_k"]
		# sel_k = model.sel_k
		# pylab.plot(xs, scipy.special.expit(-sel_k * (xs - sel_ec50)), alpha=.75)
		pylab.yticks([], [])

		pylab.title("%s - kd: %.2f - via_frac: %.2f" % (i, sel_ec50, via_frac,  ))


	# make a really nice graph for a single design
	def plot_fit_summary(model, i, fit):
		import scipy.stats
		from matplotlib import pylab
		import matplotlib.pyplot as plt


		sel_sum = model.model_selection_summary(fit)

		sel_levels = {
			k : p["selection_level"] if p["selection_level"] else 0
			for k, p in list(model.population_data.items())}

		sel_fracs = {
			k : p["P_sel"][i] / p["P_sel"].sum()
			for k, p in list(model.population_data.items())}

		pylab.xticks(
			list(sel_levels.values()), list(sel_levels.keys()))
		pylab.xlim((-1, 7))

		porder = [
			k for k, p in
			sorted(model.population_data.items(), key=lambda x: 0 if x[1]["selection_level"] is None else x[1]["selection_level"])]

		pylab.plot(
			[sel_levels[k] for k in porder],
			[sel_fracs[k] for k in porder],
			"-o",
			color="black", label="observed")

		lbl = False
		for k in sel_sum:
			n = sel_sum[k]["P_sel"].sum()
			p = sel_sum[k]["selected_fraction"][i]
			sel_level = model.population_data[k]["selection_level"]
			counts=sel_sum[k]["P_sel"][i]
			plt.text(sel_levels[k] + 0.2, sel_fracs[k], '%.0f' % counts)

			if p<=0:
				continue

			bn = scipy.stats.binom(n=n, p=p)

			parkey = model.population_data[k]["parent"]
			pylab.plot(
				[sel_levels[parkey], sel_levels[k]],
				[sel_fracs[parkey], float(bn.ppf(.5)) / n],
				"--", color="red", alpha=.25
			)




			for ci in (.68, .95, .99):
				pylab.plot(
					[sel_level] * 2, bn.ppf([ci, 1-ci]) / n,
					linewidth=10, color="red", alpha=.25,
					label="predicted" if not lbl else None
				)
				lbl=True

		pylab.legend(fontsize="large", loc="best")

		expected_f = []
		xs = []
		for pkey in porder:
			pdat = model.population_data[pkey]
			if ( pdat['parent'] is None ):
				continue
			xs.append(pkey)

			parent_pop_fraction = unorm(model.population_data[pdat['parent']]["P_sel"])[i]
			selected_fraction = model.response_impl.selection_mass(
				sel_level = pdat["selection_level"], sel_k = model.sel_k, sel_ec50 = fit['sel_ec50'][i],
				div_parent_is_expressor = model.population_data[pdat['parent']]["selection_level"] is None,
				**{param : pdat[param] for param in model.response_impl.population_params}
			)
			selected_fraction = model.min_selection_rate + (selected_fraction * (1 - model.min_selection_rate))
			sel_pop_fraction = parent_pop_fraction * (selected_fraction / model.model_populations[pkey]["Frac_sel_pop"](fit) )

			expected_f.append(sel_pop_fraction)

		pylab.plot(
			xs,
			expected_f,
			"-o",
			color="blue", label="observed")

		pylab.twinx()
		xs = numpy.linspace(-2, 8)
		sel_ec50 = fit["sel_ec50"][i]
		# sel_k = fit["sel_k"][i] if len(fit["sel_k"]) > 1 else fit["sel_k"]
		sel_k = model.sel_k
		pylab.plot(xs, scipy.special.expit(-sel_k * (xs - sel_ec50)), alpha=.75)
		pylab.yticks([], [])

		pylab.title("%s - ec50: %.2f - k: %.2f" % (i, sel_ec50, sel_k))

	def model_outlier_summary(self, params):
		selection_summary = self.model_selection_summary(params)

		for v in list(selection_summary.values()):
			logpmf = scipy.stats.binom.logpmf(
				v["P_sel"],
				n=v["P_sel"].sum(),
				p=v["P_cleave"])

			max_logpmf = scipy.stats.binom.logpmf(
				numpy.round(v["P_sel"].sum() * v["P_cleave"]),
				n=v["P_sel"].sum(),
				p=v["P_cleave"])

			sel_llh = logpmf - max_logpmf
			v["sel_log_likelihood"] = numpy.where(sel_llh != -numpy.inf, sel_llh, numpy.nan)

			sel_error = (v["P_sel"] / v["P_sel"].sum()) - v["P_cleave"]
			v["sel_log_likelihood_signed"] = -v["sel_log_likelihood"] * numpy.sign(sel_error)

		return selection_summary

	# Some of the pymc3 variables are internally stored in a transformed format and this converts
	# an entire dictionary to that format
	def to_transformed(self, val_dict):
		r = {}
		for n, val in list(val_dict.items()):
			if n in self.to_trans:
				k, f = self.to_trans[n]
				r[k] = f(val)
			else:
				r[n] = val

		return r

	# wrapper function that lets you pass theano functions
	# a dictionary of values for parameters
	def _function(self, f):
		if isinstance(f, theano.tensor.TensorVariable):
			fn = theano.function(self.model.free_RVs, f, on_unused_input="ignore")

			def call_fn(val_dict):
				val_dict = self.to_transformed(val_dict)
				# eprint([str(n) for n in self.model.free_RVs])
				return fn(*[val_dict[str(n)] for n in self.model.free_RVs])

			return call_fn
		else:
			# if it's not a theano variable, make a function that absorbs all arguments and just returns
			#  the passed value
			return lambda _: f

import unittest

class TestFractionalSelectionModel(unittest.TestCase):
	def setUp(self):
		numpy.random.seed(1663)

		test_members = 10
		num_sampled = 1e4
		self.pop_specs = dict(
			[(0, dict(parent = None, selection_level = None, selected = num_sampled))] +
			[(g, dict(parent = g - 1, selection_level = g, selected = num_sampled))
				for g in range(1, 7)]
		)
		self.test_vars = {
			"sel_k" : 1.5,
			"sel_ec50" : numpy.concatenate((
				numpy.random.normal(loc=.5, scale=.75, size=int(test_members * .9)),
				numpy.random.normal(loc=3, scale=1.5, size=int(test_members * .1)),
			)),
			"init_pop" : numpy.random.lognormal(size=test_members)
		}


	def test_basic_model(self):
		self.test_model = FractionalSelectionModel(homogenous_k=True)
		self.test_data = self.test_model.generate_data(self.pop_specs, **self.test_vars)

		test_map = self.test_model.build_model( self.test_data ).find_MAP()

		numpy.testing.assert_allclose( test_map["sel_k"], self.test_vars["sel_k"], rtol=1e-2, atol=.025)
		numpy.testing.assert_allclose( test_map["sel_ec50"], self.test_vars["sel_ec50"], rtol=.3)
