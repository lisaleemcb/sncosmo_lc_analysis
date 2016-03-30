import sncosmo
import triangle
from astropy.table import Table
import math
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import numpy as np
import collections

class LCA(object):
    """class to streamline light curve fits with SNCosmo """
    def __init__(self, data, model, vparams, truths=None):
        # super(, self).__init__()  <-- don't need this yet

        self.sn_id = 0
        self.data = data

        #source = sncosmo.get_source('salt2-extended')
        self.model = model
        #self.model = = sncosmo.Model(source=source)
        self.vparams = vparams
        self.bounds = {'c':(-0.3, 0.3), 'x1':(-3.0, 3.0)}

        self._fit_out = None
        self.fit_res = None
        self.fit_model = None

        self._mcmc_out = None
        self.mcmc_res = None
        self.mcmc_model = None

        self._nest_out = None
        self.nest_res = None
        self.nest_model = None
        # variables for time statistics


    @property
    def fit_out(self):
        return self._fit_out

    @fit_out.getter
    def fit_out(self):
        if self._fit_out is None:
            print "running chi^2 fit"

            self._fit_out = self.runMLE()
            self.fit_res = self._fit_out[0]
            self.fit_model = self._fit_out[1]

        return self._fit_out

    @fit_out.setter
    def fit_out(self, val):
        self._fit_out = val

        if val:
            self.fit_res = self._fit_out[0]
            self.fit_model = self._fit_out[1]

        else:
            self.fit_res = None
            self.fitmode = None

        return self._fit_out

    @property
    def mcmc_out(self):
        return self.mcmc_out

    @mcmc_out.getter
    def mcmc_out(self):
        if self._mcmc_out is None:
            print "running MCMC fit"

            self._mcmc_out = self.runMCMC()
            self.mcmc_res = self._mcmc_out[0]
            self.mcmc_model = self._mcmc_out[1]

        return self._mcmc_out

    @mcmc_out.setter
    def mcmc_out(self, val):
        self._mcmc_out = val

        if val:
            self.mcmc_res = self._mcmc_out[0]
            self.mcmc_model = self._mcmc_out[1]

        else:
            self.mcmc_res = None
            self.mcmc_model = None
        return self._mcmc_out

    @property
    def nest_out(self):
        return self.nest_out

    @nest_out.getter
    def nest_out(self):
        if self._nest_out is None:
            print "running nest fit"

            self._nest_out = self.run_nest()
            self.nest_res = self._nest_out[0]
            self.nest_model = self._nest_out[1]

        return self._nest_out

#methods
#--------------

# functions for changing aspects of fits
    def runMLE(self):
        MLEout = sncosmo.fit_lc(self.data, self.model, vparam_names=self.vparams,
                                    bounds=self.bounds, minsnr=3.0)

        return MLEout


    def runMCMC(self):
        mcmc_out = sncosmo.mcmc_lc(self.data, self.model, vparam_names=self.vparams,
                                    bounds=self.bounds, minsnr=3.0)

        return mcmc_out

    def run_nest(self):
        nest_out = sncosmo.nest_lc(sn, model, vparam_names=['t0', 'x0', 'x1', 'c'],
                                    bounds={'c':(-0.3, 0.3), 'x1':(-3.0, 3.0)},
                                    guess_amplitude_bound=True, minsnr=3.0, verbose=True)

    def reset_fit(self, fit):
        if fit is 'MLE':
            self._fit_out = None;
        if fit is 'MCMC':
            self._mcmc_out = None;
        if fit is 'nest':
            self._nest_out = None;

    def reset_all(self):
        self._fit_out = None;
        self._mcmc_out = None;
        self._nest_out = None;

        return

    def rerun_fits(self):
        if not self.fit_res:
            print "fit_lc() output reran"

            reset_all(self)
            self._fit_out = runMLE(self)

            self.fit_res = self._fit_out[0]
            self.fit_model = self._fit_out[1]

        if not self.mcmc_res:
            print "mcmc_lc() output reran"
            self._mcmc_out = runMCMC(self)

            self.mcmc_res = self._mcmc_out[0]
            self.mcmc_model = self._mcmc_out[1]

        if not self.nest_res:
            print "nest_lc() output reran"
            self._nest_out = run_nest(self)

            self.nest_res = self._nest_out[0]
            self.nest_model = self._nest_out[1]

# io functions
    @staticmethod
    def writeIDs(filename, ids):
        # is this bad practice?
        names='ids'
        id_table = Table([ids], names=(names,))
        id_table.write(filename, 'ids', append=True)

        return

    @staticmethod
    def readIDs(filename,):
        names='ids'
        id_read = Table.read(filename, names)
        id_list = list(id_read['ids'])

        return id_list

    def read_fits(self, filename, id):
        # fit_lc fit
        fit_read = Table.read(filename, id + '_MLEfit')

        fit_errors_dict = collections.OrderedDict()
        for colnames in fit_read.colnames:
            fit_errors_dict[colnames] = fit_read[colnames][0]

        fit_dict = fit_read.meta
        fit_dict['errors'] = fit_errors_dict

        fit_result = sncosmo.utils.Result(fit_dict)

        # mcmc_lc fit
        mcmc_read = Table.read(filename, id + '_mcmc')

        mcmc_errors_dict = collections.OrderedDict()
        for colnames in mcmc_read.colnames:
            mcmc_errors_dict[colnames] = mcmc_read[colnames][len(mcmc_read.columns[0]) - 1]

        mcmc_read.remove_row(len(mcmc_read.columns[0]) - 1)

        mcmc_samples = np.array([np.array(mcmc_read.columns[0]),
                         np.array(mcmc_read.columns[1]),
                         np.array(mcmc_read.columns[2]),
                         np.array(mcmc_read.columns[3])])

        mcmc_dict = mcmc_read.meta
        mcmc_dict['errors'] = mcmc_errors_dict
        mcmc_dict['samples'] = mcmc_samples.T

        mcmc_result = sncosmo.utils.Result(mcmc_dict)

        # nest_lc fit
        nest_read = Table.read(filename, id + '_nest')

        nest_param_dict = collections.OrderedDict()

        for colnames in nest_read.colnames:
            nest_param_dict[colnames] = nest_read[colnames][len(nest_read.columns[0]) - 1]

        nest_read.remove_row(len(nest_read.columns[0]) - 1)
        nest_read.remove_column('z')

        nest_errors_dict = collections.OrderedDict()

        for colnames in nest_read.colnames:
            nest_errors_dict[colnames] = nest_read[colnames][len(nest_read.columns[0]) - 1]

        nest_read.remove_row(len(nest_read.columns[0]) - 1)

        nest_samples = np.array([np.array(nest_read.columns[0]),
                         np.array(nest_read.columns[1]),
                         np.array(nest_read.columns[2]),
                         np.array(nest_read.columns[3])])

        nest_bounds = {}

        for colnames in nest_read.colnames:
            nest_bounds[colnames] = tuple(nest_read.meta[colnames])
            del nest_read.meta[colnames]

        nest_dict = nest_read.meta
        nest_dict['errors'] = nest_errors_dict
        nest_dict['param_dict'] = nest_param_dict
        nest_dict['samples'] = nest_samples.T
        nest_dict['bounds'] = nest_bounds

        nest_result = sncosmo.utils.Result(nest_dict)

        # now make new models instances for each fit
        count = np.arange(len(fit_result.parameters))

        fit_model_params = {}
        for number in count:
            fit_model_params[fit_result.param_names[number]] = fit_result.parameters[number]

        mcmc_model_params = {}
        for number in count:
            mcmc_model_params[mcmc_result.param_names[number]] = mcmc_result.parameters[number]

        nest_model_params = nest_result.param_dict

        fit_model = sncosmo.Model(source='salt2-extended')
        fit_model.set(**fit_model_params)

        mcmc_model = sncosmo.Model(source='salt2-extended')
        mcmc_model.set(**mcmc_model_params)

        nest_model = sncosmo.Model(source='salt2-extended')
        nest_model.set(**nest_model_params)

        fit_out = (fit_result, fit_model)
        mcmc_out = (mcmc_result, mcmc_model)
        nest_out = (nest_result, nest_model)

        return fit_out, mcmc_out, nest_out

    def write_fits(self, filename, id):
        # fit_lc fit
        if self.fit_res:
            fit_errors = self.fit_res.errors

            fit_table = Table()

            for keys in fit_errors:
                fit_table[keys] = [fit_errors[keys]]

            for key in self.fit_res.keys():
                if key == 'errors':
                    continue
                fit_table.meta[key] = self.fit_res[key]

            fit_table.write(filename, id + '_MLEfit', append=True)

        # mcmc_lc fit
        if self.mcmc_res:
            mcmc_errors = self.mcmc_res.errors
            mcmc_table = Table(self.mcmc_res.samples, names=self.mcmc_res.vparam_names)

            for key in self.mcmc_res.keys():
                if key == 'errors' or key =='samples':
                    continue
                mcmc_table.meta[key] = self.mcmc_res[key]

            mcmc_table.add_row(mcmc_errors.values())
            mcmc_table.write(filename, id + '_mcmc', append=True)

            # nest_lc fit
        if self.nest_res:
            nest_errors = self.nest_res.errors
            nest_param_dict = self.nest_res.param_dict
            nest_bounds = self.nest_res.bounds

            nest_table = Table(self.nest_res.samples, names=self.nest_res.vparam_names)

            for key in self.nest_res.keys():
                if key == 'errors' or key =='samples' or key =='param_dict' or key == 'bounds':
                    continue
                nest_table.meta[key] = self.nest_res[key]

            nest_table.add_row(nest_errors.values())

            temp_z = np.zeros((len(nest_table['t0'])))
            col_z = Table.Column(name='z', data=temp_z)
            nest_table.add_column(col_z)

            param_list = []
            for colname in nest_table.colnames:
                param_list.append(nest_param_dict[colname])

            nest_table.add_row(param_list)

            for key in nest_bounds.keys():
                nest_table.meta[key] = nest_bounds[key]

            nest_table.write(filename, id + '_nest', append=True)

        return

# visualization functions
    def plot_lightcurves(self, fits=True):
        data = self.data
        model = self.model
        fit_model = self.fit_model
        mcmc_model = self.mcmc_model
        nest_model = self.nest_model

        models = [model]
        model_names = ['model']
        if fits:
            if fit_model:
                models.append(fit_model)
                model_names.append('MLE')
            if mcmc_model:
                models.append(mcmc_model)
                model_names.append('MCMC')
            if nest_model:
                models.append(nest_model)
                model_names.append('nest')

        fig = sncosmo.plot_lc(data, model=models, model_label=model_names)

        return fig

    def plot_corner(self):
        model = self.model
        v_params = self.mcmc_res.vparam_names
        samples = self.mcmc_res.samples

        n_dim, n_samples = len(v_params), len(samples)

        # make figure

        figure = triangle.corner(samples, labels=[v_params[0], v_params[1], v_params[2], v_params[3]],
                     truths=[model.get(v_params[0]), model.get(v_params[1]),
                             model.get(v_params[2]), model.get(v_params[3])],
                     range=n_dim*[0.9999],
                     show_titles=True, title_args={"fontsize": 12}, hist_kwargs={'normed':True, 'color':'r'})

        figure.gca().annotate("mcmc sampling", xy=(0.5, 1.0), xycoords="figure fraction",
                  xytext=(0, -5), textcoords="offset points",
                  ha="center", va="top")


        mean = self.fit_res.parameters[0]
        variance = self.fit_res.covariance[0][0]
        sigma = math.sqrt(variance)
        min = self.fit_res.parameters[0] - 8 * sigma
        max = self.fit_res.parameters[0] + 8 * sigma
        x = np.linspace(min, max, 100)
        axes = figure.axes
        axes[0].plot(x,mlab.normpdf(x,mean,sigma))
        #axes[5].hist(samples[:,1], bins=20, normed=True)

    #    figure_nest = triangle.corner(nestSamples, labels=[nestVParams[0], nestVParams[1], nestVParams[2], nestVParams[3]],
    #                 truths=[model.get(nestVParams[0]), model.get(nestVParams[1]),
    #                         model.get(nestVParams[2]), model.get(nestVParams[3])],
    #                 weights=self.nest_res.weights, range=nest_ndim*[0.9999],
    #                 show_titles=True, title_args={"fontsize": 12}, hist_kwargs={'normed':True, 'color':'b'})

    #    figure.gca().annotate("nest sampling", xy=(0.5, 1.0), xycoords="figure fraction",
    #              xytext=(0, -5), textcoords="offset points",
    #              ha="center", va="top")

        return figure

    def plot_MLEGaussian(self):
            mean = 0
            variance = 1
            sigma = math.sqrt(variance)
            x = np.linspace(-3,3,100)
            plt.plot(x,mlab.normpdf(x,mean,sigma))

    def plot_trace(self):
        mcmc_vparams = self.mcmc_res.vparam_names
        mcmcSamples = self.mcmc_res.samples

        trace_fig = plt.figure(figsize=(20,8))

        mcmc1 = trace_fig.add_subplot(241)
        mcmc2 = trace_fig.add_subplot(242)
        mcmc3 = trace_fig.add_subplot(243)
        mcmc4 = trace_fig.add_subplot(244)

        mcmc1.plot(mcmcSamples[:,0])
        mcmc2.plot(mcmcSamples[:,1])
        mcmc3.plot(mcmcSamples[:,2])
        mcmc4.plot(mcmcSamples[:,3])

        mcmc1.set_title('mcmc: ' + mcmc_vparams[0])
        mcmc2.set_title('mcmc: ' + mcmc_vparams[1])
        mcmc3.set_title('mcmc: ' + mcmc_vparams[2])
        mcmc4.set_title('mcmc: ' + mcmc_vparams[3])

        trace_fig.tight_layout()

        return trace_fig


# fit statistics
    def metadata(self):
        print self.model

    def statistics(self):
        print "chi2"
        print "dof: ", self._fit_out[0].dof


    def comparefits2truth(self):
        pass

    def comparefits2fits(self):
        pass

    def calculateBias(LC):
        model = LC.model
        mcmcSamples = LC.mcmc_res.samples
        nestSamples = LC.nest_res.samples

        #fitBias = [mcmcSamples[0,x] - model.get() for x in ]
