import sncosmo
import triangle
from astropy.table import Table
import matplotlib.pyplot as plt
import numpy as np
import collections

class LC(object):
    """class to streamline light curve fits with SNCosmo """
    def __init__(self, model, data, vparams, bounds={'c':(-0.3, 0.3), 'x1':(-3.0, 3.0)}, truths=None):
        # super(, self).__init__()  <-- don't need this yet

        self.model = model
        self.data = data
        self.vparams = vparams
        self.bounds = bounds

        self._fitOut = None
        self.fitRes = None
        self.fitModel = None

        self._mcmcOut = None
        self.mcmcRes = None
        self.mcmcModel = None

        self._nestOut = None
        self.nestRes = None
        self.nestModel = None

        # variables for time statistics


    @property
    def fitOut(self):
        return self.fitOut

    @fitOut.getter
    def fitOut(self):
        if self._fitOut is None:
            print "running chi^2 fit"

            self._fitOut = sncosmo.fit_lc(self.data, self.model, vparam_names=self.vparams,
                                        bounds=self.bounds, minsnr=3.0)
            self.fitRes = self._fitOut[0]
            self.fitModel = self._fitOut[1]

        return self._fitOut

    @fitOut.setter
    def fitOut(self, val):
        self._fitOut = val

        if val:
            self.fitRes = self._fitOut[0]
            self.fitModel = self._fitOut[1]

        else:
            self.fitRes = None
            self.fitmode = None

        return self._fitOut

    @property
    def mcmcOut(self):
        return self.mcmcOut

    @mcmcOut.getter
    def mcmcOut(self):
        if self._mcmcOut is None:
            print "running MCMC fit"

            self._mcmcOut = sncosmo.mcmc_lc(self.data, self.model, vparam_names=self.vparams,
                                        bounds=self.bounds, minsnr=3.0)
            self.mcmcRes = self._mcmcOut[0]
            self.mcmcModel = self._mcmcOut[1]

        return self._mcmcOut

    @mcmcOut.setter
    def mcmcOut(self, val):
        self._mcmcOut = val

        if val:
            self.mcmcRes = self._mcmcOut[0]
            self.mcmcModel = self._mcmcOut[1]

        else:
            self.mcmcRes = None
            self.mcmcModel = None
        return self._mcmcOut

    @property
    def nestOut(self):
        return self.nestOut

    @nestOut.getter
    def nestOut(self):
        if self._nestOut is None:
            print "running nest fit"

            self._nestOut = sncosmo.nest_lc(self.data, self.model, vparam_names=self.vparams,
                                    bounds=self.bounds, guess_amplitude_bound=True,
                                     minsnr=3.0, verbose=True)
            self.nestRes = self._nestOut[0]
            self.nestModel = self._nestOut[1]

        return self._nestOut

    @nestOut.setter
    def nestOut(self, val):
        self._nestOut = val

        if val:
            self.nestRes = self._nestOut[0]
            self.nestModel = self._nestOut[1]

        else:
            self.nestRes = None
            self.nestModel = None

        return self._nestOut


#methods
#--------------

# functions for changing various parameters
    def reset(fit):
        pass

    def reset_all(self):
        self._fitOut = None;
        self._mcmcOut = None;
        self._nestOut = None;

        return

    def setBounds(self, new_bounds):
        self.bounds = new_bounds

        # resets fits now that bounds have changed
        self._fitOut = None;
        self._mcmcOut = None;
        self._nestOut = None;

        return

    def rerunFits(self):
        if not self.fitRes:
            print "fit_lc() output reran"
            self.fitRes = self._fitOut[0]
            self.fitModel = self._fitOut[1]

        if not self.mcmcRes:
            print "mcmc_lc() output reran"
            self.mcmcRes = self._mcmcOut[0]
            self.mcmcModel = self._mcmcOut[1]

        if not self.nestRes:
            print "nest_lc() output reran"
            self.nestRes = self._nestOut[0]
            self.nestModel = self._nestOut[1]

# io functions
    def readFits(self, filename, id):
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

        fitmodel_params = {}
        for number in count:
            fitmodel_params[fit_result.param_names[number]] = fit_result.parameters[number]

        mcmcmodel_params = {}
        for number in count:
            mcmcmodel_params[mcmc_result.param_names[number]] = mcmc_result.parameters[number]

        nestmodel_params = nest_result.param_dict

        fitmodel = sncosmo.Model(source='salt2-extended')
        fitmodel.set(**fitmodel_params)

        mcmcmodel = sncosmo.Model(source='salt2-extended')
        mcmcmodel.set(**mcmcmodel_params)

        nestmodel = sncosmo.Model(source='salt2-extended')
        nestmodel.set(**nestmodel_params)

        fitOut = (fit_result, fitmodel)
        mcmcOut = (mcmc_result, mcmcmodel)
        nestOut = (nest_result, nestmodel)

        return fitOut, mcmcOut, nestOut

    def writeFits(self, filename, id):
        # fit_lc fit
        if self.fitRes:
            fit_errors = self.fitRes.errors

            fit_table = Table()

            for keys in fit_errors:
                fit_table[keys] = [fit_errors[keys]]

            for key in self.fitRes.keys():
                if key == 'errors':
                    continue
                fit_table.meta[key] = self.fitRes[key]

            fit_table.write(filename, id + '_MLEfit', append=True)

        # mcmc_lc fit
        if self.mcmcRes:
            mcmc_errors = self.mcmcRes.errors
            mcmc_table = Table(self.mcmcRes.samples, names=self.mcmcRes.vparam_names)

            for key in self.mcmcRes.keys():
                if key == 'errors' or key =='samples':
                    continue
                mcmc_table.meta[key] = self.mcmcRes[key]

            mcmc_table.add_row(mcmc_errors.values())
            mcmc_table.write(filename, id + '_mcmc', append=True)

            # nest_lc fit
        if self.nestRes:
            nest_errors = self.nestRes.errors
            nest_param_dict = self.nestRes.param_dict
            nest_bounds = self.nestRes.bounds

            nest_table = Table(self.nestRes.samples, names=self.nestRes.vparam_names)

            for key in self.nestRes.keys():
                if key == 'errors' or key =='samples' or key =='param_dict' or key == 'bounds':
                    continue
                nest_table.meta[key] = self.nestRes[key]

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
    def plotLC(self, fits=True):
        data = self.data
        model = self.model
        fit_model = self.fitModel
        mcmc_model = self.mcmcModel
        nest_model = self.nestModel

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

    def plotCorner(self):
        model = self.model

        mcmcVParams = self.mcmcRes.vparam_names
        nestVParams = self.nestRes.vparam_names

        mcmcSamples = self.mcmcRes.samples
        nestSamples = self.nestRes.samples

        mcmc_ndim, mcmc_nsamples = len(mcmcVParams), len(mcmcSamples)
        nest_ndim, nest_nsamples = len(nestVParams), len(nestSamples)

        # make figure

        figure_mcmc = triangle.corner(mcmcSamples, labels=[mcmcVParams[0], mcmcVParams[1], mcmcVParams[2], mcmcVParams[3]],
                     truths=[model.get(mcmcVParams[0]), model.get(mcmcVParams[1]),
                             model.get(mcmcVParams[2]), model.get(mcmcVParams[3])],
                     range=mcmc_ndim*[0.9999],
                     show_titles=True, title_args={"fontsize": 12})

        figure_mcmc.gca().annotate("mcmc sampling", xy=(0.5, 1.0), xycoords="figure fraction",
                  xytext=(0, -5), textcoords="offset points",
                  ha="center", va="top")

        figure_nest = triangle.corner(nestSamples, labels=[nestVParams[0], nestVParams[1], nestVParams[2], nestVParams[3]],
                     truths=[model.get(nestVParams[0]), model.get(nestVParams[1]),
                             model.get(nestVParams[2]), model.get(nestVParams[3])],
                     weights=self.nestRes.weights, range=nest_ndim*[0.9999],
                     show_titles=True, title_args={"fontsize": 12})

        figure_nest.gca().annotate("nest sampling", xy=(0.5, 1.0), xycoords="figure fraction",
                  xytext=(0, -5), textcoords="offset points",
                  ha="center", va="top")

        return figure_mcmc, figure_nest

    def plotTrace(self):
        mcmc_vparams = self.mcmcRes.vparam_names
        nest_vparams = self.nestRes.vparam_names

        mcmcSamples = self.mcmcRes.samples
        nestSamples = self.nestRes.samples

        trace_fig = plt.figure(figsize=(20,8))

        mcmc1 = trace_fig.add_subplot(241)
        mcmc2 = trace_fig.add_subplot(242)
        mcmc3 = trace_fig.add_subplot(243)
        mcmc4 = trace_fig.add_subplot(244)

        nest1 = trace_fig.add_subplot(245)
        nest2 = trace_fig.add_subplot(246)
        nest3 = trace_fig.add_subplot(247)
        nest4 = trace_fig.add_subplot(248)

        mcmc1.plot(mcmcSamples[:,0])
        mcmc2.plot(mcmcSamples[:,1])
        mcmc3.plot(mcmcSamples[:,2])
        mcmc4.plot(mcmcSamples[:,3])

        mcmc1.set_title('mcmc: ' + mcmc_vparams[0])
        mcmc2.set_title('mcmc: ' + mcmc_vparams[1])
        mcmc3.set_title('mcmc: ' + mcmc_vparams[2])
        mcmc4.set_title('mcmc: ' + mcmc_vparams[3])

        nest1.plot(nestSamples[:,0])
        nest2.plot(nestSamples[:,1])
        nest3.plot(nestSamples[:,2])
        nest4.plot(nestSamples[:,3])

        nest1.set_title('nest: ' + nest_vparams[0])
        nest2.set_title('nest: ' + nest_vparams[1])
        nest3.set_title('nest: ' + nest_vparams[2])
        nest4.set_title('nest: ' + nest_vparams[3])

        trace_fig.tight_layout()

        return trace_fig


# fit statistics
    def metadata(self):
        print self.model

    def statistics(self):
        print "chi2"
        print "dof: ", self._fitOut[0].dof


    def comparefits2truth(self):
        pass

    def comparefits2fits(self):
        pass

    def calculateBias(LC):
        model = LC.model
        mcmcSamples = LC.mcmcRes.samples
        nestSamples = LC.nestRes.samples

        #fitBias = [mcmcSamples[0,x] - model.get() for x in ]
