import tkinter as tk
from tkinter import ttk
import numpy as np
from iSLAT.COMPONENTS.GUI import Data_field

class FunctionsFrame:
    def __init__(self, parent):
        self.parent = parent
        self.frame = ttk.Frame(parent)
        self.frame.pack(fill=tk.BOTH, expand=True)

        # Create a label for the frame
        self.label = ttk.Label(self.frame, text="Functions Frame")
        self.label.pack(pady=10)

        # Add a button to the frame
        self.button = ttk.Button(self.frame, text="Click Me", command=self.on_button_click)
        self.button.pack(pady=10)

    def Save(self):
        """
        Save() is connected to the "Save Line" button of the tool.
        \nThis function appends information of the the strongest line (as determined by intensity) in the spanned area graph to a csv file. 
        \nThe name of the csv file is set with the "svd_line_file" variable in the second code block above. 
        \nFor the parameters of the line that is saved, refer to the "line2save" variable in onselect().
        \nWhen starting the tool up, the "headers" variable is set to False. After apending a line to the csv for the first time, the "headers" variable is changed to False.
        """
        '''global line2save
        global headers
        global selectedline
        global linesavepath'''

        # This section is necessary for refreshing the text feed area to the left of the tool
        data_field.delete ('1.0', "end")

        try:
            linesavepath
        except NameError:
            data_field.delete ('1.0', "end")
            data_field.insert ('1.0', 'Line save file is not defined!')
        else:
            if selectedline == True:  # "selectedline" variable is determined by whether or not an area was selected in the top graph or not

                line2save.to_csv (linesavepath, mode='a', index=False, header=False)

                data_field.insert ('1.0', 'Line Saved!')
                fig.canvas.draw_idle ()
            else:
                data_field.insert ('1.0', 'No Line Selected!')
                fig.canvas.draw_idle ()
                return

        canvas.draw()

    def fit_onselect(self):
        """fit_onselect() is connected to the "Fit Line" button of the tool. \nThis function fits the line selected in the top graph using LMFIT"""
        #global selectedline

        print(' ')
        print('Fitting line with LMFIT ...')

        if selectedline == True:  # "selectedline" variable is determined by whether or not an area was selected in the top graph or not

            # using one less pixel on each side here, because of how data_region_x is defined: to include 1 more pixel on each side
            gauss_fit, gauss_fwhm, gauss_area, x_fit = fit_line(data_region_x[1], data_region_x[-2])

            dely = gauss_fit.eval_uncertainty(sigma=3)
            ax2.fill_between(x_fit, gauss_fit.best_fit - dely, gauss_fit.best_fit + dely, color=user_settings["theme"]["uncertainty_band_color"],
                            label=r'3-$\sigma$ uncertainty band')
            ax2.plot(x_fit, gauss_fit.best_fit, label='Gauss. fit', color=user_settings["theme"]["selection_color"], ls='--')

            data_field.insert(tk.END, ('\n ' + '\nGaussian fit results: ' + '\nCentroid (μm) = ' + str (
                np.round(gauss_fit.params['center'].value, decimals=5)) + ' +/- ' + str (
                np.round(gauss_fit.params['center'].stderr, decimals=5)) + '\nFWHM (km/s) = ' + str (
                np.round(gauss_fwhm[0], decimals=1)) + ' +/- ' + str (
                np.round(gauss_fwhm[1], decimals=1)) + '\nArea (erg/s/cm2) = ' + f'{gauss_area[0]:.{3}e}' +
                                        ' +/- ' + f'{gauss_area[1]:.{3}e}'))

            fig.canvas.draw_idle()
        else:
            data_field.delete('1.0', "end")
            data_field.insert('1.0', 'No Line Selected!')
            fig.canvas.draw_idle()
            return
        canvas.draw()