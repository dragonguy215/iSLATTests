import tkinter as tk
import numpy as np

def Save(data_field, linesavepath, selectedline, line2save):
    """
    Save() is connected to the "Save Line" button of the tool.
    \nThis function appends information of the strongest line (as determined by intensity) in the spanned area graph to a csv file. 
    \nThe name of the csv file is set with the "linesavepath" variable. 
    \nFor the parameters of the line that is saved, refer to the "line2save" variable.
    """
    messages = []

    if not linesavepath:
        messages.append('Line save file is not defined!')
    else:
        if selectedline:  # "selectedline" variable determines whether an area was selected in the top graph
            line2save.to_csv(linesavepath, mode='a', index=False, header=False)
            messages.append('Line Saved!')
        else:
            messages.append('No Line Selected!')
            return messages

    return messages

def fit_onselect(data_field, selectedline, data_region_x, fit_line, user_settings):
    """
    fit_onselect() is connected to the "Fit Line" button of the tool. 
    \nThis function fits the line selected in the top graph using LMFIT.
    """
    print(' ')
    print('Fitting line with LMFIT ...')

    if selectedline:  # "selectedline" variable determines whether an area was selected in the top graph
        # Fit the line using one less pixel on each side
        gauss_fit, gauss_fwhm, gauss_area, x_fit = fit_line(data_region_x[1], data_region_x[-2])

        dely = gauss_fit.eval_uncertainty(sigma=3)
        uncertainty_band = {
            "x": x_fit,
            "y_lower": gauss_fit.best_fit - dely,
            "y_upper": gauss_fit.best_fit + dely,
            "color": user_settings["theme"]["uncertainty_band_color"],
            "label": r'3-$\sigma$ uncertainty band'
        }
        fit_line_plot = {
            "x": x_fit,
            "y": gauss_fit.best_fit,
            "color": user_settings["theme"]["selection_color"],
            "label": 'Gauss. fit',
            "linestyle": '--'
        }

        fit_results = {
            "centroid": (np.round(gauss_fit.params['center'].value, decimals=5),
                         np.round(gauss_fit.params['center'].stderr, decimals=5)),
            "fwhm": (np.round(gauss_fwhm[0], decimals=1), np.round(gauss_fwhm[1], decimals=1)),
            "area": (f'{gauss_area[0]:.{3}e}', f'{gauss_area[1]:.{3}e}')
        }

        return {"uncertainty_band": uncertainty_band, "fit_line_plot": fit_line_plot, "fit_results": fit_results}
    else:
        return {"error": 'No Line Selected!'}
