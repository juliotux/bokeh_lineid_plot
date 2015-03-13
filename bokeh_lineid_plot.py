# bokeh_lineid_plot.py
#
# Copyright (C) 2015 - Julio Campagnolo <juliocampagnolo@gmail.com>
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

import numpy as np

__all__ = ['plot_line_ids']

def adjust_boxes(line_wave, box_widths, left_edge, right_edge,
                 max_iter=1000, adjust_factor=0.35,
                 factor_decrement=3.0, fd_p=0.75):
    """Ajdust given boxes so that they don't overlap.
    Parameters
    ----------
    line_wave: list or array of floats
        Line wave lengths. These are assumed to be the initial y (wave
        length) location of the boxes.
    box_widths: list or array of floats
        Width of box containing labels for each line identification.
    left_edge: float
        Left edge of valid data i.e., wave length minimum.
    right_edge: float
        Right edge of valid data i.e., wave lengths maximum.
    max_iter: int
        Maximum number of iterations to attempt.
    adjust_factor: float
        Gap between boxes are reduced or increased by this factor after
        each iteration.
    factor_decrement: float
        The `adjust_factor` itself if reduced by this factor, after
        certain number of iterations. This is useful for crowded
        regions.
    fd_p: float
        Percentage, given as a fraction between 0 and 1, after which
        adjust_factor must be reduced by a factor of
        `factor_decrement`. Default is set to 0.75.
    Returns
    -------
    wlp, niter, changed: (float, float, float)
        The new x (wave length) location of the text boxes, the number
        of iterations used and a flag to indicated whether any changes to
        the input locations were made or not.
    Notes
    -----
    This is a direct translation of the code in lineid_plot.pro file in
    NASA IDLAstro library.
    Positions are returned either when the boxes no longer overlap or
    when `max_iter` number of iterations are completed. So if there are
    many boxes, there is a possibility that the final box locations
    overlap.
    References
    ----------
    + http://idlastro.gsfc.nasa.gov/ftp/pro/plot/lineid_plot.pro
    + http://idlastro.gsfc.nasa.gov/
    """
    # Adjust positions.
    niter = 0
    changed = True
    nlines = len(line_wave)

    wlp = line_wave[:]
    while changed:
        changed = False
        for i in range(nlines):
            if i > 0:
                diff1 = wlp[i] - wlp[i - 1]
                separation1 = (box_widths[i] + box_widths[i - 1]) / 2.0
            else:
                diff1 = wlp[i] - left_edge + box_widths[i] * 1.01
                separation1 = box_widths[i]
            if i < nlines - 2:
                diff2 = wlp[i + 1] - wlp[i]
                separation2 = (box_widths[i] + box_widths[i + 1]) / 2.0
            else:
                diff2 = right_edge + box_widths[i] * 1.01 - wlp[i]
                separation2 = box_widths[i]

            if diff1 < separation1 or diff2 < separation2:
                if wlp[i] == left_edge: diff1 = 0
                if wlp[i] == right_edge: diff2 = 0
                if diff2 > diff1:
                    wlp[i] = wlp[i] + separation2 * adjust_factor
                    wlp[i] = wlp[i] if wlp[i] < right_edge else right_edge
                else:
                    wlp[i] = wlp[i] - separation1 * adjust_factor
                    wlp[i] = wlp[i] if wlp[i] > left_edge else left_edge
                changed = True
            niter += 1
        if niter == max_iter * fd_p: adjust_factor /= factor_decrement
        if niter >= max_iter: break

    return wlp, changed, niter

def plot_line_ids(plot, spec_wave, spec_flux, line_wave, line_label, **kwargs):
    '''
    Label features with automatic layout of labels.

    Parameters
    ----------
    plot: Bokeh.plotting.figure
        The plot where the features will be labeled.
    spec_wave: list or array of floats
        Wave lengths of data.
    spec_flux: list or array of floats
        Flux at each wavelength.
    line_wave: list or array of floats
        Wave length of features to be labelled.
    line_label: list of strings
        Label text for each line.

    kwargs: key value pairs
        All of these keywords are optional.
        The following keys are recognized:
            *All the keywords of the text and line properties used in Bokeh
            max_iter: int
                Maximum iterations to use. Default is set to 1000.
            label_loc: string or float
                The y location of the labels. If a string is given, it's value
                may be 'above' or 'bellow', indicating that the label will be
                located above or bellow the spectrum by a distance given by
                label_dist. If a float is given, all the labels will be placed
                in a constant y value.
            label_dist: float
                Separation between the spectrum flux and the label.
            extend: bool
                Draw a line linking the spectrum and the label.
            vert_space: float
                Space between the label and the init of the vertical line
            box_width: floats
                The separation between labels in units of the wavelenght.
            adjust_factor: float
                Gap between boxes are reduced or increased by this factor after
                each iteration.
            factor_decrement: float
                The `adjust_factor` itself if reduced by this factor, after
                certain number of iterations. This is useful for crowded
                regions.
    '''
    spec_wave = np.array(spec_wave)
    spec_flux = np.array(spec_flux)
    line_wave = np.array(line_wave)
    line_label = np.array(line_label)

    if not len(line_wave) == len(line_label):
        raise ValueError('the line_wave and line_label must have the same lengths.')

    indx = np.argsort(spec_wave)
    spec_wave[:] = spec_wave[indx]
    spec_flux[:] = spec_flux[indx]
    indx = np.argsort(line_wave)
    line_wave1 = line_wave[:] = line_wave[indx]
    line_label[:] = line_label[indx]

    label_loc = kwargs.pop('label_loc','above')
    label_dist = kwargs.pop('label_dist',1.0)
    vert_space = kwargs.pop('vert_space',0.5)

    line_flux = np.interp(line_wave, spec_wave, spec_flux)
    try:
        label_loc = float(label_loc)
        line_loc = np.array([label_loc]*len(line_wave))
    except:
        if label_loc == 'above':
            line_loc = line_flux + label_dist
        elif label_loc == 'below':
            line_loc = line_flux - label_dist
        else:
            raise ValueError('the label_loc parameter must be \'above\', \'below\' or a float.')
    line_flux1 = line_loc + vert_space*(line_flux-label_loc)/abs(line_flux-label_loc)

    max_iter = kwargs.get('max_iter', 1000)
    adjust_factor = kwargs.get('adjust_factor', 0.35)
    factor_decrement = kwargs.get('factor_decrement', 3.0)
    box_widths = kwargs.get('box_widths', 10)

    box_widths = [box_widths]*len(line_label)
    wlp, niter, changed = adjust_boxes(line_wave, box_widths,
                                       np.min(spec_wave), np.max(spec_wave),
                                       adjust_factor=adjust_factor,
                                       factor_decrement=factor_decrement,
                                       max_iter=max_iter)

    text_angle = kwargs.pop('text_angle', np.pi/2)
    text_font_size = kwargs.pop('text_font_size',None)
    text_align = kwargs.pop('text_align','left')
    text_baseline = kwargs.pop('text_baseline','bottom')

    line_color = kwargs.pop('line_color','black')
    
    plot.text(wlp, line_loc, text=line_label, text_angle=text_angle,
              text_font_size=text_font_size, text_align=text_align,
              text_baseline=text_baseline, **kwargs)
              
    if kwargs.pop('extend',True):
        for i in range(len(line_label)):
            plot.line([line_wave1[i], line_wave1[i], wlp[i]],
                      [line_flux[i], line_flux1[i], line_loc[i]],
                      line_color=line_color, **kwargs)