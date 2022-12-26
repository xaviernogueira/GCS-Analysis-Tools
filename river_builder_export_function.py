import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from typing import Tuple


def slope(
    point1: Tuple[float, float],
    point2: Tuple[float, float],
) -> float:
    """Return slope of two points."""

    if point1[0] != point2[0]:
        return float(
            (point2[1] - point1[1]) / (point2[0] - point1[0])
        )
    elif point2[1] > point1[1]:
        return math.inf
    elif point2[1] < point1[1]:
        return -math.inf
    else:
        return 0


def slope_v(x_v, y_v):
    """Return an array of slope.
        si = (y(i+1) - y(i))/(x(i+1) - x(i))
        s(-1) will be the same as s(-2) to make it the same length.
    Inputs:
    x_v - array. x values.
    y_v - array. y values.
    Output:
    s_v - array. slopes.
    """
    x1_v = x_v[:-1]
    x2_v = x_v[1:]

    y1_v = y_v[:-1]
    y2_v = y_v[1:]

    fun = np.vectorize(slope)
    s_v = fun(x1_v, y1_v, x2_v, y2_v)
    s_v = np.append(s_v, s_v[-1:])
    return s_v


def ifft_out(signal, fft, ifft_df, n, spacing):
    cos_coefs = []
    sin_coefs = []
    freq_list = []  # Frequency in cycles per reach!
    amp_list = []  # Amplitude in lnegth units
    phase_list = []

    fft_freqs = np.fft.fftfreq(fft.size, spacing)
    ifft = np.fft.ifft(fft).real
    reach_length = signal.size * spacing
    for index, i in enumerate(fft):
        if i != 0.0:
            cos_coefs.append(i.real)
            sin_coefs.append(i.imag)
            temp_fft = np.fft.fft(signal)
            np.put(temp_fft, range(index + 2, len(temp_fft)), 0.0)
            np.put(temp_fft, range(0, index + 1), 0.0)
            temp_ifft = np.fft.ifft(temp_fft).real
            if n == 1:
                ifft = temp_ifft
            amp = (np.amax(temp_ifft) - np.amin(temp_ifft)) / 2
            amp_list.append(amp)

            slope = slope_v(np.arange(0, signal.size), temp_ifft)
            slope_pre = slope[0:-1]
            slope_post = slope[1:]
            slope_change = np.multiply(slope_pre, slope_post)

            frequency = len(np.where(slope_change <= 0)[0])/2
            freq_list.append(frequency)

            ifft_df['harmonic_%s' % (index + 1)] = temp_ifft

            sub_index = 0
            # Finds when the single FFT component IFFT crossed 0, and therefore it's phase in length units
            if temp_ifft[0] < 0:
                while temp_ifft[sub_index] < 0:
                    sub_index += 1
                phase = -1 * sub_index * spacing
            elif temp_ifft[0] > 0:
                while temp_ifft[sub_index] > 0:
                    sub_index += 1
                phase = sub_index * spacing

            phase_list.append(phase)

    return [sin_coefs, cos_coefs, freq_list, amp_list, phase_list, ifft_df, ifft]


def by_fft(signal, n, spacing):
    ifft_df = pd.DataFrame()
    ifft_df['raw_series'] = signal
    fft = np.fft.fft(signal)
    if n == 0:
        ifft = np.fft.ifft(fft).real

    np.put(fft, range(n, len(fft)), 0.0)
    fft_freqs = np.fft.fftfreq(signal.size, spacing)

    out_list = ifft_out(signal, fft, ifft_df, n, spacing)
    freqs = out_list[2]
    amps = out_list[3]
    phases = out_list[4]
    ifft = out_list[-1]
    ifft_df = out_list[-2]
    ifft_df['all_%s_harmonics' % n] = ifft

    return [ifft, n, out_list[1], out_list[2], ifft_df, freqs, amps, phases]


def by_power(signal, n, spacing):
    ifft_df = pd.DataFrame()
    ifft_df['raw_series'] = signal
    fft = np.fft.fft(signal)
    if n == 0:
        ifft = np.fft.ifft(fft).real

    psd = np.abs(fft) ** 2
    indices = np.argsort(psd).tolist()
    n_indices = indices[:-n]
    np.put(fft, n_indices, 0.0)
    fft_freqs = np.fft.fftfreq(signal.size, spacing)

    out_list = ifft_out(signal, fft, ifft_df, n, spacing)
    freqs = out_list[2]
    amps = out_list[3]
    phases = out_list[4]
    ifft = out_list[-1]
    ifft_df = out_list[-2]
    ifft_df['all_%s_harmonics' % n] = ifft

    return [ifft, n, out_list[1], out_list[2], ifft_df, freqs, amps, phases]


def by_power_binned(signal, n, spacing):
    ifft_df = pd.DataFrame()
    ifft_df['raw_series'] = signal
    fft = np.fft.fft(signal)
    if n == 0:
        ifft = np.fft.ifft(fft).real

    psd = (np.abs(fft) ** 2).tolist()
    indices = []

    avg = (len(psd) / 2) / float(n)
    bins_list = []  # Stores n sub-lists of FFT components from which PSD is calculated
    last = 0.0

    while last < (len(psd) / 2):
        bins_list.append(psd[int(last):int(last + avg)])
        last += avg

    add = 0
    # Test this to make sure the add thing works out
    for i, sub_list in enumerate(bins_list):
        max_sub_index = sub_list.index(np.max(sub_list))
        indices.append(max_sub_index + add)
        add += len(sub_list)

    full_indices = np.argsort(psd).tolist()
    replace_indices = [i for i in full_indices if i not in indices]

    np.put(fft, replace_indices, 0.0)
    fft_freqs = np.fft.fftfreq(signal.size, spacing)

    out_list = ifft_out(signal, fft, ifft_df, n, spacing)
    freqs = out_list[2]
    amps = out_list[3]
    phases = out_list[4]
    ifft = out_list[-1]
    ifft_df = out_list[-2]
    ifft_df['all_%s_harmonics' % n] = ifft

    return [ifft, n, out_list[1], out_list[2], ifft_df, freqs, amps, phases]


def river_builder_harmonics(in_csv, index_field, units='', fields=[], field_names=[], r_2=0.95, n=0, methods='ALL'):
    """This function plots a N number of Fourier coefficients reconstrution of input signals. Exports coefficients to csv or text file.
    in_csv= A csv file location with evenly spaced values (string).
    sort_by (optional) allows an unsorted csv to be sorted by a input index field header (string)
    out_folder= Folder to which plots and exported coefficients are saved (string).
    index_field is the csv header corresponding to the centerline position, the units parameter can be any length unit and is used strictly for plotting (empty is default).
    fields= A list of csv headers from which signals are plotted and reconstructed (list of strings)
    field_names (optional) if specified must be a list of strings with names for plotting titles
    R_2= R^2 threshold for signal reconstruction (float). 0.95 is default.
    n (0 is default) if not 0 and an int allows a number of Fourier components to specified (int), as opposed to the standard R^2 threshold (default)
    by_power (False is default) if True takes the N highest power components first
    by_bins (False is default) if True splits the FFT components into N bins, and selects the highest power frequency from each
    to_riverbuilder (False"""

    in_df = pd.read_csv(in_csv, engine='python')
    out_folder = os.path.dirname(in_csv)
    print('CSV imported...')

    fields = [i for i in in_df.columns.values.tolist() if i != index_field]

    try:
        in_df.sort_values(index_field, inplace=True)
        index_array = in_df.loc[:, [index_field]].squeeze()
        spacing = float(index_array[1] - index_array[0])
    except:
        print('Could not sort values by  the input index field header: %s. Please either remove sort_by parameter, or correct the input field header.' % index_field)
        sys.exit()

    if methods == 'ALL':
        # Each list associated with each method stores [ifft, n, sin_coefs, cos_coefs, ifft_df]
        methods_dict = {'by_fft': [], 'by_power': [], 'by_power_binned': []}
    else:
        methods_dict = {methods: []}

    for count, field in enumerate(fields):
        try:
            field_signal = in_df.loc[:, [str(field)]].squeeze()
        except:
            print(
                'Error! Could not use csv field headers input. Please check csv headers.' % field)

        if len(field_names) == len(fields):
            field_name = field_names[count]
        elif count == 0:
            field_names = []
            field_names.append(field)
            field_name = field_names[count]
        else:
            field_names.append(field)
            field_name = field_names[count]

        if n == 0 and r_2 > 0:
            for method in methods_dict.keys():
                in_list = []
                if method == 'by_fft':
                    for i in range(1, len(field_signal)):
                        out_list = by_fft(field_signal, i, spacing)
                        temp_r2 = np.corrcoef(field_signal, out_list[0])[
                            0][1] ** 2
                        if temp_r2 >= r_2:
                            for out in out_list:
                                in_list.append(out)
                            break

                if method == 'by_power':
                    for i in range(1, len(field_signal)):
                        out_list = by_power(field_signal, i, spacing)
                        temp_r2 = np.corrcoef(field_signal, out_list[0])[
                            0][1] ** 2
                        if temp_r2 >= r_2:
                            for out in out_list:
                                in_list.append(out)
                            break

                if method == 'by_power_binned':
                    for i in range(1, len(field_signal)):
                        out_list = by_power_binned(field_signal, i, spacing)
                        temp_r2 = np.corrcoef(field_signal, out_list[0])[
                            0][1] ** 2
                        if temp_r2 >= r_2:
                            for out in out_list:
                                in_list.append(out)
                            break

                methods_dict[method].append(in_list)

        else:
            in_list = []
            for method in methods_dict.keys():
                if method == 'by_fft':
                    out_list = by_fft(field_signal, n, spacing)
                    temp_r2 = np.corrcoef(field_signal, out_list[0])[0][1] ** 2
                    for out in out_list:
                        in_list.append(out)

                if method == 'by_power':
                    out_list = by_power(field_signal, n, spacing)
                    temp_r2 = np.corrcoef(field_signal, out_list[0])[0][1] ** 2
                    for out in out_list:
                        in_list.append(out)

                if method == 'by_power':
                    out_list = by_power(field_signal, n, spacing)
                    temp_r2 = np.corrcoef(field_signal, out_list[0])[0][1] ** 2
                    for out in out_list:
                        in_list.append(out)

                methods_dict[method].append(in_list)

    for method in methods_dict.keys():
        print('Completing %s analysis...' % method)
        # Stores data for each field within one calculation method
        method_list = methods_dict[method]
        for count, field in enumerate(fields):
            text_file = open(
                out_folder + '\\%s_%s_to_riverbuilder.txt' % (field, method), 'w+')
            field_name = field_names[count]

            list = method_list[count]
            ifft_df = list[4]
            ifft_df.to_csv(out_folder + '\\%s_harmonics_%s.csv' %
                           (field, method))

            plt.plot(index_array, in_df.loc[:, str(
                field)].squeeze(), color='blue', label='Signal')
            plt.plot(index_array, list[0], color='red',
                     linestyle='--', label='Reconstructed signal')

            if units != '':
                add_units = 'in %s' % units
            else:
                add_units = ''

            plt.xlabel('Distance along centerline %s' % add_units)
            plt.ylabel('Value')
            plt.title('%s, %s method, N=%s component harmonic reconstruction' % (
                field_name, method, list[1]))
            plt.grid(b=True, which='major', color='#666666', linestyle='-')
            plt.minorticks_on()
            plt.legend(loc='lower center')

            fig_title = out_folder + '\\%s_%s_plot.png' % (field, method)
            fig = plt.gcf()
            fig.set_size_inches(12, 6)
            plt.savefig(fig_title, dpi=300, bbox_inches='tight')
            plt.cla()

            for num, amp in enumerate(list[-2]):
                if amp != 0.0:
                    # Writes in the form of COS#=(a, f, ps, MASK0) for river builder inputs
                    text_file.write('COS%s=(%s, %s, %s, MASK0)\n' % (
                        num, amp, list[-3][num], list[-1][num]))
            text_file.close()

    print('Analysis complete. Results @ %s' % out_folder)
