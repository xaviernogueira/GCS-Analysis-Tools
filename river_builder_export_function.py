import os
import math
import logging
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from typing import Union, List, Optional


def slope(
    x1: Union[float, int],
    y1: Union[float, int],
    x2: Union[float, int],
    y2: Union[float, int],
) -> float:
    """Return slope of two points."""
    if x1 != x2:
        return float((y2 - y1) / (x2 - x1))
    elif y2 > y1:
        return math.inf
    elif y2 < y1:
        return -math.inf
    else:
        return float(0)


def slope_v(
    x_v: np.array,
    y_v: np.array,
) -> np.array:
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


def ifft_out(
    signal: pd.Series,
    fft: np.ndarray,
    ifft_df: pd.DataFrame,
    n_harmonics: int,
    spacing: int,
) -> dict:

    # init lists to store output data
    cos_coefs = []
    sin_coefs = []
    freq_list = []  # Frequency in cycles per reach!
    amp_list = []  # Amplitude in lnegth units
    phase_list = []

    #fft_freqs = np.fft.fftfreq(fft.size, spacing)
    #reach_length = signal.size * spacing

    ifft = np.fft.ifft(fft).real

    for index, i in enumerate(fft):
        if i != 0.0:
            cos_coefs.append(i.real)
            sin_coefs.append(i.imag)
            temp_fft = np.fft.fft(signal)
            np.put(temp_fft, range(index + 2, len(temp_fft)), 0.0)
            np.put(temp_fft, range(0, index + 1), 0.0)
            temp_ifft = np.fft.ifft(temp_fft).real
            if n_harmonics == 1:
                ifft = temp_ifft
            amp = (np.amax(temp_ifft) - np.amin(temp_ifft)) / 2
            amp_list.append(amp)

            slope = slope_v(np.arange(0, signal.size), temp_ifft)
            slope_pre = slope[0:-1]
            slope_post = slope[1:]
            slope_change = np.multiply(slope_pre, slope_post)

            frequency = len(np.where(slope_change <= 0)[0])/2
            freq_list.append(frequency)

            ifft_df = ifft_df.copy()
            ifft_df['harmonic_%s' % (index + 1)] = temp_ifft

            sub_index = 0

            # Finds when the single FFT component IFFT crossed 0 / it's phase in length units
            if temp_ifft[0] < 0:
                while temp_ifft[sub_index] < 0:
                    sub_index += 1
                phase = -1 * sub_index * spacing
            elif temp_ifft[0] > 0:
                while temp_ifft[sub_index] > 0:
                    sub_index += 1
                phase = sub_index * spacing

            phase_list.append(phase)

    return {
        'sine_coefs': sin_coefs,
        'cosine_coefs': cos_coefs,
        'frequencies': freq_list,
        'amplitudes': amp_list,
        'phases': phase_list,
        'inverse_FFT_df': ifft_df,
        'inverse_FFT': ifft,
    }


def by_fft(
    signal: pd.Series,
    n_harmonics: int,
    spacing: int,
) -> dict:

    ifft_df = pd.DataFrame()
    ifft_df['raw_series'] = signal
    fft = np.fft.fft(signal)

    if n_harmonics == 0:
        ifft = np.fft.ifft(fft).real

    np.put(fft, range(n_harmonics, len(fft)), 0.0)
    #fft_freqs = np.fft.fftfreq(signal.size, spacing)

    out_dict = ifft_out(signal, fft, ifft_df, n_harmonics, spacing)

    # update the inverse_FFT_df
    ifft = out_dict['inverse_FFT']
    ifft_df = out_dict['inverse_FFT_df']
    ifft_df['all_%s_harmonics' % n_harmonics] = ifft

    # old version for ref: [ifft, n_harmonics, out_list[1], out_list[2], ifft_df, freqs, amps, phases]

    return {
        'inverse_FFT': out_dict['inverse_FFT'],
        'n_harmonics': n_harmonics,
        'sine_coefs': out_dict['sine_coefs'],
        'cosine_coefs': out_dict['cosine_coefs'],
        'inverse_FFT_df': ifft_df,
        'frequencies': out_dict['frequencies'],
        'amplitudes': out_dict['amplitudes'],
        'phases': out_dict['phases'],
    }


def by_power(
    signal: pd.Series,
    n_harmonics: int,
    spacing: int,
) -> dict:
    ifft_df = pd.DataFrame()
    ifft_df['raw_series'] = signal
    fft = np.fft.fft(signal)
    if n_harmonics == 0:
        ifft = np.fft.ifft(fft).real

    psd = np.abs(fft) ** 2
    indices = np.argsort(psd).tolist()
    n_indices = indices[:-n_harmonics]
    np.put(fft, n_indices, 0.0)

    #fft_freqs = np.fft.fftfreq(signal.size, spacing)

    out_dict = ifft_out(signal, fft, ifft_df, n_harmonics, spacing)

    # update the inverse_FFT_df
    ifft = out_dict['inverse_FFT']
    ifft_df = out_dict['inverse_FFT_df']
    ifft_df['all_%s_harmonics' % n_harmonics] = ifft

    return {
        'inverse_FFT': out_dict['inverse_FFT'],
        'n_harmonics': n_harmonics,
        'sine_coefs': out_dict['sine_coefs'],
        'cosine_coefs': out_dict['cosine_coefs'],
        'inverse_FFT_df': ifft_df,
        'frequencies': out_dict['frequencies'],
        'amplitudes': out_dict['amplitudes'],
        'phases': out_dict['phases'],
    }


def by_power_binned(
    signal: pd.Series,
    n_harmonics: int,
    spacing: int,
) -> dict:

    ifft_df = pd.DataFrame()
    ifft_df['raw_series'] = signal
    fft = np.fft.fft(signal)
    if n_harmonics == 0:
        ifft = np.fft.ifft(fft).real

    psd = (np.abs(fft) ** 2).tolist()
    indices = []

    avg = (len(psd) / 2) / float(n_harmonics)
    bins_list = []  # Stores n sub-lists of FFT components from which PSD is calculated
    last = 0.0

    while last < (len(psd) / 2):
        bins_list.append(psd[int(last):int(last + avg)])
        last += avg

    add = 0

    for sub_list in bins_list:
        max_sub_index = sub_list.index(np.max(sub_list))
        indices.append(max_sub_index + add)
        add += len(sub_list)

    full_indices = np.argsort(psd).tolist()
    replace_indices = [i for i in full_indices if i not in indices]

    np.put(fft, replace_indices, 0.0)

    #fft_freqs = np.fft.fftfreq(signal.size, spacing)

    out_dict = ifft_out(signal, fft, ifft_df, n_harmonics, spacing)

    # update the inverse_FFT_df
    ifft = out_dict['inverse_FFT']
    ifft_df = out_dict['inverse_FFT_df'].copy()
    ifft_df['all_%s_harmonics' % n_harmonics] = ifft

    return {
        'inverse_FFT': out_dict['inverse_FFT'],
        'n_harmonics': n_harmonics,
        'sine_coefs': out_dict['sine_coefs'],
        'cosine_coefs': out_dict['cosine_coefs'],
        'inverse_FFT_df': ifft_df,
        'frequencies': out_dict['frequencies'],
        'amplitudes': out_dict['amplitudes'],
        'phases': out_dict['phases'],
    }


def river_builder_harmonics(
    in_csv: str,
    index_field: str,
    units: str = '',
    r2_threshold: float = 0.95,
    n_harmonics: Optional[int] = None,
    methods: str = 'ALL',
    field_headers: Optional[List[str]] = None,
) -> str:
    """Generates plots and input csv/txt files to be imported into RiverBuilder.

    # of Fourier coefs reconstrution from input signals and to a csv or text file.
    This function plots a N

    Args:
        :param in_csv: A csv file location path with evenly spaced values.
        :param index_field: The csv header corresponding to the centerline position.
        :param units: Units name for labeling (i.e. m).
        :param r2_threshold: R-squared threshold to cut off signal reconstruction.
        :param n: (overrides r2_threshold) Number of harmonics to cut off signal reconstruction.
        :param methods: Which methods to use (ALL by default), one can also select
            by_fft, by_power, or by_binned.
            by_fft = adds harmonics in order of the FFT algo.
            by_power = adds the N highest power harmonic components first.
            by_bins = splits the FFT components into N bins, and selects the highest power frequency from each.

    Returns: path to the output folder.
    """

    in_df = pd.read_csv(in_csv, engine='python')

    if field_headers is None:
        fields = ['W', 'Z']
    else:
        fields = field_headers

    if len([i for i in list(in_df.columns) if i in fields]) != len(fields):
        raise KeyError(
            f'Could not find fields={fields} in {in_csv}!'
        )

    # make output directory
    out_folder = os.path.dirname(in_csv) + '//River_Builder_inputs'
    if not os.path.exists(out_folder):
        os.mkdir(out_folder)

    logging.info('CSV imported...')

    try:
        in_df.sort_values(index_field, inplace=True)
        index_array = in_df.loc[:, [index_field]].squeeze()
        spacing = float(index_array[1] - index_array[0])
    except KeyError:
        raise KeyError(
            f'Could not sort values by  the input index field header: {index_field}. '
            'Please either remove sort_by parameter, or correct the input field header.'
        )

    # init a dict to store dicts w/ field as the key, and data in a dict as values
    if methods == 'ALL':
        methods_dict = {
            'by_fft': {},
            'by_power': {},
            'by_power_binned': {},
        }
    else:
        methods_dict = {methods: {}}

    for field in fields:
        field_signal = in_df.loc[:, [str(field)]].squeeze()

        # conduct analysis using a R-squared cut off
        if n_harmonics is None and r2_threshold > 0:

            for method in methods_dict.keys():

                if method == 'by_fft':
                    logging.info(
                        f'Applying by_fft method w/ R-Squared threshold={r2_threshold}...')
                    for i in range(1, len(field_signal)):
                        out_dict = by_fft(field_signal, i, spacing)

                        temp_r2 = np.corrcoef(
                            field_signal,
                            out_dict['inverse_FFT'],
                        )[0][1] ** 2

                        if temp_r2 >= r2_threshold:
                            methods_dict[method][field] = out_dict.copy()
                            break

                if method == 'by_power':
                    logging.info(
                        f'Applying by_power method w/ R-Squared threshold={r2_threshold}...')
                    for i in range(1, len(field_signal)):
                        out_dict = by_power(field_signal, i, spacing)

                        temp_r2 = np.corrcoef(
                            field_signal,
                            out_dict['inverse_FFT'],
                        )[0][1] ** 2

                        if temp_r2 >= r2_threshold:
                            methods_dict[method][field] = out_dict.copy()
                            break

                if method == 'by_power_binned':
                    logging.info(
                        f'Applying by_power_binned method w/ R-Squared threshold={r2_threshold}...')
                    logging.info(
                        'NOTE: This method can take a while using a R2 threshold! '
                        'We recomend only using this w/ a N-Harmonics threshold.'
                    )
                    for i in range(1, len(field_signal)):
                        out_dict = by_power_binned(field_signal, i, spacing)

                        temp_r2 = np.corrcoef(
                            field_signal,
                            out_dict['inverse_FFT'],
                        )[0][1] ** 2

                        if temp_r2 >= r2_threshold:
                            methods_dict[method][field] = out_dict.copy()
                            break

        # conduct analysis using hard N-harmonics cut off
        elif n_harmonics is not None:

            for method in methods_dict.keys():
                if method == 'by_fft':
                    logging.info(
                        f'Applying by_fft method w/ a N-Harmonics threshold={n_harmonics}...')
                    out_dict = by_fft(field_signal, n_harmonics, spacing)

                    temp_r2 = np.corrcoef(
                        field_signal,
                        out_dict['inverse_FFT'],
                    )[0][1] ** 2

                if method == 'by_power':
                    logging.info(
                        f'Applying by_power method w/ a N-Harmonics threshold={n_harmonics}...')
                    out_dict = by_power(field_signal, n_harmonics, spacing)

                    temp_r2 = np.corrcoef(
                        field_signal,
                        out_dict['inverse_FFT'],
                    )[0][1] ** 2

                if method == 'by_power_binned':
                    logging.info(
                        f'Applying by_power_binned method w/ a N-Harmonics threshold={n_harmonics}...')
                    logging.info('NOTE: This method can take a while!')
                    out_dict = by_power_binned(
                        field_signal, n_harmonics, spacing)

                    temp_r2 = np.corrcoef(
                        field_signal,
                        out_dict['inverse_FFT'],
                    )[0][1] ** 2

                methods_dict[method][field] = out_dict.copy()

        else:
            raise ValueError('You must either use a R-squared threshold greater than 0 '
                             'but less than 1, or a integer N-harmonics threshold!'
                             )

    # prepare output
    for method in methods_dict.keys():
        logging.info('Completing %s analysis...' % method)

        for field in fields:
            # get the data for the given method-field combination
            data_dict = methods_dict[method][field]
            n_harmonics_used = data_dict['n_harmonics']
            ifft_df = data_dict['inverse_FFT_df']

            # init the output text file
            text_file = open(
                out_folder +
                '\\%s_%s_to_riverbuilder.txt' % (field, method), 'w+',
            )

            # save the IFFT DataFrame to csv
            ifft_df.to_csv(
                out_folder + '\\%s_harmonics_%s.csv' % (field, method)
            )

            # make a plot of the reconstructed river topo
            plt.plot(
                index_array,
                in_df.loc[:, str(field)].squeeze(),
                color='blue',
                label='Signal',
            )
            plt.plot(
                index_array,
                data_dict['inverse_FFT'],
                color='red',
                linestyle='--',
                label='Reconstructed signal',
            )

            plt.xlim(xmin=index_array.min(), xmax=index_array.max())

            if units != '':
                add_units = 'in %s' % units
            else:
                add_units = ''

            plt.xlabel('Distance along centerline %s' % add_units)
            plt.ylabel('Value')
            plt.title(
                f'{field}, {method} method, N={n_harmonics_used} '
                f'component harmonic reconstruction'
            )
            plt.grid(b=True, which='major', color='#666666', linestyle='-')
            plt.minorticks_on()
            plt.legend(loc='lower center')

            fig_title = out_folder + '\\%s_%s_plot.png' % (field, method)
            fig = plt.gcf()
            fig.set_size_inches(12, 6)
            plt.savefig(fig_title, dpi=300, bbox_inches='tight')
            plt.cla()

            # output a text file enabling import into RiverBuilder
            for num, amp in enumerate(data_dict['amplitudes']):
                if amp != 0.0:
                    freq = data_dict['frequencies'][num]
                    phase = data_dict['phases'][num]

                    # in the form of COS#=(a, f, ps, MASK0) for river builder inputs
                    text_file.write(
                        f'COS{num}=({amp}, {freq}, {phase}, MASK0)\n'
                    )
            text_file.close()

    logging.info(f'Analysis complete. Results @ {out_folder}')
    return out_folder
