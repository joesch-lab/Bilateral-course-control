import numpy as np
from scipy.signal import find_peaks, peak_widths
import pywt

def check_saccade(peaks, peak_widths, ori_smooth):
    """
    Calculates the total turn and width turn for each peak in the given data.
    :param peaks: The indices of the saccade peaks.
    :type peaks: numpy.ndarray
    :param peak_widths: The widths of the saccades.
    :type peak_widths: numpy.ndarray
    :param ori_smooth: The smoothed orientation data.
    :type ori_smooth: numpy.ndarray
    :returns: saccade_total_turn : The turns made in each saccade (in degrees) in a 12 frame window
                saccade_width_turn : The turns made in each saccade (in degrees) in a window given by peak_widths
    :rtype: (numpy.ndarray, numpy.ndarray)
    """
    saccade_total_turn = []
    for i in peaks:
        if 6 < i < len(ori_smooth) - 6:
            saccade_total_turn.append(ori_smooth[i + 5] - ori_smooth[i - 5])
        else:
            saccade_total_turn.append(0)

    saccade_width_turn = []
    for i,j in zip(peaks, peak_widths):
        if int(i + j//2) < i < len(ori_smooth) - int(i + j//2):
            saccade_width_turn.append(ori_smooth[int(i + j // 2)] - ori_smooth[int(i - j // 2)])
        else:
            saccade_width_turn.append(0)
    return np.array(saccade_total_turn), np.array(saccade_width_turn)


def swt_denoising(data, wavelet='bior2.6', level = 7):
    """
    Denoises the data using stationary wavelet transform while only keeping frequencies between 10 to 20 Hz.
    :param data: The data to be denoised.
    :type data: numpy.ndarray
    :param wavelet: The wavelet to use for denoising.
    :type wavelet: str
    :param level: The level of decomposition to use.
    :type level: int
    :returns: The denoised data.
    :rtype: numpy.ndarray
    """
    freq = pywt.scale2frequency(wavelet, 2 ** np.arange(0, level)[::-1]) / 0.016
    n = (data.shape[0]//1200) * 1200
    data = data[:n]
    data = np.append(data, np.zeros(n%128))
    swt = pywt.swt(data, wavelet, level = level)
    index= np.argwhere((freq>10) * (freq<20))
    re_data = pywt.iswt(swt[index[0,0]], wavelet)
    return re_data


def detect_saccades_cwt(angular_speed, threshold):
    """
    Detects saccades using wavelet transform and peak finding.
    :param angular_speed: The angular speed data.
    :type angular_speed: numpy.ndarray
    :param threshold: The threshold for detecting saccades.
    :type threshold: float
    :returns: The indices of the saccade peaks, the values of the peaks, 
    the widths of the peaks, and the peak widths.
    :rtype: (numpy.ndarray, numpy.ndarray, numpy.ndarray, numpy.ndarray)
    """
    signal = angular_speed
    denoised_signal = swt_denoising(signal, level = 5) ## "denoise" OR extract the frequency between 10 to 20 Hz
    cwt_peaks, cwt_peak_values, cwt_width_heights, cwt_peak_width = detect_saccades(denoised_signal, threshold= threshold)
    peaks, peak_values, width_heights, peak_width = detect_saccades(angular_speed, threshold=threshold, width=[5, 30])
    ## width for the angular speed peak has to be less than the requirement for
    ## check of the two local maxima match or not
    i = 0
    j = 0
    to_delete = []
    final_peak_values = [] ## the final peak_values that are recorded should be the peak_values of the angular_speed
    final_peaks = []
    while(i<len(cwt_peaks) and j<len(peaks)):
        if cwt_peaks[i] + 3 >= peaks[j] >= cwt_peaks[i] - 3:
            final_peak_values.append(peak_values[j])
            final_peaks.append(peaks[j])
            i = i + 1
            j = j + 1
        elif cwt_peaks[i] - 3 > peaks[j]:
            j = j + 1
        elif cwt_peaks[i] + 3 < peaks[j]:
            to_delete.append(i)
            i = i + 1
    cwt_peaks = np.delete(cwt_peaks, to_delete)
    cwt_peak_values = np.delete(cwt_peak_values, to_delete)
    cwt_width_heights = np.delete(cwt_width_heights, to_delete)
    cwt_peak_width = np.delete(cwt_peak_width, to_delete)

    return np.array(final_peaks), np.array(final_peak_values), cwt_width_heights, cwt_peak_width


def detect_saccades(angular_speed, threshold, width = [5, 30]):
    """
    Detects saccades using peak finding, called by detect_saccades_cwt.
    :param angular_speed: The angular speed data.
    :type angular_speed: numpy.ndarray
    :param threshold: The height parameter for peak finding.
    :type threshold: float
    :param width: The width parameter for peak finding, min and max.
    :type width: list
    :returns: The indices of the saccade peaks, the values of the peaks, height at width, and the peak widths.
    :rtype: (numpy.ndarray, numpy.ndarray, numpy.ndarray, numpy.ndarray)"""
    # thresh = 4 * np.median(np.divide(abs(angular_speed), 0.6745))
    height = threshold
    peaks, props = find_peaks(angular_speed, height=height, distance=5, prominence=0.8*height, wlen=30, width = width, rel_height = 1)
    ## distance = minimal horizontal distance between two peaks
    ## prominence = height from the lowest point of the peak
    ## width = with of the peak, the upper limit ensures that we do not get slow changes
    peak_values = angular_speed[peaks]
    peak_width, width_heights, _, _ = peak_widths(angular_speed, peaks, rel_height = 0.8)
    ## peak_widths is a tricky argument, it might overestimate the width, beware

    ## for the negative peaks
    peaks_neg, props_neg = find_peaks(angular_speed*-1, height=height, distance=3, prominence=0.8*height, wlen = 36,
                              width = [3, 30], rel_height = 1)
    peak_values_neg = angular_speed[peaks_neg]
    peak_width_neg, width_heights_neg, _, _ = peak_widths(angular_speed*-1, peaks_neg, rel_height = 0.8)

    ## merge positive and negative peaks
    peaks = np.concatenate((peaks, peaks_neg))
    sorted_indices = np.argsort(peaks)
    peaks = np.sort(peaks)
    peak_values = np.concatenate((peak_values, peak_values_neg))[sorted_indices]
    peak_width = np.concatenate((peak_width, peak_width_neg))[sorted_indices]
    width_heights = np.concatenate((width_heights, width_heights_neg))[sorted_indices]
    return peaks, peak_values, width_heights, peak_width
