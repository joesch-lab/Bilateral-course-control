import numpy as np

def frame_timeperiod(timestamp_data, factor=1e-3, exp=16):
    """
    Calculate the timeperiod between consecutive frames from the timestamp data. 
    This code is very messy because different cameras have different ways of recording timestamps.
    Best to write a new function for new data.
    :param timestamp_data: The timestamp data.
    :type timestamp_data: numpy.ndarray
    :param factor: The factor to multiply the timeperiod by, depending on the camera and if we have milliseconds or microseconds or nanoseconds.
    :type factor: float
    :param exp: The expected value of the timeperiod, in milliseconds.
    :type exp: int
    :returns: The timeperiod between consecutive frames.
    :rtype: numpy.ndarray
    """
    time = np.diff(timestamp_data)
    ## for the pointgrey camera, values in time will be negative everytime the clock resets to 0, that is, after every second
    ## for the basler camera, values in time will be negative everytime when a new trial starts

    if (np.argwhere(time < 0).shape[0]<100) or (np.all(time[np.argwhere(time < 0)[100, 0]+1:np.argwhere(time < 0)[100, 0]+1+10] > 0)):
        ## if there very few negative values, its a basler camera
        ## to check if there are lots of negative time values
        timeperiod = time ## the difference between the timestamps gives the frame_period
        IFI = np.nanmean(timeperiod[timeperiod>0]) ## take average of only the positive values
        timeperiod = np.insert(timeperiod, 0, IFI)
        timeperiod[np.where((timeperiod<0) | (timeperiod>IFI*10))] = IFI
        if exp-2 < IFI * factor < exp+2:
            pass
        else:
            timeperiod = timeperiod * factor
    else:
        timeperiod = timestamp_data ## the timestamps already contrained the frame_priods, no need to subtreact
        IFI = np.nanmean(timeperiod[timeperiod>0])
        timeperiod[np.where(timeperiod < 0)] = IFI
        if exp-2 < IFI * 1e-3 < exp+2:
            pass
        else:
            timeperiod = timeperiod * factor
    return timeperiod