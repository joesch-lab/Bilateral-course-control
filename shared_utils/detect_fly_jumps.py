import os
import cv2
import numpy as np
import pandas as pd
import math
from skimage import filters


def angle_from_3_points(p1, p2, p3):
    x1, y1 = p1
    x2, y2 = p2
    x3, y3 = p3

    a = np.sqrt((y3 - y2) ** 2 + (x3 - x2) ** 2)
    b = np.sqrt((y2 - y1) ** 2 + (x2 - x1) ** 2)
    c = np.sqrt((y3 - y1) ** 2 + (x3 - x1) ** 2)

    angle = (np.arccos((a**2 + b**2 - c**2)/(2*a*b))/np.pi)*180
    return angle, (a**2 + b**2 - c**2)/(2*a*b)


def find_jumps_video(Dir, first_frame=0):
    """script to detect and separate out male and female in the video"""
    length = 0
    errors = []
    ## constants for flies
    Area_max_cutoff = 1550
    Area_min_cutoff = 500
    Max_area = 3500
    body_cutoff = 5000
    divider_minor = 6
    divider_major = 3

    for file in os.listdir(Dir):
        if file.endswith('.avi'):
            filename = r'{}/{}'.format(Dir, file)
        if file.endswith('.csv') and 'corrected' not in file:
            csv_file = os.path.join(Dir, file)
            csv_data = pd.read_csv(csv_file, error_on_lines='skip')
            fly_ori = np.array(csv_data['fly_frame_ori'])
            pos_x = np.array(csv_data['pos_x'])
            pos_y = np.array(csv_data['pos_y'])

    ## find the video file
    orientation = []
    wing_angle = []
    fly_lost = []
    try:
        cap = cv2.VideoCapture(filename)
    except:
        print('The video was not found')
    cap.set(cv2.CAP_PROP_POS_MSEC, 0)

    z = -1
    real_angle = first_frame
    while ((cap.isOpened()) & (z < fly_ori.shape[0] -1)):
        ret, image = cap.read()
        # if 25 > pos_x[z + 1] or pos_x[z + 1] > image.shape[0] - 75 or 75 > pos_y[z + 1] or  pos_y[z + 1] > image.shape[0] - 75:
        #     z = z + 1
        #     wing_angle.append(0)
        #     continue
        if image is None:
            print('the end of the video {}'.format(z))
            break
        try:
            z = z + 1
            image = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)
            cv2.imshow('frame', image)
            cv2.waitKey(0)
        except:
           pass

        ## image inverted
        imagem = cv2.bitwise_not(image)

        ## median blur to remove salt and pepper noise
        imagem = cv2.medianBlur(imagem, 5)
        # cv2.imshow('frame', imagem[pos_y[z+1]-75:pos_y[z+1]+75, pos_x[z+1]-75:pos_x[z+1]+75])
        # cv2.waitKey(0)
        ## multi otsu binatization used to determine the thresolds for bg, body and wings
        try:
            thresholds = filters.threshold_multiotsu(imagem[pos_y[z+1]-75:pos_y[z+1]+75, pos_x[z+1]-75:pos_x[z+1]+75])
        except (IndexError, ValueError):
            errors.append(z)
            orientation.append(fly_ori[z])
            wing_angle.append(0)
            continue

        ## upper threshold applied to get the body only (without the wing)
        ## if the flies are close but not touching, they should appear separated here
        thresh, imagem2 = cv2.threshold(imagem, thresholds[1], 255, cv2.THRESH_BINARY)
        # cv2.imshow('frame', imagem2)
        # cv2.waitKey(0)

        ## find contours
        _, contours, _ = cv2.findContours(imagem2, cv2.RETR_LIST, cv2.CHAIN_APPROX_SIMPLE)
        area_max = 0
        fly = 0
        for i in range(len(contours)):
            Area = cv2.contourArea(contours[i])
            if len(contours[i]) > 20 and Area > 600 and Area < body_cutoff:
                # M = cv2.moments(contours[i])
                # cx = int(M['m10'] / M['m00'])
                # cy = int(M['m01'] / M['m00'])
                # if cx < 20 or cx > 40 or cy < 20 or cy > 40:
                #     cv2.imshow('frame', np.concatenate((image, imagem2)))
                #     cv2.waitKey(10)
                #     continue
                area_max = Area
                fly = i
                ellipse = cv2.fitEllipse(contours[i])
                cnt = contours
        if area_max == 0:
            fly_lost.append(z)
            # orientation.append(orientation[z-1])
            wing_angle.append(0)
            cv2.imshow('frame', np.concatenate((image, imagem2)))
            cv2.waitKey(1)
            continue
        ## centroid of fly
        M = cv2.moments(contours[fly])
        ellipse = cv2.fitEllipse(contours[fly])
        cx = int(M['m10'] / M['m00'])
        cy = int(M['m01'] / M['m00'])
        fly_angle = 180 - ellipse[2]

        center = (int(ellipse[0][0]), int(ellipse[0][1]))
        axes = (int(ellipse[1][0] / 2), int(ellipse[1][1] / 2))
        ep1 = (int(center[0] - axes[1] * np.sin((fly_angle / 180) * np.pi)),
               int(center[1] - axes[1] * np.cos((fly_angle / 180) * np.pi)))
        ep2 = (int(center[0] + axes[1] * np.sin((fly_angle / 180) * np.pi)),
               int(center[1] + axes[1] * np.cos((fly_angle / 180) * np.pi)))

        ## lower threshold applied to get body with wing
        thresh, imagem1 = cv2.threshold(imagem, thresholds[0], 255, cv2.THRESH_BINARY)
        # cv2.imshow('frame', np.concatenate((imagem1, imagem2)))
        # cv2.waitKey(0)
        imagem_new = imagem1 - imagem2
        # cv2.imshow('frame', imagem_new)
        # cv2.waitKey(0)

        ## erosion to get rid of the pixels that appear around the body
        kernel = np.ones((5, 5), np.uint8)
        erosion = cv2.erode(imagem_new, kernel, iterations=1)

        ## contour of the remaining (hopefully only wing)
        _, contours_eroded, _ = cv2.findContours(erosion, cv2.RETR_LIST, cv2.CHAIN_APPROX_SIMPLE)

        ## try to find wings here
        area_min = 200
        area_max = 1000
        wing = 0
        cx2 = []
        cy2 = []
        for i in range(len(contours_eroded)):
            Area = cv2.contourArea(contours_eroded[i])
            if area_max > Area > area_min:
                # area_max = Area
                wing = wing + 1
                M = cv2.moments(contours_eroded[i])
                ## new image created with only the male fly
                image_test = np.zeros((100, 100), np.uint8)  ## shows wings
                cv2.drawContours(image_test, contours_eroded, i, 255, -1)
                image_test = cv2.cvtColor(image_test, cv2.COLOR_GRAY2BGR)
                cx_new = int(M['m10'] / M['m00'])
                cy_new = int(M['m01'] / M['m00'])
                # cx2 = cx2 + int(M['m10'] / M['m00'])
                # cy2 = cy2 + int(M['m01'] / M['m00'])
                if cx+50 > cx_new > cx-50 and cy+50 > cy_new > cy-50:
                    cx2.append(cx_new)
                    cy2.append(cy_new)

        if len(cx2)==2:
            wing_angle.append(angle_from_3_points([cx2[0], cy2[0]], [cx, cy], [cx2[1], cy2[1]])[0])
        else:
            wing_angle.append(0)

    csv_file_new = csv_file[:-4] + 'corrected' + '.csv'
    csv_file_new = open(csv_file_new, 'w')
    df = pd.DataFrame()
    for keys in csv_data.keys():
        df[keys] = csv_data[keys][:z+1]
    df['wing_angle'] = wing_angle[:z+1]
    df['fly_lost'] = fly_lost
    df.to_csv(csv_file_new, index=False)
    csv_file_new.close()
    return fly_lost


def find_jumps(dist_pos_x, dist_pos_y, num):
    """
    Find the frames where the fly jumps based on the distance between the fly's position in consecutive frames
    :param dist_pos_x: list of x2-x1 values
    :type dist_pos_x: list
    :param dist_pos_y: list of y2-y1 values
    :type dist_pos_y: list
    :param num: distance threshold
    :type num: float
    """
    dist = math.sqrt(np.add(np.square(dist_pos_x), np.square(dist_pos_y)))

    jumps_frames = []
    for i in range(dist.shape[0]):
        if dist[i] > num:
            jumps_frames.append(i)
        else:
            pass

    jumps = []
    i = 0
    while i <= len(jumps_frames)-2:
        if jumps_frames[i+1] - jumps_frames[i] > 5:
            jumps.append(jumps_frames[i])
            i = i + 1
        else:
            jumps.append(jumps_frames[i+1])
            i = i + 2

    return jumps


def time_to_jump(jump_frames, frame_changes):
    """
    Find the time it takes for the fly to jump after the stimulus is presented
    :param jump_frames: list of frames where the fly jumps
    :type jump_frames: list
    :param frame_changes: list of frames where the stimulus changes
    :type frame_changes: list
    :return: jump_delay - list of time it takes for the fly to jump after the stimulus is presented, jump_time - list of frame number when the fly jumps
    :rtype: (list, list)
    """
    j = 0
    i = 0
    jump_delay = []
    jump_time = []
    while(i < (len(frame_changes)-1) and j < (len(jump_frames)-1)):
        if jump_frames[j] > frame_changes[i] and jump_frames[j] < frame_changes[i+1]:
            jump_delay.append(jump_frames[j]-frame_changes[i])
            jump_time.append(i)
            j += 1
            i += 1              ## if i keep this, only the first jump after loom is registered
        else:
            if jump_frames[j] < frame_changes[i]:
                j += 1
            else:
                jump_delay.append(-1)
                i += 1
    return jump_delay, jump_time


def loom_response_jump(speed_data, loom_start=124):
    """
    Find the time it takes for the fly to jump after the loom starts
    :param speed_data: list of speed values where -10 indicates the fly made a jump
    :type speed_data: list
    :param loom_start: frame number when the loom starts
    :type loom_start: int
    :return: time_to_jump - list of time it takes for the fly to jump after the loom is presented
    :rtype: list
    """
    time_to_jump = []
    # speed_data = speed_data[0]
    for i in range(len(speed_data)):
        if speed_data[i] == -10:
            if i-loom_start < -1:
                pass
            else:
                time_to_jump.append(i-loom_start)
                break

    return time_to_jump
