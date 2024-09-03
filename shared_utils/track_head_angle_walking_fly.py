import cv2
import numpy as np
import os
from skimage import filters
import skvideo.io as sk
import scipy.io
import json
import pandas as pd
from scipy import ndimage
import copy

def reorient_fly(img, angular_offset, pos_x, pos_y):
    if img is not None:
        img_small = img
        ## to remove the black panels that I originally made, if the original padding is white, remove this
        pad_1 = (0, 0)
        pad_2 = (0, 0)
        if pos_x - 75 <= 0:
            pad_2 = (75 - pos_x, 0)
            image_edge_2_before = 0
            image_edge_2_after = pos_x + 75
        if pos_x + 75 >= 736:
            pad_2 = (0, 75 + pos_x - 736)
            image_edge_2_before = pos_x - 75
            image_edge_2_after = 736
        if pos_y - 75 <= 0:
            pad_1 = (75 - pos_y, 0)
            image_edge_1_before = 0
            image_edge_1_after = pos_y + 75
        if pos_y + 75 >= 736:
            pad_1 = (0, 75 + pos_y - 736)
            image_edge_1_before = pos_y - 75
            image_edge_1_after = 736
        if np.ndim(img_small) == 2:
            img_small[:, :] = np.pad(img_small[pad_1[0]:150 - pad_1[1], pad_2[0]:150 - pad_2[1]], (pad_1, pad_2), 'constant', constant_values=(255, 255))
        else:
            img_small[:, :, 0] = np.pad(img_small[pad_1[0]:150 - pad_1[1], pad_2[0]:150 - pad_2[1], 0], (pad_1, pad_2), 'constant', constant_values=(255, 255))
            img_small[:, :, 1] = np.pad(img_small[pad_1[0]:150 - pad_1[1], pad_2[0]:150 - pad_2[1], 1], (pad_1, pad_2), 'constant', constant_values=(255, 255))
            img_small[:, :, 2] = np.pad(img_small[pad_1[0]:150 - pad_1[1], pad_2[0]:150 - pad_2[1], 2], (pad_1, pad_2), 'constant', constant_values=(255, 255))
        ##
    else:
        print('shit')

    img_small = ndimage.rotate(img_small, angular_offset, cval=255)
    if np.ndim(img_small) == 2:
        s1, s2 = img_small.shape
    else:
        s1, s2, _ = img_small.shape
    img_small = img_small[s1 // 2 - 75:s1 // 2 + 75, s2 // 2 - 75:s2 // 2 + 75]

    return img_small

def largest_cnt(contours):
    max_cnt = 0
    max_area = 0
    for i in range(0, len(contours)):
        if len(contours[i]) > 5 and len(contours[i]) < 500:
            Area = cv2.contourArea(contours[i])
            if Area > max_area:
                max_area = Area
                max_cnt = i
            else:
                continue
    return max_cnt

Dir2 = r'C:\Users\rsatapat\Documents\HS_split_SPARCD_33'
for files in os.listdir(Dir2):
    if files.endswith('.avi') and 'combined' not in files and 'resnet' not in files and 'new' in files and 'compressed' not in files:
        vidname = files
        print(vidname)
        vid = r'{}/{}'.format(Dir2, files)
        cap = cv2.VideoCapture(vid, 0)
    elif files.endswith('.csv') and 'resnet' not in files:
        csvfile = r'{}/{}'.format(Dir2, files)
    elif files.endswith('.m'):
        ori_smooth = scipy.io.loadmat(os.path.join(Dir2, files))['fly_theta']
    elif files.endswith('.json') and 'resnet' not in files:
        jsonfile = r'{}/{}'.format(Dir2, files)
        with open(jsonfile, 'r') as f:
            data = json.load(f)

csvfile = open(csvfile, 'r')
data = pd.read_csv(csvfile)
stim = list(data['stim'])
time = list(data['timestamp'])
pos_x = np.array(data['pos_x'])
pos_y = np.array(data['pos_y'])
cap = cv2.VideoCapture(vidname)

while (True):
    ## first image
    ret, img = cap.read()
    img = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)
    angular_offset = 90 - (int(ori_smooth[0]) % 360)
    img_small = reorient_fly(img, angular_offset, pos_x[0], pos_y[0])
    cv2.imshow('frame', img_small)
    if cv2.waitKey(0) & 0xFF == ord('m'):
        break
    elif cv2.waitKey(0) & 0xFF == ord(' '):
        continue
## user finds the arena, by selecting ROI
arena = cv2.selectROI(img_small, False)
print(arena)
loc = [arena[0] + arena[2] // 2, arena[1] + arena[3] // 2]
size = arena[2] // 2

template = img_small[arena[1]:arena[1] + arena[3], arena[0]:arena[0] + arena[2]]
print('template size' + str(template.shape))
cv2.imshow('frame', template)
cv2.waitKey(0)
## eroded tempalte
thresholds = filters.threshold_multiotsu(template)
ret, eroded_template = cv2.threshold(template, thresholds[0], 255, cv2.THRESH_BINARY)
eroded_template = cv2.bitwise_not(eroded_template)
# cv2.imshow('frame', img_small)
kernel = cv2.getStructuringElement(cv2.MORPH_CROSS, (10, 10))
eroded_template = cv2.erode(eroded_template, kernel, iterations=2)
im2, template_cnt, hierarchy = cv2.findContours(eroded_template, cv2.RETR_TREE, cv2.CHAIN_APPROX_SIMPLE)
template_cnt = template_cnt[largest_cnt(template_cnt)]
##
frame=0
while(True):
    ret, img = cap.read()
    img = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)
    # img = cv2.bitwise_not(img)
    # img = cv2.medianBlur(img, 5)

    angular_offset = 90 - (int(ori_smooth[frame]) % 360)
    img_small = reorient_fly(img, angular_offset, pos_x[frame], pos_y[frame])

    res = cv2.matchTemplate(img_small, template, cv2.TM_CCOEFF_NORMED)
    min_val, max_val, min_loc, max_loc = cv2.minMaxLoc(res)
    top_left = max_loc
    w, h = template.shape[::-1]
    bottom_right = (top_left[0] + w, top_left[1] + h)
    # cv2.rectangle(img_small, top_left, bottom_right, 255, 2)
    img_small = cv2.resize(img_small[top_left[1]:top_left[1] + h, top_left[0]:top_left[0] + w], (w*10, h*10))
    img_small_copy = copy.deepcopy(img_small)
    # cv2.imshow('frame', img_small)
    # cv2.waitKey(0)

    thresholds = filters.threshold_multiotsu(img_small)
    ret, img_small = cv2.threshold(img_small, thresholds[0], 255, cv2.THRESH_BINARY)
    img_small = cv2.bitwise_not(img_small)
    # cv2.imshow('frame', img_small)
    # kernel = cv2.getStructuringElement(cv2.MORPH_CROSS, (10, 10))
    erosion = cv2.erode(img_small, kernel, iterations=2) ## change the iterations and kernel size to remove legs from the image
    # cv2.imshow('frame', erosion)
    # cv2.waitKey(0)

    ## find the antennae by looking for corners in the image
    corners = cv2.goodFeaturesToTrack(erosion, 10, 0.01, 80) ## the final parameter, minimum distance is critical
    for i in corners:
        x, y = i.ravel()
        cv2.circle(erosion, (x, y), 5, 255, -1)
        cv2.circle(img_small_copy, (x, y), 5, 255, -1)
    # cv2.imshow('frame', erosion)
    # cv2.waitKey(0)
    cv2.imshow('frame', img_small_copy)
    cv2.waitKey(0)

    ## find contour
    # im2, contours, hierarchy = cv2.findContours(erosion, cv2.RETR_TREE, cv2.CHAIN_APPROX_SIMPLE)
    # cnt_img = cv2.drawContours(erosion, contours, -1, (0, 255, 0), 3)
    ## contour with the largest size, the whole fly
    # contours = contours[largest_cnt(contours)]
    # ret = cv2.matchShapes(contours, template_cnt, cv2.CONTOURS_MATCH_I1, 0.0)
    # cv2.imshow('frame', cnt_img)
    # cv2.waitKey(0)

    ## find contour
    # this, contours, hierarchy = cv2.findContours(erosion, 2, 1)
    ## contour with the largest size, the whole fly
    # cnt = contours[largest_cnt(contours)]
    ## find the convex hull, essentially a low res contour
    # hull = cv2.convexHull(cnt, returnPoints=True)
    # cnt_img = cv2.drawContours(cnt_img, [hull], 0, (255, 255, 255), 3)
    ## find the topmost point on the hull, the antennae
    # hull_sorted = hull[:, 0, :][np.argsort(hull[:, 0, 1])] # sort by the y values
    ## there might some noise, so some cleaning to make sure we only take the antennae
    # for i in range(4):
    #     for j in range(i+1, 4):
    #         pass
    #         distance = np.sqrt(np.sum(np.square(np.subtract(hull_sorted[i], hull_sorted[j]))))
    #         if 110 > distance > 100:
    #             pass
    #             head_angle = np.divide(np.arctan2(hull_sorted[j, 1]-hull_sorted[i, 1], hull_sorted[j, 0]-hull_sorted[i, 0]), np.pi)*180
    #             cv2.putText(cnt_img, str(head_angle), (200, 50), cv2.FONT_HERSHEY_SIMPLEX,  1, (255,255,255), 3)
    #             cv2.circle(cnt_img, tuple(hull_sorted[i]), 10, (255, 255, 255), -1)
    #             cv2.circle(cnt_img, tuple(hull_sorted[j]), 10, (255, 255, 255), -1)

    # defects = cv2.convexityDefects(cnt, hull)
    # for i in range(defects.shape[0]):
    #     s, e, f, d = defects[i, 0]
    #     start = tuple(cnt[s][0])
    #     end = tuple(cnt[e][0])
    #     far = tuple(cnt[f][0])
    #     cv2.line(erosion, start, end, [0, 255, 0], 2)
    #     cv2.circle(erosion, far, 5, [0, 0, 255], -1)
    # cv2.imshow('frame', cnt_img)
    # cv2.waitKey(0)

    # sobelxy = cv2.Sobel(src=imagem, ddepth=-1, dx=2, dy=2, ksize=7)
    # cv2.imshow('frame', sobelxy)
    # cv2.waitKey(0)

    # canny = cv2.Canny(image=erosion, threshold1=5, threshold2=250, apertureSize=7)
    # cv2.imshow('frame', canny)
    # cv2.waitKey(0)

    # thresholds = filters.threshold_multiotsu(imagem)
    # thresh, imagem2 = cv2.threshold(imagem, thresholds[1], 255, cv2.THRESH_BINARY)
    #
    # cv2.imshow('frame', imagem2)
    # cv2.waitKey(0)
    frame+=1
