import os
import cv2
import numpy as np

## script to check if the post-hoc fly tracking was correct
## take the video, take corrected and old orientation and plot both to see what happened
def check_tracking_results(Dir, correct_ori, heading_direction):
    for files in os.listdir(Dir):
        if files.endswith('.avi'):
            vid_file = os.path.join(Dir, files)
    cap = cv2.VideoCapture(vid_file)
    starttime = 0
    cap.set(cv2.CAP_PROP_POS_FRAMES, starttime)

    while cap.isOpened():
        ret, image = cap.read()
        if image is None:
            print('the end of the video')
            break
        else:
            image = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)
        i = int(cap.get(cv2.CAP_PROP_POS_FRAMES))

        ## white line is the orientation
        image = cv2.line(image, (30, 30), (int(100 * np.cos(((correct_ori[i] / 180) * np.pi) + np.pi / 2) + 30),
                                           int(100 * np.sin(((correct_ori[i] / 180) * np.pi) + np.pi / 2) + 30)), (255, 255, 255), 2)
        ## black line is the heading direction
        image = cv2.line(image, (30, 30), (int(100 * np.cos(((heading_direction[i] / 180) * np.pi) + np.pi / 2) + 30),
                                           int(100 * np.sin(((heading_direction[i] / 180) * np.pi) + np.pi / 2) + 30)), (0, 0, 0), 2)

        cv2.imshow('frame', cv2.resize(image, (120, 120)))
        cv2.waitKey(50)
        i += 1

    return 1