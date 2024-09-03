import psychopy
import random
import image_methods
import cv2
import numpy as np


global FRAMERATE
FRAMERATE = 60

class fly:
    def __init__(self,strain,id,sex,age):
        self.id = id
        self.sex = sex
        self.age = age
        self.strain = strain

class video:
    def __init__(self,res,fps):
        self.res = res
        self.fps = FRAMERATE

class exp:
    def __init__(self,date,time,stim,bg,fg):
        self.time = time
        self.date = date
        self.stim = stim
        self.bg = bg
        self.fg = fg

class pinwheel_stim:  ## depending on the name of the stim, make an if else statement and add different arguments ## all arguments should have a default value
    def __init__(self, direction,spatial,temp,exp_time,stim_time,on_frames,off_frames,contrast,size,closedloop_gain, stim_attributes=[]):
        self.name = 'pinwheel'
        self.direction = direction      ## 1 if counterclockwise, -1 if clockwise and 0 if no rotation ## could be a list for certain stimuli ##
        self.temp = temp                ## rotations made by the pinwheel in 1 second
        self.spatial = spatial          ## no. of bars in pinwheel
        self.contrast = contrast        ## contrast of the pattern, from 0 to 100
        self.stim_size = size
        self.stim_attributes = stim_attributes
        self.closedloop_gain = closedloop_gain

        self.stim_time = stim_time      ## w.r.t initiation of the experiment
        self.on_frames = on_frames      ## 0 = first frame or loop recorded
        self.off_frames = off_frames    ## 0 = first frame or loop recorded

        self.angle = (360 * self.temp) / (FRAMERATE * self.spatial)
        ## Source File
        self.stimulus_image = '24ringpinwheel_100_6.png'
        self.stimulus_image = r'{}\{}'.format('D:\Roshan\Project\Python_codes\Stimulus\Generate stimulus', self.stimulus_image)

    def make_stim(self, win):
        ## Generate stimulus in psychopy
        self.wheel_stim = psychopy.visual.ImageStim(
            win=win,
            image=self.stimulus_image,
            mask='circle',
            pos=(0, 0),
            ori=0,  # 0 is vertical,positive values are rotated clockwise
            size=1,
        )
        self.wheel_stim.autoDraw = False


    def show_stim(self, pos, fly_turn, closedloop_gain):
        self.wheel_stim.pos = pos
        self.wheel_stim.ori = (self.wheel_stim.ori + self.direction * self.angle - closedloop_gain * fly_turn) % 360
        self.wheel_stim.draw()

class ringpinwheel_stim:
    def __init__(self, direction,spatial,temp,exp_time,stim_time,on_frames,off_frames,contrast,size,hole_size, closedloop_gain, stim_attributes=[]):
        self.name = 'ringpinwheel'
        self.direction = direction      ## 1 if counterclockwise, -1 if clockwise and 0 if no rotation ## could be a list for certain stimuli ##
        self.temp = temp                ## rotations made by the pinwheel in 1 second
        self.spatial = spatial          ## no. of bars
        self.contrast = contrast  ## contrast of the pattern, from 0 to 100
        self.stim_size = size
        self.hole_size = hole_size
        self.closedloop_gain = closedloop_gain
        self.stim_attributes = stim_attributes

        self.stim_time = stim_time      ## w.r.t initiation of the experiment
        self.on_frames = on_frames      ## 0 = first frame or loop recorded
        self.off_frames = off_frames    ## 0 = first frame or loop recorded



        ## Source File
        self.stimulus_image = 'ringpinwheel_'+str(hole_size)+'_'+str(spatial)+'_'+str(contrast)+'.png'
        self.stimulus_image = r'{}\{}'.format('D:\Roshan\Project\Python_codes\Stimulus\Generate stimulus', self.stimulus_image)

    def make_stim(self, win):
        ## Generate stimulus in psychopy
        self.wheel_stim = psychopy.visual.ImageStim(
            win=win,
            image=self.stimulus_image,
            mask='circle',
            pos=(0, 0),
            ori=0,  # 0 is vertical,positive values are rotated clockwise
            size=1,
        )
        self.wheel_stim.autoDraw = False

    def show_stim(self, pos, fly_turn, closedloop_gain):
        self.wheel_stim.pos = pos
        self.wheel_stim.ori = (self.wheel_stim.ori + self.direction * self.angle - closedloop_gain * fly_turn) % 360
        self.wheel_stim.draw()

class looming_stim:  ## depending on the name of the stim, make an if else statement and add different arguments ## all arguments should have a default value
    def __init__(self,speed,contrast, expansion_time,stay_time,stim_time,on_frames,off_frames,stim_attributes=[]):
        self.name = 'looming'
        self.speed = speed               ## speed of expansion
        self.contrast = contrast         ## contrast of dot
        self.expansion_time = expansion_time        ## time the loom expands
        self.stay_time = stay_time          ## time dot persists after maximum expansion

        self.stim_time = stim_time      ## w.r.t initiation of the experiment
        self.on_frames = on_frames      ## frame when stimulus starts
        self.off_frames = off_frames    ## 0 = first frame or loop recorded
        self.stim_attributes = stim_attributes

    def make_stim(self,win):
        self.dot_stim = psychopy.visual.Circle(
            win=win,
            units='norm',
            radius=1/512,
            lineColor=0,
            pos=(0, 0)
        )
        self.dot_stim.autoDraw = False

    def show_stim(self,pos,ori):
        self.dot_stim.fillColor = self.contrast
        self.dot_stim.size = self.dot_stim.size + self.speed
        self.dot_stim.pos = pos
        self.dot_stim.draw()

class looming_stim_constant_speed:  ## depending on the name of the stim, make an if else statement and add different arguments ## all arguments should have a default value
    def __init__(self,speed, contrast, size, distance,toi, stay_time, stim_time, on_frames, off_frames, stim_attributes=[]):
        self.name = 'looming_constant_speed'
        self.speed = speed
        self.contrast = contrast
        self.size = size
        self.distance = distance
        self.toi = toi

        def angle_subtended(self, frames):
            return (self.size) / ((self.distance - self.speed * (frames / (self.toi * 60))))

        def radius_in_px(self):
            result = []
            for i in range(self.toi * 60):
                result.append((self.angle_subtended(i) * 2 * (1024 / 110)))
            return result

        self.stay_time = stay_time
        self.on_frames = on_frames      ## frame when stimulus starts
        self.off_frames = off_frames    ## 0 = first frame or loop recorded
        self.stim_attributes = stim_attributes

class bout:
    def generate_random_values(self, dict):
        output_list = []
        for keys in dict:
            output_list.appen(random.choice(dict['keys']))
        return output_list

    ##function for generating mask
    def ROI(image, x, y):
        height, width = image.shape
        loc1 = x
        size1 = y
        circle_img = np.zeros((height, width), np.uint8)
        cv2.circle(circle_img, tuple(loc1), size1, 1, thickness=-1)
        masked_data = cv2.bitwise_and(image, image, mask=circle_img)
        return masked_data

    ## create bout object
    def __init__(self, bout_duration, stimuli=[], stimulus_duration = []): ## stimulus_duration is a list of lists with start_frame and end_frame
        ## Pinwheel stimulus parameters
        pinwheel_stim_parameters = {
            'DIRECTION': [1, -1],
            'NUM_BARS': [4, 6, 12],
            'FREQ': [2, 5, 10, 15, 20],
            'PINWHEEL_CONTRAST': [100, 90, 75, 50],
            'SIZE': [2, 1, 0.5],
            'CLOSEDLOOP_GAIN': [1, -0, -1]}

        ## Looming stimulus parameters
        looming_stim_parameters = {
            'SPEED': [1, 2, 4, 8, 12],
            'LOOM_CONTRAST': [-1, -0.9, -0.75, -0.5, 1],
            'EXPANSION_TIME': [20, 30, 40],
            'STAY_TIME': [10, 20, 30]}
        self.duration = bout_duration
        self.stimuli = stimuli
        self.parameters = []
        self.stimuli_to_show = []
        i = 0
        for stimulus in stimuli:
            if stimulus == 'pinwheel':
                self.parameters.append(self.generate_random_values(pinwheel_stim_parameters))
                self.stimuli_to_show.append(pinwheel_stim(self.parameters[i]))
            elif stimulus == 'looming':
                self.parameters.append(self.generate_random_values(looming_stim_parameters))
                self.stimuli_to_show.append(looming_stim(self.parameters[i]))
            i = i + 1

    ##

    def start(self, win, writer_csv, writer_small, camera, ):
        W_RATIO = 512
        W_REM = 0
        H_RATIO = 512
        H_REM = 0

        stimuli_to_show = []
        i = 0
        for stimulus in self.stimuli_to_show:
            stimulus.make_stim(win)
            i = i + 1


        while i<self.current_frame:
            image, timestamp = image_methods.grab_image(camera, 'ptgrey')
            cv_image = ROI(image, loc, size)
            ret, diff_img = cv2.threshold(cv_image, 60, 255, cv2.THRESH_BINARY)
            this, contours, hierarchy = cv2.findContours(diff_img, cv2.RETR_LIST, cv2.CHAIN_APPROX_SIMPLE)
            pos, ellipse, cnt = fly_track.fly_court_pos(contours)
            if len(pos) != 1:
                frames_lost += 1  # counts number of frames in which the location of fly could not be determined
                print('fly lost')
                pos.append(position)
                if frames_lost > 150:  # if the number of lost frames is more than 15, abort experiment
                    print('Too many frames are being lost')
                    break
            else:
                frames_lost = 0
                fly_ori = ellipse[2]
                fly_turn = fly_ori - fly_ori_last

                if fly_turn > 90:
                    fly_turn = -(180 - fly_turn)
                elif fly_turn < -90:
                    fly_turn = 180 + fly_turn

                fly_frame_ori = fly_frame_ori_last + fly_turn
            x = (pos[0][0] - 512) / W_RATIO + W_REM
            y = (pos[0][1] - 512) / H_RATIO + H_REM


            for stimulus in stimuli_to_show:
                if self.current_frame >= stimulus.start_frame and self.current_frame < stimulus.end_frame:
                    stimulus.show_stim()

            win.flip()

            writer_csv.writerow([pos[0][0], pos[0][1], fly_frame_ori, timestamp, dir, loom_status[0]])
            cropped_image = cv_image[int(pos[0][1]) - 30:int(pos[0][1]) + 30, int(pos[0][0]) - 30:int(pos[0][0]) + 30]
            if cropped_image.shape == (60, 60):
                writer_small.writeFrame(cropped_image)
            else:
                print("cropping failure...")
                writer_small.writeFrame(filler_frame)

            position = pos[0]
            fly_ori_last = fly_ori
            fly_frame_ori_last = fly_frame_ori
            frame_number = frame_number + 1
            video_frame = video_frame + 1
        return 1

