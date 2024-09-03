class fly:
    def __init__(self,strain,id,sex,age):
        self.id = id
        self.sex = sex
        self.age = age
        self.strain = strain

class video:
    def __init__(self,res,fps):
        self.res = res
        self.fps = fps

class exp:
    def __init__(self,date,time,stim,bg,fg):
        self.time = time
        self.date = date
        self.stim = stim
        self.bg = bg
        self.fg = fg

class stim_inst:  ## depending on the name of the stim, make an if else statement and add different arguments ## all arguments should have a default value
    def __init__(self,stim, direction,spatial,temp,exp_time,stim_time,on_frames,off_frames,contrast,size,stim_attributes=[]):
        self.stim = stim                ## name of stimulus e.g. 'pinwheel','double_pinwheel'  ## list with name of stimulus and name of file
        self.direction = direction      ## 1 if counterclockwise, -1 if clockwise and 0 if no rotation ## could be a list for certain stimuli ##
        self.temp = temp                ## no. of bars in the pinwheel
        self.spatial = spatial          ## rotations made by the pinwheel in 1 second
        self.exp_time = exp_time        ## w.r.t initiation of the experiment
        self.stim_time = stim_time      ## w.r.t initiation of the experiment
        self.on_frames = on_frames      ## 0 = first frame or loop recorded
        self.off_frames = off_frames    ## 0 = first frame or loop recorded
        self.contrast = contrast        ## contrast of the pattern, from 0 to 100
        self.stim_size = size
        self.stim_attributes = stim_attributes