function [stim_dir,nrep]=reconstruct_scanning_rect_with_pause_fullres(logfile)
    %read log file and recover parameters
    params=read_scanning_rect_with_pause(logfile);
    nrep=params.nrep;
    freezeF = params.freezefr;
    dx=params.dx;
    rect_width=params.rwidth;
    rect_height=params.rheight;
    nfrst=params.nfr;
    
    
    
    debug=0;
    %reconstruct the stimulus
    % the size of the on screen window in pixels
    screenXpixels=608;
    screenYpixels = 342;
    hor_screen_span_deg=140;

    wall_thickness =5;
    wall_offset=params.wall_offset;
    % inner walls [left right top bottom]
    innerwalls = [0+wall_offset(1)+wall_thickness screenXpixels-wall_offset(1)-wall_thickness...
                  0+wall_offset(2)+wall_thickness screenYpixels-wall_offset(2)-wall_thickness];  
    
    xoffset=innerwalls(1);
    yoffset=innerwalls(3);
    
    stim_size = [innerwalls(2)-innerwalls(1)+1, innerwalls(4)-innerwalls(3)+1];
    stim_dir = zeros(nfrst,1,4); %only one rect, 4 vals for rect: (x,y,dx,dy)   
    
    
    w_2=rect_width/2;
    h_2=rect_height/2;
    initPos_hor = [innerwalls(1)+w_2 innerwalls(3)+h_2];
    initPos_ver = [innerwalls(1)+h_2 innerwalls(4)-w_2];
    rectPos=initPos_hor;
    dotDirection = [1 0];
    initDirection= [1 0];

    frame_count=0;        
    repi=0;
    freeseframe=0;
    
    while repi<nrep  
        if freeseframe>=freezeF  %save direction only for the moving rect
            stim_dir(frame_count+1,1,:)=[rectPos(1), rectPos(2),dotDirection(1),dotDirection(2)];
        else %if freezing save position of the rect but no direction
            stim_dir(frame_count+1,1,:)=[rectPos(1), rectPos(2),0, 0];
        end
%         if freeseframe>=freezeF %save direction of motion
%             if dotDirection(1)~=0 %moving horizontally
%                 
%                 hordir(frame_count+1,rectPos(1), rectPos(2))=dotDirection(1);
%             else
%                 verdir(frame_count+1,rectPos(1), rectPos(2))=dotDirection(2);
%             end 
%             stim_dir(frame_count+1,1,:)=[rectPos(1), rectPos(2),dotDirection(1),dotDirection(2)];
%         end

        if freeseframe>=freezeF  
            rectPos=rectPos+dotDirection.*dx;

            %     check for collision with walls
            %     vertical walls   
            if dotDirection(1)~=0 %moving horizontally
                if rectPos(1)+w_2>innerwalls(2) %end of one hor direction        
                        rectPos(1)=rectPos(1)-dotDirection(1).*dx;
                        dotDirection=-dotDirection;                       
                end  
                if rectPos(1)-w_2<innerwalls(1) %end of both directons of hor line       
                        rectPos(1)=rectPos(1)-dotDirection(1).*dx;
                        rectPos(2)=rectPos(2)+rect_height;
                        freeseframe=0;
                        dotDirection=-dotDirection;                       
                        if rectPos(2)+h_2>innerwalls(4) %next hor line would hit ceiling, start moving vertically 
                            rectPos= initPos_ver;                        
                            dotDirection=[0, -1];
                        end
                end    
            else %moving vertically
                if rectPos(2)+w_2>innerwalls(4)                     
                    rectPos(2)=rectPos(2)-dotDirection(2).*dx;
                    dotDirection=-dotDirection;   
                    rectPos(1)=rectPos(1)+rect_height;            
                    freeseframe=0;
                    if rectPos(1)+h_2>innerwalls(2) %end of both ver directions  
                       rectPos= initPos_hor;
                       dotDirection=[1, 0];
                    end
                end        
                if rectPos(2)-w_2<innerwalls(3) %end of one vertical direction       
                    rectPos(2)=rectPos(2)-dotDirection(2).*dx;
                    dotDirection=-dotDirection;            
                end
            end 
        else %after the jump, freese for a few frames
            freeseframe=freeseframe+1;
        end

        frame_count=frame_count+1;   
        if sum(abs(rectPos-initPos_hor))==0 && sum(abs(dotDirection-initDirection))==0 && freeseframe==0
            repi=repi+1;
        end      
    end
    stim_dir(frame_count+1:end,:,:)=[];
    %     to verify make a low res video
%     if debug
%         v=VideoWriter('scanning_pause_rec.avi');
%         open(v);
%         for i=1:nfrst
%             writeVideo(v,uint8(frame_array_lr(:,:,i))*255);
%         end
%         close(v);
%     end
end

function params = read_scanning_rect_with_pause(logfile)

    % read the file with stimuli parameters
    f=fopen(logfile);
    tline=fgetl(f);
    while  ~feof(f) && isempty(tline)
        tline=fgetl(f);
    end

    if isempty(strfind(tline,'scanning_rect'))
        disp('The file does not contain information about scanning rect stimulus');
        fclose(f);
        return;
    end

    tline=fgetl(f); %Start time:
    tline=fgetl(f); %End time:
    params = struct('ifi',-1,'wall_offset',-1,'nrep',-1,'rwidth',-1,'rheight',-1,'dx',-1,'freezefr',-1,'tbefore',-1,'tafter',-1,'nfr',-1);
         
    while ~feof(f)
        tline=fgetl(f);
        params = init_param(tline, params);
    end  
    
    fclose(f);
end

function params = init_param(strline, params)
    semicolon = strfind(strline,':')+1;
    varname=strline(1:semicolon);
    varval=strline(semicolon+1:end);
    
    switch varname
        case 'Flip rate: '
            params.ifi=str2double(varval); %flip rate : frames per second
        case 'Wall offset: '
            nn = strsplit(varval, {' '});
            nn= nn(~cellfun('isempty',nn));  
            params.wall_offset =cellfun(@str2double, nn)'; %wall offset
        case 'Number of repetitions: '
            params.nrep =str2double(varval); 
        case 'Rectangle width: '
            params.rwidth = str2double(varval);
        case 'Rectangle height: '
            params.rheight = str2double(varval);
        case 'Speed: '
            params.dx = str2double(varval);
        case 'Number of frames to freeze: '
            params.freezefr = uint32(str2double(varval));
        case 'Background time before: '
            params.tbefore = str2double(varval);
        case 'Background time after: '
            params.tafter = str2double(varval);
        case 'Frames presented: '
            params.nfr = str2double(varval); 
    end     
end