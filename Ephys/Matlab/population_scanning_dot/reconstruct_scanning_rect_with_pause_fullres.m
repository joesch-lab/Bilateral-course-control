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
    
    while repi<1%nrep  
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
