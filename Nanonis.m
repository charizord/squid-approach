
classdef Nanonis
    
	properties (SetAccess = public)
        Computer=General.Check_Computer();
        instr=General.Find_instr('TCPIP-192.168.236.1'); 
        Zscanner_Length=22;  % Length of extention in um
        XYscanner_Length=40;
        Axis=['0','x','y','z'];
    end
    methods(Static)

        %%   Nanonis comunication
        function ss = Header(mes,BS)
           % This function convert the command to a binar value 
           % Body size (int) (4 bytes) is the size of the message body in bytes.
           commName=pad(sprintf('%X',mes),64,'0');
           bodySizVal=BS;
           bodySize=pad(sprintf('%X',bodySizVal),8,'left','0');
           Res=pad(sprintf('%X',1),4,'left','0');
           Ress=pad(sprintf('%X',0),4,'left','0');
           ss = strcat(commName,bodySize, Res,Ress);
       end
       
        function Send(varargin)
           % Sends the message to Nanonis
           % message format sending: handle to Nanonis connection (nis) ;
           % String with the command to send ; Number of bytes of Body
           % message ; For each argument in Body we send 2 variables (seperated by comma),
           % string with the data type (example 'int32') and the value 
           % in the appropriate format 
           % Example : Nanonis.Send(nis,'Signals.ValGet', 8, 'int', DCchannel, 'uint32',0);
           Nis=Nanonis;
           flushinput(Nis.instr);
           por = Nanonis.Header(varargin{2}, varargin{3});
           for i=2:(nargin-1)/2
               switch varargin{2*i}
                   case 'int'
                       scase = dec2hex(varargin{2*i+1});
                       scaseS = pad(scase,8, 'left', '0');
                   case 'uint16'
                       scase = dec2hex(varargin{2*i+1});
                       scaseS = pad(scase,4, 'left', '0');
                   case 'uint32' 
                       scase = dec2hex(varargin{2*i+1});
                       scaseS = pad(scase,8, 'left', '0');
                   case 'float32'
                       scase = num2hex(single(varargin{2*i+1}));
                       scaseS = pad(scase,8, 'left', '0');
                   case 'float64'
                       scase = num2hex(varargin{2*i+1});
                       scaseS = pad(scase,16, 'left', '0');
                   case 'string'
                       scaseS= [];
                       siz = varargin{2*i-1};
                      % scase= dec2hex(uint32(varargin{2*i+1}));
                      %updated 2021.10.28 in order to use new approach app:
                       scase = dec2hex(uint32(char(varargin{2*i+1})));
                       for j=1:siz
                           scaseS = strcat(scaseS,scase(j,1), scase(j,2));
                       end
               end
               por = strcat(por, scaseS);      
           end
           for i=0:(size(por,2)/2-1)
               cs = strcat(por(2*i+1), por(2*i+2));
               commM(i+1) = uint8(hex2dec(cs));
           end
           fwrite(Nis.instr, commM);     % Send the final message to nanonis    
        end
       
       

        function varargout = Receive(varargin)
            % Receives the message from Nanonis
            % message format sending: handle to Nanonis connection (nis) ; String
            % (or number of string seperated by comma) of the Data type we are
            % receiving ; number of bytes receiving
            % Example: Nanonis.Send(nis,'Motor.PosGet', 0);
            % [aa, bb, cc]  = Nanonis.Receive(nis, 'float64', 'float64', 'float64',24)
            Nis=Nanonis;

            DataReceived = fread(Nis.instr,40+varargin{nargin});
            k = 1;
            dk = 41;
            vc = 1;
            % Convert the data from nanonis to unic
            for i=2:(nargin-1)
               switch varargin{i}
                       case 'uint8' 
                           varargout{vc} = DataReceived(dk);
                           vc = vc+1; dk = dk+1;
                       case 'string'
                           siz = varargout{vc-1};
                           varargout{vc} = char(DataReceived(dk:(dk+siz-1)));
                           vc = vc+1; dk = dk+siz;
                       case 'int'
                           varargout{vc} = DataReceived(dk+3) + DataReceived(dk+2)*2^8 + DataReceived(dk+1)*2^16;
                           if (DataReceived(dk) > 127);
                                varargout{vc} = varargout{vc} + (DataReceived(dk)-256)*2^24;
                           else
                               varargout{vc} = varargout{vc} + DataReceived(dk)*2^24;
                           end
                           vc = vc+1; dk = dk+4;
                       case 'uint16'
                           varargout{vc} = DataReceived(dk+1) + DataReceived(dk)*2^8;
                           vc = vc+1; dk = dk+2;
                       case 'uint32' 
                           varargout{vc} = DataReceived(dk+3) + DataReceived(dk+2)*2^8 + DataReceived(dk+1)*2^16 + DataReceived(dk)*2^24;
                           vc = vc+1; dk = dk+4;
                       case 'float32'
                           mes = dec2hex(DataReceived(dk));
                           for j=1:3
                               mes1 = dec2hex(DataReceived(dk+j));
                               if (mes1 == '0')
                                  mes1 = strcat(mes1,'0');
                               end
                               if (size(mes1,2) == 1)
                                  mes1 = strcat('0',mes1);
                               end
                               mes = strcat(mes,mes1);                     
                           end
                           varargout{vc} = typecast(uint32(hex2dec(mes)),'single');
                           clear mm; clear mes; clear mm; 
                           vc = vc+1; dk = dk+4;
                       case 'float64'
                           mes = dec2hex(DataReceived(dk));
                           for j=1:7
                               mes1 = dec2hex(DataReceived(dk+j));
                               if (mes1 == '0')
                                  mes1 = strcat(mes1,'0');
                               end
                               if (size(mes1,2) == 1)
                                  mes1 = strcat('0',mes1);
                               end
                               mes = strcat(mes,mes1);
                           end
                           varargout{vc} = hex2num(mes);
                           clear mm; clear mes; clear mm; vc = vc+1; dk = dk+8; 
                        case '1D array string'   
                            % had hoc
                           siz = varargout{vc-2};
                           varargout{vc} = char(DataReceived(dk:(dk+siz-1)));
                           vc = vc+1; dk = dk+siz;
                           
                        case '1D array float64'
                           Data_Length= varargout{vc-1};    % allways?
                           Vec=zeros(1,Data_Length);
                           for f=1:Data_Length
                               mes = dec2hex(DataReceived(dk));
                               for j=1:7
                                   mes1 = dec2hex(DataReceived(dk+j));
                                   if (mes1 == '0')
                                      mes1 = strcat(mes1,'0');
                                   end
                                   if (size(mes1,2) == 1)
                                      mes1 = strcat('0',mes1);
                                   end
                                   mes = strcat(mes,mes1);
                               end
                               Vek(f)=hex2num(mes);
                               dk = dk+8;   
                           end
                           varargout{vc} = Vek;
                           vc = vc+1; 
                           
                        case '2D array float32'
                            Data_Size=[varargout{vc-2},varargout{vc-1}];  % Rows and columns
                            M=zeros(Data_Size);
                            for R=1:Data_Size(1)
                                for C=1:Data_Size(2)
                                   mes = dec2hex(DataReceived(dk));
                                   for j=1:3
                                       mes1 = dec2hex(DataReceived(dk+j));
                                       if (mes1 == '0')
                                          mes1 = strcat(mes1,'0');
                                       end
                                       if (size(mes1,2) == 1)
                                          mes1 = strcat('0',mes1);
                                       end
                                       mes = strcat(mes,mes1);                     
                                   end
                                   M(R,C) = typecast(uint32(hex2dec(mes)),'single');
                                   dk = dk+4;
                                end
                            end
                            varargout{vc} = M;
                            vc = vc+1; 
                   end 
                end 
        end    

        %% 
        function [V]= Get(NanonisChannel)
            % Get the voltage on a Nanonis channel
            Nis=Nanonis;
            Nanonis.Send(Nis.instr, 'Signals.ValGet',8,'int',NanonisChannel,'uint32', 1);
            V=Nanonis.Receive(Nis.instr,'float32', 4);
            
            if abs(V)<1e-5
                for i=1:5
                   pause(0.01)
                   Nanonis.Send(Nis.instr,'Pll.FreqShiftGet',4,'int',1);
                   V=Nanonis.Receive(Nis.instr,'float32',4);
                   if abs(V)>1e-5
                       break
                   end
                end
            end
        end
        
        function [Average,s]= Get_Average(NanonisChannel,Avg)
        % Ger an avrage voltage on a Nanonis channel
            Nis=Nanonis;
            V=zeros(1,Avg);
            for i = 1:length(V)
                Nanonis.Send(Nis.instr, 'Signals.ValGet',8,'int',NanonisChannel,'uint32', 1);
                V(i) = Nanonis.Receive(Nis.instr,'float32', 4);
            end
            % in order to clean the read from false values every read
            % bigger than 2 std does not get in to the average
            s = std(V);
            condition = find(abs(V-mean(V))<=2*s);
            Average = mean(V(condition));
        end
        
        function [Time_Of_Scan,Scan_Name]=Scan(Nanonis_Def)
            Nis=Nanonis;
            % Nanonis_Def{Filke_Name,Xcenter[um],Ycenter[um],Width[um],Hight[um],Pixels,Lines,SpeedPerLine,Angle}
            BaseNameNanonis=Nanonis_Def{1};
            Xcenter=Nanonis_Def{2}*1e-6;
            Ycenter=Nanonis_Def{3}*1e-6;
            Width=Nanonis_Def{4}*1e-6;
            Height=Nanonis_Def{5}*1e-6;
            Pixels=Nanonis_Def{6};
            Lines=Nanonis_Def{7};
            Speed=Nanonis_Def{8};
            Angle=Nanonis_Def{9};
            
            %instrument connections
            % Define frame size
            Nanonis.Send(Nis.instr, 'Scan.FrameSet', 20, 'float32',Xcenter,'float32',Ycenter,'float32',Width, 'float32',Height, 'float32',Angle);            
            pause(0.2);
            % Define speed
            Nanonis.Send(Nis.instr,'Scan.SpeedSet', 22, 'float32', 50*1e-6, 'float32', 50*1e-6,'float32', Speed,'float32', Speed, 'uint16', 2, 'float32', 1);        
            pause(0.5);
            %   Set name of scan
            Nanonis.Send(Nis.instr,'Scan.PropsGet',0);
            [qw,qe,qr,qt,mes] = Nanonis.Receive(Nis.instr,'uint32','uint32', 'uint32', 'int', 'string', 30);
            sb = size(BaseNameNanonis,2); 
            if (qt == sb)
                BaseName = strcat(mes','00001');
            else
                num = str2num(mes(sb+1:end)');
                num = num+1;
                snum = pad(num2str(num),5,'left','0');
                BaseName = strcat(mes(1:sb)',snum);
            end
            Nanonis.Send(Nis.instr,'Scan.PropsSet', 20+size(BaseName,2),'uint32', qw,'uint32', qe, 'uint32', qr, 'int' ,size(BaseName,2), 'string', BaseName, 'int', 0);
            
            pause(0.5);
            % Lines and Pixels
            Nanonis.Send(Nis.instr,'Scan.BufferGet',0);
            ncan = Nanonis.Receive(Nis.instr,'int',4);
            switch ncan
                case 1
                    Nanonis.Send(Nis.instr,'Scan.BufferGet',0 );
                    [vb,vf,vrtt, vrr] = Nanonis.Receive(Nis.instr,'int','int', 'int', 'int',16);
                    pause(0.2);
                    Nanonis.Send(Nis.instr,'Scan.BufferSet',16,'int',vb,'int',vf, 'int',Pixels, 'int', Lines);
                case 2
                    Nanonis.Send(Nis.instr,'Scan.BufferGet',0);
                    [vb,vff,vf,Pixss, Linss] = Nanonis.Receive(Nis.instr,'int','int','int', 'int', 'int',20);
                    pause(0.2);
                    Nanonis.Send(Nis.instr,'Scan.BufferSet',20,'int',vb,'int',vff,'int',vf, 'int',Pixels, 'int', Lines);
                case 3
                    Nanonis.Send(Nis.instr,'Scan.BufferGet',0);
                    [vb,vff,vf,vfa,Pixss, Linss] = Nanonis.Receive(Nis.instr,'int','int','int','int', 'int', 'int',24);
                    pause(0.2);
                    Nanonis.Send(Nis.instr,'Scan.BufferSet',24,'int',vb,'int',vff,'int',vf,'int',vfa, 'int',Pixels, 'int', Lines);
            end
            pause(0.5);
            Nanonis.Send(Nis.instr,'Scan.Action',6,'uint16', 0, 'uint32', 1);
            pause(Speed*4)
            Nanonis.Send(Nis.instr,'Scan.Action',6,'uint16', 0, 'uint32', 1);
            Time_Of_Scan=seconds(2*Lines*Speed);    % Not include the time to reach the sratr of the scan
            
            Scan_Name=BaseName;
            
            % delete here?
            
            
        end 

        function Val=If_Scan()
            Nis=Nanonis;
            % Check scan status
            Nanonis.Send(Nis.instr,'Scan.StatusGet',0);
            Val = Nanonis.Receive(Nis.instr, 'uint32',4);
        end
        
        function Check_Scan()
        % The function check if scan and if Vfb jumped
        % If so, it's blink
             try
                while Nanonis.If_Scan()   
                    try
                        Vfb=Nanonis.Get(1); 
                        if abs(Vfb)>9.85
                            DAC.Blink(2)
                        end
                        pause(0.5)
                    catch err
                        disp(err.message)
                    end
                end
                
            catch err
                
                while Nanonis.If_Scan()   
                    try
                        Vfb=Nanonis.Get(1); 
                        if abs(Vfb)>9.85
                            DAC.Blink(2)
                        end
                        pause(0.5)
                    catch err
                        disp(err.message)
                    end
                end
             end
            
        end
        
        
        function Set(NanonisChannel,V)
           % Sets voltage V on chanels AO(2,3,4,8)
            Nis=Nanonis;
            Nanonis.Send(Nis.instr,'UserOut.ValSet',8,'int',NanonisChannel,'float32',V);
        end
        
        function Blink(channel) 
            Nanonis.Set(channel,0)
            pause(3)
            Nanonis.Set(channel,5)
        end
        
        function Noise=Get_Noise(channel)
            
            Length=100;
            N=zeros(1,Length);
            for i=1:Length
               N(i)=Nanonis.Get(channel);
            end
            Max=max(N);
            Min=min(N);
            Noise=abs(Max)-abs(Min);
        end
        
        
        % Error
        %????
        function [X,Y]=BiasSweep(Vi, Vo, steps, period, varargin)
           % Reads voltage form Vi to Vo in steps on chanell chan
%            Vi=0
%            Vo=1
%            steps=201
%            period=30
           
            Nis=Nanonis;            
            V=linspace(Vi,Vo,steps);
            chan = 1;
            Nanonis.Send(Nis.instr,'BiasSwp.Open',0); pause(0.1);
            Nanonis.Send(Nis.instr,'BiasSwp.PropsSet',8, 'uint16', cast(steps, 'uint16'), 'uint16', cast(period, 'uint16'), 'uint16', 0, 'uint16', 0 ); pause(0.1);
            Nanonis.Send(Nis.instr,'BiasSwp.LimitsSet',8, 'float32', Vi, 'float32', Vo); pause(0.1);
            
            Timer=tic;
            Nanonis.Send(Nis.instr,'BiasSwp.Start',24,'uint32',1,'uint32',1,'uint32',0,'int',1,'string','T','uint32',0);
            
            if toc(Timer)<2
            Nanonis.Send(Nis.instr,'BiasSwp.Start',24,'uint32',1,'uint32',1,'uint32',0,'int',1,'string','T','uint32',0);
            end
            
            
            [a,b,c,d,e,Data]=Nanonis.Receive(Nis.instr,'int','int','1D array string','int','int','2D array float32',40+8+28*1+8+4*1*2*(steps+1));
            
            X=Data(1,:);
            Y=Data(2,:);
        end

        function Position = Get_Scanner_Z()
           % Gets Z scanner position
           Nis=Nanonis;
           Nanonis.Send(Nis.instr,'ZCtrl.ZPosGet', 0);
           Position = Nanonis.Receive(Nis.instr, 'float32',4);
           Position=Position.*-10^6;
           
           if Position<1e-3
                for i=1:5
                   pause(0.01)
                   Nanonis.Send(Nis.instr,'ZCtrl.ZPosGet', 0);
                   Position = Nanonis.Receive(Nis.instr, 'float32',4);
                   Position=Position*-10^6;
                   if Position>1e-3
                       break
                   end
                end
           end
           
        end

        function Set_Scanner_Z(Dis)
           % Sets Z scanner position in um
           Nis=Nanonis;
           V=-1e-6*Dis;
           Nanonis.Send(Nis.instr,'ZCtrl.ZPosSet', 4, 'float32', V);
        end
       
        function [Data_f0,Data_df,Data_Y_size,Data_Y]=Spectrum_Get()
            Nis=Nanonis;
            Nanonis.Send(Nis.instr,'SpectrumAnlzr.DataGet',4, 'int',1);
            [Data_f0,Data_df,Data_Y_size]=Nanonis.Receive(Nis.instr,'float64','float64','int',20);
            Nanonis.Send(Nis.instr,'SpectrumAnlzr.DataGet',4, 'int',1);
            [Data_f0,Data_df,Data_Y_size,Data_Y]=Nanonis.Receive(Nis.instr,'float64','float64','int','1D array float64',20+8*Data_Y_size);
        end
        
        function Reset_Spectrum()
            Nis=Nanonis;
            Nanonis.Send(Nis.instr,'SpectrumAnlzr.FFTWindowGet',4, 'int',1);
            S=Nanonis.Receive(Nis.instr, 'uint16',2);
            Nanonis.Send(Nis.instr,'SpectrumAnlzr.FFTWindowSet',6, 'int',1,'uint16',0);
            pause(1)
            Nanonis.Send(Nis.instr,'SpectrumAnlzr.FFTWindowSet',6, 'int',1,'uint16',S);
            pause(1) 
        end
        %% Motors control

        function Motors_Set_Frequency_Voltage(Frequency,Voltage,axis)
           % axis = '0' (all axes), 'x', 'y' or 'z'
            Nis=Nanonis;
            [~,Loc]=find(Nis.Axis==axis);
            Axis = cast(Loc-1,'uint16');  % Index start at zero
            Nanonis.Send(Nis.instr,'Motor.FreqAmpSet',10,'float32',Frequency,'float32',Voltage, 'uint16', Axis);
        end

        function Motors_Move_Steps(N_Step,axis)
           % Motor moves n step in direction s: s = 'x+', 'x-', 'y+', 'y-',
           % 'z+', 'z-'
            Nis=Nanonis;
            switch axis
               case 'x+'
                   a = 0;
               case 'x-'
                   a = 1;
               case 'y+'   % For 4 K system y+ and y- are opposite
                   a = 2;
               case 'y-'
                   a = 3;
               case 'z+'
                   a = 4;
               case 'z-'
                   a = 5;
            end
            Nanonis.Send(Nis.instr,'Motor.StartMove', 14,'uint32', a, 'uint16', cast(N_Step, 'uint16') ,'uint32', 0, 'uint32', 1 );
        end
       
        function Z_enc = Get_Encoder_Z()  % Returns Z_Pos in um
            % Gets Z encoder value
            Nis=Nanonis;
            Nanonis.Send(Nis.instr,'Motor.PosGet', 0);
            [X_pos,Y_pos, Z_enc]=Nanonis.Receive(Nis.instr,'float64','float64','float64',24);
            Z_enc = Z_enc*1e6;
       end
       
       % Fix
       function Position = Get_Scanner_XY(axis)
           % Gets X Y scanner position in FollowMe Motor; Motor =  'x', 'y' 
            Nis=Nanonis;
            Nanonis.Send(Nis.instr,'FolMe.XYPosGet', 4, 'uint32',0 );
            [X,Y] = Nanonis.Receive(Nis.instr,'float64','float64', 16);
            [~,Loc]=find(Nis.Axis==axis);
            
            if Loc==2
                   Position=X.*10^6;
            else
                   Position=Y.*10^6;
            end
            
            if Position<0
               Position=0;
            elseif Position>Nis.XYscanner_Length
               Position=Nis.XYscanner_Length;
            end
           
       end
       
       
       %????
       function DataLog(BaseName, AveragePoints)
            Nis=Nanonis;
            %Opens Data Loger
            Nanonis.Send(Nis.instr,'DataLog.Open',0);
            %Sets chanel(s) to log, channels are AI2 i AI3 
            Nanonis.Send(Nis.instr,'DataLog.ChsSet',12, 'int', 2, 'int', 1, 'int', 2 );
            % Proporties set, 
            Nanonis.Send(Nis.instr,'DataLog.PropsSet',34+size(BaseName,2),'uint16', 1, 'int' , -1, 'int', -1, 'float32', -1, 'int', AveragePoints, 'int', size(BaseName,2), 'string', BaseName, 'int', 0, 'int', 0, 'int', 0  );
            %Starts data Loging
            Nanonis.Send(Nis.instr,'DataLog.Start',0);
       end
       %????
       function DataLogStop()
            Nis=Nanonis;
            %Stops Data Loger
            Nanonis.Send(Nis.instr,'DataLog.Stop',0);
       end    
        %% OscillationControl
        function OutputOn()   
             Nis=Nanonis;
             Nanonis.Send(Nis.instr,'PLL.OutOnOffSet',8,'int', 1,'uint32',1);
        end
        
        function OutputOff()   
             Nis=Nanonis;
             Nanonis.Send(Nis.instr,'PLL.OutOnOffSet',8,'int', 1,'uint32',0);
        end
        
        function PllOn()
            %set Pll controller on
            Nis=Nanonis;
            Nanonis.Send(Nis.instr,'PLL.AmpCtrlOnOffSet',8,'int', 1,'uint32',1);
            Nanonis.Send(Nis.instr,'PLL.PhasCtrlOnOffSet',8,'int', 1,'uint32',1);
        end
        
        function PllOff()
            %set Pll controller off
            Nis=Nanonis;
            Nanonis.Send(Nis.instr,'PLL.AmpCtrlOnOffSet',8,'int', 1,'uint32',0);
            Nanonis.Send(Nis.instr,'PLL.PhasCtrlOnOffSet',8,'int', 1,'uint32',0);
        end
        
        function fshift=Getfshift()
            Nis=Nanonis;
            Nanonis.Send(Nis.instr,'Pll.FreqShiftGet',4,'int',1);
            fshift=Nanonis.Receive(Nis.instr,'float32',4);

            if abs(fshift)<2e-5
                for i=1:5
                   pause(0.01)
                   Nanonis.Send(Nis.instr,'Pll.FreqShiftGet',4,'int',1);
                   fshift=Nanonis.Receive(Nis.instr,'float32',4);
                   if abs(fshift)>2e-5
                       break
                   end
                end
            end
        end 
        
        
        function touch_point=Poke(retract,zlimit,avg,step,Sigma_Multiplay,Mu_avg)

%             Poke approach until freq. shift crosses threshold. return
%             touch point and retract
            if ~exist('zlimit','var')
                zlimit=22;
            end
            if ~exist('avg','var')
                avg=3;
            end
            if ~exist('step','var')
                step=0.001;
            end
            
            if ~exist('Sigma_Multiplay','var')
                Sigma_Multiplay=6;
            end
            
            if ~exist('Mu_avg','var')
                Mu_avg=100;
            end
            
%             turn on OC
            Nanonis.OutputOn();
            Nanonis.PllOn();
            pause(3);
%             Zscanner step:
%             
%             average 1000 points and obtain mean and SD
            fshift=zeros(1,Mu_avg);
            for i=1:length(fshift)
                fshift(i)=Nanonis.Getfshift();
                pause(0.1);
            end
            
            mu=mean(fshift);
            sigma=std(fshift);
            
%             set limits
            upp_lim=mu+Sigma_Multiplay*sigma;
            low_lim=mu-Sigma_Multiplay*sigma;
            
%             open figure and set limits
            figure;
            title('TF poke approach')
            xlabel('time(s)')
            ylabel('freq. shift')
            u=yline(upp_lim,'-.r');
            l=yline(low_lim,'-.r');
            h=animatedline;
            tic;

            while 1
                Zs=Nanonis.Get_Scanner_Z();
                
%                 perform averaging on freq. shift value
                f_avg=zeros(1,avg);
                for i=1:avg
                    f_avg(i)=Nanonis.Getfshift();
                    pause(0.1)
                end
                f=mean(f_avg);
                addpoints(h,double(toc),double(f));
                
                condition1=or(f>upp_lim, f<low_lim);
                condition2=Zs>=zlimit;
                if condition1
                    Nanonis.Set_Scanner_Z(Zs-retract);
                    disp('crossed freq. shift threshold, retracted '+string(retract)+'um , touch-point='+string(Zs));
                    touch_point=Zs;
                    
                    break
                end
                if condition2
                    Nanonis.Set_Scanner_Z(Zs-retract);
                    disp('reached zlimit');
                    touch_point=0;
                    break
                end

%                 update mean
                fshift=circshift(fshift,-1);
                fshift(end)=f;
                mu=mean(fshift);
                upp_lim=mu+Sigma_Multiplay*sigma;
                
                delete(u);
                u=yline(upp_lim,'-.r');
                low_lim=mu-Sigma_Multiplay*sigma;
                delete(l);
                l=yline(low_lim,'-.r');
                 
                Nanonis.Set_Scanner_Z(Zs+step);
                pause(0.01)
            end
            Nanonis.PllOff();
            Nanonis.OutputOff();
             
        end

      function touch_point = Alt_Poke(retract_nm, speed_nm_s, scanner_limit, STD, Sigma_no, mu_avg)
            
            set_ZCtrl(10,1e-6,speed_nm_s,retract_nm) % Reseting Z-Controller
            Nanonis.SafeTip_SetThreshold(Sigma_no*STD*1e-3);
            Nanonis.OutputOn();
            Nanonis.PllOn();
            pause(1)
            fshift = construct_mu(mu_avg,3); % Initial mu calculation

            pause(3);
            Nanonis.SafeTip_SetOnOff(1) % SafeTip On

            Nanonis.ZCtrl_SetOnOff(1) % Turning Z-Controller On
            disp('Z-Controller On')
            pause(0.1)

            % Main Loop

            while 1
                drawnow

                if Nanonis.Get_Scanner_Z() >= scanner_limit % Scanner Limit Reached
                    Nanonis.ZCtrl_SetOnOff(0)
                    touch_point = 0;
                        disp('Scanner Limit reached. Retracting Tip')
                        Nanonis.Set_Scanner_Z(Nanonis.Get_Scanner_Z - retract_nm*10^-3);
                    % The following if condition is for the different modes

                    break
                end

                if Nanonis.ZCtrl_GetOnOff()==0 % ZCtrl is Off (Assuming SafeTip Triggered)
                    pause(0.2)
                    touch_point = double(Nanonis.Get_Scanner_Z()) + retract_nm*10^-3; % in um;
                    disp(char(strcat('SafeTip Triggered, touch point=',string(touch_point))))
                    General.Beep();
                    break
                end
                fshift = update_mu(fshift,3); % Updating Average
        %         ZScannerValue = round(double(Nanonis.Get_Scanner_Z()),3); % Field

            end

        %     ZScannerValue = round(double(Nanonis.Get_Scanner_Z()),3);
            Nanonis.SafeTip_SetOnOff(2); %turn off safe tip
            Nanonis.PllOff();
            Nanonis.OutputOff();

            function set_ZCtrl(Ctrl_index,PI_Const,tip_speed,retract)

                Nanonis.ZCtrl_SetCtrl(Ctrl_index) % Choosing specific controller
                Nanonis.ZCtrl_SetOnOff(0) % Making sure controller is off

                Nanonis.ZCtrl_SetPoint(PI_Const) % Sets Setpoint 
                % (should be same value as time constant and very small, as it may cause a jump at the start of the approach)
                Nanonis.ZCtrl_SetGain(tip_speed*10^-9,PI_Const,0) % Tip Speed and time constant (same as setpoint)
                Nanonis.SafeTip_SetOnOff(2) % Turning SafeTip Off
                Nanonis.ZCtrl_SetHome(retract*10^-9) % Setting home relative and retract amount

            end
            
            function fshift = construct_mu(mu_avg,mu_ch)

                fshift = zeros(mu_avg,1);
                for i=1:mu_avg
                    fshift(i) = Nanonis.Getfshift();
                    pause(0.01)
                end
                mu = mean(fshift);
                Nanonis.Set(mu_ch,mu)

            end

            function new_fshift = update_mu(fshift,mu_ch)

                new_fshift = circshift(fshift,-3);
                for i=1:3
                    new_fshift(end + i - 3) = Nanonis.Getfshift();
                    pause(0.01)
                end
                mu = mean(new_fshift);
                Nanonis.Set(mu_ch,mu)

            end


        end
        
        function TF_Guard(retract)
            
            %This function help to scan save, if there is a jump in the TF
            %it retracts
            Nis=Nanonis;
            %turn on OC
            Nanonis.OutputOn();
            Nanonis.PllOn();
            pause(3);
            
            %average 100 points and obtain mean and SD
            fshift=zeros(1,100);
            for i=1:length(fshift)
                fshift(i)=Nanonis.Getfshift();
                pause(0.1);
            end
            mu=mean(fshift);
            sigma=std(fshift);
            
            %set limits
            upp_lim=mu+6*sigma;
            low_lim=mu-6*sigma;
            
            %open figure ant set limits
            figure;
            title('TF Guard')
            xlabel('time(s)')
            ylabel('freq. shift')
            u=yline(upp_lim,'-.r');
            l=yline(low_lim,'-.r');
            h=animatedline;
            tic;

            while 1
                %perform averaging on freq. shift value
                f=Nanonis.Getfshift();
                
                
                if or(f>upp_lim, f<low_lim)
                    Nanonis.Send(Nis.instr,'Scan.Action',6,'uint16', 2, 'uint32', 1);  % Stop the scan
                    Zs=Nanonis.Get_Scanner_Z();
                    Nanonis.Set_Scanner_Z(Zs-retract);
                    General.Beep()
                    disp('TF alarm');
                    addpoints(h,double(toc),double(f));
                    break
                end
                addpoints(h,double(toc),double(f));
                %update mean
                fshift=circshift(fshift,-1);
                fshift(end)=f;
                mu=mean(fshift);
                upp_lim=mu+6*sigma;
                low_lim=mu-6*sigma;
                delete(u);
                delete(l);
                u=yline(upp_lim,'-.r');
                l=yline(low_lim,'-.r');
                 
                pause(0.01)
            end
            Nanonis.PllOff();
            Nanonis.OutputOff();
             
        end
        
        function DemodGet()
            Nis=Nanonis;
            Nanonis.Send(Nis.instr,'Pll.DemodInputGet',2,'uint16',1);
            [i,f]=Nanonis.Receive(Nis.instr,'uint16',4,'uint16',4);
        end
        function DemodFilterGet()
            Nis=Nanonis;
            Nanonis.Send(Nis.instr,'LockIn.DemodLPFilterGet',2,'uint16',1);
            [order,cutoff]=Nanonis.Receive(Nis.instr,'int',4,'float32',8);
        end
        
         %% Z-Controller
        
        function ZCtrl_SetPos(Pos)
            % Sets the Z position of the tip
            Nis=Nanonis;
            Nanonis.Send(Nis.instr,'ZCtrl.ZPosSet',4,'float32',Pos);
        end
        
        function Pos = ZCtrl_GetPos()
            % Returns the current Z position of the tip
            Nis=Nanonis;
            Nanonis.Send(Nis.instr,'ZCtrl.ZPosGet',4,'uint32',1);
            Pos = Nanonis.Receive(Nis.instr,'float32', 4);
        end
        
        function ZCtrl_SetOnOff(bool)
            % Turns Z-Controller on or off
            Nis=Nanonis;
            Nanonis.Send(Nis.instr, 'ZCtrl.OnOffSet',4,'uint32',bool);
        end
        
        function bool = ZCtrl_GetOnOff()
            % Returns the status of the Z-Controller
            Nis=Nanonis;
            Nanonis.Send(Nis.instr,'ZCtrl.OnOffGet',4,'uint32',1);
            bool = Nanonis.Receive(Nis.instr,'uint32', 4);
        end
        
        function ZCtrl_SetPoint(setpoint)
            % Sets the setpoint of the Z-Controller
            Nis=Nanonis;
            Nanonis.Send(Nis.instr, 'ZCtrl.SetpntSet',4,'float32',setpoint);
        end
        
        function ZCtrl_SetHome_Raw(a, Home)
            % WARNING! 
            % This function for some reason moves the Z Scanner
            % the 'Home' amount. Do not use this function, instead
            % use the one below
            
            % Sets the current status of the Z-Controller Home switch
            % (Absolute or Relative) and its corresponding position.
            % a = 1 : Abslute    a = 2 : Relative
            Nis=Nanonis;
            Nanonis.Send(Nis.instr, 'ZCtrl.HomePropsSet',6,'uint16',a,'float32', Home);
        end
        
        function ZCtrl_SetHome(Home)
            app.Z_Pos = -Nanonis.Get_Scanner_Z(); % Getting Position
            Nanonis.ZCtrl_LimitsEnabledSet(1) % Enabling Limits
            Nanonis.ZCtrl_LimitsSet(app.Z_Pos*10^-6,app.Z_Pos*10^-6) % Fixing scanner in place
            Nanonis.ZCtrl_SetHome_Raw(2, Home) % Setting home relative and retract amount
            Nanonis.ZCtrl_LimitsSet(0, -22*10^-6) % Zeroing Limits
            Nanonis.ZCtrl_LimitsEnabledSet(0) % Disabling Limits
        end
        
        function ZCtrl_SetGain(P,T,I)
            % Sets the setpoint of the Z-Controller
            Nis=Nanonis;
            Nanonis.Send(Nis.instr, 'ZCtrl.GainSet',12,'float32',P,'float32',T,'float32',I);
        end
        
        function ZCtrl_SetCtrl(n)
            % Sets the active Z-Controller.
            Nis=Nanonis;
            Nanonis.Send(Nis.instr, 'ZCtrl.ActiveCtrlSet',4,'int',n);
        end
        
        function ZCtrl_Withdraw
            % Switches off the Z-Controller and then fully withdraws the tip
            Nis=Nanonis;
            Nanonis.Send(Nis.instr, 'ZCtrl.Withdraw',8,'uint32',0,'int',0);
        end
        
        function ZCtrl_LimitsEnabledSet(bool)
            % Enables or disables the Z position limits.
            Nis=Nanonis;
            Nanonis.Send(Nis.instr, 'ZCtrl.LimitsEnabledSet',4,'uint32',bool);
        end
        
        function ZCtrl_LimitsSet(high, low)
            % Sets the Z position high and low limits in meters.
            Nis=Nanonis;
            Nanonis.Send(Nis.instr, 'ZCtrl.LimitsSet',8,'float32',high, 'float32', low);
        end
        
%         function Ctrl_List = ZCtrl_GetCtrl()
%             % Returns the list of Z-Controllers and the index of the active controller
%             Nis=Nanonis;
%             Nanonis.Send(Nis.instr, 'ZCtrl.CtrlListGet',4,'uint32',1);
%             Ctrl_List = Nanonis.Receive(Nis.instr,'int',4,'int',4,'1D array string',4*44,'int',4);
%         end
        %% Safe Tip
        function m = SafeTip_SetOnOff(a)
            % Switches the Safe Tip feature on or off.
            % a = 1 : On   a = 2 : Off
            Nis=Nanonis;
            Nanonis.Send(Nis.instr, 'SafeTip.OnOffSet',2,'uint16',a);
        end
        
        function [V]= SafeTip_GetOnOff()
            % Get the Saftip Status
            Nis=Nanonis;
            Nanonis.Send(Nis.instr, 'SafeTip.OnOffGet',4,'uint32', 1);
            V=Nanonis.Receive(Nis.instr,'uint16', 2);
        end
        
        function SafeTip_SetThreshold(Threshold)
            % Sets the Safe Tip configuration.
            % Two first inputs are irrelevant, last input is for threshold
            Nis=Nanonis;
            Nanonis.Send(Nis.instr, 'SafeTip.PropsSet',8,'uint16',0,'uint16',0,'float32',Threshold);
        end
                %% Data Logger
        function DataLog_Open()
            % Opens the Data Logger module.
            Nis=Nanonis;
            Nanonis.Send(Nis.instr,'DataLog.Open',0);
        end
        
        function DataLog_Start()
            % Starts the acquisition in the Data Logger module.
            Nis=Nanonis;
            Nanonis.Send(Nis.instr,'DataLog.Start',0);
        end
        
        function DataLog_Stop()
            % Stops the acquisition in the Data Logger module.
            Nis=Nanonis;
            Nanonis.Send(Nis.instr,'DataLog.Stop',0);
        end
        
        function DataLog_ChsSet(ch_index)
            % Sets the list of recorded channels in the Data Logger module.
            % This function will record only one channel because I dind't
            % put the effort to code it properly (Maor)
            % The channel Index is as assigned in the Signals Manager
            Nis=Nanonis;
            Nanonis.Send(Nis.instr,'DataLog.ChsSet', 8, 'int', 1, 'int', ch_index);
        end
        
        function DataLog_PropsSet(BaseName, Acq_dur_s)
            % Sets the acquisition configuration and the save options in the Data Logger module.
            % I've set alot of the parameters in advance. for a full explanation, check
            % TCPProtocl: DataLog.PropsSet
            Nis=Nanonis;
            Nanonis.Send(Nis.instr,'DataLog.PropsSet',34+strlength(BaseName),'uint16', 2, 'int' , 0, 'int', 0, 'float32', Acq_dur_s, 'int', 1, 'int', strlength(BaseName), 'string', BaseName, 'int', 0, 'int', 0, 'int', 0  );
        end
        
        %% Utilities
        
        function Uti_RTOversamplSet(n)
            % Sets the Real-time oversampling in the TCP Receiver.
            Nis=Nanonis;
            Nanonis.Send(Nis.instr,'Util.RTOversamplSet', 4, 'int', n);
        end
%         
        function SF=Get_Session_Folder()
            Nis=Nanonis;
            Nanonis.Send(Nis.instr,'Util.SessionPathGet',0);
            path_size=Nanonis.Receive(Nis.instr,'int',4);
            Nanonis.Send(Nis.instr,'Util.SessionPathGet',0);
            [~,path]=Nanonis.Receive(Nis.instr,'int','string',4+path_size);
            SF=strcat(path');
        end
        
        
        function Check_Vfb_Lock()
            % Check if Vfb jump if so blink
            try
                while Nanonis.If_Scan()   
                    try
                        Vfb=Nanonis.Get(1); 
                        if abs(Vfb)>9.85
                            DAC.Blink(2)
                        end
                        pause(0.5)
                    catch err
                        disp(err.message)
                    end
                end
            catch err
                while Nanonis.If_Scan()   
                    try
                        Vfb=Nanonis.Get(1); 
                        if abs(Vfb)>9.85
                            DAC.Blink(2)
                        end
                        pause(0.5)
                    catch err
                        disp(err.message)
                    end
                end
            end
            
        end
        
   end
end