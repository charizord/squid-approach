classdef Squid_Functions
    
	properties (SetAccess = public)
        Computer=General.Check_Computer();
    end
    
    methods(Static)
        
        function [Rmb]=R_Magic_Box(V)
            R=V/1e-5;
            if R>1531 || R<1302
                disp('Wrong input, please enter a value between 1362-1531')
            else      
                Rmb=((1020-R)*1531+523602)/(R-1531);
            end
        end
        
        function [SOT_Diameter] = Tip_Diameter(Field_Period)
            SOT_Diameter = 2000*sqrt(20/(3.14*Field_Period));
        end
        
        function [Field_Period] = Field_Period(TipDiameter)
            Field_Period = 80*10^6/(3.14*(TipDiameter^2));
        end
        
        function [Rp,Isotmax]=Iv(Instrument,Output_Chan,Input_Chan,Vf,Rbias,Rfb,Steps)
            SOT=Squid_Functions;
            switch SOT.Computer
                case '1.5K'
                    Turn_Ratio=3.74;   % Update when know or take from file
                    Rshunt=1; %[Ohm] % Per system
                case '4K'
                    Turn_Ratio=2;
                    Rshunt=1; %[Ohm] % Per system
            end

            if ~ exist('Rbias','var')
                Rbias=2550; %[Ohm]
            end
            if ~ exist('Rfb','var')
                Rfb=11000; %[Ohm]
            end
            if ~ exist('Steps','var')
                Steps=201; 
            end
            
            V_Bias=linspace(0,Vf,Steps);
            switch Instrument
                case 'DAC'
                    DAC.Set(0,0);
                    V_fb = DAC.Ramp(Output_Chan, Input_Chan, 0, Vf,Steps, 1,30);
                    DAC.Set(0,0);
                case 'Nanonis'
                    [~,V_fb] = Nanonis.BiasSweep(0, Vf,Steps-1,10);
                case 'Keithley'
                    % Future
            end
            I_SOT=(V_fb-V_fb(1))./(Turn_Ratio*Rfb);
%             % Plot

%             
            IV=figure(666);
            set(IV,'Name','IV','NumberTitle','off')
            hold on
            grid on
            plot(V_Bias,I_SOT.*1e6)
            xlabel('V_b_i_a_s [V]')
            ylabel('I_{SOT} [\muA]')
            title('Iv')
            
            I_Bias=V_Bias./Rbias; %[A]
            I_Shunt=I_Bias-I_SOT; %[A]
            Point=floor(Steps/4);
            Rp=I_Shunt(Point)/I_SOT(Point)*Rshunt;
            Isotmax=max(I_SOT.*1e6);
        end
        
        function [Excaitation]=TF_Coupling(Bdc,Dis,Bac,Sens)
            
            % This function calculate the excaitation of the tip in nm
            Bac_Tag=(Bac/10)*Sens;
            Excaitation=(Bac_Tag/Bdc)*Dis*1e3; 

        end
        %% 04/02 Tip calibration:
        
        function [Fields,Vfb] = Tip_Calibration_Keithley(Initial_Field,Final_Field,Number_OF_Fields,keit_mag,keitVbias_num)
            
            % If not bias from keithley?
            Vbias = Keithley.Get(keitVbias_num);
                
            % Save
            Mag_Vfb_file=General.Add_Time2File(['Response',num2str(Vbias)])
            Mag_Vfb_Title={'Field[G]','Vfb[V]'};
            General.File_make(Mag_Vfb_file,Mag_Vfb_Title);
            
            NanonisChannel=1;
            Avg=100;
            
            Fields=linspace(Initial_Field,Final_Field,Number_OF_Fields);
            Vfb=Fields.*0;
            DAC.Blink(2)
            
            figure()
            hold on
            xlim([min(Fields)-1 max(Fields)+1])
            xlabel('Field [G]')
            ylabel('Vfb [V]')
            title(strcat('Vfb Vs Field at Vbias=',num2str(Vbias),'[V]'))
            grid on
            
            for i=1:length(Fields)
                    Magnet.Set_Keithley(Fields(i),keit_mag)
                    Vfb(i)= Nanonis.Get_Average(NanonisChannel,Avg);
                    if abs(Vfb(i))>9 %if lost lock
                        DAC.Blink(instruments,2)
                        Vfb(i)= Nanonis.Get_Average(NanonisChannel,Avg);
                    end
                    plot(Fields(i),Vfb(i),'.')
                    General.Save(Mag_Vfb_file,{Fields(i),Vfb(i)})
                    percentage_of_run = round(100*i./Number_OF_Fields)
            end
            Magnet.Control_Keithley(instruments,0,keit_mag)
            hold off
            cftool(Fields,Vfb)
        end

        
        function Noise_Characterization(Vb_src,Vb_Chan,Vfb_src,Vfb_Chan,Vfb_Average,Mag_instr,Fi,Ff,Fn,Vi,Vf,Vn,Blink)
            % Vb_Src = 'DAC','Nanonis','Keithley'
            % Vfb_Src = 'DAC','Nanonis','Keithley'
            % Mag_instr= 1 - PS, 2-5 Keithley
            
            if ~exist('Blink','Var')
               Blink=1; 
            end
            % Create file
            d=datestr(datetime);
            File_Name=[d(13:14),'-',d(16:17),'-',d(19:20),'_Noise_Characterization'];
            Noise_Characterization_file=General.Add_Time2File(File_Name);
            Mag_Vfb_Title={'Vb[V]','Field[G]','Vfb[V]','SpecDataF0','SpecDatadF','SpecDataStep','SpecDataY'};  % Think
            General.File_Make(Noise_Characterization_file,Mag_Vfb_Title);
            
            Fields=linspace(Fi,Ff,Fn);
            Voltages=linspace(Vi,Vf,Vn);
            
            if Blink
                DAC.Blink(2);
            end
            
            for j=1:Fn     % Go over all Fields
                % Set the field
                Magnet.Set(Mag_instr,Fields(j),0); 
                pause(1)
                if Blink
                    DAC.Blink(2);
                end
                for i=1:Vn % Go over all Vbias
                    switch Vb_src
                        case 'DAC'
                            DAC.Set(Vb_Chan,Voltages(i));
                        case 'Nanonis'
                            Nanonis.Set(Vb_Chan,Voltages(i))
                        case 'Keithley'
                            Keithley.Set(Vb_Chan,Voltages(i))
                    end
                    if Blink
                        DAC.Blink(2);
                    end
                    Nanonis.Reset_Spectrum()  % Reset the Nanonis spectrom
                    pause(20)
                % Measure Vfb depends in instrument
                    switch Vfb_src
                        case 'DAC'
                            Vfb=DAC.Get_Average(Vfb_Chan, Vfb_Average);
                        case 'Nanonis'
                            Vfb=Nanonis.Get_Average(Vfb_Chan,Vfb_Average);
                        case 'Keithley'
                            Vfb=Keithley.Get_Average(Vfb_Chan,Vfb_Average);
                    end
                % Take spectrum from Nanonis - future from 
                    [Data_f0,Data_df,Data_Y_size,Data_Y]=Nanonis.Spectrum_Get();
                    General.Save(Noise_Characterization_file,{Voltages(i),Fields(j),Vfb,Data_f0,Data_df,Data_Y_size,Data_Y});   %save
                end
                disp(['percentage of run is: ', num2str(round(100*(j./Fn))),' %'])
            end
            
%             
%             
%             for i=1:Vn  % go over all V bias
%                 % Set Vb depends in instrument
%                 switch Vb_src
%                     case 'DAC'
%                         DAC.Set(Vb_Chan,Voltages(i));
%                     case 'Nanonis'
%                         Nanonis.Set(Vb_Chan,Voltages(i))
%                     case 'Keithley'
%                         Keithley.Set(Vb_Chan,Voltages(i))
%                 end
%                 for j=1:Fn     % Go over all Fields
%                     % Set the field
%                     Magnet.Set(Mag_instr,Fields(j),0); 
%                     pause(1)
%                     if Blink
%                         DAC.Blink(2);
%                     end
%                     Nanonis.Reset_Spectrum()  % Reset the Nanonis spectrom
% 
%                     % Measure Vfb depends in instrument
%                     switch Vfb_src
%                         case 'DAC'
%                             Vfb=DAC.Get_Average(Vfb_Chan, Vfb_Average);
% %                             Vfb=DAC.Get_Average(Vb_Chan,Vfb_Chan);
%                         case 'Nanonis'
%                             Vfb=Nanonis.Get_Average(Vfb_Chan,Vfb_Average);
%                         case 'Keithley'
%                             Vfb=Keithley.Get_Average(Vb_Chan,Vfb_Average);
%                     end
%                     % Take spectrum from Nanonis - future from 
%                     [Data_f0,Data_df,Data_Y_size,Data_Y]=Nanonis.Spectrum_Get();
%                     General.Save(Noise_Characterization_file,{Voltages(i),Fields(j),Vfb,Data_f0,Data_df,Data_Y_size,Data_Y});   %save
%                 end
%                 disp(['percentage of run is: ', num2str(round(100*i./Vn)),' %'])
%             end
            
            % Set bias to zero
            switch Vb_src
                case 'DAC'
                    DAC.Set(Vb_Chan,0);
                case 'Nanonis'
                    Nanonis.Set(Vb_Chan,0);
                case 'Keithley'
                    Keithley.Set(Vb_Chan,0);
            end
            Magnet.Set(Mag_instr,0,1); % Set mag to zero

        end

        
        function [Fields,Vfb] = Noise_Characterization_Nanonis(Initial_Field,Final_Field,Number_OF_Fields,Vi,Vf,Vn,Mag_instr,V_src,Vfb_Read)


            % Create file
            Noise_Characterization_file=General.Add_Time2File('Noise_Characterization')
            Mag_Vfb_Title={'Field[G]','Vb[V]','Vfb[V]','Noise_Amp'};  % Think
            General.File_make(Noise_Characterization_file,Mag_Vfb_Title);
            
            
            Avg=100;
            NanonisChannel=1;
            
            Fields=linspace(Initial_Field,Final_Field,Number_OF_Fields);
            Voltages=linspace(Vi,Vf,Vn);
            
            Vfb_Avg=zeros(Number_OF_Fields,Vn);
            Vfb_Min=zeros(Number_OF_Fields,Vn);
            Vfb_Max=zeros(Number_OF_Fields,Vn);
            Vfb_Noise=zeros(Number_OF_Fields,Vn);
            
            DAC.Blink(2)
            

            for i=1:Number_OF_Fields    % go over all fields
                
                % Set the field

                Magnet.Set(Mag_instr,Fields(i),0);
                
                for j=1:Vn     % Go over all voltages
                    % Set the voltage
                    Keithley.Set(V_src,Voltages(j));    % Case Keithelly
                    
                    V_Meas=zeros(Avg,1);
                    for i=1:Avg
                        V_Meas(i)=Nanonis.Get(NanonisChannel)
                        if abs(V_Meas(i))>9 %if lost lock
                            DAC.Blink(2)
                            V_Meas(i)=Nanonis.Get(NanonisChannel)   % Case Nanonnis
                        end
                    end
                    Vfb_Avg(i,j)=mean(V_Meas);
                    Vfb_Min(i,j)=min(V_Meas);
                    Vfb_Max(i,j)=max(V_Meas);
                    Vfb_Noise(i,j)=Vfb_Max-Vfb_Min(i,j);  % abs???
                    
                    General.Save(Noise_Characterization_file,{Fields(i),Voltages(j),Vfb_Avg(i,j),Vfb_Min(i,j),Vfb_Max(i,j),Vfb_Noise(i,j)})
                end
                disp(['percentage of run is', num2str(round(100*i./Number_OF_Fields))])
            end
            
            [R_Min,C_Min]=find(min(min(Vfb_Noise))==Vfb_Noise);
            
            disp('Lowest noise value for')
            disp('')
            disp(['Field: ',num2str(Fields(R_Min)),' [G]'])
            disp('')
            disp(['Voltage: ',num2str(Voltages(C_Min)),' [V]'])
            
            % Set mag to zero
            
            Magnet.Set(Mag_instr,0,1);
            
        end
        % Response
        function Response = GetResponse(Mag_Inst,Center_Field, Delta_Field, NanonisChannel,Mode)
            % Gets the response at field Center_Field with steps Delta_Field
            
            if ~exist('Mode','Var')
               Mode=0; 
            end
            Length=3;
            Vfb = zeros(1,Length);
            Magn=zeros(1,Length);
            Field=[Center_Field+Delta_Field,Center_Field,Center_Field-Delta_Field];
            
            for k=1:Length
                Magnet.Set(Mag_Inst,Field(k),0);
                Magn(k) = Magnet.Get(Mag_Inst);
                pause(5)
                Vfb(k) = Nanonis.Get_Average(NanonisChannel,10);
            end
            
            % Set the currents in the leads to zero
            Magnet.Set(Mag_Inst,Center_Field,Mode);
            
            [a,b] =  polyfit(Magn,Vfb,1);   
            Response = a(1);
        end
        
        function [Vbmax, FieldMax] = FindMaxResponse (VbI, VbF, VbS, BI, BF, BS,dh, NanonisChannel)
            
            % Finds the point [Vbias, field] with the maximum response and plots the response figure        
            Vbias=linspace(VbI,VbF,VbS)
            FieldB=linspace(BI,BF, BS)
            Vbmax = Vbias(1);
            FieldMax = FieldB(1);
            % data matrix consists of dV/dB data
            data = zeros(size(Vbias,2), size(FieldB,2));
            %instrument connections
            t = instrfind('Type', 'tcpip');
            Ket = instrfind('Type', 'gpib', 'BoardIndex', 7, 'PrimaryAddress', 27, 'Tag', '');
            for j = 1:size(FieldB,2)
                sF = FieldB(j);
                m = Magnet.Control(sF,0);
                Nanonis.Send(t, 'Signals.ValGet',8,'int',NanonisChannel,'uint32', 1);
                pomV = Nanonis.Receive(t,'float32', 4);
                for i=1:size(Vbias,2)
                    Vb = [];   
                    % 
                    dataSB(i,j) = sF;
                    command = strcat(":sour:volt:lev ",num2str(Vbias(i)));
                    fprintf(Ket,command);            
                    m = Magnet.Control(sF+1*dh,0);
                    Nanonis.Send(t, 'Signals.ValGet',8,'int',NanonisChannel,'uint32', 1);
                    Vb(1) = Nanonis.Receive(t,'float32', 4);                    
                    Vb(2) = pomV;
                    m = Magnet.Control(sF-dh,0);
                    Nanonis.Send(t, 'Signals.ValGet',8,'int',NanonisChannel,'uint32', 1);
                    Vb(3) = Nanonis.Receive(t,'float32', 4);              
                    data(i,j) = (Vb(3)-Vb(1))/(2*dh);
                    if (abs(data(i,j)) == max(abs(data)))
                        Vbmax = Vbias(i);
                        FieldMax = FieldB(j);
                    end
                    
                end
            end
        figure (23)  
        surf(Vbias,FieldB,data)
        view(2)
        shading interp;
        colormap jet
        fsize = 20;
        colorbar EastOutside
        h=colorbar; 
        end
        
        function [Rp,Isotmax,I_SOT,V_Bias]=IvSave(Instrument,Output_Chan,Input_Chan,Vf,Rbias,Rfb,Steps,Folder)
            SOT=Squid_Functions;
            switch SOT.Computer
                case '1.5K'
                    Turn_Ratio=3.74;   % Update when know or take from file
                    Rshunt=1; %[Ohm] % Per system
                case '4K'
                    Turn_Ratio=2;
                    Rshunt=1; %[Ohm] % Per system
            end

            if ~ exist('Rbias','var')
                Rbias=2550; %[Ohm]
            end
            if ~ exist('Rfb','var')
                Rfb=11000; %[Ohm]
            end
            if ~ exist('Steps','var')
                Steps=201; 
            end
            
            V_Bias=linspace(0,Vf,Steps);
            switch Instrument
                case 'DAC'
                    DAC.Set(0,0);
                    V_fb = DAC.Ramp(Output_Chan, Input_Chan, 0, Vf,Steps, 1,30);
                    DAC.Set(0,0);
                case 'Nanonis'
                    [~,V_fb] = Nanonis.BiasSweep(0, Vf,Steps-1,10);
                case 'Keithley'
                    % Future
            end
            I_SOT=(V_fb-V_fb(1))./(Turn_Ratio*Rfb);
%             % Plot

%             
            IV=figure(666);
            set(IV,'Name','IV','NumberTitle','off')
            hold on
            grid on
            plot(V_Bias,I_SOT.*1e6)
            xlabel('V_b_i_a_s [V]')
            ylabel('I_{SOT} [\muA]')
            title('Iv')
            
            save(strcat(Folder,'IV up to',num2str(Vf),'V','.mat'),'V_Bias','I_SOT')
            
            I_Bias=V_Bias./Rbias; %[A]
            I_Shunt=I_Bias-I_SOT; %[A]
            Point=floor(Steps/4);
            Rp=I_Shunt(Point)/I_SOT(Point)*Rshunt;
            Isotmax=max(I_SOT.*1e6);
        end
    end
end