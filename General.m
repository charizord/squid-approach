classdef General
    
    properties (SetAccess = public)
        Computer=General.Check_Computer();
    end
    
    methods(Static)
 

        %%     Instruments
        
        function Connect()
            Gen=General;
            A=readtable(fullfile(cd(),'\Metadata\Instruments.xlsx'),'sheet',Gen.Computer);
            for i=1:length(A.SerialNumber)
                if ~isempty(A.SerialNumber{i})
                    try
                        instr_name=A.SerialNumber{i};
                        Device=A.Device{i};
                        General.Close_open_Instrument(instr_name,Device)
                    catch er
                        disp([A.Device{i},A.SerialNumber{i}])
                        disp(er.message)
                    end
                end
            end

        end
        
        function  Close_open_Instrument(instr_name,Device,close_open)
            % open or delete n instruments.
            % input: 1 X n cell array with strings as instruments names.
            % and an optional boolean for for opening or closing mode. 
            
            Gen=General;
            if ~ exist('close_open','var') % set opening mode as default
                close_open = 1;
            end
            
            instruments = instrfind;
            
            % delete instrument if open 
            for j=1:length(instruments)
                if instruments(j).Name(end-4:end) == instr_name(end-4:end)
                   delete(instruments(j))
                end
            end
            

            % if in opening mode - open instrument
            if close_open
                switch instr_name(1:4)
                    case 'GPIB'    % For GPIB connection
                        gpib_num = abs(str2double(instr_name(end-1:end)));  
                        switch Gen.Computer
                            case '1.5K'
                                Instrument = gpib('ni',0,gpib_num);
                            case '4K'
                                Instrument = gpib('Agilent',7,gpib_num);
                        end
                        
                        switch Device
                            case 'Magnet_PS'
                                fopen(Instrument);
                                Magnet.Setup(); 
                            case 'Keithley'
                                set(Instrument, 'InputBufferSize', 12000);
                                fopen(Instrument); 
                            case 'SRS'
                                set(Instrument, 'InputBufferSize', 360000);
                                fopen(Instrument);
                            otherwise
                                fopen(Instrument);
                        end

                    case 'Seri' % For Serial connection (DAC & Arduino)
                        switch length(instr_name)
                            case 11
                                Instrument = serial(instr_name(end-3:end)); 
                            case 12
                                Instrument = serial(instr_name(end-4:end)); 
                        end
                        fopen(Instrument);
                        if strcmp(Device,'DAC') % If DAC - connect DAC
                            DAC.Setup();
                        end

                    case 'TCPI' %   For Nanonis 
                            switch Device
                                case 'Nanonis'
                                    Instrument = tcpip('192.168.236.1', 6501,'NetworkRole', 'client');
                                    set(Instrument, 'InputBufferSize', 2^15);
                                    fopen(Instrument);

                                case 'Magnet'
                                    Ip='192.168.236.2';
                                    Port=4444;
                                    Instrument =tcpip(Ip, Port,'NetworkRole', 'client');
                                    fopen(Instrument);
                                    Magnet.Setup(); 
                            end
                        
                    case 'VISA' %   For Lock-in & Keithley 2450

                        switch Device
                            case 'Keithley'
                                switch Gen.Computer
                                    case '1.5K'
                                        % Add more Kiet
                                        switch instr_name
                                            case 'VISA-04499542'
                                                model_code = convertCharsToStrings(instr_name(6:end));
                                                Instrument = visa('NI',strcat('USB0::0x05E6::0x6500::',model_code,'::0::INSTR'));
                                            case 'VISA-04509845'
                                                model_code = convertCharsToStrings(instr_name(6:end));
                                                Instrument = visa('NI',strcat('USB0::0x05E6::0x2460::',model_code,'::0::INSTR'));
                                            case 'VISA-04499369'
                                                model_code = convertCharsToStrings(instr_name(6:end));
                                                Instrument = visa('NI',strcat('USB0::0x05E6::0x6500::',model_code,'::0::INSTR'));                                               
                                        end


                                    case '4K' % Think
                                        
                                        
                                        switch instr_name
                                            case 'VISA-04499369'
                                                model_code = convertCharsToStrings(instr_name(6:end));
                                                Instrument = visa('NI',strcat('USB0::0x05E6::0x6500::',model_code,'::0::INSTR'));
                                        end
 
                                    case '0.3K'
                                        model_code = convertCharsToStrings(instr_name(end-7:end));
                                        Instrument = visa('NI', strcat('USB0::0x05E6::0x2450::',model_code,'::0::INSTR'));
                                
                                end
                            case 'Lock_in'
                                Instrument = visa('NI', 'USB0::0xB506::0x2000::002982::0::INSTR');
%                             case 'Oschillator'
%                                  model_code = convertCharsToStrings(instr_name(6:end));
%                                  Instrument = visa('NI',strcat('USB0::0x2A8D::0x0396::CN',model_code,'::0::INSTR'));
                        end
                        fopen(Instrument);
                end
            end
        end 
        function  instruments = Close_open_instr(instr_names,close_open)
            % open or delete n instruments.
            % input: 1 X n cell array with strings as instruments names.
            % and an optional boolean for for opening or closing mode. 
            % output: 
            % opening mode: 1 X n intrument array.
            % closing mode: empty array
            % instruments = General.Close_open_instr({'Serial-COM10','GPIB7-1','TCPIP-192.168.236.1','GPIB7-23','GPIB7-24','GPIB7-26','GPIB7-27'});
            % instruments = General.Close_open_instr({'Serial-COM10','GPIB7-1','TCPIP-192.168.236.1','GPIB7-23','GPIB7-24','GPIB7-26','GPIB7-27','GPIB7-22','GPIB7-23','GPIB7-8','VISA982-0'});
            % set opening mode as default
            Gen=General;
            if ~ exist('close_open','var')
                close_open = 1;
            end
            % delete all open instruments 
            for i = 1:length(instr_names)
                instr_name = instr_names{i};
                instruments = instrfind;
                for j=1:length(instruments)
                    if instruments(j).Name(end-4:end) == instr_name(end-4:end)
                       delete(instruments(j))
                    end
                end
            % only on opening mode - open all instruments
                if close_open
                    
                    switch instr_name(1:4)
                        case 'GPIB'    % For GPIB connection
                            gpib_num = abs(str2double(instr_name(end-1:end))); 
                            switch Gen.Computer
                                case '1.5K'
                                    Instrument = gpib('ni',0,gpib_num);
                                case '4K'
                                	Instrument = gpib('Agilent',7,gpib_num);
                            end
                            
                            if gpib_num>1
                                set(Instrument, 'InputBufferSize', 12000);
                                fopen(Instrument);
                            else
                                fopen(Instrument);
                            	Magnet.Setup(); 
                            end
                            
                        case 'Seri' % For Serial connection (DAC)
                            switch length(instr_name)
                                case 11
                                    Instrument = serial(instr_name(end-3:end)); 
                                case 12
                                	Instrument = serial(instr_name(end-4:end)); 
                            end
                            fopen(Instrument);

                            if or(instr_name(end-1:end)=='10',instr_name(end:end)=='4')
                                switch Gen.Computer
                                    case '1.5K'
                                        DAC.Setup();
                                    case '4K'
                                        DAC.Setup();
                                end
                            end
                            
                        case 'TCPI' %   For Nanonis 
                            switch instr_name
                                case 'Nanonis'
                                    Instrument = tcpip('192.168.236.1', 6501,'NetworkRole', 'client');
                                    set(Instrument, 'InputBufferSize', 2^15);
                                case 'Magnet'
                                    Ip='192.168.236.2';
                                    Port=4444;
                                    Instrument =tcpip(Ip, Port,'NetworkRole', 'client');
                            end
                            fopen(Instrument);
            

                        case 'VISA'
                            Instrument = visa('NI', 'USB0::0xB506::0x2000::002982::0::INSTR');
                            fopen(Instrument);
                    end
                end
            end
            % set the output
            instruments = instrfind;
        end 
        
        function  instrument = Find_instr(instr_name)
            Instr_Find=instrfind;
            for i = 1:length(Instr_Find)
                if  strcmp(Instr_Find(i).Name(end-4:end),instr_name(end-4:end))  % Smarter way
                    instrument = Instr_Find(i);
                end
            end
        end 
        
        %%  Time
        
        function [time_cell,timedata] = Time_cell(tdat)
            
            if ~ exist('tdat','var')
                 tdat = datetime;
            end
            str_time = datestr(tdat);
            year = str_time(8:11);
            mon = str_time(4:6);
            day = str_time(1:2);
            hourmin = [str_time(13:14),'.',str_time(16:17)];
            num_mon = num2str(month(tdat));
            time_cell = {year,mon,day,hourmin,num_mon};
            if nargout == 2
                 timedata = tdat;
            end
        end
        
        %%  Save
        
        function File_Make(file_name,titles)  
        % Create file with titels    
            if ~isfile(file_name)
                [file_path,~]=fileparts(file_name); %sdgdsf
                % If the folder doesn't exist it creates it 
                if ~ exist(file_path,'dir')
                    mkdir(file_path)
                end
                writecell(titles,file_name,'Delimiter',' ');
            end
        end
        
        function Save(file_name,data)
        % Add vector for existing file    
            writecell(data,file_name,'WriteMode','append','Delimiter',' ');
        end
        
        function Save_SOT_Characterization(keit_num,Time,V_l_lim,V_u_lim,V_steps,H_l_lim,H_u_lim,H_steps) 
            % extracting the file path 
      
            %% saving the data
            inf_filename = General.Add_Time2File([Time,' Info'])
            inf_filename(end-3:end) = '.inf';
            % Make file
            
            fileID = fopen(inf_filename,'w');
            fprintf(fileID,'%12s\r\n','This is an I-V measurement at various magnetic fields.');
            fprintf(fileID,'%6.6f %12s %12s\r\n',V_l_lim,'V','to');
            fprintf(fileID,'%6.6f %12s\r\n',V_u_lim,'V');
            fprintf(fileID,'%6d %12s\r\n',V_steps,'steps');
            fprintf(fileID,'%6.6f %12s %12s\r\n',H_l_lim,'T','to');
            fprintf(fileID,'%6.6f %12s\r\n',H_u_lim,'T');
            fprintf(fileID,'%6d %12s\r\n',H_steps,'steps');
            fprintf(fileID,'%6d %12s\r\n',keit_num,'keit_num');
            
            
            fclose(fileID);
            
        end

        function Save_SJ_Thermal_Characterization(keit_num,Time,V_l_lim,V_u_lim,V_steps,V_Heater_l_lim,V_Heater_u_lim,V_Heater_steps) 
            % extracting the file path 
      
            %% saving the data
            inf_filename = General.Add_Time2File([Time,' Info'])
            inf_filename(end-3:end) = '.inf';
            % Make file
            
            fileID = fopen(inf_filename,'w');
            fprintf(fileID,'%12s\r\n','This is an I-V measurement at various Temperatures.');
            fprintf(fileID,'%6.6f %12s %12s\r\n',V_l_lim,'V','to');
            fprintf(fileID,'%6.6f %12s\r\n',V_u_lim,'V');
            fprintf(fileID,'%6d %12s\r\n',V_steps,'steps');
            fprintf(fileID,'%6.6f %12s %12s\r\n',V_Heater_l_lim,'V','to');
            fprintf(fileID,'%6.6f %12s\r\n',V_Heater_u_lim,'V');
            fprintf(fileID,'%6d %12s\r\n',V_Heater_steps,'steps');
            fprintf(fileID,'%6d %12s\r\n',keit_num,'keit_num');
            
            
            fclose(fileID);
            
        end
        
        function [File_Name] = Add_Time2File(fileName_finale)
            % Enter a full dir of a file name with txt end
            Gen=General;
            time_cell=General.Time_cell;
            switch Gen.Computer
                case '0.3K'
                    File_path=['C:\Users\owner\Google Drive\0.3K microscope\data','\',time_cell{1},'\',time_cell{2},'\',time_cell{3}];
                
                case '1.5K'
                    File_path=['D:\Google Drive\1.5K microscope\data','\',time_cell{1},'\',time_cell{2},'\',time_cell{3}];
                case '4K'
                    File_path=['C:\Users\Owner\Google Drive\4K microscope\data','\',time_cell{1},'\',time_cell{2},'\',time_cell{3}];
                case 'VUP'
                    File_path=['C:\Users\Danziger\Google Drive\Vup\data','\',time_cell{1},'\',time_cell{2},'\',time_cell{3}];
            end
            File_Name=fullfile(File_path,[fileName_finale,'.txt']);   
        end
        
        %%
        % Response needs to be in SOT
        
        function SetMail()
            %sendmail( mail address , title , contnent, file location); % Send mail and file example 
            mail = 'qil@mail.huji.ac.il';
            psswd = 'Zaq12wsx';
            host = 'smtp.gmail.com';
            port  = '465';
            m_subject = 'subject';
            m_text = 'test';
            setpref( 'Internet','E_mail', mail );
            setpref( 'Internet', 'SMTP_Server', host );
            setpref( 'Internet', 'SMTP_Username', mail );
            setpref( 'Internet', 'SMTP_Password', psswd );
            props = java.lang.System.getProperties;
            props.setProperty( 'mail.smtp.user', mail );
            props.setProperty( 'mail.smtp.host', host );
            props.setProperty( 'mail.smtp.port', port );
            props.setProperty( 'mail.smtp.starttls.enable', 'true' );
            props.setProperty( 'mail.smtp.debug', 'true' );
            props.setProperty( 'mail.smtp.auth', 'true' );
            props.setProperty( 'mail.smtp.socketFactory.port', port );
            props.setProperty( 'mail.smtp.socketFactory.class', 'javax.net.ssl.SSLSocketFactory' );
            props.setProperty( 'mail.smtp.socketFactory.fallback', 'false' );
        end

        function [Computer]=Check_Computer()
            [status, result] = dos('getmac'); % result should contain the MAC address of the current computer
%             result = '18-60-24-84-03-89';
            if contains(result,'18-60-24-84-03-89') == true
                Computer = '4K';
            elseif contains(result,"84-A9-3E-70-D4-3A") == true
                Computer = '1.5K';
            elseif contains(result,"84-A9-3E-70-D5-59") == true
                Computer = '0.3K';
            elseif contains(result,"34-17-EB-98-3D-17") == true
                Computer = 'VUP';
            else
                disp('not able to recognize the computer')
            end
%             switch cd
%                 case 'C:\Users\owner\Dropbox\Scripts'       % Find better conditions, maybe IP
%                     Computer='1.5K';
%                 case 'C:\Users\Owner\Dropbox\Scripts'                           % Find better conditions
%                     Computer='4K';
%                 case 'C:\Users\owner\Google Drive\Scripts'
%                     Computer='0.3K';
%                 case 'C:\Users\Danziger\Dropbox\Scripts'
%                     Computer='VUP';
%             end
        end
        %%    Temperature
        
        function [T] = Volt2Kelvin(V)
            disp('a')
            switch General.Check_Computer()
                case '1.5K'% Load the 1.5K calibration
                    load(fullfile(cd(),'Metadata','Temp_calibration_1_5K.mat'));   % New directory
                    Temp=Temp_Calibration.Temp;
                    Resistance=Temp_Calibration.Res;
                    R_Meas=V/(100e-6);
                    T=interp1(Resistance,Temp,R_Meas);
                    
                case '4K' % Load the 4K calibration
                    
                    load(fullfile(cd(),'Metadata','Temp_calebration_4K.mat'));   % New directory
                    Temp=Temp_calebration{1};
                    Volt=Temp_calebration{2};
                    C=Temp_calebration{3};

                    n=find(Volt>V,1); 
                    if isempty(n) || n==1
                        T=-1;
                        return
                    end
            
                    dV=Volt(n)-Volt(n-1);
                    dT=Temp(n)-Temp(n-1);
                    dX=V-Volt(n-1);
                    S=[Temp(n-1),dT/dV-dV*(2*C(n-1)+C(n))/6,C(n-1)/2,(C(n)-C(n-1))/(6*dV)];
                    T=S(1)+S(2).*dX+S(3)*(dX.^2)+S(4)*(dX.^3);
            end
        end
        
        function [Temperature] = Temp_Meas(Instrument,IP_Chan,OffSet)
            %% measure temperature
            %   Instrument='DAC','Nanonis','Keithley'
            switch Instrument
                case 'DAC'
                    Voltage = DAC.Get(IP_Chan);
                case 'Nanonis'
                    Voltage = Nanonis.Get(IP_Chan);
                case 'Keithley'
                    Voltage = Keithley.Get(IP_Chan)-OffSet;
            end
            Temperature=General.Volt2Kelvin(Voltage);
        end
        
        
        function Beep()
            
            x=1:1000;
            y=sin(x);
            sound(y/10,2000)
            pause(1)
            sound(y/10,2000)
            pause(1)
            sound(y/10,2000)   
            pause(3)
        end
        
    end
end



