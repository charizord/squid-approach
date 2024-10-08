classdef Overnight
    
	properties (SetAccess = public)
        Computer=General.Check_Computer();
    end
   methods(Static)

function Approach(Vbias,scan_step, Sigma_Multiplay,Mu_avg,Vfinal,Rfb,VSteps)
%Overnight Approach conducts a controlled, automated approach, tests tip crash (by comparing responses and Ivs after every stepper change),
%tests for false alarm and noise. Uses Nanonis.Poke Retracts by 1um when fshift triggered (TF).
%
%Input arguments:
%Vbias - working Vbias; 
%scan step - length of scanner step for Poke; 
%Sigma_Multiplay - Sigma threshold for Poke; Mu_avg - Number of samples for Poke Mu; 
%Vfinal - endpoint for Iv range; 
%Rfb- Rfeedback for Iv ;
%Vsteps- Voltage step increments for IV test;

%Algorthm:
% Approach_loop 


%%
%standard values for error tolerance
Max_Vfb_delta=0.01; %biggest error between Vfb and baseline Vfb permitted before a crash is indicated.
IV_RMSE_threshold=0.01; %max RMSE between baseline IV and IV conducted at tip crash test
effective_zero=0.01; %effective zero for IV liearity test
%Header values to compare to:

%Vfeedback:
Base_V_feedback=Nanonis.Get_Average(2,10); %Base V feedback to compare to throughout Poke iterations.
%Baseline I-V
[Base_Vb,Base_I]=Overnight.Iv('DAC',0,0,Vfinal,Rfb,VSteps);

%%
%resetting DAC, Scanner
DAC.Set(0,Vbias)
Nanonis.Set_Scanner_Z(0);
t=0;

%Poke While loop. t can be 0 (limit of scanner reached), -100 (potential crash detected using Vfb comparison),
%or the value of the Zscanner when fshift triggers
while t==0
    t=Overnight.Poke(0.2,15,3,scan_step,Sigma_Multiplay,Mu_avg,Base_V_feedback,Max_Vfb_delta);
    close all
    if t==0
        Nanonis.Set_Scanner_Z(0);
        Nanonis.Motors_Move_Steps(60,'z-'); %Notice axis direction
    
    elseif t==-100 %testing for crash
        DAC.Blink(2)
        %double checking Vfb
        v2fb=Nanonis.Get(2);
        [New_Vb,New_I]=Overnight.Iv('DAC',0,0,Vfinal,Rfb,VSteps);
        RMSE=sqrt(mean((New_I-Base_I).^2));
            %compare if IV are identical
            if RMSE>IV_RMSE_threshold %not identical
                    %check if new Iv is linear
                    if any(abs(diff(diff(y)))>effective_zero) %linear test to see if any 2nd derivative is non zero 
                    disp("IV linear - Tip Crashed!");
                    break
                    else%if not linear - soft crash
                    disp("IV not identical. Iv non linear. Possible soft crash");
                    break
                    end
            else %IV base and IV new are identical. Carry on.
                 Base_V_feedback=v2fb; %newbaseline for V feedback
                 t=0;%returns to start of while loop
            end
                
    
    else%stopped but no crash suspected - maybe you are done!
    %false alarm test - poke again
        Nanonis.Set_Scanner_Z(0);
        j=Overnight.Poke(1,15,3,scan_step,Sigma_Multiplay,Mu_avg,Base_V_feedback,Max_Vfb_delta);
        close all
        if (j<t+0.1 && j>t-0.1) %second poke gives same result. Congratulations, you have arrived at your destination!
            stp=Nanonis.Get_Encoder_Z();
            scn=t;
            disp("Loop complete. Location:  Stepper: " +string(stp) + " Scanner: " +string(scn) +". Retracted by 1 um.");  
            break
        elseif j==0 %false alarm. Back to the start
                t=0;
        else  %stopped twice, but not at same distance - error
            stp=Nanonis.Get_Encoder_Z();
            disp("Error: fshift threshold crossed twice at 2 different Scanner positions. Stepper: "+string(stp)+" Scanner positions: " +string(t) +" and " +string(j));
            break
        end
    end
end
end


%%
%creating a new poke function
function touch_point=Poke(retract,zlimit,avg,step,Sigma_Multiplay,Mu_avg, Base_V_feedback,delta_err)

%             Poke approach until freq. shift crosses threshold. return
%             touch point and retract
            if ~exist('zlimit','var')
                zlimit=15;
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
            
            BV=Base_V_feedback;
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
%             figure;
%             title('TF poke approach')
%             xlabel('time(s)')
%             ylabel('freq. shift')
%             u=yline(upp_lim,'-.r');
%             l=yline(low_lim,'-.r');
%             h=animatedline;
%             tic;

            while 1
                Zs=Nanonis.Get_Scanner_Z();
                
%                 perform averaging on freq. shift value
                f_avg=zeros(1,avg);
                for i=1:avg
                    f_avg(i)=Nanonis.Getfshift();
                    pause(0.1)
                end
                f=mean(f_avg);
%                 addpoints(h,double(toc),double(f));
                
                %this part for the overnight: compare Vfb to baseline BV
                %taken at start of overnight loop
                vfb=Nanonis.Get(2);
                delta_vfb=abs(BV-vfb);
                
                condition1=or(f>upp_lim, f<low_lim);
                condition2=Zs>=zlimit;
                condition3=delta_vfb>delta_err;  %difference in Vfb between base and measured bigger than err
                if condition1
                    Nanonis.Set_Scanner_Z(Zs-retract);
                    disp('crossed freq. shift threshold, retracted '+string(retract*5)+'um , touch-point='+string(Zs));
                   for i=1:4
                    a=Nanonis.Get_Scanner_Z();
                    Nanonis.Set.Scanner_Z(a-retract)
                    pause(0.5)
                   end
                    touch_point=Zs;
                    break
                end
                if condition2
                    Nanonis.Set_Scanner_Z(Zs-retract);
                    touch_point=0;
                    break
                end
                if condition3  %potential broken tip, sends
                    Nanonis.Set_Scanner_Z(Zs-retract);
                    touch_point=-100;
                    disp('Potenial Crash, testing. Scanner position: ' +string(Zs))
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
                 
                Nanonis.Set_Scanner_Z(Zs+ step);
                pause(0.01)
            end
            Nanonis.PllOff();
            Nanonis.OutputOff();
            end


    
%%
%modified Iv
function [Vb,Isot]=Iv(Instrument,Output_Chan,Input_Chan,Vf,Rfb,Steps)
            SOT=Squid_Functions;
            switch SOT.Computer
                case '1.5K'
                    Turn_Ratio=3.74;   % Update when know or take from file
                    Rshunt=1; %[Ohm] % Per system
                case '4K'
                    Turn_Ratio=2;
                    Rshunt=1; %[Ohm] % Per system
            end

%             if ~ exist('Rbias','var')
%                 Rbias=2550; %[Ohm]
%             end
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
            
            Vb=V_Bias;
            Isot=I_SOT;
 
%             % Plot

%             
%             IV=figure(666);
%             set(IV,'Name','IV','NumberTitle','off')
%             hold on
%             grid on
%             plot(V_Bias,I_SOT.*1e6)
%             xlabel('V_b_i_a_s [V]')
%             ylabel('I_{SOT} [\muA]')
%             title('Iv')
%             
%             I_Bias=V_Bias./Rbias; %[A]
%             I_Shunt=I_Bias-I_SOT; %[A]
%             Point=floor(Steps/4);
%             Rp=I_Shunt(Point)/I_SOT(Point)*Rshunt;
%             Isotmax=max(I_SOT.*1e6);
        end
   end
end
