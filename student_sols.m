%NO_PFILE
function [funs, student_id] = student_sols()
%STUDENT_SOLS Contains all student solutions to problems.

% ----------------------------------------
%               STEP 1
% ----------------------------------------
% Set to your birthdate / the birthdate of one member in the group.
% Should a numeric value of format YYYYMMDD, e.g.
% student_id = 19900101;
% This value must be correct in order to generate a valid secret key.
student_id = 19981205;


% ----------------------------------------
%               STEP 2
% ----------------------------------------
% Your task is to implement the following skeleton functions.
% You are free to use any of the utility functions located in the same
% directory as this file as well as any of the standard matlab functions.


pos=load('hip2.mat');


    function h = gen_filter()
        f_max=0.05; %Hz
        f_block=0.1; %Hz

        dt = 1;% s
        f_s = 1/dt; %Hz

        N = 60; %filter order

        f = [0,f_max,f_block,f_s/2]/(f_s/2); %*pi rad/sample
        a = [0,    1,      0,    0].*f*pi*f_s; %filter amplitude

        h = firpm(N,f,a,'differentiator'); %TODO: This line is missing some code!
    
        figure
        stem(0:N, h);
        
        figure
        [H,w]=freqz(h,1,512);
        plot(f,a,w/pi,abs(H));
        legend('ideal','FIR design');
        
        
        
        Nsamples=500;
        h_euler=[1/dt, -1/dt];
        vel.observation=conv(pos.noisy_position, h_euler)*3.6;
        vel.observation=vel.observation(1:Nsamples);
        vel.signal=conv(pos.true_position, h_euler)*3.6;
        vel.signal=vel.signal(1:Nsamples);
        t_euler=(0:Nsamples-1)/dt;
        
        vel.fir_obs=conv(pos.noisy_position,h)*3.6;
        vel.fir_sig=conv(pos.true_position,h)*3.6;
        t_filtered=(0:length(vel.fir_obs)-1)/dt;
        
        lc=lines(6);
        display(lc)
        
        figure
        xlabel 'Time (s)', ylabel 'Velocity (km/h)', hold on, grid on
        axis([0 600 0 220]);
        plot(t_euler,vel.signal,'-','Color',lc(1,:));
        plot(t_filtered,vel.fir_obs,'-','Color',lc(2,:));
        plot(t_filtered,vel.fir_sig,'--','Color',lc(3,:));
        legend('Euler filter (true)','FIR filter (measured)','FIR filter (true)');
        
        figure
        xlabel 'Time (s)', ylabel 'Velocity (km/h)', hold on, grid on
        axis([0 600 0 220]);
        plot(t_euler,vel.signal,'-','Color',lc(1,:));
        plot(t_filtered(N/2:end)-N/2,vel.fir_obs(N/2:end),'-','Color',lc(2,:));
        plot(t_filtered(N/2:end)-N/2,vel.fir_sig(N/2:end),'--','Color',lc(3,:));
        legend('Euler filter (true)','FIR filter (measured)','FIR filter (true)');
        
        figure
        xlabel 'Time (s)', ylabel 'Velocity (km/h)', hold on, grid on;
        plot(t_euler,vel.observation,'Color',lc(1,:));
        plot(t_euler,vel.signal,'--','Color',lc(2,:));
        legend('Measured', 'True');
    end

funs.gen_filter = @gen_filter;

% This file will return a structure with handles to the functions you have
% implemented. You can call them if you wish, for example:
% funs = student_sols();
% some_output = funs.some_function(some_input);

end

