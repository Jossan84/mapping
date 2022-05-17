function [ t, e, nIter ] = newtonMethod(t0,err ,f,f1,f2)

        t = t0;
        nIter = 0;
        e = 100;
            while e>err                                    
                ti    = t - (f1(t)/f2(t));   % Iteration

                e     = abs(((ti-t)/t)*100); % Error [%]
                t     = ti;
                nIter = nIter + 1;           % Nº iterations
            end
        
        %% Plot 
        steps = 1000;
        t_plot = linspace(0,1,steps);
        
        for i=1:steps
            resF(i)  = f(t_plot(i));
            resF1(i) = f1(t_plot(i));
            resF2(i) = f2(t_plot(i));
        end
        
        figure;
        ax(1) = subplot(3,1,1);
        plot(t_plot,resF,'b.-');
        hold on;
        plot(t,f(t),'*r');
        title('f');
        xlabel('t');
        grid on
        hold off
        
        
        ax(2) = subplot(3,1,2);
        plot(t_plot,resF1,'b.-');
        hold on;
        plot(t,f1(t),'*r');
        title('f1');
        xlabel('t');
        grid on
        hold off
        
        ax(3) = subplot(3,1,3);
        plot(t_plot,resF2,'b.-');
        hold on
        plot(t,f2(t),'*r');
        title('f2');
        xlabel('t');
        grid on
        hold off
        
        linkaxes(ax,'x');
        

end

